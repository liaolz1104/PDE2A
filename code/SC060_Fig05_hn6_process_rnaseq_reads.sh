#!/usr/bin/env bash
set -euo pipefail

ENV_DIR=/data/liaolz/envs/hn6_rnaseq
BASE=/data/liaolz/PDE2A
RAW="$BASE/data/HN6_KD_OE_raw_fastq/CleanData"
REF="$BASE/ref/gencode_v44"
OUT="$BASE/results/HN6_KD_OE_FASTQ_bridge_20260805"
MANIFEST="$OUT/sample_manifest.tsv"
THREADS=20

export PATH="$ENV_DIR/bin:$PATH"
mkdir -p "$OUT"/{fastqc,star_index,star,featurecounts,logs}

mapfile -t FASTQS < <(tail -n +2 "$MANIFEST" | awk -F '\t' -v root="$RAW" '{print root "/" $4; print root "/" $5}')
if [[ $(find "$OUT/fastqc" -maxdepth 1 -name '*_fastqc.zip' | wc -l) -ne 24 ]]; then
  fastqc --threads 12 --outdir "$OUT/fastqc" "${FASTQS[@]}" 2>&1 | tee "$OUT/logs/fastqc.log"
fi
if [[ ! -s "$OUT/fastqc/multiqc_report.html" ]]; then
  multiqc --force --outdir "$OUT/fastqc" "$OUT/fastqc" 2>&1 | tee "$OUT/logs/multiqc.log"
fi

if [[ ! -s "$OUT/QC_APPROVED_NO_TRIM" ]]; then
  echo "QC reports generated. Inspect them before alignment; create QC_APPROVED_NO_TRIM only if trimming is not indicated." >&2
  exit 2
fi

if [[ ! -s "$REF/GRCh38.primary_assembly.genome.fa" ]]; then
  gzip -dc "$REF/GRCh38.primary_assembly.genome.fa.gz" > "$REF/GRCh38.primary_assembly.genome.fa"
fi
if [[ ! -s "$REF/gencode.v44.annotation.gtf" ]]; then
  gzip -dc "$REF/gencode.v44.annotation.gtf.gz" > "$REF/gencode.v44.annotation.gtf"
fi
if [[ ! -s "$REF/gencode.v44.gene_id_to_name.tsv" ]]; then
  awk 'BEGIN{FS="\t"; OFS="\t"; print "gene_id","gene_name"}
    $3=="gene" {
      if (match($9, /gene_id "[^"]+"/)) {gid=substr($9,RSTART+9,RLENGTH-10)} else {next}
      if (match($9, /gene_name "[^"]+"/)) {gname=substr($9,RSTART+11,RLENGTH-12)} else {next}
      print gid,gname
    }' "$REF/gencode.v44.annotation.gtf" > "$REF/gencode.v44.gene_id_to_name.tsv"
fi

if [[ ! -s "$OUT/star_index/Genome" ]]; then
  STAR --runThreadN "$THREADS" --runMode genomeGenerate \
    --genomeDir "$OUT/star_index" \
    --genomeFastaFiles "$REF/GRCh38.primary_assembly.genome.fa" \
    --sjdbGTFfile "$REF/gencode.v44.annotation.gtf" \
    --sjdbOverhang 149 2>&1 | tee "$OUT/logs/star_index.log"
fi

align_one() {
  local sample="$1" read1="$2" read2="$3"
  sample_out="$OUT/star/$sample"
  mkdir -p "$sample_out"
  if [[ ! -s "$sample_out/Aligned.sortedByCoord.out.bam" ]]; then
    STAR --runThreadN "$THREADS" \
      --genomeDir "$OUT/star_index" \
      --readFilesIn "$RAW/$read1" "$RAW/$read2" \
      --readFilesCommand zcat \
      --outFileNamePrefix "$sample_out/" \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outSAMattributes Standard > "$OUT/logs/star_${sample}.log" 2>&1
    samtools index -@ 8 "$sample_out/Aligned.sortedByCoord.out.bam"
  fi
  echo "STAR complete: $sample"
}

while IFS=$'\t' read -r sample condition comparison read1 read2; do
  [[ "$sample" == "sample" ]] && continue
  align_one "$sample" "$read1" "$read2" &
  while [[ $(jobs -pr | wc -l) -ge 4 ]]; do
    wait -n
  done
done < "$MANIFEST"
wait

mapfile -t BAMS < <(tail -n +2 "$MANIFEST" | awk -F '\t' -v out="$OUT" '{print out "/star/" $1 "/Aligned.sortedByCoord.out.bam"}')
for strand in 0 1 2; do
  featureCounts -T "$THREADS" -p --countReadPairs -B -C -s "$strand" \
    -a "$REF/gencode.v44.annotation.gtf" -g gene_id \
    -o "$OUT/featurecounts/gene_counts_s${strand}.txt" "${BAMS[@]}" \
    2>&1 | tee "$OUT/logs/featureCounts_s${strand}.log"
done

echo "Primary processing complete: $(date -Is)" | tee "$OUT/logs/primary_processing.complete"
