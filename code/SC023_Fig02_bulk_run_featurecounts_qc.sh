#!/usr/bin/env bash
set -euo pipefail

# Checkpoint 4 only: featureCounts and unsupervised count-level QC.
source ./config.env

OUT4="${OUTDIR}/04_featurecounts"
FIG_FC="${OUTDIR}/07_figures/QC/featureCounts"
FIG_COUNT="${OUTDIR}/07_figures/QC/counts"
mkdir -p "${OUT4}" "${FIG_FC}" "${FIG_COUNT}" logs

[[ "${FEATURECOUNTS_STRANDNESS}" == "2" ]] || { echo "ERROR: FEATURECOUNTS_STRANDNESS must be 2"; exit 2; }
[[ -f "${GTF}" && -s "${GTF}" ]] || { echo "ERROR: frozen GTF missing: ${GTF}"; exit 2; }
[[ -f "${MANIFEST}" && -s "${MANIFEST}" ]] || { echo "ERROR: manifest missing: ${MANIFEST}"; exit 2; }
[[ -f config.env.before_featurecounts ]] || { echo "ERROR: config backup missing"; exit 2; }
[[ $(awk 'END{print NR-1}' "${MANIFEST}") -eq 12 ]] || { echo "ERROR: manifest must contain 12 samples"; exit 2; }

bam_count=0
while IFS=$'\t' read -r sample _; do
  [[ "${sample}" == "sample_id" ]] && continue
  bam="${OUTDIR}/02_star/${sample}/${sample}.Aligned.sortedByCoord.out.bam"
  bai="${bam}.bai"
  [[ -s "${bam}" && -s "${bai}" ]] || { echo "ERROR: BAM/BAI missing for ${sample}"; exit 3; }
  bam_count=$((bam_count + 1))
done < "${MANIFEST}"
[[ "${bam_count}" -eq 12 ]] || { echo "ERROR: expected 12 BAMs"; exit 3; }

if [[ -e "${OUT4}/FEATURECOUNTS_CONFIG_AUDIT.md" ]]; then
  echo "ERROR: refusing to overwrite existing Checkpoint 4 audit outputs"
  exit 4
fi
fc_version=$(featureCounts -v 2>&1 | awk '/featureCounts v/{print $2}')
cat > "${OUT4}/FEATURECOUNTS_CONFIG_AUDIT.md" <<EOF
# featureCounts configuration audit

- RSeQC classification: reverse/antisense stranded
- featureCounts parameter: -s 2
- Reason: 12/12 RSeQC concordant
- featureCounts version: ${fc_version}
- GTF: ${GTF}
- Samples/BAMs: 12/12, ordered by the frozen manifest
- Paired-fragment parameters: -p --countReadPairs -B -C
- Feature aggregation: -t exon -g gene_id
EOF

if [[ ! -s "${OUT4}/featureCounts_gene_id.txt" || ! -s "${OUT4}/featureCounts_gene_id.txt.summary" ]]; then
  bash scripts/04_featurecounts.sh
fi

/usr/bin/python3 scripts/11_checkpoint4_prepare.py \
  --featurecounts "${OUT4}/featureCounts_gene_id.txt" \
  --summary "${OUT4}/featureCounts_gene_id.txt.summary" \
  --manifest "${MANIFEST}" --gtf "${GTF}" --outdir "${OUT4}"

Rscript scripts/12_checkpoint4_qc.R \
  "${OUT4}/count_matrix_gene_id.tsv" \
  "${OUT4}/metadata_all12.tsv" \
  "${OUT4}/FEATURECOUNTS_RUN_AUDIT.tsv" \
  "${OUT4}" "${OUTDIR}/07_figures/QC"

echo "Checkpoint 4 calculations finished. STOP: no differential expression was run."
