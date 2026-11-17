#!/usr/bin/env bash
set -euo pipefail
source ./config.env

if [[ "${FEATURECOUNTS_STRANDNESS}" != "0" &&
      "${FEATURECOUNTS_STRANDNESS}" != "1" &&
      "${FEATURECOUNTS_STRANDNESS}" != "2" ]]; then
    echo "ERROR: strandedness has not been confirmed."
    echo "Run RSeQC infer_experiment first and set FEATURECOUNTS_STRANDNESS to 0, 1, or 2."
    exit 1
fi

mkdir -p "${OUTDIR}/04_featurecounts"

bams=()
while IFS=$'\t' read -r sample original patient region region_label batch include r1w r2w r1 r2; do
  bams+=("${OUTDIR}/02_star/${sample}/${sample}.Aligned.sortedByCoord.out.bam")
done < <(tail -n +2 "${MANIFEST}")

featureCounts \
  -T "${THREADS}" \
  -p --countReadPairs \
  -B -C \
  -s "${FEATURECOUNTS_STRANDNESS}" \
  -a "${GTF}" \
  -g gene_id \
  -t exon \
  -o "${OUTDIR}/04_featurecounts/featureCounts_gene_id.txt" \
  "${bams[@]}"

echo "featureCounts finished."
