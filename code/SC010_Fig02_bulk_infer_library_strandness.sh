#!/usr/bin/env bash
set -euo pipefail
source ./config.env

if [[ ! -f "${BED12}" ]]; then
  echo "BED12 not found: ${BED12}"
  echo "Set BED12 in config.env before running RSeQC."
  exit 2
fi

mkdir -p "${OUTDIR}/03_rseqc"

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample original patient region region_label batch include r1w r2w r1 r2; do
    bam="${OUTDIR}/02_star/${sample}/${sample}.Aligned.sortedByCoord.out.bam"
    echo "[RSeQC infer_experiment] ${sample}"
    infer_experiment.py -r "${BED12}" -i "${bam}" \
      > "${OUTDIR}/03_rseqc/${sample}.infer_experiment.txt"
done

echo "Inspect results/03_rseqc/*.infer_experiment.txt."
echo "Then set FEATURECOUNTS_STRANDNESS=0/1/2 in config.env."
