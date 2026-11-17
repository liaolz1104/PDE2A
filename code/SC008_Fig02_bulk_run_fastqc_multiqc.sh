#!/usr/bin/env bash
set -euo pipefail
source ./config.env

mkdir -p "${OUTDIR}/01_fastqc/raw_reports" "${OUTDIR}/01_fastqc/multiqc"

python3 scripts/00_manifest_preflight.py --manifest "${MANIFEST}" --path-column wsl

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample original patient region region_label batch include r1w r2w r1 r2; do
    echo "[FastQC] ${sample}"
    fastqc -t "${THREADS}" -o "${OUTDIR}/01_fastqc/raw_reports" "${r1}" "${r2}"
done

multiqc -o "${OUTDIR}/01_fastqc/multiqc" "${OUTDIR}/01_fastqc/raw_reports"

echo "FastQC/MultiQC finished."
echo "NOTE: all supplied FASTQs are under CleanData directories."
echo "Do not perform a second trimming step unless MultiQC shows a concrete adapter/quality problem."
