#!/usr/bin/env bash
set -euo pipefail
source ./config.env

mkdir -p "${OUTDIR}/02_star" "${OUTDIR}/02_star/logs"

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample original patient region region_label batch include r1w r2w r1 r2; do
    echo "[STAR] ${sample}"
    mkdir -p "${OUTDIR}/02_star/${sample}"
    STAR \
      --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX}" \
      --readFilesIn "${r1}" "${r2}" \
      --readFilesCommand zcat \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM MD \
      --quantMode GeneCounts \
      --outFileNamePrefix "${OUTDIR}/02_star/${sample}/${sample}."

    samtools index -@ "${THREADS}" \
      "${OUTDIR}/02_star/${sample}/${sample}.Aligned.sortedByCoord.out.bam"
done

echo "STAR alignment finished."
