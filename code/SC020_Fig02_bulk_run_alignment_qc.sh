#!/usr/bin/env bash
set -euo pipefail

# Continue only the frozen technical QC gates. This script never calls featureCounts.
source ./config.env

PYTHON_BIN="/usr/bin/python3"
QC_SCRIPT="scripts/09_checkpoint23_qc.py"
STAR_DIR="${OUTDIR}/02_star"
RSEQC_DIR="${OUTDIR}/03_rseqc"
LOG_DIR="logs"
CONDA_EXE="/data/liaolz/miniconda3/bin/conda"
RSEQC_ENV="/data/liaolz/envs/pde2a_rseqc_env"

[[ -f "${MANIFEST}" ]] || { echo "ERROR: manifest missing: ${MANIFEST}"; exit 2; }
[[ -f "${QC_SCRIPT}" ]] || { echo "ERROR: QC script missing: ${QC_SCRIPT}"; exit 2; }
[[ -d "${STAR_DIR}" ]] || { echo "ERROR: STAR output directory missing: ${STAR_DIR}"; exit 2; }
[[ -f "${GTF}" ]] || { echo "ERROR: matching GTF missing: ${GTF}"; exit 2; }
mkdir -p "${STAR_DIR}" "${RSEQC_DIR}" "${OUTDIR}/07_figures/QC/STAR" "${OUTDIR}/07_figures/QC/RSeQC" "${LOG_DIR}"

echo "[$(date --iso-8601=seconds)] Monitoring the existing STAR queue."
while true; do
  "${PYTHON_BIN}" "${QC_SCRIPT}" status --manifest "${MANIFEST}" --star-dir "${STAR_DIR}"
  completed=$(awk -F '\t' 'NR>1 && $5=="COMPLETED" {n++} END {print n+0}' "${STAR_DIR}/STAR_RUN_STATUS.tsv")
  running=$(awk -F '\t' 'NR>1 && $5=="RUNNING" {n++} END {print n+0}' "${STAR_DIR}/STAR_RUN_STATUS.tsv")
  failed=$(awk -F '\t' 'NR>1 && $5=="FAILED" {n++} END {print n+0}' "${STAR_DIR}/STAR_RUN_STATUS.tsv")
  echo "[$(date --iso-8601=seconds)] STAR completed=${completed}/12 running=${running} failed=${failed}"
  if [[ "${failed}" -gt 0 ]]; then
    echo "ERROR: partial/failed STAR output detected; no automatic rerun is permitted."
    exit 3
  fi
  [[ "${completed}" -eq 12 ]] && break
  if [[ "${running}" -eq 0 ]] && ! pgrep -f '02_star_align.sh' >/dev/null; then
    echo "ERROR: STAR queue stopped before 12/12 completion."
    exit 4
  fi
  sleep 60
done

echo "[$(date --iso-8601=seconds)] Starting Checkpoint 2 BAM/STAR audit."
"${PYTHON_BIN}" "${QC_SCRIPT}" checkpoint2 \
  --manifest "${MANIFEST}" \
  --star-dir "${STAR_DIR}" \
  --figure-dir "${OUTDIR}/07_figures/QC/STAR"

grep -q '^\*\*Checkpoint 2 = PASS\.\*\*$' "${STAR_DIR}/QC_CHECKPOINT2_AUDIT.md" || {
  echo "ERROR: Checkpoint 2 did not pass; RSeQC remains blocked."
  exit 5
}

[[ -x "${CONDA_EXE}" ]] || { echo "ERROR: conda executable missing: ${CONDA_EXE}"; exit 6; }
if [[ ! -x "${RSEQC_ENV}/bin/infer_experiment.py" || ! -x "${RSEQC_ENV}/bin/gtfToGenePred" || ! -x "${RSEQC_ENV}/bin/genePredToBed" ]]; then
  echo "[$(date --iso-8601=seconds)] Creating isolated RSeQC/UCSC environment: ${RSEQC_ENV}"
  "${CONDA_EXE}" create -y -p "${RSEQC_ENV}" -c conda-forge -c bioconda \
    "rseqc=5.0.4" ucsc-gtftogenepred ucsc-genepredtobed
fi
INFER_EXPERIMENT="${RSEQC_ENV}/bin/infer_experiment.py"
GTF_TO_GENEPRED="${RSEQC_ENV}/bin/gtfToGenePred"
GENEPRED_TO_BED="${RSEQC_ENV}/bin/genePredToBed"

if [[ ! -s "${BED12}" ]]; then
  mkdir -p "$(dirname "${BED12}")"
  gene_pred="${BED12%.bed12}.genePred"
  echo "[$(date --iso-8601=seconds)] UCSC conversion: frozen GTF -> genePred -> BED12"
  "${GTF_TO_GENEPRED}" -genePredExt -geneNameAsName2 "${GTF}" "${gene_pred}.tmp"
  "${GENEPRED_TO_BED}" "${gene_pred}.tmp" "${BED12}.tmp"
  mv "${gene_pred}.tmp" "${gene_pred}"
  mv "${BED12}.tmp" "${BED12}"
fi
[[ -s "${BED12}" ]] || { echo "ERROR: BED12 unavailable after UCSC conversion: ${BED12}"; exit 7; }

audit_bam="${STAR_DIR}/P08_C/P08_C.Aligned.sortedByCoord.out.bam"
"${PYTHON_BIN}" "${QC_SCRIPT}" audit-bed12 \
  --bed12 "${BED12}" --bam "${audit_bam}" --output-dir "${RSEQC_DIR}/bed12_audit"
grep -q 'Audit status: PASS' "${RSEQC_DIR}/bed12_audit/BED12_AUDIT.md" || {
  echo "ERROR: BED12 structural/coordinate/chromosome audit failed."
  exit 8
}

# Small read-only smoke run proves that RSeQC can parse the BED12 before the 12-sample run.
"${INFER_EXPERIMENT}" -r "${BED12}" -i "${audit_bam}" -s 1000 \
  > "${RSEQC_DIR}/bed12_audit/RSeQC_BED12_READABILITY_SMOKE.txt.tmp"
mv "${RSEQC_DIR}/bed12_audit/RSeQC_BED12_READABILITY_SMOKE.txt.tmp" \
  "${RSEQC_DIR}/bed12_audit/RSeQC_BED12_READABILITY_SMOKE.txt"

echo "[$(date --iso-8601=seconds)] Running RSeQC infer_experiment for 12 BAMs."
while IFS=$'\t' read -r sample _; do
  [[ "${sample}" == "sample_id" ]] && continue
  bam="${STAR_DIR}/${sample}/${sample}.Aligned.sortedByCoord.out.bam"
  out="${RSEQC_DIR}/${sample}.infer_experiment.txt"
  [[ -s "${bam}" ]] || { echo "ERROR: BAM missing: ${bam}"; exit 8; }
  if [[ ! -s "${out}" ]]; then
    "${INFER_EXPERIMENT}" -r "${BED12}" -i "${bam}" > "${out}.tmp"
    mv "${out}.tmp" "${out}"
  fi
done < "${MANIFEST}"

"${PYTHON_BIN}" "${QC_SCRIPT}" checkpoint3 \
  --manifest "${MANIFEST}" \
  --rseqc-dir "${RSEQC_DIR}" \
  --figure-dir "${OUTDIR}/07_figures/QC/RSeQC"

echo "[$(date --iso-8601=seconds)] Checkpoint 2-3 workflow stopped before featureCounts as required."
