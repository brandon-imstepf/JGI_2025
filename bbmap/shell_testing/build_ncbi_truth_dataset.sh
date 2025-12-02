#!/bin/bash
set -euo pipefail

usage() {
cat <<'EOF'
build_ncbi_truth_dataset.sh
----------------------------------
Generate training TSVs using "NCBI reference is truth".
For each genome in the input folder, all CallGenes nofilter candidates are labeled:
  - Positive if they intersect the NCBI reference GFF
  - Negative otherwise
Every candidate is preserved (no subsampling). The script writes:
  1) A combined TSV (positives + negatives)
  2) Aggregate positives-only / negatives-only TSVs
  3) Per-genome TSVs under positives/ and negatives/ folders

Required args:
  in=<folder>         Folder containing *.fna.gz + *.gff.gz pairs

Optional args:
  outdir=<path>       Root output directory (default: ../../training_datasets)
  label=<name>        Dataset label used in filenames (default: basename of input)
  run_tag=<tag>       Custom run tag (default: timestamp)
  Any other tokens are passed directly to callgenes.sh
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CALLGENES="${REPO_DIR}/callgenes.sh"
GFFSETOP="${REPO_DIR}/gffsetop.sh"
GFF2TSV="${REPO_DIR}/gff2tsv.sh"

INFOLDER=""
OUTDIR="${REPO_DIR}/../training_datasets"
LABEL=""
RUN_TAG=""
PASSTHRU_ARGS=()

for arg in "$@"; do
    case "$arg" in
        in=*) INFOLDER="${arg#in=}";;
        outdir=*) OUTDIR="${arg#outdir=}";;
        label=*) LABEL="${arg#label=}";;
        run_tag=*) RUN_TAG="${arg#run_tag=}";;
        -h|--help) usage; exit 0;;
        *) PASSTHRU_ARGS+=("$arg");;
    esac
done

if [[ -z "$INFOLDER" ]]; then
    echo "ERROR: missing required arg in=<folder>" >&2
    usage
    exit 1
fi

if [[ ! -d "$INFOLDER" ]]; then
    echo "ERROR: input folder $INFOLDER does not exist" >&2
    exit 1
fi

if [[ -z "$LABEL" ]]; then
    LABEL="$(basename "$INFOLDER")"
fi
if [[ -z "$RUN_TAG" ]]; then
    RUN_TAG="$(date +%Y%m%d_%H%M%S)"
fi

TRUTH_MODE="ncbi_truth"
RUN_DIR="${OUTDIR}/${LABEL}_${TRUTH_MODE}_${RUN_TAG}"
COMBINED_TSV="${RUN_DIR}/${LABEL}.${TRUTH_MODE}.all.tsv.gz"
AGG_POS_TSV="${RUN_DIR}/${LABEL}.${TRUTH_MODE}.positives.tsv.gz"
AGG_NEG_TSV="${RUN_DIR}/${LABEL}.${TRUTH_MODE}.negatives.tsv.gz"
POS_DIR="${RUN_DIR}/positives"
NEG_DIR="${RUN_DIR}/negatives"

mkdir -p "${RUN_DIR}" "${POS_DIR}" "${NEG_DIR}"
rm -f "${COMBINED_TSV}" "${AGG_POS_TSV}" "${AGG_NEG_TSV}"

TMPDIR="$(mktemp -d)"
echo "Output directory: ${RUN_DIR}"
echo "Temp directory: ${TMPDIR}"

cleanup() {
    rm -rf "${TMPDIR}"
}
trap cleanup EXIT

for ref_fasta in "${INFOLDER}"/*.fna.gz; do
    [[ -e "$ref_fasta" ]] || continue
    base="$(basename "$ref_fasta" .fna.gz)"
    ref_gff="${INFOLDER}/${base}.gff.gz"
    if [[ ! -f "$ref_gff" ]]; then
        echo "[WARN] Missing reference GFF for ${base}, skipping"
        continue
    fi

    echo "--- Processing ${base} ---"
    ALL_GFF="${TMPDIR}/${base}.all.gff"
    POS_GFF="${TMPDIR}/${base}.positives.gff"
    NEG_GFF="${TMPDIR}/${base}.negatives.gff"

    if ! "${CALLGENES}" in="${ref_fasta}" outgff="${ALL_GFF}" truegenes="${ref_gff}" cds nofilter "${PASSTHRU_ARGS[@]}" >/dev/null; then
        echo "[ERROR] CallGenes nofilter run failed for ${base}; skipping genome" >&2
        continue
    fi

    if ! "${GFFSETOP}" in_a="${ALL_GFF}" in_b="${ref_gff}" out="${POS_GFF}" op=intersect >/dev/null; then
        echo "[ERROR] GffSetOperator intersect failed for ${base}; skipping genome" >&2
        continue
    fi
    if ! "${GFFSETOP}" in_a="${ALL_GFF}" in_b="${ref_gff}" out="${NEG_GFF}" op=subtract >/dev/null; then
        echo "[ERROR] GffSetOperator subtract failed for ${base}; skipping genome" >&2
        continue
    fi

    POS_TSV_TMP="${TMPDIR}/${base}.positives.tsv"
    NEG_TSV_TMP="${TMPDIR}/${base}.negatives.tsv"

    if ! "${GFF2TSV}" in="${POS_GFF}" out="${POS_TSV_TMP}" label=1 >/dev/null; then
        echo "[ERROR] gff2tsv failed for positives on ${base}; skipping genome" >&2
        continue
    fi
    if ! "${GFF2TSV}" in="${NEG_GFF}" out="${NEG_TSV_TMP}" label=0 >/dev/null; then
        echo "[ERROR] gff2tsv failed for negatives on ${base}; skipping genome" >&2
        continue
    fi

    if [[ -s "${POS_TSV_TMP}" ]]; then
        gzip -c "${POS_TSV_TMP}" >> "${AGG_POS_TSV}"
        gzip -c "${POS_TSV_TMP}" >> "${COMBINED_TSV}"
        gzip -c "${POS_TSV_TMP}" > "${POS_DIR}/${base}.positives.tsv.gz"
    else
        rm -f "${POS_TSV_TMP}"
    fi

    if [[ -s "${NEG_TSV_TMP}" ]]; then
        gzip -c "${NEG_TSV_TMP}" >> "${AGG_NEG_TSV}"
        gzip -c "${NEG_TSV_TMP}" >> "${COMBINED_TSV}"
        gzip -c "${NEG_TSV_TMP}" > "${NEG_DIR}/${base}.negatives.tsv.gz"
    else
        rm -f "${NEG_TSV_TMP}"
    fi

    rm -f "${ALL_GFF}" "${POS_GFF}" "${NEG_GFF}"
done

echo "--- Dataset build complete ---"
echo "Combined TSV : ${COMBINED_TSV}"
echo "Positives TSV: ${AGG_POS_TSV}"
echo "Negatives TSV: ${AGG_NEG_TSV}"
echo "Per-genome TSVs stored under ${POS_DIR} and ${NEG_DIR}"
