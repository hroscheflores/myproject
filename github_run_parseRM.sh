#!/usr/bin/env bash

set -euo pipefail

SCRIPTNAME="path/to/parseRM.pl"
REPEAT_LIBRARY="path/to/repeat_library.fa.classified"
SAMPLE_ID="sample_id"

RM_DIR="path/to/repeatmasker_output_directory"
RM_OUTFILE="sample.fa.out"
FASTA_SRC="path/to/sample.fa"
FASTA_NAME="sample.fa"

OUTDIR="path/to/output_directory"

mkdir -p "${OUTDIR}"
cd "${RM_DIR}"

log() {
  echo "[parseRM][${SAMPLE_ID}] $*" >&2
}

log "Starting parseRM run"

if [[ ! -s "${FASTA_SRC}" ]]; then
  log "ERROR: FASTA missing or empty: ${FASTA_SRC}"
  exit 2
fi

if [[ ! -s "${RM_OUTFILE}" ]]; then
  log "ERROR: RepeatMasker .out file missing or empty: ${RM_OUTFILE}"
  exit 2
fi

if [[ ! -s "${REPEAT_LIBRARY}" ]]; then
  log "ERROR: repeat library missing or empty: ${REPEAT_LIBRARY}"
  exit 2
fi

if [[ ! -e "${FASTA_NAME}" ]]; then
  cp -p "${FASTA_SRC}" "${FASTA_NAME}"
fi

perl "${SCRIPTNAME}" -i "${RM_OUTFILE}" -p -f "${FASTA_NAME}" -r "${REPEAT_LIBRARY}" -v
perl "${SCRIPTNAME}" -i "${RM_OUTFILE}" -a 5,15 -k -v
perl "${SCRIPTNAME}" -i "${RM_OUTFILE}" -l 50,1 -k -v

for f in parseRM_* "${RM_OUTFILE}".*summary* "${RM_OUTFILE}".*landscape* "${RM_OUTFILE}".*bins*; do
  [[ -e "${f}" ]] && cp -p "${f}" "${OUTDIR}/"
done

log "Job complete"