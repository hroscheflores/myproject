#!/usr/bin/env bash

set -euo pipefail

GRAFFITE_DIR="path/to/GraffiTE"
SINGULARITY_IMAGE="path/to/graffite.sif"
NEXTFLOW_BIN="nextflow"

ASSEMBLIES_CSV="path/to/assemblies.csv"
TE_LIBRARY="path/to/repeat_library.fa.classified"
REFERENCE_GENOME="path/to/reference_genome.fasta"
READS_CSV="path/to/input_reads.csv"

GRAPH_METHOD="pangenie"
WORKDIR="graffite_run"

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "Starting GraffiTE run"
echo "Working directory: ${WORKDIR}"

"${NEXTFLOW_BIN}" run "${GRAFFITE_DIR}/main.nf" \
    -profile cluster \
    --assemblies "${ASSEMBLIES_CSV}" \
    --TE_library "${TE_LIBRARY}" \
    --reference "${REFERENCE_GENOME}" \
    --reads "${READS_CSV}" \
    --graph_method "${GRAPH_METHOD}" \
    -with-singularity "${SINGULARITY_IMAGE}" \
    -work-dir "${WORKDIR}"

echo "GraffiTE run complete"