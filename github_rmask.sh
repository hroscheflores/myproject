#!/usr/bin/env bash
#SBATCH -A r00110
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH -N 1
#SBATCH -t 0-16:00
#SBATCH -p general
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hrosche@iu.edu

set -euo pipefail

echo "[${SLURM_JOB_ID:-NOJOB}] Running RepeatMasker on brem1_de_novo_prefixed (pangenome) at $(date)"

RepeatMasker \
  -gff \
  -pa 32 \
  -lib "/N/project/bsf/hector/delivery_from_sean/repeat_library_4_samples_with_unknowns/ref_INhi_895_Italian_consensi.fa.classified" \
  -dir "/N/scratch/hrosche/pangenome/TE/he_et_al/brem1_de_novo_prefixed_pangenome_library" \
  "/N/project/bsf/hector/pangenome/input_genomes/sanitized/brem1_de_novo_prefixed.fa"

echo "[${SLURM_JOB_ID:-NOJOB}] Finished brem1_de_novo_prefixed (pangenome) at $(date)"
