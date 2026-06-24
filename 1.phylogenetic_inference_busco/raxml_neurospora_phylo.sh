#!/bin/bash
#SBATCH --job-name=raxml_neurospora_phylo
#SBATCH --output=raxml_neurospora_phylo_%j.out
#SBATCH --error=raxml_neurospora_phylo_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G

module load raxml-ng/1.2.2


WORK_DIR="/shared/home/toufiktakhi/neurosporagbs"
MSA="${WORK_DIR}/busco_single_copy_core_sequences_aligned_clean_best.fasta"
PARTITIONS="${WORK_DIR}/busco_single_copy_core_sequences_aligned_clean_best.partitions.txt"
PREFIX="${WORK_DIR}/neurospora_phylo"

echo "Starting the construction of the RAxML-NG tree..."


raxml-ng \
  --all \
  --msa $MSA \
  --model $PARTITIONS \
  --bs-trees autoMRE \
  --threads auto{16} \
  --prefix $PREFIX

echo "Done calculating RAxML-NG."
