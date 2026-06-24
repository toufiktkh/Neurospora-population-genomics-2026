#!/bin/bash
#SBATCH --job-name=filter_aln
#SBATCH --output=filter_aln_%j.out
#SBATCH --error=filter_aln_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G

module load trimal/1.4.1 
module load mafft/7.525

mypath=/shared/home/toufiktakhi/neurosporagbs/
path=${mypath}busco_single_copy_core_sequences_aligned_clean/
pathout=${mypath}busco_single_copy_core_sequences_aligned_clean_best/
stats=${mypath}busco_single_copy_core_sequences_aligned_clean_best.txt
fasta=${mypath}busco_single_copy_core_sequences_aligned_clean_best.fasta

rm -rf $pathout
mkdir $pathout

python3 ${mypath}filter_aln_qual.py \
    --path $path \
    --pathout $pathout \
    --stats $stats \
    --fasta $fasta
