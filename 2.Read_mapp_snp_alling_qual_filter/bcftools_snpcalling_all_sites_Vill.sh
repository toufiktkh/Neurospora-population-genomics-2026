#!/bin/bash
#SBATCH --job-name=Vill_snpcalling_all_sites
#SBATCH --output=snpcalling_Vill_all_sites_%j.out
#SBATCH --error=snpcalling_Vill_all_sites_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=40G

module load bcftools
module load htslib

reference=$1
bams=$2
output=$3

bcftools mpileup \
    -f $reference \
    -d 1000 \
    -Q 20 \
    -a AD,DP \
    --threads 30 \
    ${bams}*.bam | \
bcftools call \
    -m \
    -f GQ,GP \
    --ploidy 1 \
    --threads 30 \
    -Oz \
    -o ${output}.gz

tabix -p vcf ${output}.gz
