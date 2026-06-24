#!/bin/bash
#SBATCH --job-name=extract_cleaned_dp_mq_qual
#SBATCH --output=extract_cleaned_dp_mq_%j.out
#SBATCH --error=extract_cleaned_dp_mq_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G

module load bcftools
module load htslib

mypath="/shared/home/toufiktakhi/neurosporagbs/"
out_base="${mypath}dp_mq_all_sites_domefire/"

echo "Extracting stats for E-DF1..."

# 1- E-DF1 
# Generating a text file listing the individuals in the exact order
bcftools query -l ${mypath}E-DF1_all_sites_filtered.vcf.gz > ${out_base}readmapping_stats_DF1/all_sites_sample_columns.txt
echo "Saved sample column order to sample_columns.txt"

#extraction
bcftools query -f '%CHROM\t%POS[\t%DP]\n' \
    ${mypath}E-DF1_all_sites_filtered.vcf.gz | \
    gzip > ${out_base}readmapping_stats_DF1/all_sites_DP.txt.gz

bcftools query -f '%CHROM\t%POS\t%MQ\n' \
    ${mypath}E-DF1_all_sites_filtered.vcf.gz | \
    gzip > ${out_base}readmapping_stats_DF1/all_sites_MQ.txt.gz

bcftools query -f '%QUAL\n' \
    ${mypath}E-DF1_all_sites_filtered.vcf.gz | \
    gzip > ${out_base}readmapping_stats_DF1/all_sites_QUAL.txt.gz


echo "Extracting stats for E-DF3..."

# 2 E-DF3 Extraction
bcftools query -l ${mypath}E-DF3_all_sites_filtered.vcf.gz > ${out_base}readmapping_stats_DF3/all_sites_sample_columns.txt
echo "Saved sample column order to sample_columns.txt"

#extraction
bcftools query -f '%CHROM\t%POS[\t%DP]\n' \
    ${mypath}E-DF3_all_sites_filtered.vcf.gz | \
    gzip > ${out_base}readmapping_stats_DF3/all_sites_DP.txt.gz

bcftools query -f '%CHROM\t%POS\t%MQ\n' \
    ${mypath}E-DF3_all_sites_filtered.vcf.gz | \
    gzip > ${out_base}readmapping_stats_DF3/all_sites_MQ.txt.gz

bcftools query -f '%QUAL\n' \
    ${mypath}E-DF3_all_sites_filtered.vcf.gz | \
    gzip > ${out_base}readmapping_stats_DF3/all_sites_QUAL.txt.gz


echo "Extraction complete for both DF1 and DF3"
