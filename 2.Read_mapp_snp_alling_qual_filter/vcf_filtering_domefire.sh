#!/bin/bash
#SBATCH --job-name=vcf_filtering_all_sites_DF
#SBATCH --output=vcf_filtering_all_sites_DF_%j.out
#SBATCH --error=vcf_filtering_all_sites_DF_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

module load bcftools
module load htslib

mypath="/shared/home/toufiktakhi/neurosporagbs/"

# ----- For E-DF1 ------
bcftools filter -S . \
    -e 'FORMAT/DP < 3' \
    ${mypath}E-DF1_all_sites_raw.vcf.gz | \
bcftools view \
    -e 'TYPE!="ref" && (QUAL < 30 || MQ < 30)' \
    -Oz -o ${mypath}E-DF1_all_sites_filtered.vcf.gz

tabix -p vcf ${mypath}E-DF1_all_sites_filtered.vcf.gz
echo "E-DF1 done"

# ----- For E-DF3 ------
bcftools filter -S . \
    -e 'FORMAT/DP < 3' \
    ${mypath}E-DF3_all_sites_raw.vcf.gz | \
bcftools view \
    -e 'TYPE!="ref" && (QUAL < 30 || MQ < 30)' \
    -Oz -o ${mypath}E-DF3_all_sites_filtered.vcf.gz

tabix -p vcf ${mypath}E-DF3_all_sites_filtered.vcf.gz
echo "E-DF3 done"

echo "All filtering complete!"
