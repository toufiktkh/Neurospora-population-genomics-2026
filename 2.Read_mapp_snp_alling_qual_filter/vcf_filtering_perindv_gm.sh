#!/bin/bash
#SBATCH --job-name=vcf_filtering_gm
#SBATCH --output=vcf_filtering_gm_%j.out
#SBATCH --error=vcf_filtering_gm_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G

module load bcftools
module load htslib

mypath=/shared/home/toufiktakhi/neurosporagbs/

# 1. Drop the FORMAT/GP annotation globally
# 2. Set individual genotypes to missing (./.) if DP < 3
# 3. Filter variant sites but keep invariant sites (TYPE!="ref")

bcftools annotate -x FORMAT/GP ${mypath}gm_neutre-22_all_sites_raw.vcf.gz | \
bcftools filter -S . -e 'FORMAT/DP < 3' | \
bcftools view -e 'TYPE!="ref" && (QUAL < 30 || MQ < 30)' -Oz -o ${mypath}gm_neutre-22_all_sites_filtered.vcf.gz

tabix -p vcf ${mypath}gm_neutre-22_all_sites_filtered.vcf.gz

echo "Filtering done for La Grande Motte"
echo "Individuals remaining:"
bcftools query -l ${mypath}gm_neutre-22_all_sites_filtered.vcf.gz | wc -l
echo "Sites remaining:"
bcftools view -H ${mypath}gm_neutre-22_all_sites_filtered.vcf.gz | wc -l
