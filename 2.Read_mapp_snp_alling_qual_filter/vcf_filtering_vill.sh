#!/bin/bash
#SBATCH --job-name=vcf_filtering_all_sites_vill
#SBATCH --output=vcf_filtering_all_sites_vill_%j.out
#SBATCH --error=vcf_filtering_all_sites_vill_%j.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

module load bcftools
module load htslib

mypath=/shared/home/toufiktakhi/neurosporagbs/

# 1. Remove the 4 low-coverage individuals (^ removes them)
# 2. Drop the FORMAT/GP annotation globally
# 3. Set individual genotypes to missing (./.) if DP < 3
# 4. Filter variant sites but PROTECT invariant sites
# 5. Remove any site missing in more than 20% of the population
bcftools view -s ^Vill.B5BIS,Vill.C5.3,Vill.C4.1,Vill.C7.1 ${mypath}Vill-A1-3_all_sites_raw.vcf.gz | \
bcftools annotate -x FORMAT/GP | \
bcftools filter -S . -e 'FORMAT/DP < 3' | \
bcftools view -e 'TYPE!="ref" && (QUAL < 30 || MQ < 30)' -Oz -o ${mypath}Vill-A1-3_all_sites_filtered.vcf.gz

# Index the final output
tabix -p vcf ${mypath}Vill-A1-3_all_sites_filtered.vcf.gz

echo "Filtering done."
echo "Remaining sites :"
bcftools view -H ${mypath}Vill-A1-3_all_sites_filtered.vcf.gz | wc -l
