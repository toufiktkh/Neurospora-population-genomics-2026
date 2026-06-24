#!/bin/bash
#SBATCH --job-name=variant_del
#SBATCH --output=variant_del_%j.out
#SBATCH --error=variant_del_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=64G

module load bcftools
module load htslib

mypath=~/neurosporagbs/test_12_06_26/
THREADS=30

# DOMEFIRE — use the best-mapped / subset outputs

# --- E-DF1 all sites (SNPs + invariants) ---
bcftools view --threads $THREADS --exclude-types indels,mnps,other \
    ${mypath}E-DF1_best_all_sites_filtered.vcf.gz | \
bcftools norm --threads $THREADS -d all \
    -Oz -o ${mypath}E-DF1_snps_invariants_dedup.vcf.gz
tabix -p vcf ${mypath}E-DF1_snps_invariants_dedup.vcf.gz

# --- E-DF1 variants only (SNPs only) ---
bcftools view --threads $THREADS --exclude-types indels,mnps,other \
    ${mypath}E-DF1_variants_filtered.vcf.gz | \
bcftools norm --threads $THREADS -d all \
    -Oz -o ${mypath}E-DF1_snps_only_dedup.vcf.gz
tabix -p vcf ${mypath}E-DF1_snps_only_dedup.vcf.gz

# --- E-DF3 all sites ---
bcftools view --threads $THREADS --exclude-types indels,mnps,other \
    ${mypath}E-DF3_best_all_sites_filtered.vcf.gz | \
bcftools norm --threads $THREADS -d all \
    -Oz -o ${mypath}E-DF3_snps_invariants_dedup.vcf.gz
tabix -p vcf ${mypath}E-DF3_snps_invariants_dedup.vcf.gz

# --- E-DF3 variants only ---
bcftools view --threads $THREADS --exclude-types indels,mnps,other \
    ${mypath}E-DF3_variants_filtered.vcf.gz | \
bcftools norm --threads $THREADS -d all \
    -Oz -o ${mypath}E-DF3_snps_only_dedup.vcf.gz
tabix -p vcf ${mypath}E-DF3_snps_only_dedup.vcf.gz




# VILLEVEYRAC

# --- all sites ---
bcftools view --threads $THREADS --exclude-types indels,mnps,other \
    ${mypath}Vill-A1-3_all_sites_filtered.vcf.gz | \
bcftools norm --threads $THREADS -d all \
    -Oz -o ${mypath}Vill-A1-3_snps_invariants_dedup.vcf.gz
tabix -p vcf ${mypath}Vill-A1-3_snps_invariants_dedup.vcf.gz

# --- variants only ---
bcftools view --threads $THREADS --exclude-types indels,mnps,other \
    ${mypath}Vill-A1-3_variants_only_filtered.vcf.gz | \
bcftools norm --threads $THREADS -d all \
    -Oz -o ${mypath}Vill-A1-3_snps_only_dedup.vcf.gz
tabix -p vcf ${mypath}Vill-A1-3_snps_only_dedup.vcf.gz



# LA GRANDE MOTTE

# --- all sites ---
bcftools view --threads $THREADS --exclude-types indels,mnps,other \
    ${mypath}LGM_all_sites_filtered.vcf.gz | \
bcftools norm --threads $THREADS -d all \
    -Oz -o ${mypath}LGM_snps_invariants_dedup.vcf.gz
tabix -p vcf ${mypath}LGM_snps_invariants_dedup.vcf.gz

# --- variants only ---
bcftools view --threads $THREADS --exclude-types indels,mnps,other \
    ${mypath}LGM_variants_only_filtered.vcf.gz | \
bcftools norm --threads $THREADS -d all \
    -Oz -o ${mypath}LGM_snps_only_dedup.vcf.gz
tabix -p vcf ${mypath}LGM_snps_only_dedup.vcf.gz

echo "Done. Site counts:"
for f in E-DF1 E-DF3 Vill-A1-3 LGM; do
    echo "$f all-sites dedup: $(bcftools view -H ${mypath}${f}_snps_invariants_dedup.vcf.gz | wc -l)"
    echo "$f SNPs-only dedup: $(bcftools view -H ${mypath}${f}_snps_only_dedup.vcf.gz | wc -l)"
done
