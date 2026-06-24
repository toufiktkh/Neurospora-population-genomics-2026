#!/bin/bash
#SBATCH --job-name=plink_pca_filtering_DF
#SBATCH --output=plink_pca_filtering_DF_%j.out
#SBATCH --error=plink_pca_filtering_DF_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module load plink/1.90b6.18

module load plink/1.90b6.18
# E-DF1 ----
# Filtering on missing genotypes (Missing > 60%)
plink --vcf ~/neurosporagbs/E-DF1_snps_invariants_dedup.vcf.gz \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:#:\$1:\$2 \
    --mind 0.6 \
    --make-bed \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_samples_cleaned

# Filtering SNPs from the rest of the isoaltes
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_samples_cleaned \
    --allow-extra-chr \
    --geno 0.4 \
    --maf 0.05 \
    --make-bed \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_variants_cleaned

# Linkage Disequilibrium pruning
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_variants_cleaned \
    --allow-extra-chr \
    --indep-pairwise 50 5 0.5 \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_pruned_data

# Run the PCA on the independent SNPs
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_variants_cleaned \
    --allow-extra-chr \
    --extract ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_pruned_data.prune.in \
    --pca 10 \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_pca_final

#E-DF3 ----
## Filtering on missing genotypes (Missing > 80%) (Reason : only 9 isolates, and since it's an ALL sites VCF file, we have all invariant sites, which can contain a lot of missing data caused by GBS sequencing)
plink --vcf ~/neurosporagbs/E-DF3_snps_invariants_dedup.vcf.gz \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:#:\$1:\$2 \
    --mind 0.8 \
    --make-bed \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_samples_cleaned
# Filtering SNPs from the rest of the isoaltes
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_samples_cleaned \
    --allow-extra-chr \
    --geno 0.4 \
    --maf 0.05 \
    --make-bed \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_filtered_snps

# Linkage Disequilibrium pruning
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_filtered_snps \
    --allow-extra-chr \
    --indep-pairwise 50 5 0.5 \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_pruned_data

# Run the PCA on the independent SNPs
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_filtered_snps \
    --allow-extra-chr \
    --extract ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_pruned_data.prune.in \
    --pca 10 \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_pca




# ---------- SplitsTree ----------
# E-DF1
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_variants_cleaned  \
      --allow-extra-chr \
      --recode vcf-iid \
      --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_snps_splitstree

# Convert to Phylip / Fasta format
python3 vcf2phylip.py -i ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_snps_splitstree.vcf -n


# E-DF3
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_filtered_snps  \
      --allow-extra-chr \
      --recode vcf-iid \
      --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_snps_splitstree

# Convert to Phylip / Fasta format
python3 vcf2phylip.py -i ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_snps_splitstree.vcf -n


