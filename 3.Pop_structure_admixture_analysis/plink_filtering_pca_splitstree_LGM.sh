#!/bin/bash
#SBATCH --job-name=plink_pca_filtering_LGM
#SBATCH --output=plink_pca_filtering_LGM_%j.out
#SBATCH --error=plink_pca_filtering_LGM_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module load plink/1.90b6.18

# Filtering on missing genotypes (Missing > 60%)
plink --vcf ~/neurosporagbs/LGM_snps_invariants_dedup.vcf.gz \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:#:\$1:\$2 \
    --mind 0.6 \
    --make-bed \
    --out ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_samples_cleaned

# Filtering SNPs from the rest of the isolates
plink --bfile ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_samples_cleaned \
    --allow-extra-chr \
    --geno 0.4 \
    --maf 0.05 \
    --make-bed \
    --out ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_variants_cleaned

# Linkage Disequilibrium pruning
plink --bfile ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_variants_cleaned \
    --allow-extra-chr \
    --indep-pairwise 50 5 0.5 \
    --out ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_pruned_data

# Run the PCA on the independent SNPs
plink --bfile ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_variants_cleaned \
    --allow-extra-chr \
    --extract ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_pruned_data.prune.in \
    --pca 10 \
    --out ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_pca

# ---------- Splitstree ----------
# Convert PCA data back into a VCF
plink --bfile ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_variants_cleaned \
      --allow-extra-chr \
      --recode vcf-iid \
      --out ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_snps_splitstree

# Convert to Phylip / Fasta format
python3 vcf2phylip.py -i ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_snps_splitstree.vcf -n
