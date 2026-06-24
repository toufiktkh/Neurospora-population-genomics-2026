#!/bin/bash
#SBATCH --job-name=genetic_distances
#SBATCH --output=calculate_genetic_dist_%j.out
#SBATCH --error=calculate_genetic_dist_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G


module load plink/1.90b6.18

# 1. DF1
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_variants_cleaned --distance square 1-ibs --allow-extra-chr --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_genetic_distances

# 2. DF3
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_filtered_snps --distance square 1-ibs --allow-extra-chr --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_genetic_distances

# 3. Villeveyrac
plink --bfile ~/neurosporagbs/pca_splitstree_vill/all_sites/Vill_variants_cleaned --distance square 1-ibs --allow-extra-chr --out ~/neurosporagbs/pca_splitstree_vill/all_sites/Vill_genetic_distances

# 4. LGM
plink --bfile ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_variants_cleaned --distance square 1-ibs --allow-extra-chr --out ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_genetic_distances
