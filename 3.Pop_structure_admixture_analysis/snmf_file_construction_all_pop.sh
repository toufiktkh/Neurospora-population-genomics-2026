#!/bin/bash
#SBATCH --job-name=vcf_construction_snmf
#SBATCH --partition=fast
#SBATCH --cpus-per-task=8
#SBATCH --mem=30GB
#SBATCH --time=02:00:00
#SBATCH --output=snmf_vcf_construction_%j.log

module load plink/1.90b6.18

# Using the generated plink pruned files to extract the vcf 
# ------------------------------------For DomeFire------------------------------------
#EDF1
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_variants_cleaned \
    --allow-extra-chr \
    --extract ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_pruned_data.prune.in \
    --recode vcf-iid \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_snps_only_snmf
#EDF3
plink --bfile ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_filtered_snps \
    --allow-extra-chr \
    --extract ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_pruned_data.prune.in \
    --recode vcf-iid \
    --out ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_snps_only_snmf

# ------------------------------------For Vill------------------------------------
plink --bfile ~/neurosporagbs/pca_splitstree_vill/all_sites/Vill_variants_cleaned \
    --allow-extra-chr \
    --extract ~/neurosporagbs/pca_splitstree_vill/all_sites/Vill_pruned_data.prune.in \
    --recode vcf-iid \
    --out ~/neurosporagbs/pca_splitstree_vill/all_sites/Vill_snps_only_snmf

# ------------------------------------For LGM------------------------------------
plink --bfile ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_variants_cleaned \
    --allow-extra-chr \
    --extract ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_pruned_data.prune.in \
    --recode vcf-iid \
    --out ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_snps_only_snmf
