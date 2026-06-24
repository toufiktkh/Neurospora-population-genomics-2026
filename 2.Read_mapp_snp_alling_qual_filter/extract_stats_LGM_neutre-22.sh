#!/bin/bash
#SBATCH --job-name=extract_dp_mq_gm
#SBATCH --output=extract_dp_gm_%j.out
#SBATCH --error=extract_dp_mq_gm_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G

module load bcftools
module load htslib

VCF_IN="$HOME/neurosporagbs/gm_neutre-22_all_sites_filtered.vcf.gz"
OUT_DIR="$HOME/neurosporagbs/dp_mq_all_sites_gm"

mkdir -p "${OUT_DIR}"

# Save sample order
bcftools query -l "${VCF_IN}" > "${OUT_DIR}/sample_columns.txt"

# QUAL per site
zcat "${VCF_IN}" | bcftools query -f "%QUAL\n" | gzip > "${OUT_DIR}/gm_neutre-22_QUAL.txt.gz"
echo "Qual extraction complete"

# MQ per site
bcftools query -f '%CHROM\t%POS\t%MQ\n' "${VCF_IN}" | gzip > "${OUT_DIR}/all_sites_MQ.txt.gz"
echo "MQ extraction complete."

# DP matrix (CHROM, POS, Indiv1_DP, Indiv2_DP, ...)
bcftools query -f '%CHROM\t%POS[\t%DP]\n' "${VCF_IN}" | gzip > "${OUT_DIR}/all_samples_DP_matrix.txt.gz"
echo "DP extraction complete."

# AD matrix
bcftools query -f '%CHROM\t%POS[\t%AD]\n' "${VCF_IN}" | gzip > "${OUT_DIR}/all_samples_AD_matrix.txt.gz"
echo "AD extraction complete."

echo "All extractions done."
