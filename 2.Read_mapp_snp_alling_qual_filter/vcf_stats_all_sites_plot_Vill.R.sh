#!/bin/bash
#SBATCH --job-name=plot_VILL
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=plot_VILL_%j.out
#SBATCH --error=plot_VILL_%j.err

module load r
module load htslib
module load bcftools

Rscript vcf_stats_all_sites_plot_Vill.R
