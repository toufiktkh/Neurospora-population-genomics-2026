#!/bin/bash
#SBATCH --job-name=plot_DF3
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=plot_DF3_%j.out
#SBATCH --error=plot_DF3_%j.err

module load r
module load bcftools
module load htslib

Rscript vcf_stats_all_sites_plot_DF3.R
