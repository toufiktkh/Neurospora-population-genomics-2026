#!/bin/bash
#SBATCH --job-name=plot_LGM
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=plot_LGM_%j.out
#SBATCH --error=plot_LGM_%j.err

module load r
module load bcftools
module load htslib

Rscript vcf_stats_all_sites_plot_LGM.R
