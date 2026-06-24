library(poppr)
library(vcfR)
library(ggplot2)

setwd("/shared/home/toufiktakhi/neurosporagbs/")

pdf("histograms_clone_correction_snps_only.pdf", width = 8, height = 5)

process_site <- function(vcf_file, site_name, regex_pattern, threshold = 0.05) {
  cat("\n--- Processing site:", site_name, "with threshold:", threshold, "---\n")
  
  #VCF data to "genind" format
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  vcf_genind <- vcfR2genind(vcf)
  
  # Convert to genclone
  obj <- as.genclone(vcf_genind)
  saveRDS(obj, paste0(site_name, "_genclone.rds"))
  cat("Genclone object saved to:", paste0(site_name, "_genclone.rds"), "\n")
  
  #Assign strata
  samples <- indNames(obj)
  strata_labels <- gsub(regex_pattern, "\\1", samples)
  strata_labels[strata_labels == samples] <- "Other"
  strata(obj) <- data.frame(group = strata_labels)
  
  # Calculate genetic distance and NORMALIZE
  dist_matrix <- diss.dist(obj) / nLoc(obj)
  
  # Generate the histogram using ggplot2
  dist_df <- data.frame(distance = as.vector(dist_matrix))
  
  p <- ggplot(dist_df, aes(x = distance)) +
    geom_histogram(bins = 80, fill = "#4A6984", color = "#FFFFFF", size = 0.2, alpha = 0.9) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "#C15C3D", size = 0.7) +
    annotate("text", x = threshold, y = Inf, label = paste0("Threshold (", threshold*100, "%)"), 
             hjust = -0.1, vjust = 2, color = "#C15C3D", fontface = "italic", size = 3.5) +
    labs(
      title = paste("Pairwise genetic distances of the", site_name, "population"),
      x = "Proportion of different SNPs (0 = identical, 1 = different)",
      y = "Count of pairwise comparisons"
    ) +
    theme_classic(base_size = 11) + 
    theme(
      text = element_text(family = "sans", color = "#222222"),
      plot.title = element_text(face = "bold", size = 12, color = "#111111", margin = margin(b=10)),
      axis.line = element_line(color = "#444444", size = 0.4),
      axis.ticks = element_line(color = "#444444", size = 0.4),
      panel.grid.major.y = element_line(color = "#F0F0F0", size = 0.5) 
    )
    
  print(p)
  
  # Identify and Filter MLGs
  cat("Total individuals:", nInd(obj), "\n")
  cat("Strict MLGs (100% identity):", mlg(obj, quiet = TRUE), "\n")
  
  mlg.filter(obj, distance = dist_matrix) <- threshold
  cat("MLLs after filtering (threshold =", threshold, "):", nmll(obj), "\n")
  
  # Clone Correction
  obj_cc <- clonecorrect(obj, strata = ~group)
  kept_samples <- indNames(obj_cc)
  cat(site_name, ":", nInd(obj), "initial ->", length(kept_samples), "after correction\n")
  
  #Export the list of kept individuals
  output_txt <- paste0(site_name, "_snps_only_samples_to_keep.txt")
  write.table(kept_samples, output_txt, row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("Sample list saved to:", output_txt, "\n")
  
  return(kept_samples)
}

#Data Execution

# DF1
process_site("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_snps_only_snmf.vcf.gz", 
             "E-DF1", ".*_(T\\d+)S.*", threshold = 0.03)

# DF3
process_site("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_snps_only_snmf.vcf.gz", 
             "E-DF3", ".*_(T\\d+)S.*", threshold = 0.06)

# Villeveyrac
process_site("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_vill/all_sites/Vill_snps_only_snmf.vcf.gz", 
             "Vill", "^Vill[.-]([A-Z][0-9]+).*", threshold = 0.01)

# La Grande Motte
process_site("/shared/home/toufiktakhi/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_snps_only_snmf.vcf.gz", 
             "LGM", "^(LGM\\d+)\\..*", threshold = 0.02)

dev.off()
cat("\nDone. Check histograms_clone_correction_snps_only.pdf\n")
