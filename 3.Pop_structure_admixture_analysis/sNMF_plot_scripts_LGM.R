library(LEA)
library(ggplot2)
library(tidyr)

setwd("/shared/home/toufiktakhi/neurosporagbs/")

# 1. Getting best K for NEUROSPORA SP. GM (K = 3)
cat("Extraction and alignment of data for GM (K=3)...\n")
project_gm  <- load.snmfProject("/shared/home/toufiktakhi/neurosporagbs/admixture_analysis/LGM_snps_only_clone_corr.snmfProject")
gm_samples  <- system("bcftools query -l ~/neurosporagbs/pca_splitstree_LGM/all_sites/LGM_snps_only_clone_corr.vcf.gz", intern = TRUE)

# Extract site names (e.g., LGM2.1 -> LGM2)
gm_transect <- gsub("^(LGM\\d+)\\..*", "\\1", gm_samples)
gm_transect[!grepl("^LGM\\d+$", gm_transect)] <- "Ref"

k <- 3
best_run <- which.min(cross.entropy(project_gm, K = k))
Q_matrix <- Q(project_gm, K = k, run = best_run)
colnames(Q_matrix) <- paste0("Pop", 1:k)

gm_k <- as.data.frame(Q_matrix)
gm_k$sample <- gm_samples
gm_k$transect <- gm_transect

gm_multi_K_df <- pivot_longer(gm_k,
                              cols = starts_with("Pop"),
                              names_to = "population",
                              values_to = "proportion")

gm_multi_K_df <- gm_multi_K_df[!is.na(gm_multi_K_df$proportion), ]

# Fixate properly the levels of factors
gm_multi_K_df$population <- factor(gm_multi_K_df$population, levels = paste0("Pop", 1:k))
sample_order_gm <- gm_samples[order(gm_transect)]
gm_multi_K_df$sample <- factor(gm_multi_K_df$sample, levels = unique(sample_order_gm))

# Graph for GM (Flipped Horizontally)
cat("Generating the aligned graph for GM...\n")
p_gm <- ggplot(gm_multi_K_df, aes(x = sample, y = proportion, fill = population)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(. ~ transect, scales = "free_x", space = "free_x") + # Flipped facet layout to columns
    scale_fill_brewer(palette = "Set1", name = "Ancestral\nClusters", na.value = "transparent") +
    labs(title = "sNMF Ancestry Proportions of Neurospora sp. from La Grande Motte site at K = 3",
         x = "Isolates (grouped by sampling tree)", y = "Ancestry Proportion") + 
    theme_minimal() +
    theme(
        axis.text.x = element_text(size = 9, face = "bold", angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(face = "bold", size = 9),
        panel.spacing.x = unit(0.1, "lines"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )


ggsave("~/neurosporagbs/snmf/LGM_sNMF_bestK.pdf", p_gm, width = 16, height = 6, dpi = 300)
cat("Done! Output saved to ~/neurosporagbs/snmf/LGM_sNMF_bestK.pdf\n")
