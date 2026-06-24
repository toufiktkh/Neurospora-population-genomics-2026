library(LEA)
library(ggplot2)
library(tidyr)

setwd("/shared/home/toufiktakhi/neurosporagbs/")

# 1. Getting best K for NEUROSPORA SP. DF1 (K = 3)
cat("Extraction of data for DF1 (K = 3)...\n")
project_df1  <- load.snmfProject("/shared/home/toufiktakhi/neurosporagbs/admixture_analysis/DF1_snps_only_clone_corr.snmfProject")
df1_samples  <- system("bcftools query -l ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF1_snps_only_clone_corr_snmf.vcf.gz", intern = TRUE)

df1_transect <- gsub(".*_(T\\d+)S.*", "\\1", df1_samples)
df1_transect[!grepl("^T\\d+$", df1_transect)] <- "Ref"

k1 <- 3
best_run1 <- which.min(cross.entropy(project_df1, K = k1))
Q_matrix1 <- Q(project_df1, K = k1, run = best_run1)
colnames(Q_matrix1) <- paste0("Pop", 1:k1)

df1_k <- as.data.frame(Q_matrix1)
df1_k$sample <- df1_samples
df1_k$transect <- df1_transect

df1_multi_K_df <- pivot_longer(df1_k,
                               cols = starts_with("Pop"),
                               names_to = "population",
                               values_to = "proportion")

df1_multi_K_df <- df1_multi_K_df[!is.na(df1_multi_K_df$proportion), ]

df1_multi_K_df$population <- factor(df1_multi_K_df$population, levels = paste0("Pop", 1:k1))
sample_order_df1 <- df1_samples[order(df1_transect)]
df1_multi_K_df$sample <- factor(df1_multi_K_df$sample, levels = unique(sample_order_df1))

# Graph for DF1
cat("Generating the graph for DF1...\n")
p_df1 <- ggplot(df1_multi_K_df, aes(x = sample, y = proportion, fill = population)) + 
    geom_bar(stat = "identity", width = 1) +
    facet_grid(. ~ transect, scales = "free_x", space = "free_x") + 
    scale_fill_brewer(palette = "Set1", name = "Ancestral\nClusters", na.value = "transparent") +
    labs(title = "sNMF Ancestry Proportions of Neurospora sp. (E-DF1) at K = 3",
         x = "Isolates (grouped by Sample tree)", y = "Ancestry Proportion") + 
    theme_minimal() +
    theme(
        axis.text.x = element_text(size = 8, face = "bold", angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(face = "bold", size = 9), # Top labels will now be T1, T2... and Ref
        panel.spacing.x = unit(0.1, "lines"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )

ggsave("~/neurosporagbs/snmf/DF1_sNMF_bestK.pdf", p_df1, width = 16, height = 6, dpi = 300) 

# 2. Getting best K for N. DISCRETA PS4B DF3 (K = 1)
cat("Extraction of data for DF3 (K = 1)...\n")
project_df3  <- load.snmfProject("/shared/home/toufiktakhi/neurosporagbs/admixture_analysis/DF3_snps_only_clone_corr.snmfProject")
df3_samples  <- system("bcftools query -l ~/neurosporagbs/pca_splitstree_domefire/all_sites/DF3_snps_only_clone_corr.vcf.gz", intern = TRUE)

df3_transect <- gsub(".*_(T\\d+)S.*", "\\1", df3_samples)
df3_transect[!grepl("^T\\d+$", df3_transect)] <- "Ref"

k3 <- 1
best_run3 <- which.min(cross.entropy(project_df3, K = k3))
Q_matrix3 <- Q(project_df3, K = k3, run = best_run3)
colnames(Q_matrix3) <- "Pop1"

df3_k <- as.data.frame(Q_matrix3)
df3_k$sample <- df3_samples
df3_k$transect <- df3_transect

df3_multi_K_df <- pivot_longer(df3_k,
                               cols = starts_with("Pop"),
                               names_to = "population",
                               values_to = "proportion")

df3_multi_K_df <- df3_multi_K_df[!is.na(df3_multi_K_df$proportion), ]

df3_multi_K_df$population <- factor(df3_multi_K_df$population, levels = "Pop1")
sample_order_df3 <- df3_samples[order(df3_transect)]
df3_multi_K_df$sample <- factor(df3_multi_K_df$sample, levels = unique(sample_order_df3))

# Graph for DF3
cat("Generating the graph for DF3...\n")
p_df3 <- ggplot(df3_multi_K_df, aes(x = sample, y = proportion, fill = population)) + 
    geom_bar(stat = "identity", width = 1) +
    facet_grid(. ~ transect, scales = "free_x", space = "free_x") + 
    scale_fill_brewer(palette = "Set1", name = "Ancestral\nClusters", na.value = "transparent") +
    labs(title = "sNMF Ancestry Proportions of N. discreta PS4B (E-DF3) at K = 1",
         x = "Isolates (grouped by Sample tree)", y = "Ancestry Proportion") + 
    theme_minimal() +
    theme(
        axis.text.x = element_text(size = 8, face = "bold", angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(face = "bold", size = 9),
        panel.spacing.x = unit(0.1, "lines"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )

ggsave("~/neurosporagbs/snmf/DF3_sNMF_bestK.pdf", p_df3, width = 10, height = 6, dpi = 300) 

cat("Finished successfuly ! The coherent graphs are saved in snmf/\n")
