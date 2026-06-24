library(ggplot2)
library(data.table)

stats_dir    <- "/shared/home/toufiktakhi/neurosporagbs/test_12_06_26/dp_mq_all_sites/domefire/readmapping_stats_DF1/"
output_pdf   <- "/shared/home/toufiktakhi/neurosporagbs/test_12_06_26/dp_mq_all_sites/domefire/E-DF1_all_sites_stats_figures.pdf"
output_table <- "/shared/home/toufiktakhi/neurosporagbs/test_12_06_26/dp_mq_all_sites/domefire/E-DF1_all_sites_stats_summary.csv"
vcf_file     <- "/shared/home/toufiktakhi/neurosporagbs/test_12_06_26/E-DF1_all_sites_filtered.vcf.gz"

dp_file   <- file.path(stats_dir, "all_sites_DP.txt.gz")
mq_file   <- file.path(stats_dir, "all_sites_MQ.txt.gz")
qual_file <- file.path(stats_dir, "all_sites_QUAL.txt.gz")

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

cat("Extracting sample names...\n")
clean_samples <- system(paste("bcftools query -l", vcf_file), intern = TRUE)

cat("Reading Depth (DP) matrix...\n")
dp_mat <- fread(dp_file, header = FALSE, fill = TRUE, na.strings = ".")
colnames(dp_mat) <- c("CHROM", "POS", clean_samples)

# Mask DP < 3 as missing
for (col in clean_samples)
    set(dp_mat, j = col, value = replace(dp_mat[[col]], dp_mat[[col]] < 3, NA))

# SUMMARY TABLE - from FULL data
cat("Calculating summary stats from full data...\n")
summary_list <- list()
for (s in clean_samples) {
    summary_list[[s]] <- data.table(
        sample    = s,
        mean_DP   = round(mean(dp_mat[[s]], na.rm = TRUE), 2),
        median_DP = round(median(dp_mat[[s]], na.rm = TRUE), 2),
        sd_DP     = round(sd(dp_mat[[s]], na.rm = TRUE), 2)
    )
}
summary_table <- rbindlist(summary_list)
write.csv(summary_table, output_table, row.names = FALSE)

# PER-INDIVIDUAL DP HISTOGRAMS — from FULL data, BEFORE downsampling
# Store each individual's full DP vector
cat("Extracting full per-individual DP vectors...\n")
indiv_dp_full <- list()
for (s in clean_samples) {
    v <- dp_mat[[s]]
    indiv_dp_full[[s]] <- v[!is.na(v)]
}

# GLOBAL BOXPLOT — with downsampling to 100000 to avoid crash while readin the all_sites file
cat("Downsampling for global boxplot...\n")
set.seed(42)
if (nrow(dp_mat) > 100000) {
    dp_mat_sample <- dp_mat[sample(1:nrow(dp_mat), 100000), ]
} else {
    dp_mat_sample <- dp_mat
}
dp_long_plot <- melt(dp_mat_sample, id.vars = c("CHROM", "POS"),
                     variable.name = "sample", value.name = "DP")

rm(dp_mat, dp_mat_sample); gc()

# MQ and QUAL
cat("Reading MQ matrix...\n")
mq_mat <- fread(mq_file, header = FALSE, fill = TRUE, na.strings = ".")
colnames(mq_mat) <- c("CHROM", "POS", "MQ")
if (nrow(mq_mat) > 500000) mq_mat <- mq_mat[sample(1:.N, 500000)]

cat("Reading QUAL matrix...\n")
qual_df <- fread(qual_file, header = FALSE, col.names = "QUAL", na.strings = ".")
qual_vals <- qual_df$QUAL[!is.na(qual_df$QUAL)]
if (length(qual_vals) > 500000) qual_vals <- sample(qual_vals, 500000)
rm(qual_df); gc()

# FIGURES
cat("Generating plots...\n")
pdf(output_pdf, width = 16, height = 8)

# 1. Global DP boxplot (downsampled)
print(
    ggplot(dp_long_plot[!is.na(DP)], aes(x = sample, y = DP)) +
    geom_boxplot(fill = "aquamarine4", outlier.size = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    labs(title = "Sequencing depth (DP) per individual (subsample of 100k sites) - E-DF1",
         x = "Individual", y = "DP")
)

# 2. Zoomed boxplot
print(
    ggplot(dp_long_plot[!is.na(DP) & DP <= 375], aes(x = sample, y = DP)) +
    geom_boxplot(fill = "aquamarine4", outlier.size = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    labs(title = "Sequencing depth (DP) per individual (zoom 0-375, 100k subsample) - E-DF1",
         x = "Individual", y = "DP")
)

# 3. MQ histogram
print(
    ggplot(mq_mat[!is.na(MQ)], aes(x = MQ)) +
    geom_histogram(bins = 50, fill = "lightskyblue4", color = "white") +
    labs(title = "Mapping quality (MQ) distribution - E-DF1", x = "MQ", y = "Count")
)

# 4. QUAL histogram
print(
    ggplot(data.frame(QUAL = qual_vals), aes(x = QUAL)) +
    geom_histogram(bins = 50, fill = "lightpink4", color = "white") +
    labs(title = "Variant quality (QUAL) distribution - E-DF1", x = "QUAL", y = "Count")
)

# 5. Per-individual DP histograms from FULL data
for (indiv in clean_samples) {
    v <- indiv_dp_full[[indiv]]
    print(
        ggplot(data.frame(DP = v[v <= 375]), aes(x = DP)) +
        geom_histogram(bins = 50, fill = "aquamarine4", color = "white") +
        labs(title = paste("DP distribution (full data) -", indiv),
             x = "Depth", y = "Count")
    )
}

dev.off()
cat("\nDone.\n")
