library(ggplot2)
library(data.table)

stats_dir    <- "/shared/home/toufiktakhi/neurosporagbs/dp_mq_EDF1_stats/"
output_pdf   <- "/shared/home/toufiktakhi/neurosporagbs/E-DF1_cleaned_figures.pdf"
output_table <- "/shared/home/toufiktakhi/neurosporagbs/E-DF1_cleaned_summary.csv"
dp_file   <- file.path(stats_dir, "final_DP_matrix.txt.gz")
mq_file   <- file.path(stats_dir, "final_MQ.txt.gz")
qual_file <- file.path(stats_dir, "final_QUAL.txt.gz")

samples <- readLines(file.path(stats_dir, "sample_columns.txt"))
clean_samples <- gsub(".*/", "", samples)
clean_samples <- gsub("\\..*", "", clean_samples)

cat("Reading cleaned data...\n")
dp_mat <- fread(dp_file, header = FALSE, fill = TRUE)
colnames(dp_mat) <- c("CHROM", "POS", clean_samples)

mq_mat <- fread(mq_file, header = FALSE, fill = TRUE)
colnames(mq_mat) <- c("CHROM", "POS", "MQ")

qual_vals <- as.numeric(readLines(gzcon(file(qual_file, "rb"))))
qual_vals <- qual_vals[!is.na(qual_vals)]

# rewrite DP to long format
cat("Reshaping data...\n")
dp_long <- melt(dp_mat, id.vars = c("CHROM", "POS"), variable.name = "sample", value.name = "DP")
dp_long[, DP := as.numeric(DP)]

summary_rows <- list()

cat("Generating FINAL plots...\n")
pdf(output_pdf, width = 16, height = 8)

# box plots
# DP boxplot
print(
    ggplot(dp_long[!is.na(DP)], aes(x = sample, y = DP)) +
    geom_boxplot(fill = "steelblue", outlier.size = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) + # Angle 45 is easier to read
    labs(title = "RAW sequencing depth (DP) per individual - E-DF1", x = "Individual", y = "DP")
)
# DP Boxplot (Zoomed 0-375)
print(
    ggplot(dp_long[!is.na(DP) & DP <= 375], aes(x = sample, y = DP)) +
    geom_boxplot(fill = "steelblue", outlier.size = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    labs(title = "CLEANED sequencing depth (DP) per individual (Min DP 3) - E-DF1", x = "Individual", y = "DP")
)

# MQ Histogram

print(
    ggplot(mq_mat[!is.na(MQ)], aes(x = MQ)) +
    geom_histogram(bins = 50, fill = "darkorange", color = "white") +
    labs(title = "Cleaned Mapping Quality (MQ) distribution (Min MQ 30) - E-DF1", x = "MQ", y = "Count")
)

#QUAL Histogram
print(
    ggplot(data.frame(QUAL = qual_vals), aes(x = QUAL)) +
    geom_histogram(bins = 50, fill = "tomato", color = "white") +
    labs(title = "CLEANED QUAL distribution (Min QUAL 30) - E-DF1", x = "QUAL", y = "Count")
)

# loop per-Individual DP Histograms
for (indiv in clean_samples) {
    indiv_dp <- dp_long[sample == indiv, DP]
    indiv_dp <- indiv_dp[!is.na(indiv_dp)]

    print(
        ggplot(data.frame(DP = indiv_dp[indiv_dp <= 375]), aes(x = DP)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white") +
        labs(title = paste("DP distribution -", indiv), x = "Depth", y = "Count")
    )
}

dev.off()

summary_table <- dp_long[, .(
    mean_DP = round(mean(DP, na.rm = TRUE), 2),
    median_DP = median(DP, na.rm = TRUE)
), by = sample]

write.csv(summary_table, output_table, row.names = FALSE)

cat("\nSuccess! Final graphs and table generated.\n")
