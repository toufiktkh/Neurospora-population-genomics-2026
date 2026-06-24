library(poppr)

# Function to test thresholds for a specific site
test_site <- function(rds_file, site_name) {
  cat("\n--- Testing Thresholds for:", site_name, "---\n")
  obj <- readRDS(rds_file)
  
  # Calculate normalized distance (This is the fast part)
  dist_matrix <- diss.dist(obj) / nLoc(obj)
  
  # Test a range of thresholds
  thresholds <- c(0.001, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05)
  
  for (t in thresholds) {
    mlg.filter(obj, distance = dist_matrix) <- t
    cat("Threshold:", t, " -> Lineages (MLLs):", nmll(obj), "\n")
  }
}

# Run it on the .rds files if they exist
if (file.exists("/shared/home/toufiktakhi/neurosporagbs/DF1_genclone.rds")) test_site("/shared/home/toufiktakhi/neurosporagbs/DF1_genclone.rds", "DF1")
if (file.exists("/shared/home/toufiktakhi/neurosporagbs/DF3_genclone.rds")) test_site("/shared/home/toufiktakhi/neurosporagbs/DF3_genclone.rds", "DF3")
if (file.exists("/shared/home/toufiktakhi/neurosporagbs/Vill_genclone.rds")) test_site("/shared/home/toufiktakhi/neurosporagbs/Vill_genclone.rds", "Vill")
if (file.exists("/shared/home/toufiktakhi/neurosporagbs/GM_genclone.rds")) test_site("/shared/home/toufiktakhi/neurosporagbs/GM_genclone.rds", "LGM")
