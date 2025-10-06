library(dplyr)
library(ggplot2)
library(gridExtra)

source("Utils.R")
# read in patents

sheets <- openxlsx::getSheetNames("~/Documents/Projects/MAGEA3/patent.xlsx")
patents_list = lapply(sheets, function(sheet){
  patent = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/patent.xlsx", sheet = sheet) %>%
    mutate(mismatch = stringdist::stringdist("KVAELVHFL", peptide, method = "hamming"),
           source = sheet) %>%
    distinct(gene, peptide, source, .keep_all = T) %>%
    filter(!is.na(peptide))
})

names(patents_list) = sheets

patents_df = do.call("rbind", patents_list) %>%
  distinct(gene, peptide, source, .keep_all = T) %>%
  #mutate(EpiTox = ifelse(peptide %in% cutoff_4$peptide, "Yes", "No")) %>%
  filter(!is.na(peptide)) %>%
  select(-source, -Note) %>%
  distinct (.keep_all = T) %>%
  arrange(mismatch)
rownames(patents_df) = 1:nrow(patents_df)

# read in epitox cutoff 4
cutoff_4 = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  select(chr, uniprot, Gene.Names, peptide, mismatch, blosum_similarity, Wildtype,
         affinity, processing_score, presentation_percentile,
         Ranking_score, Rank) %>%
  distinct(peptide, .keep_all = T) %>%
  mutate(is_known_offtarget = ifelse(peptide %in% patents_df$peptide, T, F)) %>%
  relocate(is_known_offtarget, .after = "Ranking_score") %>%
  arrange(-Ranking_score)



# calculate enrichment
N <- nrow(cutoff_4)
K <- sum(cutoff_4$is_known_offtarget)

cat("\n=== DATA SUMMARY ===\n")
cat("Total OT (N):", N, "\n")
cat("Known found OTs (K):", K, "\n")


thresholds <- c(50, 100, 200, 300, 500, 750, 1000, nrow(cutoff_4))

# Calculate enrichment for each threshold
results_list <- lapply(thresholds, function(n) {
  calculate_enrichment(n, cutoff_4, "is_known_offtarget", K, N)
})

# Convert to data frame
results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(
    threshold = x$n,
    found = x$k,
    expected = round(x$expected, 2),
    enrichment = round(x$enrichment, 2),
    fisher_pval = formatC(x$fisher_pval, format = "e", digits = 2),
    hypergeo_pval = formatC(x$hypergeo_pval, format = "e", digits = 2),
    odds_ratio = round(x$odds_ratio, 2)
  )
}))

# Print results
cat("\n=== ENRICHMENT ANALYSIS RESULTS ===\n\n")
cat("Total predictions (N):", N, "\n")
cat("Known off-targets (K):", K, "\n\n")
print(results_df, row.names = FALSE)

cat("\n=== INTERPRETATION ===\n")
cat("Enrichment > 1: Your ranking is better than random\n")
cat("p-value < 0.05: Enrichment is statistically significant\n")
cat("Odds ratio: How much more likely to find off-target in top vs bottom\n\n")

offtarget_column = "is_known_offtarget"
#=========================================
# Create a finer grid for smooth curve
fine_thresholds <- seq(10, N, by = 10)
fine_enrichments <- sapply(fine_thresholds, function(n) {
  k <- sum(cutoff_4[[offtarget_column]][1:n])
  if (k == 0) return(0)
  enrichment <- (k / n) / (K / N)
  return(enrichment)
})

# Prepare data for plotting
enrichment_df <- data.frame(
  threshold = fine_thresholds,
  enrichment = fine_enrichments
)

# Add key threshold points
key_points_df <- results_df[, c("threshold", "enrichment")]
key_points_df$enrichment <- as.numeric(key_points_df$enrichment)

# Plot 1: Enrichment curve
p1 <- ggplot(enrichment_df, aes(x = threshold, y = enrichment)) +
  geom_line(color = "#A2C510", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#E63946", size = 1) +
  geom_point(data = key_points_df, aes(x = threshold, y = enrichment),
             color = "#A2C510", size = 3) +
  labs(x = "Top N Off-targets",
       y = "Enrichment Factor",
       title = "Enrichment vs Threshold") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank())

# Plot 2: Cumulative recovery curve
cumulative_k <- sapply(fine_thresholds, function(n) {
  sum(cutoff_4[[offtarget_column]][1:n])
})
cumulative_percent <- (cumulative_k / K) * 100
random_recovery <- (fine_thresholds / N) * 100

recovery_df <- data.frame(
  threshold = rep(fine_thresholds, 2),
  percent_recovered = c(cumulative_percent, random_recovery),
  type = rep(c("EpiTox", "Random"), each = length(fine_thresholds))
)

# "#A2C510", "#FBB800", "#99CFE9", "#2C4255"
p2 <- ggplot(recovery_df, aes(x = threshold, y = percent_recovered, color = type)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("EpiTox" = "#A2C510", "Random" = "#2C4255")) +
  labs(x = "Top N Off-targets",
       y = "% Known Off-Targets Recovered",
       title = "Cumulative Recovery of Known Off-Targets",
       color = "") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = c(0.8, 0.3),
        legend.background = element_rect(fill = "white", color = "gray80"),
        panel.grid.minor = element_blank())

# Plot 3: -log10(p-value) across thresholds
log_pvals <- -log10(as.numeric(results_df$fisher_pval))
pval_df <- data.frame(
  threshold = results_df$threshold,
  log_pval = log_pvals
)

p3 <- ggplot(pval_df, aes(x = threshold, y = log_pval)) +
  geom_line(color = "#BED658", size = 1.2) +
  geom_point(color = "#A2C510", size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "#2C4255", size = 1) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed",
             color = "#FBB800", size = 1) +
  annotate("text", x = max(pval_df$threshold) * 0.7, y = -log10(0.05) + 0.3,
           label = "p = 0.05", color = "#2C4255", size = 3.5) +
  annotate("text", x = max(pval_df$threshold) * 0.7, y = -log10(0.01) + 0.3,
           label = "p = 0.01", color = "#FBB800", size = 3.5) +
  labs(x = "Top N Off-targets",
       y = "-log10(p-value)",
       title = "Statistical Significance across Thresholds") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank())

# Plot 4: Distribution of known off-target positions
known_positions_for_plot <- which(cutoff_4[[offtarget_column]])
positions_df <- data.frame(position = known_positions_for_plot)
median_pos <- median(known_positions_for_plot)

p4 <- ggplot(positions_df, aes(x = position)) +
  geom_histogram(bins = 20, fill = "#E3EEB7", color = "#A2C510", alpha = 0.8) +
  geom_vline(xintercept = median_pos, linetype = "dashed",
             color = "#2C4255", size = 1.2) +
  annotate("text", x = median_pos + N * 0.1, y = Inf,
           label = paste("Median:", median_pos),
           color = "#2C4255", vjust = 2, size = 4) +
  labs(x = "Rank Position",
       y = "Frequency",
       title = "Distribution of Known Off-Targets\nin Ranked List") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank())

# Arrange all plots in a 2x2 grid
grid.arrange(p1, p2, p3, p4, ncol = 2)



# ============================================
# ADDITIONAL USEFUL METRICS
# ============================================

cat("\n=== ADDITIONAL METRICS ===\n\n")

# Get positions of known off-targets for these metrics
known_positions_all <- which(cutoff_4[[offtarget_column]])

# Median rank of known off-targets
cat("Median rank of known off-targets:", median(known_positions_all), "\n")

# What % of known off-targets in top 10%, 20%, etc?
percentiles <- c(0.1, 0.2, 0.3, 0.5)
for (p in percentiles) {
  threshold_p <- ceiling(N * p)
  k_in_percentile <- sum(cutoff_4[[offtarget_column]][1:threshold_p])
  percent_recovered <- (k_in_percentile / K) * 100
  cat(sprintf("Top %.0f%% (%d predictions): %d/%d (%.1f%%) known off-targets\n",
              p*100, threshold_p, k_in_percentile, K, percent_recovered))
}

# Area under enrichment curve (rough estimate)
auc_enrichment <- sum(fine_enrichments * 10) / N
cat("\nArea Under Enrichment Curve:", round(auc_enrichment, 2), "\n")
cat("(Higher is better; >1 means better than random overall)\n")








