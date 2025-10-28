# Precision and Positive Predictive Value (PPV) Analysis
# =========================================================

# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)

biocopy_colors = c("#A2C510", "#2C4255", "#99CFE9", "#FBB800", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#939597", "#3373A1", "#D14B47", "#98D0BC", "#4F3D7F", "#F08000", "#64686A", "#338232")

# SETUP
# ===============
# - is_known_offtarget: TRUE/FALSE for known off-targets
# - is_experimentally_tested: TRUE/FALSE for predictions you tested in the lab
# - is_validated_offtarget: TRUE/FALSE for experimentally confirmed off-targets
#   (this includes BOTH known and novel off-targets found in experiments)

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

# read in SCORE
score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  arrange(lev_dist)

# read in epitox cutoff 4
predictions_df = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  select(chr, uniprot, Gene.Names, peptide, mismatch, blosum_similarity, Wildtype,
         affinity, processing_score, presentation_percentile,
         Ranking_score, Rank) %>%
  distinct(peptide, .keep_all = T) %>%
  mutate(is_known_offtarget = ifelse(peptide %in% patents_df$peptide, T, F),
         is_experimentally_tested = ifelse(peptide %in% score_reference$peptide, T, F),
         is_validated_offtarget = ifelse(peptide %in% score_reference$peptide[score_reference$Outcome == "Binder"], T, F)
         ) %>%
  relocate(is_known_offtarget, .after = "Ranking_score") %>%
  arrange(-Ranking_score, affinity)

predictions_df$rank = 1:nrow(predictions_df)
# ==========================================
# DATA SUMMARY
# ==========================================

cat("\n=== DATA SUMMARY ===\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat(sprintf("Total candidates: %d\n", nrow(predictions_df)))
cat(sprintf("Known off-targets: %d\n", sum(predictions_df$is_known_offtarget)))
cat(sprintf("Experimentally tested: %d (%.1f%%)\n",
            sum(predictions_df$is_experimentally_tested),
            sum(predictions_df$is_experimentally_tested)/nrow(predictions_df)*100))

tested_df <- predictions_df[predictions_df$is_experimentally_tested, ]
n_tested <- nrow(tested_df)
n_validated <- sum(tested_df$is_validated_offtarget)
n_known_in_tested <- sum(tested_df$is_known_offtarget)
n_known_validated <- sum(tested_df$is_known_offtarget & tested_df$is_validated_offtarget)
n_novel_validated <- sum(tested_df$is_validated_offtarget & !tested_df$is_known_offtarget)

cat(sprintf("\nExperimental Results:\n"))
cat(sprintf("  Validated off-targets: %d/%d (%.1f%%)\n",
            n_validated, n_tested, n_validated/n_tested*100))
cat(sprintf("    ├─ Known (in database): %d\n", n_known_validated))
cat(sprintf("    └─ Novel (discovered): %d\n", n_novel_validated))
cat(sprintf("  False positives: %d (%.1f%%)\n",
            n_tested - n_validated, (n_tested - n_validated)/n_tested*100))
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# ==========================================
# ANALYSIS 1: Binned Discovery Rates
# ==========================================

cat("=== BINNED ANALYSIS ===\n\n")

# Create bins (you can adjust bin size)
bin_size <- 100
max_tested_rank <- max(tested_df$rank)
bins <- seq(0, max_tested_rank, by = bin_size)
if (max(bins) < max_tested_rank) {
  bins <- c(bins, max_tested_rank)
}

# Calculate discovery rate per bin
bin_analysis <- data.frame()
for (i in 1:(length(bins)-1)) {
  bin_start <- bins[i] + 1
  bin_end <- bins[i+1]

  bin_data <- tested_df[tested_df$rank >= bin_start & tested_df$rank <= bin_end, ]
  n_in_bin <- nrow(bin_data)

  if (n_in_bin > 0) {
    n_validated_in_bin <- sum(bin_data$is_validated_offtarget)
    discovery_rate <- n_validated_in_bin / n_in_bin * 100

    bin_analysis <- rbind(bin_analysis, data.frame(
      bin_label = sprintf("%d-%d", bin_start, bin_end),
      bin_midpoint = (bin_start + bin_end) / 2,
      bin_start = bin_start,
      bin_end = bin_end,
      n_tested = n_in_bin,
      n_validated = n_validated_in_bin,
      discovery_rate = discovery_rate
    ))
  }
}

cat("Discovery Rate by Rank Bin:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
print(bin_analysis[, c("bin_label", "n_tested", "n_validated", "discovery_rate")],
      row.names = FALSE, digits = 1)
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# ==========================================
# ANALYSIS 2: Statistical Test for Non-Random Distribution
# ==========================================

cat("=== STATISTICAL VALIDATION ===\n\n")

# Chi-square test: Are off-targets uniformly distributed across bins?
expected_per_bin <- n_validated / nrow(bin_analysis)
observed <- bin_analysis$n_validated
expected <- rep(expected_per_bin, nrow(bin_analysis))

# Only perform chi-square if we have enough data
if (sum(observed) >= 5 && all(expected >= 1)) {
  chi_result <- chisq.test(observed, p = expected/sum(expected))

  cat("Chi-Square Test for Uniform Distribution:\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat(sprintf("Null hypothesis: Off-targets are uniformly distributed across rank bins\n"))
  cat(sprintf("Chi-square statistic: %.2f\n", chi_result$statistic))
  cat(sprintf("P-value: %s\n", formatC(chi_result$p.value, format = "e", digits = 2)))

  if (chi_result$p.value < 0.05) {
    cat("\n✓ SIGNIFICANT: Off-targets are NOT uniformly distributed (p < 0.05)\n")
    cat("  → EpiTox's ranking successfully concentrates off-targets!\n")
  } else {
    cat("\n✗ NOT SIGNIFICANT: Cannot reject uniform distribution (p >= 0.05)\n")
    cat("  → Ranking may not be better than random\n")
  }
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
}

# Permutation test: Compare top vs bottom
top_quartile_n <- ceiling(n_tested * 0.25)
bottom_quartile_n <- ceiling(n_tested * 0.25)

top_quartile <- tested_df[1:top_quartile_n, ]
bottom_quartile <- tested_df[(n_tested - bottom_quartile_n + 1):n_tested, ]

top_validated <- sum(top_quartile$is_validated_offtarget)
bottom_validated <- sum(bottom_quartile$is_validated_offtarget)

cat("Top vs Bottom Quartile Comparison:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat(sprintf("Top 25%% (ranks 1-%d):\n", top_quartile_n))
cat(sprintf("  Validated: %d/%d (%.1f%%)\n",
            top_validated, top_quartile_n, top_validated/top_quartile_n*100))
cat(sprintf("\nBottom 25%% (ranks %d-%d):\n", n_tested - bottom_quartile_n + 1, n_tested))
cat(sprintf("  Validated: %d/%d (%.1f%%)\n",
            bottom_validated, bottom_quartile_n, bottom_validated/bottom_quartile_n*100))

# Fisher's exact test for top vs bottom
contingency <- matrix(c(top_validated, top_quartile_n - top_validated,
                        bottom_validated, bottom_quartile_n - bottom_validated),
                      nrow = 2, byrow = TRUE)
fisher_result <- fisher.test(contingency, alternative = "greater")

cat(sprintf("\nFisher's Exact Test: p = %s\n",
            formatC(fisher_result$p.value, format = "e", digits = 2)))
cat(sprintf("Odds Ratio: %.2f\n", fisher_result$estimate))

if (fisher_result$p.value < 0.05) {
  cat("\n✓ Top quartile is significantly enriched vs bottom quartile!\n")
} else {
  cat("\n✗ No significant difference between top and bottom quartiles\n")
}
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# ==========================================
# ANALYSIS 3: Saturation Analysis
# ==========================================

cat("=== SATURATION ANALYSIS ===\n\n")

# Cumulative discoveries
tested_df_sorted <- tested_df[order(tested_df$rank), ]
tested_df_sorted$cumulative_validated <- cumsum(tested_df_sorted$is_validated_offtarget)
tested_df_sorted$cumulative_tested <- 1:nrow(tested_df_sorted)

# Calculate discovery rate in last 20% of tested predictions
last_20_pct_start <- ceiling(n_tested * 0.8)
last_20_pct <- tested_df_sorted[last_20_pct_start:n_tested, ]
discoveries_in_last_20 <- sum(last_20_pct$is_validated_offtarget)

cat("Saturation Assessment:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat(sprintf("Total validated off-targets: %d\n", n_validated))
cat(sprintf("Discoveries in last 20%% of tested: %d (%.1f%%)\n",
            discoveries_in_last_20, discoveries_in_last_20/n_validated*100))

if (discoveries_in_last_20 <= n_validated * 0.1) {
  cat("\n✓ GOOD SATURATION: <10% of discoveries in last 20% tested\n")
  cat("  → You've likely found most off-targets in your tested range\n")
} else if (discoveries_in_last_20 <= n_validated * 0.2) {
  cat("\n⚠ MODERATE SATURATION: 10-20% of discoveries in last 20% tested\n")
  cat("  → You may have found most off-targets, but some may remain\n")
} else {
  cat("\n✗ POOR SATURATION: >20% of discoveries in last 20% tested\n")
  cat("  → Still discovering at high rate; more off-targets likely remain\n")
}
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# ==========================================
# ANALYSIS 4: Estimate Remaining Off-Targets
# ==========================================

cat("=== ESTIMATION OF REMAINING OFF-TARGETS ===\n\n")

# Use discovery rate in last tested bins to extrapolate
last_bins <- tail(bin_analysis, 2)
avg_discovery_rate_late <- mean(last_bins$discovery_rate) / 100

n_untested <- nrow(predictions_df) - n_tested
estimated_remaining <- round(n_untested * avg_discovery_rate_late)

cat("Extrapolation to Untested Region:\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat(sprintf("Untested predictions: %d\n", n_untested))
cat(sprintf("Avg discovery rate in last 2 bins: %.2f%%\n", avg_discovery_rate_late * 100))
cat(sprintf("Estimated remaining off-targets: ~%d\n", estimated_remaining))
cat(sprintf("Estimated total off-targets: %d (found) + %d (estimated) = %d\n",
            n_validated, estimated_remaining, n_validated + estimated_remaining))

if (estimated_remaining <= n_validated * 0.2) {
  cat("\n✓ Likely found >80% of all off-targets\n")
} else if (estimated_remaining <= n_validated * 0.5) {
  cat("\n⚠ Likely found 67-80% of all off-targets\n")
} else {
  cat("\n⚠ May have found <67% of all off-targets\n")
}
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# ==========================================
# VISUALIZATIONS
# ==========================================

# Plot 1: Discovery rate by rank bin
p1 <- ggplot(bin_analysis, aes(x = bin_midpoint, y = discovery_rate)) +
  geom_bar(stat = "identity", fill = "#A2C510", alpha = 0.8) +
  geom_line(color = "#2C4255", linewidth = 1.2) +
  geom_point(color = "#2C4255", size = 3) +
  labs(x = "Rank Position (Midpoint of Bin)",
       y = "Discovery Rate (%)",
       title = "Discovery Rate Across Ranked List",
       subtitle = "Are off-targets concentrated at the top?") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        panel.grid.minor = element_blank())

# Plot 2: Cumulative discoveries (saturation curve)
p2 <- ggplot(tested_df_sorted, aes(x = rank, y = cumulative_validated)) +
  geom_line(color = "#A2C510", size = 1.5) +
  geom_point(data = tested_df_sorted[tested_df_sorted$is_validated_offtarget, ],
             aes(x = rank, y = cumulative_validated),
             color = "#2C4255", size = 2, alpha = 0.6) +
  labs(x = "Rank Position",
       y = "Cumulative Off-Targets Discovered",
       title = "Saturation Curve",
       subtitle = "Does discovery rate plateau?") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        panel.grid.minor = element_blank())

# Plot 3: Distribution of validated off-targets across rank
validated_ranks <- tested_df$rank[tested_df$is_validated_offtarget]

p3 <- ggplot(data.frame(rank = validated_ranks), aes(x = rank)) +
  geom_histogram(bins = 20, fill = "#A2C510", alpha = 0.7, color = "white") +
  geom_vline(xintercept = median(validated_ranks),
             color = "#2C4255", linetype = "dashed", size = 1.2) +
  annotate("text", x = median(validated_ranks) * 1.2, y = Inf, vjust = 2,
           label = sprintf("Median: %d", round(median(validated_ranks))),
           color = "#2C4255", size = 4) +
  labs(x = "Rank Position",
       y = "Number of Off-Targets",
       title = "Distribution of Validated Off-Targets",
       subtitle = "Where do true off-targets appear in ranking?") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        panel.grid.minor = element_blank())

# Plot 4: Comparison of quartiles
quartile_data <- data.frame(
  quartile = c("Top 25%", "25-50%", "50-75%", "Bottom 25%"),
  discovery_rate = numeric(4),
  n_validated = numeric(4),
  n_tested = numeric(4)
)

quartile_breaks <- quantile(tested_df$rank, probs = c(0, 0.25, 0.5, 0.75, 1))
for (i in 1:4) {
  q_data <- tested_df[tested_df$rank > quartile_breaks[i] &
                        tested_df$rank <= quartile_breaks[i+1], ]
  quartile_data$n_tested[i] <- nrow(q_data)
  quartile_data$n_validated[i] <- sum(q_data$is_validated_offtarget)
  quartile_data$discovery_rate[i] <- sum(q_data$is_validated_offtarget) / nrow(q_data) * 100
}

quartile_data$quartile <- factor(quartile_data$quartile,
                                 levels = c("Top 25%", "25-50%", "50-75%", "Bottom 25%"))

p4 <- ggplot(quartile_data, aes(x = quartile, y = discovery_rate, fill = quartile)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%d/%d\n(%.1f%%)", n_validated, n_tested, discovery_rate)),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = biocopy_colors#c("Top 25%" = "#06A77D",
                              # "25-50%" = "#2E86AB",
                              # "50-75%" = "#F77F00",
                              # "Bottom 25%" = "#E63946")
                    ) +
  labs(x = "Rank Quartile",
       y = "Discovery Rate (%)",
       title = "Discovery Rate by Quartile",
       subtitle = "Is top quartile enriched?") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "none",
        panel.grid.minor = element_blank())

# Display all plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

# ==========================================
# KEY TAKEAWAYS
# ==========================================

cat("\n=== KEY TAKEAWAYS FOR YOUR PAPER ===\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat(sprintf("1. Tested %d/%d candidates (%.1f%%)\n",
            n_tested, nrow(predictions_df), n_tested/nrow(predictions_df)*100))
cat(sprintf("2. Discovered %d validated off-targets\n", n_validated))
cat(sprintf("3. Top 25%% discovery rate: %.1f%% vs Bottom 25%%: %.1f%%\n",
            top_validated/top_quartile_n*100, bottom_validated/bottom_quartile_n*100))
cat(sprintf("4. Statistical validation: %s\n",
            ifelse(exists("fisher_result") && fisher_result$p.value < 0.05,
                   "SIGNIFICANT enrichment at top",
                   "Check statistical tests above")))
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
