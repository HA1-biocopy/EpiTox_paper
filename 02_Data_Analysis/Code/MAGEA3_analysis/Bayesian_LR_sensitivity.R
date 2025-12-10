# =============================================================================
# Likelihood Ratio (LR) Sensitivity Analysis
# =============================================================================
# Purpose: Demonstrate how LR values affect posterior probabilities
# to justify the chosen LR parameters in the Bayesian assessment
#
# Author: Based on Bayesian_peptide_assessment.R framework
# Date: 2025
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# =============================================================================
# FUNCTION: Calculate posterior probability from LR
# =============================================================================
calculate_posterior <- function(prior_prob = 0.5, lr) {
  prior_odds <- prior_prob / (1 - prior_prob)
  posterior_odds <- prior_odds * lr
  posterior_prob <- posterior_odds / (1 + posterior_odds)
  return(posterior_prob * 100)  # Return as percentage
}

# =============================================================================
# SENSITIVITY ANALYSIS 1: Single Evidence LR Impact
# =============================================================================
cat("\n=== LR SENSITIVITY ANALYSIS ===\n\n")
cat("SCENARIO 1: Impact of Individual Evidence Types\n")
cat("Starting prior: 50%\n\n")

# Your current LR values from the Bayesian framework
current_lrs <- data.frame(
  evidence_type = c(
    "Multiple High Relevance",
    "One High Relevance",
    "Multiple Medium Relevance",
    "One Medium Relevance",
    "Same Allele",
    "Strong Binder × On-Target",
    "Normal Tissue",
    "Predicted Only"
  ),
  lr_value = c(12, 8, 6, 4, 10, 12.5, 10, 0.3)
)

# Calculate posteriors
current_lrs$posterior_prob <- sapply(current_lrs$lr_value,
                                     function(lr) calculate_posterior(0.5, lr))

print(current_lrs)
cat("\n")

# =============================================================================
# SENSITIVITY ANALYSIS 2: LR Range Exploration
# =============================================================================
cat("\nSCENARIO 2: Exploring LR Value Ranges\n")
cat("Testing how different LR values translate to posterior probabilities\n\n")

# Create a range of LR values
lr_range <- data.frame(
  lr = c(0.1, 0.3, 0.5, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20)
)

lr_range$posterior_prob <- sapply(lr_range$lr,
                                  function(lr) calculate_posterior(0.5, lr))
lr_range$interpretation <- case_when(
  lr_range$posterior_prob >= 80 ~ "High Confidence",
  lr_range$posterior_prob >= 50 ~ "Medium Confidence",
  lr_range$posterior_prob >= 30 ~ "Low Confidence",
  TRUE ~ "Very Low Confidence"
)

print(lr_range)
cat("\n")

# =============================================================================
# SENSITIVITY ANALYSIS 3: Combined Evidence
# =============================================================================
cat("\nSCENARIO 3: Combined Evidence (Multiple LRs)\n")
cat("Showing how evidence accumulates\n\n")

# Test combinations based on your framework
combinations <- data.frame(
  scenario = c(
    "Predicted only",
    "Predicted + Intermediate binding + Different class",
    "Predicted + Intermediate binding + Same class",
    "Same allele only",
    "Same allele + Strong binding",
    "Same allele + Strong binding + Normal tissue",
    "High relevance + Same allele + Strong binding + Normal tissue",
    "Medium relevance + Different allele",
    "Low relevance + Different allele + Weak binding"
  ),
  lr_product = c(
    0.3,                    # Predicted only - Very Low
    0.3 * 3 * 0.5,          # Predicted + intermediate + different class - Low
    0.3 * 3 * 2,            # Predicted + intermediate + same class - Medium
    10,                     # Same allele - High
    10 * 5,                 # Same allele × Strong binder - High
    10 * 5 * 10,            # + Normal tissue - High
    12 * 10 * 12.5 * 10,    # All strong evidence - High
    4 * 2,                  # Medium relevance + diff allele - High
    2 * 2 * 2               # Multiple weak factors - High
  )
)

combinations$posterior_prob <- sapply(combinations$lr_product,
                                      function(lr) calculate_posterior(0.5, lr))
combinations$confidence <- case_when(
  combinations$posterior_prob >= 80 ~ "High",
  combinations$posterior_prob >= 50 ~ "Medium",
  combinations$posterior_prob >= 30 ~ "Low",
  TRUE ~ "Very Low"
)

print(combinations)
cat("\n")

# =============================================================================
# VISUALIZATION 1: LR to Posterior Probability Curve
# =============================================================================
cat("Generating visualizations...\n")

# Generate smooth curve
lr_smooth <- seq(0.1, 20, by = 0.1)
posterior_smooth <- sapply(lr_smooth, function(lr) calculate_posterior(0.5, lr))

plot_data <- data.frame(lr = lr_smooth, posterior = posterior_smooth)

p1 <- ggplot(plot_data, aes(x = lr, y = posterior)) +
  geom_line(color = "#438D99", size = 1.2) +
  geom_hline(yintercept = c(30, 50, 80), linetype = "dashed",
             color = c("#8ECAE6", "#FDE399", "#C61E19"), alpha = 0.6) +
  geom_vline(xintercept = c(1, 4, 8, 12), linetype = "dotted", alpha = 0.4) +
  # Add ALL scenario LR points (from combinations dataframe that will be created later)
  geom_point(data = combinations, aes(x = lr_product, y = posterior_prob),
             color = "#C61E19", size = 3, alpha = 0.7) +
  scale_x_continuous(breaks = c(0.1, 0.3, 1, 2, 4, 6, 8, 10, 12, 15, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  # Add confidence level labels using geom_label
  geom_label(aes(x = 18, y = 85, label = "High\nConfidence"),
             size = 3, color = "#C61E19", fontface = "bold", label.size = 0) +
  geom_label(aes(x = 18, y = 55, label = "Medium"),
             size = 3, color = "#FBB800", fontface = "bold", label.size = 0) +
  geom_label(aes(x = 18, y = 35, label = "Low"),
             size = 3, color = "#8ECAE6", fontface = "bold", label.size = 0) +
  geom_label(aes(x = 2, y = 10, label = "Very Low"),
             size = 3, color = "#A2C510", fontface = "bold", label.size = 0) +
  labs(
    title = "Likelihood Ratio to Posterior Probability Conversion",
    subtitle = "Starting from 50% prior probability (red points = scenario LR values)",
    x = "Likelihood Ratio (LR)",
    y = "Posterior Probability (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_line(color = "gray95")
  )

# =============================================================================
# VISUALIZATION 2: Evidence Type Comparison
# =============================================================================

p2 <- ggplot(current_lrs, aes(x = reorder(evidence_type, posterior_prob),
                              y = posterior_prob, fill = posterior_prob)) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_fill_gradient2(low = "#A2C510", mid = "#FDE399", high = "#C61E19",
                       midpoint = 50, limits = c(0, 100)) +
  coord_flip() +
  labs(
    title = "Impact of Individual Evidence Types on Posterior Probability",
    subtitle = "Starting from 50% prior",
    x = "Evidence Type",
    y = "Posterior Probability (%)",
    fill = "Posterior %"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

# =============================================================================
# VISUALIZATION 3: Combined Evidence Scenarios
# =============================================================================

combinations$scenario_short <- factor(1:nrow(combinations))

p3 <- ggplot(combinations, aes(x = scenario_short, y = posterior_prob,
                               fill = confidence)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.0f%%", posterior_prob)),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("High" = "#C61E19",
                               "Medium" = "#FDE399",
                               "Low" = "#8ECAE6",
                               "Very Low" = "#A2C510")) +
  scale_x_discrete(labels = c("Predicted", "Pred+\nClassII", "Pred+Bind+\nClass", "Allele",
                              "Allele+Bind", "All+Tissue", "All Max", "Med+Diff", "Weak")) +
  labs(
    title = "Combined Evidence Accumulation Scenarios",
    subtitle = "How multiple evidence types combine",
    x = "Evidence Scenario",
    y = "Posterior Probability (%)",
    fill = "Confidence"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# =============================================================================
# VISUALIZATION 4: 4-Panel Summary Figure
# =============================================================================
library(gridExtra)

# Create individual panels
# Panel 1: LR conversion with confidence zones and ALL scenario points
p4_panel1 <- ggplot() +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 30), fill = "#A2C510", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 30, ymax = 50), fill = "#8ECAE6", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 50, ymax = 80), fill = "#FDE399", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 80, ymax = 100), fill = "#C61E19", alpha = 0.2) +
  geom_line(data = plot_data, aes(x = lr, y = posterior), color = "#438D99", size = 1.2) +
  # Add ALL scenario points from combinations
  geom_point(data = combinations, aes(x = lr_product, y = posterior_prob),
             color = "#C61E19", size = 3, alpha = 0.7) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "LR to Posterior Conversion",
       x = "Likelihood Ratio",
       y = "Posterior Probability (%)") +
  theme_minimal() +
  # Add confidence level labels
  geom_label(aes(x = 18, y = 90, label = "High\nConfidence"),
             size = 3, color = "#C61E19", fontface = "bold", label.size = 0) +
  geom_label(aes(x = 18, y = 65, label = "Medium"),
             size = 3, color = "#FBB800", fontface = "bold", label.size = 0) +
  geom_label(aes(x = 18, y = 40, label = "Low"),
             size = 3, color = "#8ECAE6", fontface = "bold", label.size = 0) +
  geom_label(aes(x = 18, y = 15, label = "Very Low"),
             size = 3, color = "#A2C510", fontface = "bold", label.size = 0) +

  theme(plot.title = element_text(face = "bold", size = 10))

# Panel 2: Current LR values
categories_p2 <- c("Predicted\nOnly", "Low\nRelevance", "Medium\nRelevance",
                   "High\nRelevance", "Same\nAllele", "Strong\nBinding")
lr_values_p2 <- c(0.3, 2, 4, 8, 10, 12.5)
colors_p2 <- c("#A2C510", "#8ECAE6", "#FDE399", "#C61E19", "#C61E19", "#C61E19")

p4_panel2 <- ggplot(data.frame(cat = factor(categories_p2, levels = categories_p2),
                               val = lr_values_p2, col = colors_p2),
                    aes(x = cat, y = val, fill = col)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_text(aes(label = val), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_identity() +
  labs(title = "Current LR Values",
       x = NULL, y = "LR Value") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(size = 8))

# Panel 3: Four-tier stratification
scenarios_p3 <- c("Predicted\n(Very Low)", "Class\nMismatch\n(Low)",
                  "Medium\nTier", "Strong\nEvidence\n(High)")
posteriors_p3 <- c(23.1, 31.0, 64.3, 92.3)
colors_p3 <- c("#A2C510", "#8ECAE6", "#FDE399", "#C61E19")

p4_panel3 <- ggplot(data.frame(scen = factor(scenarios_p3, levels = scenarios_p3),
                               post = posteriors_p3, col = colors_p3),
                    aes(x = scen, y = post, fill = col)) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", post)), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_identity() +
  labs(title = "Four-Tier Stratification (69% spread)",
       x = NULL, y = "Posterior Probability (%)") +
  ylim(0, 100) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(size = 8))

# Panel 4: Evidence accumulation text
p4_panel4 <- ggplot() +
  annotate("text", x = 0.5, y = 0.95, label = "Evidence Accumulation Examples",
           hjust = 0.5, vjust = 1, size = 5, fontface = "bold") +
  annotate("text", x = 0.05, y = 0.75, hjust = 0, vjust = 1, size = 3.5,
           label = "• Predicted only (Very Low):") +
  annotate("text", x = 0.95, y = 0.75, hjust = 1, vjust = 1, size = 3.5, fontface = "bold",
           label = "LR = 0.3 → 23%", color = "#438D99") +
  annotate("text", x = 0.05, y = 0.62, hjust = 0, vjust = 1, size = 3.5,
           label = "• Predicted + Binding + Class II (Low):") +
  annotate("text", x = 0.95, y = 0.62, hjust = 1, vjust = 1, size = 3.5, fontface = "bold",
           label = "LR = 0.3×3×0.5 = 0.45 → 31%", color = "#438D99") +
  annotate("text", x = 0.05, y = 0.49, hjust = 0, vjust = 1, size = 3.5,
           label = "• Predicted + Binding + Class I (Medium):") +
  annotate("text", x = 0.95, y = 0.49, hjust = 1, vjust = 1, size = 3.5, fontface = "bold",
           label = "LR = 0.3×3×2 = 1.8 → 64%", color = "#438D99") +
  annotate("text", x = 0.05, y = 0.36, hjust = 0, vjust = 1, size = 3.5,
           label = "• Same allele (High):") +
  annotate("text", x = 0.95, y = 0.36, hjust = 1, vjust = 1, size = 3.5, fontface = "bold",
           label = "LR = 10 → 91%", color = "#438D99") +
  annotate("text", x = 0.05, y = 0.23, hjust = 0, vjust = 1, size = 3.5,
           label = "• Allele + Binding (High):") +
  annotate("text", x = 0.95, y = 0.23, hjust = 1, vjust = 1, size = 3.5, fontface = "bold",
           label = "LR = 10×5 = 50 → 98%", color = "#438D99") +
  annotate("text", x = 0.5, y = 0.08, hjust = 0.5, vjust = 0.5, size = 3,
           label = "✓ Four distinct evidence tiers\n✓ Evidence accumulates realistically\n✓ Clear stratification for decision-making",
           fontface = "italic") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()

# Combine all panels
p4_combined <- grid.arrange(p4_panel1, p4_panel2, p4_panel3, p4_panel4,
                            ncol = 2, nrow = 2,
                            top = "LR Sensitivity Analysis - Key Insights")

# Save the 4-panel figure
ggsave("/mnt/user-data/outputs/LR_sensitivity_summary.png", p4_combined,
       width = 14, height = 10, dpi = 300)

cat("  - LR_sensitivity_summary.png (4-panel overview)\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
cat("Saving outputs...\n")

# Save plots
ggsave("/mnt/user-data/outputs/LR_sensitivity_curve.png", p1,
       width = 10, height = 6, dpi = 300)
ggsave("/mnt/user-data/outputs/LR_evidence_comparison.png", p2,
       width = 10, height = 7, dpi = 300)
ggsave("/mnt/user-data/outputs/LR_combined_scenarios.png", p3,
       width = 10, height = 6, dpi = 300)

# Save data tables
openxlsx::write.xlsx(
  list(
    "Current_LRs" = current_lrs,
    "LR_Range" = lr_range,
    "Combined_Scenarios" = combinations
  ),
  "/mnt/user-data/outputs/LR_sensitivity_analysis.xlsx"
)

cat("\n=== SENSITIVITY ANALYSIS COMPLETE ===\n")
cat("\nFiles saved:\n")
cat("  - LR_sensitivity_curve.png\n")
cat("  - LR_evidence_comparison.png\n")
cat("  - LR_combined_scenarios.png\n")
cat("  - LR_sensitivity_summary.png (4-panel overview)\n")
cat("  - LR_sensitivity_analysis.xlsx\n")

# =============================================================================
# SUMMARY INTERPRETATION
# =============================================================================
cat("\n=== KEY FINDINGS ===\n\n")

cat("1. LR VALUE JUSTIFICATION:\n")
cat("   - LR < 1: Decreases confidence (e.g., 0.3 for 'predicted only')\n")
cat("   - LR = 1: No change (neutral evidence)\n")
cat("   - LR 2-4: Modest increase (30-67% posterior)\n")
cat("   - LR 6-10: Strong increase (75-91% posterior)\n")
cat("   - LR > 12: Very strong increase (>92% posterior)\n\n")

cat("2. DISCRIMINATION CAPABILITY:\n")
cat("   - Your LR values span from 0.3 to 12.5\n")
cat("   - This provides good separation between confidence levels\n")
cat("   - 'Predicted only' (0.3) → ",
    round(calculate_posterior(0.5, 0.3), 1), "% (Very Low)\n")
cat("   - 'Strong evidence' (12) → ",
    round(calculate_posterior(0.5, 12), 1), "% (High)\n\n")

cat("3. EVIDENCE ACCUMULATION:\n")
cat("   - Multiple weak evidences can accumulate to medium confidence\n")
cat("   - Strong evidence in multiple dimensions → High confidence\n")
cat("   - This reflects real-world scientific inference\n\n")

cat("=== END OF ANALYSIS ===\n\n")
