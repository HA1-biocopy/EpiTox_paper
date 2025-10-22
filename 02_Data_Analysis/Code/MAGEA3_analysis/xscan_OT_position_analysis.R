# Peptide Dissimilarity Position Analysis
# Analyzes where off-target peptides differ from target, classified by position type

library(tidyverse)
source("xscan_OT_utils.R")

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

# Example target sequence and classifications
target_sequence <- "KVAELVHFL"

target_classifications <- anchor_results %>%
  mutate(WT_AA = strsplit(target_sequence, "")[[1]]) %>%
  select(Position, WT_AA, Category)

offtargets_data = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  select(peptide, Wildtype, blosum_similarity, mismatch,
         affinity, presentation_score) %>%
  distinct(.keep_all = T)

# Run analysis
cat("=== ANALYZING OFF-TARGET PEPTIDES ===\n\n")

# Perform complete analysis
analysis_results <- analyze_offtargets(target_sequence, target_classifications, offtargets_data)

# Create comprehensive summary table
summary_table <- create_summary_table(analysis_results)


categorized <- summary_table %>%
  mutate(
    Category = case_when(
      Anchor_Mismatches > 0 & Exposed_Mismatches > 0 ~ "Both Anchor & Backbone",
      Anchor_Mismatches > 0 & Exposed_Mismatches == 0 ~ "Anchor Only",
      Anchor_Mismatches == 0 & Exposed_Mismatches > 0 ~ "Backbone Only",
      TRUE ~ "Perfect Match"
    )
  )

biocopy_colors = c("#A2C510", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")

# just for the plots
PieChart(Category, data = categorized, hole = 0,
         fill = biocopy_colors,
         color="white",
         main = "Off-target substitution profile",
         labels_cex = 0.6)

ggplot(categorized, aes(Total_Mismatches, fill = Category)) +
  geom_bar(position = "dodge")

# Grouped bar plot
p_grouped <- plot_mismatch_by_type(summary_table)
print(p_grouped)

# Stacked bar plot
p_stacked <- plot_mismatch_stacked(summary_table)
print(p_stacked)

# Heatmap
p_heatmap <- plot_mismatch_heatmap(summary_table)
print(p_heatmap)

# For large datasets
summary_plots <- plot_dataset_summary(summary_table)

# View individual plots
print(summary_plots$pie_chart)
print(summary_plots$risk_distribution)

# Or view all at once
dashboard <- plot_summary_dashboard(summary_table)

