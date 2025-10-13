# Peptide Anchor Position Analysis for HLA Binding
# This script identifies anchor positions based on substitution matrix data

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# ============================================================================
# 1. DATA LOADING AND TRANSFORMATION
# ============================================================================

# Function to convert long format (peptide, affinity) to position x substitution matrix
convert_to_matrix <- function(df, wt_peptide) {
  # df should have columns: peptide, affinity
  # wt_peptide: wild-type peptide sequence (string, e.g., "ALYNKVFLK")

  # Convert wild-type peptide to character vector
  wt_sequence <- strsplit(wt_peptide, "")[[1]]
  L <- nchar(wt_peptide)

  # All possible amino acids
  aa_list <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

  # Initialize matrix to store results
  result_list <- list()

  # Process each mutant peptide
  for (i in 1:nrow(df)) {
    mut_peptide <- df$peptide[i]
    mut_affinity <- df$affinity[i]
    mut_sequence <- strsplit(mut_peptide, "")[[1]]

    # Find which position was mutated
    diff_positions <- which(wt_sequence != mut_sequence)

    # Handle single substitutions
    if (length(diff_positions) == 1) {
      pos <- diff_positions[1]
      mut_aa <- mut_sequence[pos]

      # Store in list format: position_AA = affinity
      key <- paste0("P", pos, "_", mut_aa)
      result_list[[key]] <- mut_affinity
    }
  }

  # Build matrix from list
  substitution_matrix <- matrix(NA, nrow = L, ncol = length(aa_list),
                                dimnames = list(paste0("P", 1:L), aa_list))

  for (key in names(result_list)) {
    parts <- strsplit(key, "_")[[1]]
    pos <- as.numeric(gsub("P", "", parts[1]))
    aa <- parts[2]

    if (aa %in% aa_list) {
      substitution_matrix[pos, aa] <- result_list[[key]]
    }
  }

  return(substitution_matrix)
}


# input
target = data.frame(peptide = "KVAELVHFL", affinity = 14.39613)
substitution_matrix <- openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/Xscan_MAGEA3.xlsx") %>%
  select(peptide, affinity) %>%
  convert_to_matrix(df = ., wt_peptide = as.character(target$peptide[1]))


# ============================================================================
# 2A. CALCULATE METRICS FOR ANCHOR POSITION IDENTIFICATION
# ============================================================================
analyze_anchor_positions <- function(substitution_matrix, wt_kd, peptide_sequence) {

  L <- nrow(substitution_matrix)

  # Initialize results dataframe
  results <- data.frame(
    Position = 1:L,
    WT_AA = peptide_sequence,
    WT_KD = wt_kd,
    Mean_Mut_KD = NA,
    Median_Mut_KD = NA,
    Min_Mut_KD = NA,
    Max_Mut_KD = NA,
    Fold_Change_Mean = NA,
    Fold_Change_Max = NA,
    Fold_Change_Min = NA,  # Added: best mutation
    Delta_KD_Mean = NA,
    Delta_KD_Max = NA,
    Sensitivity_Score = NA,
    CV = NA,
    Percent_Better_Than_WT = NA,  # NEW: % of mutations that improve binding
    Best_Improvement = NA,  # NEW: best improvement over WT
    Is_Anchor_Classical = FALSE,  # Classical definition
    Is_Permissive = FALSE  # NEW: allows many substitutions
  )

  # Calculate metrics for each position
  for (i in 1:L) {
    mut_kds <- substitution_matrix[i, ]

    results$Mean_Mut_KD[i] <- mean(mut_kds, na.rm = TRUE)
    results$Median_Mut_KD[i] <- median(mut_kds, na.rm = TRUE)
    results$Min_Mut_KD[i] <- min(mut_kds, na.rm = TRUE)
    results$Max_Mut_KD[i] <- max(mut_kds, na.rm = TRUE)

    # Fold changes
    results$Fold_Change_Mean[i] <- results$Mean_Mut_KD[i] / wt_kd
    results$Fold_Change_Max[i] <- results$Max_Mut_KD[i] / wt_kd
    results$Fold_Change_Min[i] <- results$Min_Mut_KD[i] / wt_kd

    # Delta KD (absolute change)
    results$Delta_KD_Mean[i] <- results$Mean_Mut_KD[i] - wt_kd
    results$Delta_KD_Max[i] <- results$Max_Mut_KD[i] - wt_kd

    # Sensitivity score: normalized range of KD values
    results$Sensitivity_Score[i] <- (max(mut_kds, na.rm = TRUE) - min(mut_kds, na.rm = TRUE)) / wt_kd

    # Coefficient of variation (CV)
    results$CV[i] <- sd(mut_kds, na.rm = TRUE) / mean(mut_kds, na.rm = TRUE)

    # METRICS for permissiveness analysis
    # What % of mutations maintain or improve binding?
    better_than_wt <- sum(mut_kds <= wt_kd, na.rm = TRUE)
    total_valid <- sum(!is.na(mut_kds))
    results$Percent_Better_Than_WT[i] <- (better_than_wt / total_valid) * 100

    # Best improvement (if any mutation is better than WT)
    best_mut_kd <- min(mut_kds, na.rm = TRUE)
    results$Best_Improvement[i] <- wt_kd / best_mut_kd  # >1 means improvement
  }

  # Classical anchor identification: high fold-change mean and high sensitivity
  fc_threshold <- quantile(results$Fold_Change_Mean, 0.7, na.rm = TRUE)
  sens_threshold <- quantile(results$Sensitivity_Score, 0.7, na.rm = TRUE)

  results$Is_Anchor_Classical <- (results$Fold_Change_Mean > fc_threshold &
                                    results$Sensitivity_Score > sens_threshold)

  # Permissive positions - allow many substitutions without losing binding
  # High % better than WT OR low fold change mean
  perm_threshold <- quantile(results$Percent_Better_Than_WT, 0.7, na.rm = TRUE)
  results$Is_Permissive <- (results$Percent_Better_Than_WT > perm_threshold)

  return(results)
}

# ============================================================================
# 2B. ADDITIONAL ANALYSIS: POSITION INTERACTIONS
# ============================================================================

analyze_position_compensation <- function(anchor_results) {

  cat("\n=== POSITION CLASSIFICATION ===\n\n")

  # Classify positions into categories
  anchor_results$Category <- "Neutral"
  anchor_results$Category[anchor_results$Is_Anchor_Classical] <- "Classical Anchor"
  anchor_results$Category[anchor_results$Is_Permissive] <- "Permissive"
  anchor_results$Category[anchor_results$Is_Anchor_Classical &
                            anchor_results$Is_Permissive] <- "Flexible Anchor"

  cat("Position categories:\n")
  print(table(anchor_results$Category))

  cat("\n\nClassical Anchors (high sensitivity, intolerant to mutations):\n")
  classical <- anchor_results %>%
    filter(Category == "Classical Anchor") %>%
    select(Position, WT_AA, Fold_Change_Mean, Percent_Better_Than_WT)
  print(classical)

  cat("\n\nPermissive Positions (tolerate many substitutions):\n")
  permissive <- anchor_results %>%
    filter(Category == "Permissive") %>%
    select(Position, WT_AA, Fold_Change_Mean, Percent_Better_Than_WT)
  print(permissive)

  cat("\n\nFlexible Anchors (high sensitivity BUT also allow good alternatives):\n")
  flexible <- anchor_results %>%
    filter(Category == "Flexible Anchor") %>%
    select(Position, WT_AA, Fold_Change_Mean, Percent_Better_Than_WT, Best_Improvement)
  print(flexible)

  return(anchor_results)
}

# Run analysis
anchor_results <- analyze_anchor_positions(substitution_matrix, target$affinity, target$peptide)
anchor_results <- analyze_position_compensation(anchor_results)


# ============================================================================
# 3. VISUALIZATION
# ============================================================================
# Plot 1: Update fill aesthetic
L = nchar(target$peptide)

p1 = ggplot(anchor_results, aes(x = Position, y = Fold_Change_Mean)) +
  geom_col(aes(fill = Category), color = "#2C4255") +  # Changed from Is_Anchor to Category
  geom_hline(yintercept = 1, linetype = "dashed", color = "#A2C510") +
  scale_fill_manual(values = c("Classical Anchor" = "#A2C510",
                               "Permissive" = "steelblue",
                               "Flexible Anchor" = "purple",
                               "Neutral" = "#F0F0F0"),
                    name = "Position Type") +
  labs(title = "Mean Fold Change in KD by Position",
       x = "Peptide Position",
       y = "Mean Fold Change\n(Mutant KD / WT KD)") +
  theme_minimal() +
  theme(text = element_text(size = 10),
        legend.position = "none") +
  scale_x_continuous(breaks = 1:L)


  # rest stays the same...

  # Plot 2: COMPLETELY REPLACED - now shows permissiveness
  p2 = ggplot(anchor_results, aes(x = Position, y = Percent_Better_Than_WT)) +
  geom_col(aes(fill = Category), color = "#2C4255") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#A2C510") +
  scale_fill_manual(values = c("Classical Anchor" = "#A2C510",
                               "Permissive" = "steelblue",
                               "Flexible Anchor" = "purple",
                               "Neutral" = "#F0F0F0"),
                    name = "Position Type") +
  labs(title = "Position Permissiveness",
       subtitle = "% of mutations that maintain/improve binding",
       x = "Peptide Position",
       y = "% Better or Equal to WT") +
  theme_minimal() + theme(text = element_text(size = 10),
                          legend.position = "none") +
    scale_x_continuous(breaks = 1:L)

# Plot 3: COMPLETELY REPLACED - now Fold Change vs Permissiveness
  # Replace the entire p3 plot with this improved version:
  p3 <- ggplot(anchor_results, aes(x = Fold_Change_Mean, y = Percent_Better_Than_WT)) +
    geom_point(aes(color = Category, size = Sensitivity_Score)) +
    geom_text(aes(label = Position), vjust = -0.8, size = 3) +
    geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c("Classical Anchor" = "#A2C510",
                                  "Permissive" = "steelblue",
                                  "Flexible Anchor" = "purple",
                                  "Neutral" = "gray70"),
                       name = "Position Type") +
    scale_size_continuous(name = "Sensitivity", range = c(2, 6)) +
    labs(title = "Position Classification Map",
         subtitle = "Classical anchors (high FC, low %), Permissive (low FC, high %)",
         x = "Mean Fold Change",
         y = "% Better or Equal to WT") +
    theme_minimal() +
    theme(text = element_text(size = 10),
          legend.position = "right",
          legend.box = "vertical",
          legend.box.just = "left",
          legend.key.size = unit(0.8, "cm"),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10, face = "bold"))


# Plot 4: Heatmap of substitution effects
# Convert to fold-change matrix for visualization
fc_matrix <- substitution_matrix / target$affinity

# Transpose so rows are amino acids and columns are positions
fc_matrix_t <- t(fc_matrix)
colnames(fc_matrix_t) <- paste0("P", 1:nchar(target$peptide))

# Convert to long format for ggplot
heatmap_data <- as.data.frame(fc_matrix_t) %>%
  rownames_to_column("AA") %>%
  pivot_longer(-AA, names_to = "Position", values_to = "FoldChange")

p4 <- ggplot(heatmap_data, aes(x = Position, y = AA, fill = log10(FoldChange))) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2(low = "white", high = "#A2C510",
                       midpoint = 0, na.value = "#F0F0F0",
                       name = "Log10\nFold Change") +
  labs(title = "Log10 Fold Change Heatmap (Mutant KD / WT KD)",
       x = "Peptide Position",
       y = "Amino Acid") +
  theme_minimal() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Replace the combined_plot section with this:
combined_plot <- grid.arrange(
  arrangeGrob(p1, p2, ncol = 2),
  arrangeGrob(p3, p4, ncol = 2),
  nrow = 2,
  heights = c(0.7, 1.1),  # Changed: less space for overview (p1,p2), more for details (p3,p4)
  top = textGrob("",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Display combined plot
print(combined_plot)
# ============================================================================
# 4. SUMMARY OUTPUT
# ============================================================================

# Replace the entire summary section with:
cat("\n=== ANCHOR POSITION ANALYSIS SUMMARY ===\n\n")
cat("Peptide Sequence:", paste(peptide_sequence, collapse = ""), "\n")
cat("Peptide Length:", L, "\n")
cat("Wild-type KD:", wt_kd, "nM\n\n")

cat("=== KEY INSIGHT FOR YOUR EXPERIMENTAL PARADOX ===\n")
cat("Positions classified as 'Permissive' or 'Flexible Anchor' can explain\n")
cat("why peptides bind despite having 'wrong' anchor residues:\n")
cat("- They may tolerate multiple amino acid types\n")
cat("- Alternative residues can maintain binding\n")
cat("- Classical anchor rules may be too strict\n\n")

cat("Identified Classical Anchors (intolerant to mutations):\n")
classical_pos <- anchor_results %>%
  filter(Category == "Classical Anchor") %>%
  select(Position, WT_AA, Fold_Change_Mean, Percent_Better_Than_WT)
print(classical_pos)

cat("\n\nIdentified Permissive Positions (tolerate many substitutions):\n")
permissive_pos <- anchor_results %>%
  filter(Category == "Permissive") %>%
  select(Position, WT_AA, Fold_Change_Mean, Percent_Better_Than_WT)
print(permissive_pos)

cat("\n\nComplete Results Table:\n")
print(anchor_results %>% select(Position, WT_AA, Category, Fold_Change_Mean,
                                Percent_Better_Than_WT, Sensitivity_Score))
