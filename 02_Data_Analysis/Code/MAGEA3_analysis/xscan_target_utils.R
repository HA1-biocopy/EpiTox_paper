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

  cat("Position categories (counts):\n")
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
