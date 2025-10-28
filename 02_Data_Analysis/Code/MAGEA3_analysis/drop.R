# ============================================================================
# Multi-Feature Biophysical Scoring for TCRm Off-Target Assessment
# Lean, Interpretable Approach
# ============================================================================

library(Peptides)
library(tidyverse)
library(corrplot)
library(ggrepel)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Top N
N_DISAGREEMENTS <- 20
N_VALIDATION <- 30

# Define your target peptide
TARGET_PEPTIDE <- "SLYNTVATL"  # REPLACE with your actual target

# Define backbone positions (TCR-facing)
BACKBONE_POSITIONS <- 3:8  # P3-P8 for 9-mers

# Feature weights (biology-informed)
# Backbone features get higher weight as they face TCR
WEIGHTS <- c(
  gravy = 1.0,
  pI = 1.0,
  charge = 1.0,
  aromaticity = 1.0,
  aliphatic_index = 1.0,
  backbone_gravy = 2.0,      # 2x weight
  backbone_charge = 2.0,      # 2x weight
  backbone_aromaticity = 2.0  # 2x weight
)

# Normalize weights to sum to 1
WEIGHTS <- WEIGHTS / sum(WEIGHTS)

# ============================================================================
# FUNCTIONS
# ============================================================================

compute_biophysical_features <- function(peptide_seq, backbone_pos = BACKBONE_POSITIONS) {
  # Global features
  gravy_global <- hydrophobicity(peptide_seq, scale = "KyteDoolittle")
  pI_global <- pI(peptide_seq)
  charge_global <- charge(peptide_seq, pH = 7.4)
  aromat_global <- aIndex(peptide_seq)
  aliph_global <- aliphatic(peptide_seq)

  # Backbone-specific features (TCR-facing positions)
  backbone_seq <- substr(peptide_seq, min(backbone_pos), max(backbone_pos))
  gravy_backbone <- hydrophobicity(backbone_seq, scale = "KyteDoolittle")
  charge_backbone <- charge(backbone_seq, pH = 7.4)
  aromat_backbone <- aIndex(backbone_seq)

  tibble(
    gravy = gravy_global,
    pI = pI_global,
    charge = charge_global,
    aromaticity = aromat_global,
    aliphatic_index = aliph_global,
    backbone_gravy = gravy_backbone,
    backbone_charge = charge_backbone,
    backbone_aromaticity = aromat_backbone
  )
}

compute_weighted_distance <- function(peptide_features, target_features, weights = WEIGHTS) {
  # Scale features (z-score normalization)
  feature_names <- names(weights)

  differences <- sapply(feature_names, function(feat) {
    (peptide_features[[feat]] - target_features[[feat]])
  })

  # Weighted Euclidean distance
  weighted_sq_diff <- weights * (differences^2)
  distance <- sqrt(sum(weighted_sq_diff))

  return(distance)
}

compute_feature_contributions <- function(peptide_features, target_features, weights = WEIGHTS) {
  # Return contribution of each feature to total distance
  feature_names <- names(weights)

  contributions <- sapply(feature_names, function(feat) {
    diff_sq <- (peptide_features[[feat]] - target_features[[feat]])^2
    weighted_contribution <- weights[feat] * diff_sq
    return(weighted_contribution)
  })

  # Return as percentage of total distance
  contributions / sum(contributions) * 100
}

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

# Load your data
# peptide_data <- read_csv("your_peptide_data.csv")
# Expected columns: peptide_id, sequence, anchor_intact, affinity, blosum62

# For demonstration, assuming you have:
# peptide_data <- tibble(
#   peptide_id = ...,
#   sequence = ...,
#   anchor_intact = ...,  # factor with 3 levels
#   affinity = ...,
#   blosum62 = ...
# )

# Compute affinity+BLOSUM combined score
peptide_data <- peptide_data %>%
  mutate(affinityBLOSUM_score = scale(affinity)[,1] + scale(blosum62)[,1])

# ============================================================================
# COMPUTE BIOPHYSICAL FEATURES
# ============================================================================

# Target peptide features
target_features <- compute_biophysical_features(TARGET_PEPTIDE)

# All peptides
peptide_data <- peptide_data %>%
  rowwise() %>%
  mutate(features = list(compute_biophysical_features(sequence))) %>%
  unnest(features) %>%
  ungroup()

# Feature columns
feature_cols <- names(WEIGHTS)

# Scale all features for distance calculation
peptide_data_scaled <- peptide_data %>%
  mutate(across(all_of(feature_cols), ~scale(.)[,1], .names = "{.col}_scaled"))

target_features_scaled <- target_features %>%
  mutate(across(all_of(feature_cols), ~scale(.)[,1], .names = "{.col}_scaled"))

# ============================================================================
# COMPUTE BIOPHYSICAL DISTANCE FROM TARGET
# ============================================================================

peptide_data <- peptide_data %>%
  rowwise() %>%
  mutate(
    biophys_distance = compute_weighted_distance(
      cur_data() %>% select(ends_with("_scaled")),
      target_features_scaled,
      WEIGHTS
    )
  ) %>%
  ungroup()

# ============================================================================
# CREATE RANKINGS
# ============================================================================

peptide_data <- peptide_data %>%
  mutate(
    # Rank by affinity+BLOSUM (higher score = higher rank)
    rank_affinityBLOSUM = rank(-affinityBLOSUM_score, ties.method = "min"),

    # Rank by biophysical distance (smaller distance = more similar = higher rank)
    rank_biophysical = rank(biophys_distance, ties.method = "min"),

    # Compute disagreement
    rank_disagreement = abs(rank_affinityBLOSUM - rank_biophysical),

    # Normalize disagreement to 0-1 scale
    disagreement_score = rank_disagreement / max(rank_disagreement)
  )

# ============================================================================
# EXPLORATORY ANALYSIS
# ============================================================================

cat("=== Summary Statistics ===\n")
cat("Total peptides:", nrow(peptide_data), "\n")
cat("Correlation between rankings:",
    cor(peptide_data$rank_affinityBLOSUM, peptide_data$rank_biophysical), "\n\n")

# Feature correlations
cat("=== Feature Correlations ===\n")
cor_matrix <- peptide_data %>%
  select(all_of(feature_cols)) %>%
  cor(use = "complete.obs")

print(round(cor_matrix, 2))

# Visualize correlations
corrplot(cor_matrix, method = "circle", type = "upper",
         title = "Biophysical Feature Correlations",
         mar = c(0, 0, 2, 0))

# Feature distributions by anchor status
peptide_data %>%
  pivot_longer(cols = all_of(feature_cols),
               names_to = "feature",
               values_to = "value") %>%
  ggplot(aes(x = anchor_intact, y = value, fill = anchor_intact)) +
  geom_boxplot() +
  facet_wrap(~feature, scales = "free_y") +
  labs(title = "Biophysical Features by Anchor Status",
       x = "Anchor Status", y = "Feature Value") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ============================================================================
# RANKING COMPARISON
# ============================================================================

# Scatter plot: two rankings
p1 <- ggplot(peptide_data, aes(x = rank_affinityBLOSUM, y = rank_biophysical)) +
  geom_point(aes(color = disagreement_score), alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  scale_color_viridis_c(option = "plasma", name = "Disagreement") +
  labs(
    title = "Comparison of Ranking Methods",
    subtitle = "Points near diagonal = agreement; far from diagonal = disagreement",
    x = "Rank by Affinity + BLOSUM62",
    y = "Rank by Biophysical Distance"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p1)

# Distribution of disagreement scores
p2 <- ggplot(peptide_data, aes(x = disagreement_score)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of Ranking Disagreements",
    x = "Disagreement Score (0 = perfect agreement, 1 = maximum disagreement)",
    y = "Count"
  ) +
  theme_minimal()

print(p2)

# ============================================================================
# IDENTIFY TOP DISAGREEMENTS
# ============================================================================


top_disagreements <- peptide_data %>%
  arrange(desc(rank_disagreement)) %>%
  slice_head(n = N_DISAGREEMENTS) %>%
  select(peptide_id, sequence, rank_affinityBLOSUM, rank_biophysical,
         rank_disagreement, affinity, blosum62, biophys_distance,
         anchor_intact, everything())

cat("\n=== Top", N_DISAGREEMENTS, "Ranking Disagreements ===\n")
print(top_disagreements %>%
        select(peptide_id, sequence, rank_affinityBLOSUM, rank_biophysical, rank_disagreement))

# ============================================================================
# FEATURE ATTRIBUTION FOR DISAGREEMENTS
# ============================================================================

# Compute feature contributions for top disagreements
top_disagreements_contrib <- top_disagreements %>%
  rowwise() %>%
  mutate(
    contributions = list(compute_feature_contributions(
      cur_data() %>% select(all_of(feature_cols)),
      target_features,
      WEIGHTS
    ))
  ) %>%
  unnest_wider(contributions, names_sep = "_") %>%
  ungroup()

# Reshape for visualization
contrib_long <- top_disagreements_contrib %>%
  select(peptide_id, sequence, rank_affinityBLOSUM, rank_biophysical,
         starts_with("contributions_")) %>%
  pivot_longer(cols = starts_with("contributions_"),
               names_to = "feature",
               values_to = "contribution_pct") %>%
  mutate(feature = str_remove(feature, "contributions_"))

# Visualize feature contributions
p3 <- ggplot(contrib_long,
             aes(x = reorder(peptide_id, -rank_biophysical),
                 y = contribution_pct,
                 fill = feature)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_brewer(palette = "Set3", name = "Feature") +
  labs(
    title = "Feature Contributions to Biophysical Distance",
    subtitle = paste("Top", N_DISAGREEMENTS, "ranking disagreements"),
    x = "Peptide ID",
    y = "Contribution to Distance (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p3)

# ============================================================================
# SELECT CANDIDATES FOR EXPERIMENTAL VALIDATION
# ============================================================================

# Strategy: Combine top candidates from both rankings + high disagreements

top_affinity <- peptide_data %>%
  arrange(rank_affinityBLOSUM) %>%
  slice_head(n = 15) %>%
  mutate(selection_reason = "Top_AffinityBLOSUM")

top_biophysical <- peptide_data %>%
  arrange(rank_biophysical) %>%
  slice_head(n = 15) %>%
  mutate(selection_reason = "Top_Biophysical")

top_disagreement_candidates <- peptide_data %>%
  arrange(desc(rank_disagreement)) %>%
  slice_head(n = 10) %>%
  mutate(selection_reason = "High_Disagreement")

# Combine and remove duplicates
validation_candidates <- bind_rows(
  top_affinity,
  top_biophysical,
  top_disagreement_candidates
) %>%
  distinct(peptide_id, .keep_all = TRUE) %>%
  arrange(rank_affinityBLOSUM)

cat("\n=== Candidates for Experimental Validation ===\n")
cat("Total unique candidates:", nrow(validation_candidates), "\n\n")

# Summary by selection reason
validation_candidates %>%
  count(selection_reason) %>%
  print()

# Visualize validation candidates in ranking space
p4 <- ggplot(peptide_data, aes(x = rank_affinityBLOSUM, y = rank_biophysical)) +
  geom_point(alpha = 0.3, color = "gray70") +
  geom_point(data = validation_candidates,
             aes(color = selection_reason),
             size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_brewer(palette = "Set1", name = "Selection Reason") +
  labs(
    title = "Candidates Selected for Experimental Validation",
    x = "Rank by Affinity + BLOSUM62",
    y = "Rank by Biophysical Distance"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p4)

# ============================================================================
# EXPORT RESULTS
# ============================================================================

# Full dataset with all rankings
write_csv(peptide_data %>%
            select(peptide_id, sequence, anchor_intact,
                   affinity, blosum62, affinityBLOSUM_score,
                   all_of(feature_cols), biophys_distance,
                   rank_affinityBLOSUM, rank_biophysical,
                   rank_disagreement, disagreement_score),
          "full_peptide_rankings.csv")

# Validation candidates
write_csv(validation_candidates %>%
            select(peptide_id, sequence, selection_reason,
                   rank_affinityBLOSUM, rank_biophysical,
                   affinity, blosum62, biophys_distance,
                   anchor_intact),
          "validation_candidates.csv")

# Top disagreements with feature contributions
write_csv(top_disagreements_contrib,
          "top_disagreements_feature_attribution.csv")

cat("\n=== Analysis Complete ===\n")
cat("Results exported to:\n")
cat("  - full_peptide_rankings.csv\n")
cat("  - validation_candidates.csv\n")
cat("  - top_disagreements_feature_attribution.csv\n")
