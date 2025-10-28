# ============================================================================
# Multi-Feature Biophysical Scoring for TCRm Off-Target Assessment
# Lean, Interpretable Approach
# ============================================================================

library(Peptides)
library(tidyverse)
library(corrplot)
library(ggrepel)

source("Multi_features_ranking_Utils.R")
# ============================================================================
# CONFIGURATION
# ============================================================================

# Top N
N_DISAGREEMENTS <- 20
N_VALIDATION <- 30

# Define your target peptide
TARGET_PEPTIDE <- "KVAELVHFL"  # REPLACE with your actual target

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
# LOAD AND PREPARE DATA
# ============================================================================

# Load data
known_ot = read.csv("../../data/known_patents_OTs.csv") %>%
  mutate(peptide = ifelse(peptide == "TVAELVQFV", "TVAELVQFL", peptide))
bayesian = openxlsx::read.xlsx(("../../data/Bayesian/MAGEA3_full_results.xlsx")) %>%
  dplyr::select(peptide_id, confidence_level, evidence_chain) %>%
  dplyr::rename(id = peptide_id)

positions_categories = openxlsx::read.xlsx("../../data/OT_positions_classification.xlsx") %>%
  mutate(anchor_status = case_when(
    Category == "Anchor Only" ~ "anchor_mismatch",
    Category == "Perfect Match" ~ "Intact",
    Category == "Backbone Only" ~ "Intact",
    TRUE ~ Category),
    id = paste0(uniprot, "_", peptide)) %>%
  dplyr::select(id, anchor_status)

epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  dplyr::rename(GeneName = Gene.Names) %>%
  dplyr::select(uniprot, GeneName, peptide, Wildtype, blosum_similarity, mismatch,
         affinity, presentation_score, Ranking_score, Rank) %>%
  mutate(id = paste0(uniprot, "_", peptide),
         known_peptide = ifelse(peptide %in% known_ot$peptide, "Known", "unkwon"),
         Rank = ifelse(Rank == "Random", "Low", Rank)) %>%
  distinct(id, .keep_all = T) %>%
  dplyr::rename(Bi_features_score = Ranking_score) %>%
  merge(., positions_categories, by ="id") %>%
  merge(., bayesian, by ="id") %>%
  relocate(anchor_status, .after = "mismatch")

#============================================================================
# COMPUTE MULTI-BIOPHYSICAL FEATURES
# ============================================================================

# Target peptide features
target_features <- compute_biophysical_features(TARGET_PEPTIDE)

# All peptides - compute features
# Use group_by(id) with group_modify to handle duplicates properly
peptide_data <- epitox %>%
  group_by(id) %>%
  group_modify(~ {
    features <- compute_biophysical_features(.x$peptide[1])
    bind_cols(.x, features)
  }) %>%
  ungroup()

# Multi-biophysical feature columns
multibiophys_cols <- c("gravy", "pI", "charge", "aromaticity", "aliphatic_index",
                       "backbone_gravy", "backbone_charge", "backbone_aromaticity")

# Calculate scaling parameters from all peptides
scaling_params <- peptide_data %>%
  dplyr::select(all_of(multibiophys_cols)) %>%
  summarise(across(everything(),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))))

# Scale peptide features using computed parameters
peptide_data <- peptide_data %>%
  mutate(
    gravy_scaled = apply_scaling(gravy,
                                 scaling_params$gravy_mean,
                                 scaling_params$gravy_sd),
    pI_scaled = apply_scaling(pI,
                              scaling_params$pI_mean,
                              scaling_params$pI_sd),
    charge_scaled = apply_scaling(charge,
                                  scaling_params$charge_mean,
                                  scaling_params$charge_sd),
    aromaticity_scaled = apply_scaling(aromaticity,
                                       scaling_params$aromaticity_mean,
                                       scaling_params$aromaticity_sd),
    aliphatic_index_scaled = apply_scaling(aliphatic_index,
                                           scaling_params$aliphatic_index_mean,
                                           scaling_params$aliphatic_index_sd),
    backbone_gravy_scaled = apply_scaling(backbone_gravy,
                                          scaling_params$backbone_gravy_mean,
                                          scaling_params$backbone_gravy_sd),
    backbone_charge_scaled = apply_scaling(backbone_charge,
                                           scaling_params$backbone_charge_mean,
                                           scaling_params$backbone_charge_sd),
    backbone_aromaticity_scaled = apply_scaling(backbone_aromaticity,
                                                scaling_params$backbone_aromaticity_mean,
                                                scaling_params$backbone_aromaticity_sd)
  )

# Scale target features using the SAME parameters
target_features_scaled <- target_features %>%
  mutate(
    gravy_scaled = apply_scaling(gravy,
                                 scaling_params$gravy_mean,
                                 scaling_params$gravy_sd),
    pI_scaled = apply_scaling(pI,
                              scaling_params$pI_mean,
                              scaling_params$pI_sd),
    charge_scaled = apply_scaling(charge,
                                  scaling_params$charge_mean,
                                  scaling_params$charge_sd),
    aromaticity_scaled = apply_scaling(aromaticity,
                                       scaling_params$aromaticity_mean,
                                       scaling_params$aromaticity_sd),
    aliphatic_index_scaled = apply_scaling(aliphatic_index,
                                           scaling_params$aliphatic_index_mean,
                                           scaling_params$aliphatic_index_sd),
    backbone_gravy_scaled = apply_scaling(backbone_gravy,
                                          scaling_params$backbone_gravy_mean,
                                          scaling_params$backbone_gravy_sd),
    backbone_charge_scaled = apply_scaling(backbone_charge,
                                           scaling_params$backbone_charge_mean,
                                           scaling_params$backbone_charge_sd),
    backbone_aromaticity_scaled = apply_scaling(backbone_aromaticity,
                                                scaling_params$backbone_aromaticity_mean,
                                                scaling_params$backbone_aromaticity_sd)
  )

# ============================================================================
# COMPUTE MULTI-BIOPHYSICAL DISTANCE (ANCHOR STATUS-AWARE)
# ============================================================================

peptide_data <- peptide_data %>%
  rowwise() %>%
  mutate(
    # Get conditional weights for this peptide based on anchor status
    conditional_weights = list(compute_conditional_weights(anchor_status)),

    # Compute multi-biophysical distance using conditional weights
    multibiophys_distance = {
      feature_names <- names(unlist(conditional_weights))
      weights_vec <- unlist(conditional_weights)

      # Calculate weighted distance
      squared_diffs <- sapply(feature_names, function(feat) {
        feat_scaled <- paste0(feat, "_scaled")
        diff <- cur_data()[[feat_scaled]] - target_features_scaled[[feat_scaled]]
        return(diff^2)
      })

      sqrt(sum(weights_vec * squared_diffs))
    }
  ) %>%
  ungroup()

# ============================================================================
# NORMALIZE SCORES TO COMPARABLE 0-1 SCALE
# ============================================================================

peptide_data <- peptide_data %>%
  mutate(
    # Bi_features is already 0-1 (keep as is)
    #Bi_features_score = Bi_features,

    # Transform multibiophys_distance to 0-1 similarity score
    # (inverted so higher = more similar = higher priority, matching Bi_features direction)
    multibiophys_similarity = 1 - (
      (multibiophys_distance - min(multibiophys_distance)) /
        (max(multibiophys_distance) - min(multibiophys_distance))
    ),

    # Compute disagreement on 0-1 scale
    score_disagreement = abs(Bi_features_score - multibiophys_similarity),

    # Also create ranks for selection/visualization
    rank_Bi_features = rank(-Bi_features_score, ties.method = "min"),
    rank_multibiophys = rank(-multibiophys_similarity, ties.method = "min"),
    rank_disagreement = abs(rank_Bi_features - rank_multibiophys)
  )

# ============================================================================
# EXPLORATORY ANALYSIS
# ============================================================================
biocopy_colors = c("#A2C510", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")

cat("=== Summary Statistics ===\n")
cat("Total peptides:", nrow(peptide_data), "\n")
cat("Correlation between Bi_features and Multi-biophysical rankings:",
    cor(peptide_data$rank_Bi_features, peptide_data$rank_multibiophys, method = "spearman"), "\n\n")

cat("=== Distribution by Anchor Status ===\n")
peptide_data %>%
  count(anchor_status) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  print()

ggplot(peptide_data, aes(Bi_features_score, multibiophys_distance)) +
  geom_point() +
  geom_smooth() +
  theme_light()

cor(peptide_data$Bi_features_score, peptide_data$multibiophys_distance, method = "spearman")
#========================================================================
# Paper summary
#========================================================================
peptide_data %>%
  dplyr::select(Rank) %>%
  table %>%
  prop.table()

peptide_data %>%
  dplyr::select(Bi_features_score) %>%
  summary()

peptide_data %>%
  group_by(Rank) %>%
  summarise(mu = median(Bi_features_score))

peptide_data %>%
  dplyr::select(Rank, Wildtype) %>%
  table %>%
  prop.table()

peptide_data %>%
  dplyr::select(Rank, confidence_level) %>%
  table %>%
  prop.table()

peptide_data %>%
  filter(Wildtype == "No") %>%
  group_by(Rank) %>%
  summarise(mu = median(Bi_features_score))

ggplot(peptide_data, aes(confidence_level, Bi_features_score, fill = Rank)) +
  geom_boxplot(position = "dodge") +
  theme_light() +
  labs(x = "") +
  ggtitle("Bayesian classification vs. Bi-features") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = biocopy_colors)

ggplot(peptide_data) +
  geom_histogram(aes(multibiophys_similarity), fill = biocopy_colors[1], alpha = 0.7, col = "white") +
  geom_histogram(aes(Bi_features_score), fill = biocopy_colors[2], alpha = 0.7, col = "white") +
  labs(x = "Ranking scores distribution", caption = "Green = Bi features\nBlue = Multibiophys") +
  theme_light()

# Compare target to dataset mean for each feature
target_comparison <- data.frame(
  feature = multibiophys_cols,
  target_value = unlist(target_features[multibiophys_cols]),
  mean_value = colMeans(peptide_data[multibiophys_cols], na.rm = TRUE),
  sd_value = sapply(peptide_data[multibiophys_cols], sd, na.rm = TRUE)
) %>%
  mutate(
    z_score = (target_value - mean_value) / sd_value,
    interpretation = case_when(
      abs(z_score) > 2 ~ "Extreme",
      abs(z_score) > 1 ~ "Unusual",
      TRUE ~ "Typical"
    )
  )

# What are the multibiophys_similarity scores for top Bi_features peptides?
top_bi <- peptide_data %>%
  arrange(desc(Bi_features_score)) %>%
  slice_head(n = 50)

summary(top_bi$multibiophys_similarity)

# Compare to full dataset
cat("Top 50 Bi_features peptides:\n")
cat("  Mean multibiophys_similarity:", mean(top_bi$multibiophys_similarity), "\n")
cat("Full dataset:\n")
cat("  Mean multibiophys_similarity:", mean(peptide_data$multibiophys_similarity), "\n")

print(target_comparison)
# Feature correlations
cat("\n=== Multi-Biophysical Feature Correlations ===\n")
cor_matrix <- peptide_data %>%
  dplyr::select(all_of(multibiophys_cols), affinity, blosum_similarity) %>%
  cor(use = "complete.obs")

print(round(cor_matrix, 2))

# Visualize correlations
corrplot(cor_matrix, method = "circle", type = "upper",
         col = colorRampPalette(c("#99CFE9", "white", "#A2C510"))(200),
         title = "Multi-Biophysical Feature Correlations",
         mar = c(0, 0, 2, 0))

# Feature distributions by anchor status
peptide_data %>%
  pivot_longer(cols = all_of(multibiophys_cols),
               names_to = "feature",
               values_to = "value") %>%
  ggplot(aes(x = anchor_status, y = value, fill = anchor_status)) +
  geom_boxplot() +
  facet_wrap(~feature, scales = "free_y", ncol = 4) +
  labs(title = "Multi-Biophysical Features by Anchor Status",
       x = "Anchor Status", y = "Feature Value") +
  theme_light() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("Both Anchor & Backbone" = "#FBB800",
                               "Intact" = "#99CFE9", "anchor_mismatch" = "#C61E19"))

# Compare rankings within each anchor status group
anchor_comparison <- peptide_data %>%
  group_by(anchor_status) %>%
  summarise(
    n = n(),
    mean_affinity = mean(affinity, na.rm = TRUE),
    mean_blosum = mean(blosum_similarity, na.rm = TRUE),
    mean_Bi_features = mean(Bi_features_score, na.rm = TRUE),
    mean_multibiophys_distance = mean(multibiophys_distance, na.rm = TRUE),
    correlation_ranks = cor(rank_Bi_features, rank_multibiophys)
  )

cat("\n=== Ranking Correlation by Anchor Status ===\n")
print(anchor_comparison)

# Visualize distance distributions by anchor status
ggplot(peptide_data, aes(x = anchor_status, y = multibiophys_distance, fill = anchor_status)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.9) +
  labs(
    title = "Multi-Biophysical Distance Distribution by Anchor Status",
    subtitle = "Conditional weighting schemes applied based on mismatch location",
    x = "Anchor Status",
    y = "Multi-Biophysical Distance from Target"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values =c("Both Anchor & Backbone" = "#FBB800",
                              "Intact" = "#99CFE9", "anchor_mismatch" = "#C61E19"))

# ============================================================================
# RANKING COMPARISON: BI_FEATURES VS MULTI-BIOPHYSICAL
# ============================================================================

# Scatter plot: two ranking methods
p1 <- ggplot(peptide_data, aes(x = rank_Bi_features, y = rank_multibiophys)) +
  geom_point(aes(color = score_disagreement), alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  #scale_color_gradient(low = "#A2C510", high = "#EAECEE", name = "Disagreement") +
  scale_color_gradient(low = "#EAECEE", high = "#A2C510", name = "Disagreement") +
  labs(
    title = "Bi_features vs Multi-Biophysical Ranking",
    subtitle = "Points near diagonal = agreement; far from diagonal = disagreement",
    x = "Rank by Bi_features (Affinity + BLOSUM62)",
    y = "Rank by Multi-Biophysical Distance"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p1)

# Distribution of disagreement scores
p2 <- ggplot(peptide_data, aes(x = score_disagreement)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of Ranking Disagreements",
    subtitle = "Between Bi_features and Multi-Biophysical methods",
    x = "Disagreement Score (0 = perfect agreement, 1 = maximum disagreement)",
    y = "Count"
  ) +
  theme_minimal()

print(p2)

# Disagreement by anchor status
p3 <- ggplot(peptide_data, aes(x = anchor_status, y = score_disagreement, fill = anchor_status)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.9) +
  labs(
    title = "Ranking Disagreement by Anchor Status",
    subtitle = "Do certain mismatch types show more disagreement?",
    caption = "Disagreement = abs(Bi_features - Multi_features)",
    x = "Anchor Status",
    y = "Disagreement Score"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values =c("Both Anchor & Backbone" = "#FBB800",
                              "Intact" = "#99CFE9", "anchor_mismatch" = "#C61E19"))

print(p3)

# ============================================================================
# IDENTIFY TOP DISAGREEMENTS
# ============================================================================

top_disagreements <- peptide_data %>%
  arrange(desc(rank_disagreement)) %>%
  slice_head(n = N_DISAGREEMENTS) %>%
  dplyr::select(uniprot, peptide, anchor_status,
         rank_Bi_features, rank_multibiophys, rank_disagreement,
         affinity, blosum_similarity, Bi_features_score, multibiophys_distance,
         everything())

cat("\n=== Top", N_DISAGREEMENTS, "Ranking Disagreements ===\n")
print(top_disagreements %>%
        dplyr::select(uniprot, peptide, anchor_status,
               rank_Bi_features, rank_multibiophys, rank_disagreement))

# Top disagreements within each anchor status group
top_disagreements_by_anchor <- peptide_data %>%
  group_by(anchor_status) %>%
  arrange(desc(rank_disagreement)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  dplyr::select(uniprot, peptide, anchor_status,
         rank_Bi_features, rank_multibiophys, rank_disagreement)

cat("\n=== Top Disagreements by Anchor Status ===\n")
print(top_disagreements_by_anchor)

# ============================================================================
# FEATURE ATTRIBUTION FOR DISAGREEMENTS
# ============================================================================

# Compute feature contributions for top score disagreements
top_disagreements_contrib <- top_disagreements %>%
  rowwise() %>%
  mutate(
    contributions = list(compute_feature_contributions(
      cur_data() %>%
        dplyr::select(all_of(multibiophys_cols)),
      target_features,
      unlist(conditional_weights)
    ))
  ) %>%
  unnest_wider(contributions, names_sep = "_") %>%
  ungroup()

# Reshape for visualization
contrib_long <- top_disagreements_contrib %>%
  dplyr::select(id, peptide, Bi_features_score, multibiophys_similarity,
         starts_with("contributions_")) %>%
  pivot_longer(cols = starts_with("contributions_"),
               names_to = "feature",
               values_to = "contribution_pct") %>%
  mutate(feature = str_remove(feature, "contributions_"))

# Visualize feature contributions
p4 <- ggplot(contrib_long,
             aes(x = reorder(id, -multibiophys_similarity),
                 y = contribution_pct,
                 fill = feature)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_brewer(palette = "Set3", name = "Feature") +
  labs(
    title = "Feature Contributions to Multi-Biophysical Distance",
    subtitle = paste("Top", N_DISAGREEMENTS, "score disagreements"),
    x = "Peptide ID",
    y = "Contribution to Distance (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p4)

# ============================================================================
# SELECT CANDIDATES FOR EXPERIMENTAL VALIDATION
# ============================================================================


# Strategy: Combine top candidates from both rankings + high disagreements

top_Bi_features <- peptide_data %>%
  arrange(rank_Bi_features) %>%
  slice_head(n = 15) %>%
  mutate(selection_reason = "Top_Bi_features")

top_multibiophys <- peptide_data %>%
  arrange(rank_multibiophys) %>%
  slice_head(n = 15) %>%
  mutate(selection_reason = "Top_MultiBiophysical")

top_disagreement_candidates <- peptide_data %>%
  arrange(desc(rank_disagreement)) %>%
  slice_head(n = 10) %>%
  mutate(selection_reason = "High_Disagreement")

# Combine and remove duplicates
validation_candidates <- bind_rows(
  top_Bi_features,
  top_multibiophys,
  top_disagreement_candidates
) %>%
  distinct(id, .keep_all = TRUE) %>%
  arrange(rank_Bi_features)

cat("\n=== Candidates for Experimental Validation ===\n")
cat("Total unique candidates:", nrow(validation_candidates), "\n\n")

# Summary by selection reason
validation_candidates %>%
  count(selection_reason) %>%
  print()

# Visualize validation candidates in ranking space
p5 <- ggplot(peptide_data, aes(x = rank_Bi_features, y = rank_multibiophys)) +
  geom_point(alpha = 0.3, color = "gray70") +
  geom_point(data = validation_candidates,
             aes(color = selection_reason),
             size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  #scale_color_brewer(palette = "Set1", name = "Selection Reason") +
  labs(
    title = "Top Candidates for Experimental Validation",
    subtitle = "Comparing Bi_features and Multi-Biophysical rankings",
    x = "Rank by Bi_features (Affinity + BLOSUM62)",
    y = "Rank by Multi-Biophysical Distance"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = biocopy_colors)


ggplot(peptide_data, aes(x = rank_Bi_features, y = rank_multibiophys)) +
  geom_point(alpha = 0.3, color = "gray70") +
  geom_point(data = validation_candidates,
             aes(color = selection_reason),
             size = 3, alpha = 0.8) +
  geom_text_repel(data = validation_candidates,
                  aes(label = id),
                  size = 3,
                  max.overlaps = 11,
                  box.padding = 0.5,
                  point.padding = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Top Candidates for Experimental Validation",
    subtitle = "Comparing Bi_features and Multi-Biophysical rankings",
    x = "Rank by Bi_features (Affinity + BLOSUM62)",
    y = "Rank by Multi-Biophysical Distance"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = biocopy_colors)
print(p5)

# ============================================================================
# EXPORT RESULTS
# ============================================================================

# Full dataset with all rankings
openxlsx::write.xlsx(peptide_data %>%
            select(id, mismatch, Wildtype,
                   anchor_status, known_peptide,
                   confidence_level, evidence_chain,
                   affinity, blosum_similarity, Bi_features_score,
                   all_of(multibiophys_cols), multibiophys_distance,
                   rank_Bi_features, rank_multibiophys,
                   rank_disagreement, score_disagreement),
          "../../data/full_peptide_rankings.xlsx")

# Validation candidates
openxlsx::write.xlsx(validation_candidates %>%
            select(id, selection_reason, anchor_status,
                   rank_Bi_features, rank_multibiophys,
                   affinity, blosum_similarity, Bi_features_score,
                   multibiophys_distance),
          "../../data/validation_candidates.xlsx")

# Top disagreements with feature contributions
openxlsx::write.xlsx(top_disagreements_contrib,
          "../../data/top_disagreements_feature_attribution.xlsx")

cat("\n=== Analysis Complete ===\n")
cat("Results exported to:\n")
cat("  - full_peptide_rankings.xlsx\n")
cat("  - validation_candidates.xlsx\n")
cat("  - top_disagreements_feature_attribution.xlsx\n")
