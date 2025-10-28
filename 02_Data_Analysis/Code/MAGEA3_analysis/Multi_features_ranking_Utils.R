# ============================================================================
# FUNCTIONS
# ============================================================================

compute_biophysical_features <- function(peptide_seq, backbone_pos = BACKBONE_POSITIONS) {
  # Global features
  gravy_global <- hydrophobicity(peptide_seq, scale = "KyteDoolittle")
  pI_global <- pI(peptide_seq)
  charge_global <- charge(peptide_seq, pH = 7.4)

  # Aromaticity: fraction of Phe, Trp, Tyr
  aromatic_residues <- c("F", "W", "Y")
  seq_split <- strsplit(peptide_seq, "")[[1]]
  aromat_global <- sum(seq_split %in% aromatic_residues) / length(seq_split)

  # Aliphatic index
  aliph_global <- aIndex(peptide_seq)

  # Backbone-specific features (TCR-facing positions)
  backbone_seq <- substr(peptide_seq, min(backbone_pos), max(backbone_pos))
  gravy_backbone <- hydrophobicity(backbone_seq, scale = "KyteDoolittle")
  charge_backbone <- charge(backbone_seq, pH = 7.4)

  # Backbone aromaticity
  backbone_split <- strsplit(backbone_seq, "")[[1]]
  aromat_backbone <- sum(backbone_split %in% aromatic_residues) / length(backbone_split)

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

compute_conditional_weights <- function(anchor_status) {
  # Define weight schemes based on biological context

  if (anchor_status == "intact") {
    # Anchors intact, only backbone mismatches
    # → HLA binding maintained, TCR recognition is key
    # → Backbone features get highest weight
    weights <- c(
      gravy = 1.0,
      pI = 1.0,
      charge = 1.0,
      aromaticity = 1.0,
      aliphatic_index = 1.0,
      backbone_gravy = 3.0,
      backbone_charge = 3.0,
      backbone_aromaticity = 3.0
    )
  } else if (anchor_status == "anchor_mismatch") {
    # Only anchor mismatches, backbone intact
    # → HLA binding affected, presentation altered
    # → Global features more relevant
    weights <- c(
      gravy = 2.0,
      pI = 2.0,
      charge = 2.0,
      aromaticity = 2.0,
      aliphatic_index = 2.0,
      backbone_gravy = 1.0,
      backbone_charge = 1.0,
      backbone_aromaticity = 1.0
    )
  } else {  # "both" - mismatches everywhere (~80% of data)
    # Mismatches at both anchor and backbone
    # → Complex interplay, no clear dominant factor
    # → Moderate backbone emphasis (conservative)
    weights <- c(
      gravy = 1.0,
      pI = 1.0,
      charge = 1.0,
      aromaticity = 1.0,
      aliphatic_index = 1.0,
      backbone_gravy = 1.5,  # Slight emphasis
      backbone_charge = 1.5,
      backbone_aromaticity = 1.5
    )
  }

  # Normalize to sum to 1
  weights / sum(weights)
}

# Function to apply scaling using pre-computed parameters
apply_scaling <- function(value, mean_val, sd_val) {
  (value - mean_val) / sd_val
}

# ADD THIS FUNCTION:
compute_feature_contributions <- function(peptide_features, target_features, weights) {
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
