# =============================================================================
# preprocess_experimental_data.R - Context-Weighted Version
# =============================================================================
# Preprocessing Pipeline for Experimental Peptide Data
#
# This version implements:
# 1. Context-weighted evidence (allele + normal tissue priority)
# 2. Flexible disease matching (optional target_disease)
# 3. Study relevance scoring
# =============================================================================

library(dplyr)
library(tidyr)

# Source utility functions
#source("utils_exp.R")

#' Aggregate experimental evidence with context weighting
#'
#' CRITICAL: Weights studies by relevance to target context
#' Priority: Target allele > Normal tissue > Disease match
#'
#' @param experimental_df Raw experimental data frame
#' @param target_allele Target HLA allele (e.g., "HLA-A*02:01")
#' @param target_disease Target disease (optional, can be NULL)
#' @param target_class Target HLA class (default: "I")
#'
#' @return Aggregated data frame with context-weighted evidence

aggregate_experimental_evidence <- function(experimental_df,
                                            target_allele,
                                            target_disease = NULL,
                                            target_class = "I") {

  cat("\n=== AGGREGATING EXPERIMENTAL EVIDENCE (CONTEXT-WEIGHTED) ===\n")
  cat("Input rows:", nrow(experimental_df), "\n")
  cat("Unique peptides:", n_distinct(experimental_df$id), "\n")
  cat("Target allele:", target_allele, "\n")
  cat("Target disease:", ifelse(is.null(target_disease), "Not specified (optional)", target_disease), "\n")
  cat("Target class:", target_class, "\n\n")

  # Calculate context relevance for each study entry
  experimental_df <- experimental_df %>%
    mutate(
      context_relevance = mapply(
        calculate_context_relevance,
        study_allele = hla_allele,
        study_is_normal = is_normal,
        study_disease = disease,
        study_class = hla_class,
        MoreArgs = list(
          target_allele = target_allele,
          target_disease = target_disease,
          target_class = target_class
        )
      )
    )

  cat("Context relevance calculated for all entries\n")
  cat("Relevance range:", round(min(experimental_df$context_relevance), 2), "-",
      round(max(experimental_df$context_relevance), 2), "\n\n")

  # Aggregate by peptide
  aggregated <- experimental_df %>%
    group_by(id) %>%
    summarise(
      # -----------------------------------------------------------------------
      # CONTEXT-WEIGHTED EVIDENCE
      # -----------------------------------------------------------------------
      n_studies_total = n_distinct(study_ref, na.rm = TRUE),

      # Count studies by relevance tier
      n_high_relevance = sum(context_relevance >= 1.2, na.rm = TRUE),     # Target allele + context
      n_medium_relevance = sum(context_relevance >= 0.8 & context_relevance < 1.2, na.rm = TRUE),  # Target allele
      n_low_relevance = sum(context_relevance >= 0.4 & context_relevance < 0.8, na.rm = TRUE),     # Partial match
      n_minimal_relevance = sum(context_relevance < 0.4, na.rm = TRUE),   # Weak match

      # Maximum relevance (best study)
      max_relevance = max(context_relevance, na.rm = TRUE),

      # Categorize evidence based on context relevance
      evidence = case_when(
        n_high_relevance >= 2 ~ "multiple_high_relevance",
        n_high_relevance >= 1 ~ "one_high_relevance",
        n_medium_relevance >= 2 ~ "multiple_medium_relevance",
        n_medium_relevance >= 1 ~ "one_medium_relevance",
        n_low_relevance >= 2 ~ "multiple_low_relevance",
        n_low_relevance >= 1 ~ "one_low_relevance",
        TRUE ~ "predicted"
      ),

      study_refs = paste(unique(study_ref), collapse = ";"),

      # -----------------------------------------------------------------------
      # HLA ALLELE-SPECIFIC BINDING
      # -----------------------------------------------------------------------
      binding_info = list(get_allele_specific_binding(hla_allele, Binding, target_allele)),

      has_target_allele = binding_info[[1]]$has_target_allele,
      target_allele_binding = binding_info[[1]]$target_binding,
      any_allele_binding = binding_info[[1]]$any_binding,
      all_bindings = binding_info[[1]]$all_bindings,

      binding_affinity = ifelse(!is.na(target_allele_binding),
                                target_allele_binding,
                                any_allele_binding),

      # -----------------------------------------------------------------------
      # HLA ALLELE
      # -----------------------------------------------------------------------
      hla_allele = get_best_allele(hla_allele, target_allele),

      # -----------------------------------------------------------------------
      # HLA CLASS
      # -----------------------------------------------------------------------
      hla_class = {
        classes <- hla_class[!is.na(hla_class) & hla_class != ""]
        if (length(classes) == 0) {
          NA_character_
        } else if (target_class %in% classes) {
          target_class
        } else {
          first(classes)
        }
      },

      # -----------------------------------------------------------------------
      # DISEASE/TISSUE CONTEXT
      # -----------------------------------------------------------------------
      diseases = paste(unique(disease[disease != "" & !is.na(disease)]), collapse = ";"),

      is_normal_any = any(!is.na(is_normal) &
                            (is_normal == TRUE |
                               tolower(as.character(is_normal)) == "true" |
                               tolower(as.character(is_normal)) == "yes" |
                               is_normal == 1),
                          na.rm = TRUE),

      disease_tissue = categorize_disease(
        disease_string = diseases[1],
        is_normal_flag = is_normal_any,
        target_disease = target_disease
      ),

      # -----------------------------------------------------------------------
      # METADATA
      # -----------------------------------------------------------------------
      parent_protein = first(parent_source_antigen_name),
      methods = paste(unique(Experimental_method), collapse = ";"),

      .groups = "drop"
    ) %>%
    select(-binding_info)

  cat("\n--- Aggregation Summary ---\n")
  cat("Output rows (unique peptides):", nrow(aggregated), "\n")
  cat("Total studies per peptide - Min:", min(aggregated$n_studies_total),
      "Max:", max(aggregated$n_studies_total),
      "Median:", median(aggregated$n_studies_total), "\n")
  cat("Peptides with high relevance studies:", sum(aggregated$n_high_relevance > 0), "\n")
  cat("Peptides with target allele:", sum(aggregated$has_target_allele, na.rm = TRUE), "\n")
  cat("Peptides from normal tissue:", sum(aggregated$is_normal_any, na.rm = TRUE), "\n\n")

  cat("Evidence distribution (context-weighted):\n")
  print(table(aggregated$evidence, useNA = "ifany"))
  cat("\nBinding distribution:\n")
  print(table(aggregated$binding_affinity, useNA = "ifany"))
  cat("\nTissue distribution:\n")
  print(table(aggregated$disease_tissue, useNA = "ifany"))

  return(aggregated)
}

#' Merge experimental evidence to peptide prediction list
#'
#' @param peptides_df Your full peptide list (predictions)
#' @param experimental_summary Aggregated experimental evidence
#' @return Merged data frame ready for Bayesian assessment

merge_to_peptide_list <- function(peptides_df, experimental_summary) {

  cat("\n=== MERGING TO PEPTIDE LIST ===\n")
  cat("Peptide list size:", nrow(peptides_df), "\n")
  cat("Experimental evidence for:", nrow(experimental_summary), "peptides\n")

  merged <- peptides_df %>%
    left_join(experimental_summary, by = "id") %>%
    mutate(
      # Flag whether we have experimental data
      has_experimental_data = !is.na(n_studies_total) & n_studies_total > 0,

      # Fill NAs for peptides without experimental data
      evidence = if_else(is.na(evidence), "predicted", evidence),
      n_studies_total = replace_na(n_studies_total, 0),
      n_high_relevance = replace_na(n_high_relevance, 0),
      n_medium_relevance = replace_na(n_medium_relevance, 0),
      n_low_relevance = replace_na(n_low_relevance, 0),
      has_target_allele = replace_na(has_target_allele, FALSE),

      # Add context relevance category for interpretation
      context_relevance_category = case_when(
        max_relevance >= 1.2 ~ "High (target allele + context)",
        max_relevance >= 0.8 ~ "Medium (target allele)",
        max_relevance >= 0.4 ~ "Low (partial match)",
        max_relevance > 0 ~ "Minimal (weak match)",
        TRUE ~ "None (predicted only)"
      )
    )

  cat("\nMerge results:\n")
  cat("  With experimental data:", sum(merged$has_experimental_data), "\n")
  cat("  Without experimental data (predicted only):", sum(!merged$has_experimental_data), "\n")
  cat("  Match rate:",
      sprintf("%.1f%%", sum(merged$has_experimental_data)/nrow(merged)*100), "\n")

  if (sum(merged$has_experimental_data) > 0) {
    cat("\nContext relevance distribution:\n")
    print(table(merged$context_relevance_category[merged$has_experimental_data], useNA = "ifany"))
  }

  return(merged)
}

#' Complete preprocessing pipeline with context weighting
#'
#' @param peptides_df Your peptide list (must have 'id' column)
#' @param experimental_df Raw experimental data from IEDB
#' @param target_allele Target HLA allele
#' @param target_disease Target disease (OPTIONAL - can be NULL)
#' @param target_class Target HLA class
#' @return Preprocessed data ready for Bayesian assessment

preprocess_for_bayesian <- function(peptides_df,
                                    experimental_df,
                                    target_allele,
                                    target_disease = NULL,
                                    target_class = "I") {

  cat("\n")
  cat("======================================================================\n")
  cat("CONTEXT-WEIGHTED EXPERIMENTAL EVIDENCE PREPROCESSING\n")
  cat("======================================================================\n")
  cat("Priority: Target Allele > Normal Tissue > Disease Match\n")
  cat("======================================================================\n")

  # Step 1: Aggregate experimental evidence with context weighting
  experimental_summary <- aggregate_experimental_evidence(
    experimental_df = experimental_df,
    target_allele = target_allele,
    target_disease = target_disease,
    target_class = target_class
  )

  # Step 2: Merge to peptide list
  preprocessed <- merge_to_peptide_list(
    peptides_df = peptides_df,
    experimental_summary = experimental_summary
  )

  # Step 3: Format for Bayesian assessment
  cat("\n=== FORMATTING FOR BAYESIAN ASSESSMENT ===\n")

  final_df <- preprocessed %>%
    mutate(
      peptide_id = id,
      sequence = id,
      binding_metric = NA_character_
    ) %>%
    select(
      # Required columns
      peptide_id,
      sequence,
      evidence,
      hla_allele,
      binding_affinity,
      binding_metric,
      disease_tissue,
      hla_class,

      # Context relevance info
      has_experimental_data,
      has_target_allele,
      n_studies_total,
      n_high_relevance,
      n_medium_relevance,
      n_low_relevance,
      max_relevance,
      context_relevance_category,

      # Metadata
      parent_protein,
      original_diseases = diseases,
      study_references = study_refs,

      # Keep all other original columns
      everything(),
      -id
    )

  cat("Final dataset ready!\n")
  cat("Rows:", nrow(final_df), "\n")
  cat("Columns:", ncol(final_df), "\n\n")

  # Run quality checks
  check_preprocessing_quality(final_df)

  cat("\n")
  cat("======================================================================\n")
  cat("PREPROCESSING COMPLETE\n")
  cat("======================================================================\n\n")

  return(final_df)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Load your data
# peptides_df <- read.csv("predicted_peptides.csv")
# experimental_df <- read.csv("iedb_data.csv")  # Must have 'is_normal' column!
#
# # Run preprocessing (disease is optional!)
# preprocessed_data <- preprocess_for_bayesian(
#   peptides_df = peptides_df,
#   experimental_df = experimental_df,
#   target_allele = "HLA-A*02:01",
#   target_disease = NULL,  # Optional! Can specify "melanoma" or leave NULL
#   target_class = "I"
# )
