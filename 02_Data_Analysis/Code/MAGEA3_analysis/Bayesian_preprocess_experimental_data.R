# =============================================================================
# Bayesian_preprocess_experimental_data.R
# =============================================================================
# Preprocessing Pipeline for Experimental Peptide Data
#
# This file handles:
# 1. Aggregating experimental data (multiple rows per peptide)
# 2. Maintaining allele-binding relationships
# 3. Merging to peptide prediction list
# 4. Formatting for Bayesian assessment
# =============================================================================

library(dplyr)
library(tidyr)

# Source utility functions
source("utils_exp.R")

#' Aggregate experimental evidence by peptide
#'
#' CRITICAL: This function properly handles peptides with multiple experimental
#' entries (different studies, alleles, etc.) and maintains allele-binding relationships
#'
#' @param experimental_df Raw experimental data frame from IEDB with columns:
#'   - id: peptide identifier
#'   - study_ref: study reference ID
#'   - Binding: binding result (Positive-High, Positive-Intermediate, Positive-Low, Negative)
#'   - hla_class: HLA class (I, II, non classical)
#'   - hla_allele: HLA allele (e.g., HLA-A*02:01 or "human")
#'   - disease: disease context
#'   - parent_source_antigen_name: source protein
#'   - Experimental_method: experimental method used
#' @param target_allele Target HLA allele (e.g., "HLA-A*02:01")
#' @param target_disease Target disease for relevance scoring (optional)
#' @param target_class Target HLA class (default: "I")
#'
#' @return Aggregated data frame with one row per peptide

aggregate_experimental_evidence <- function(experimental_df,
                                            target_allele,
                                            target_disease = NULL,
                                            target_class = "I") {

  cat("\n=== AGGREGATING EXPERIMENTAL EVIDENCE ===\n")
  cat("Input rows:", nrow(experimental_df), "\n")
  cat("Unique peptides:", n_distinct(experimental_df$id), "\n")
  cat("Target allele:", target_allele, "\n")
  cat("Target disease:", ifelse(is.null(target_disease), "Not specified", target_disease), "\n")
  cat("Target class:", target_class, "\n\n")

  aggregated <- experimental_df %>%
    group_by(id) %>%
    summarise(
      # -----------------------------------------------------------------------
      # EVIDENCE LEVEL: Count unique studies
      # -----------------------------------------------------------------------
      n_studies = n_distinct(study_ref, na.rm = TRUE),
      study_refs = paste(unique(study_ref), collapse = ";"),

      evidence = case_when(
        n_studies >= 3 ~ "multiple_studies",
        n_studies == 2 ~ "one_study",
        n_studies == 1 ~ "one_database",
        TRUE ~ NA_character_
      ),

      # -----------------------------------------------------------------------
      # HLA ALLELE-SPECIFIC BINDING (CRITICAL: Maintains allele-binding relationship)
      # -----------------------------------------------------------------------
      # Get binding info specific to target allele
      binding_info = list(get_allele_specific_binding(hla_allele, Binding, target_allele)),

      has_target_allele = binding_info[[1]]$has_target_allele,
      target_allele_binding = binding_info[[1]]$target_binding,
      any_allele_binding = binding_info[[1]]$any_binding,
      all_bindings = binding_info[[1]]$all_bindings,

      # Final binding affinity: Use target-specific if available, otherwise any
      binding_affinity = ifelse(!is.na(target_allele_binding),
                                target_allele_binding,
                                any_allele_binding),

      # -----------------------------------------------------------------------
      # HLA ALLELE: Prioritize target allele if present
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
          target_class  # Prioritize target class
        } else {
          first(classes)
        }
      },

      # -----------------------------------------------------------------------
      # DISEASE/TISSUE CONTEXT
      # -----------------------------------------------------------------------
      diseases = paste(unique(disease[disease != "" & !is.na(disease)]), collapse = ";"),
      disease_tissue = categorize_disease(diseases[1], target_disease),

      # -----------------------------------------------------------------------
      # METADATA (for reference)
      # -----------------------------------------------------------------------
      parent_protein = first(parent_source_antigen_name),
      methods = paste(unique(Experimental_method), collapse = ";"),

      .groups = "drop"
    ) %>%
    # Remove the temporary list column
    select(-binding_info)

  cat("\n--- Aggregation Summary ---\n")
  cat("Output rows (unique peptides):", nrow(aggregated), "\n")
  cat("Studies per peptide - Min:", min(aggregated$n_studies),
      "Max:", max(aggregated$n_studies),
      "Median:", median(aggregated$n_studies), "\n")
  cat("Peptides with target allele:", sum(aggregated$has_target_allele, na.rm = TRUE), "\n")
  cat("Peptides with specific alleles:", sum(!is.na(aggregated$hla_allele)), "\n\n")

  cat("Evidence distribution:\n")
  print(table(aggregated$evidence, useNA = "ifany"))
  cat("\nBinding distribution:\n")
  print(table(aggregated$binding_affinity, useNA = "ifany"))

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
      has_experimental_data = !is.na(n_studies) & n_studies > 0,

      # Fill NAs for peptides without experimental data
      evidence = if_else(is.na(evidence), "predicted", evidence),
      n_studies = replace_na(n_studies, 0),
      has_target_allele = replace_na(has_target_allele, FALSE)
    )

  cat("\nMerge results:\n")
  cat("  With experimental data:", sum(merged$has_experimental_data), "\n")
  cat("  Without experimental data (predicted only):", sum(!merged$has_experimental_data), "\n")
  cat("  Match rate:",
      sprintf("%.1f%%", sum(merged$has_experimental_data)/nrow(merged)*100), "\n")

  return(merged)
}

#' Complete preprocessing pipeline
#'
#' @param peptides_df Your peptide list (must have 'id' column)
#' @param experimental_df Raw experimental data from IEDB
#' @param target_allele Target HLA allele
#' @param target_disease Target disease (optional)
#' @param target_class Target HLA class
#' @return Preprocessed data ready for Bayesian assessment

preprocess_for_bayesian <- function(peptides_df,
                                    experimental_df,
                                    target_allele,
                                    target_disease = NULL,
                                    target_class = "I") {

  cat("\n")
  cat("======================================================================\n")
  cat("PEPTIDE EXPERIMENTAL EVIDENCE PREPROCESSING PIPELINE\n")
  cat("======================================================================\n")

  # Step 1: Aggregate experimental evidence
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
      # Ensure required columns exist with proper format
      peptide_id = id,
      sequence = id,  # Assuming id contains sequence; adjust if needed
      binding_metric = NA_character_  # Not using numeric binding in this dataset
    ) %>%
    select(
      # Required columns for Bayesian assessment
      peptide_id,
      sequence,
      evidence,
      hla_allele,
      binding_affinity,
      binding_metric,
      disease_tissue,
      hla_class,

      # Additional useful columns
      has_experimental_data,
      has_target_allele,
      n_studies,
      parent_protein,
      original_diseases = diseases,
      study_references = study_refs,

      # Keep all other original columns from peptides_df
      everything(),
      -id  # Remove duplicate id column
    )

  cat("Final dataset ready for Bayesian assessment!\n")
  cat("Rows:", nrow(final_df), "\n")
  cat("Columns:", ncol(final_df), "\n\n")

  # Run quality checks
  check_preprocessing_quality(final_df)

  cat("\n")
  cat("======================================================================\n")
  cat("PREPROCESSING COMPLETE - Ready for Bayesian assessment\n")
  cat("======================================================================\n\n")

  return(final_df)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Step 1: Load your data
# peptides_df <- read.csv("predicted_peptides.csv")  # Your peptide list
# experimental_df <- read.csv("iedb_experimental_data.csv")  # IEDB data
#
# # Step 2: Run preprocessing
# preprocessed_data <- preprocess_for_bayesian(
#   peptides_df = peptides_df,
#   experimental_df = experimental_df,
#   target_allele = "HLA-A*02:01",
#   target_disease = "melanoma",
#   target_class = "I"
# )
#
# # preprocessed_data is now ready for Bayesian assessment!
