# =============================================================================
# utils_exp.R
# =============================================================================
# Utility Functions for Experimental Evidence Processing
#
# This file contains seven helper functions used by the preprocessing and
# Bayesian assessment scripts.
# =============================================================================

library(dplyr)

#' Extract tissue type from disease name
#'
#' @param disease Disease name string
#' @return Tissue type or NA
extract_tissue_type <- function(disease) {
  tissue_map <- list(
    "melanoma" = "skin",
    "breast" = "breast",
    "lung" = "lung",
    "colon" = "colon|colorectal",
    "ovarian" = "ovarian",
    "leukemia" = "blood|leukemia|lymphoma|myeloma",
    "glioblastoma" = "brain|glioblastoma",
    "pancreatic" = "pancreas",
    "bladder" = "bladder|urinary",
    "cervical" = "cervix|cervical",
    "kidney" = "kidney|renal",
    "pleural" = "pleural|mesothelioma",
    "uterine" = "uterine|uterus",
    "diabetes" = "pancreas",
    "ankylosing" = "spine|joint",
    "uveitis" = "eye",
    "neuroblastoma" = "neural|nervous",
    "meningioma" = "brain|meninges",
    "osteosarcoma" = "bone",
    "neoplasm" = NA_character_  # Generic, can't map to specific tissue
  )

  for (disease_key in names(tissue_map)) {
    if (grepl(disease_key, disease, ignore.case = TRUE)) {
      return(tissue_map[[disease_key]])
    }
  }
  return(NA_character_)
}

#' Categorize disease relevance relative to target
#'
#' @param disease_string Semicolon-separated disease names
#' @param target_disease Target disease for comparison
#' @return Category: "normal", "same_disease", "similar", or NA
categorize_disease <- function(disease_string, target_disease = NULL) {
  if (is.na(disease_string) || disease_string == "") {
    return(NA_character_)
  }

  # CRITICAL: Check explicit normal flag first
  if (!is.na(is_normal_flag) && (is_normal_flag == TRUE ||
                                 tolower(as.character(is_normal_flag)) == "true" ||
                                 tolower(as.character(is_normal_flag)) == "yes" ||
                                 is_normal_flag == 1)) {
    return("normal")
  }

  # If disease string is missing/empty, return NA (not "normal"!)
  if (is.na(disease_string) || disease_string == "") {
    return(NA_character_)
  }

  # Split multiple diseases
  diseases_list <- strsplit(disease_string, ";")[[1]]
  diseases_list <- trimws(diseases_list)

  # Check for normal/healthy tissue (highest priority)
  if (any(grepl("healthy|normal", diseases_list, ignore.case = TRUE))) {
    return("normal")
  }

  # If target disease specified, check for match
  if (!is.null(target_disease) && target_disease != "") {
    # Exact match
    if (any(grepl(target_disease, diseases_list, ignore.case = TRUE))) {
      return("same_disease")
    }

    # Similar disease (same organ/tissue type)
    target_tissue <- extract_tissue_type(target_disease)
    if (!is.na(target_tissue)) {
      if (any(grepl(target_tissue, diseases_list, ignore.case = TRUE))) {
        return("similar")
      }
    }
  }

  # Has disease info but not matching target
  return("similar")
}

#' Get allele-specific binding information
#' CRITICAL: Maintains the relationship between HLA allele and binding result
#'
#' @param alleles Vector of HLA alleles for a peptide
#' @param bindings Vector of binding values for a peptide
#' @param target_allele Target HLA allele to prioritize
#' @return Named list with binding info for target allele
get_allele_specific_binding <- function(alleles, bindings, target_allele) {

  # Initialize result
  result <- list(
    has_target_allele = FALSE,
    target_binding = NA_character_,
    any_binding = NA_character_,
    all_bindings = paste(unique(bindings), collapse = ";")
  )

  # Clean alleles - remove "human" generic entries
  specific_mask <- !is.na(alleles) & tolower(alleles) != "human" & alleles != ""

  if (!any(specific_mask)) {
    # No specific allele information
    # Use any binding info available
    if (length(bindings) > 0 && any(!is.na(bindings))) {
      result$any_binding <- categorize_binding(bindings[!is.na(bindings)])
    }
    return(result)
  }

  # Check if target allele is present
  target_mask <- specific_mask & toupper(alleles) == toupper(target_allele)

  if (any(target_mask)) {
    result$has_target_allele <- TRUE
    # Get binding specifically for target allele
    target_bindings <- bindings[target_mask]
    result$target_binding <- categorize_binding(target_bindings)
  }

  # Also get best binding across all alleles (fallback)
  all_specific_bindings <- bindings[specific_mask]
  if (length(all_specific_bindings) > 0) {
    result$any_binding <- categorize_binding(all_specific_bindings)
  }

  return(result)
}

#' Categorize binding values (handles multiple values, prioritizes strongest)
#'
#' @param binding_values Vector of binding values
#' @return Single categorized binding: "strong", "intermediate", "weak", "negative", or NA
categorize_binding <- function(binding_values) {
  if (length(binding_values) == 0 || all(is.na(binding_values))) {
    return(NA_character_)
  }

  binding_values <- binding_values[!is.na(binding_values)]

  # Prioritize: Negative > Strong > Intermediate > Weak
  # (Negative is explicit experimental result, very informative!)
  if (any(grepl("Negative", binding_values, ignore.case = TRUE))) {
    return("negative")
  }
  if (any(grepl("Positive-High", binding_values, ignore.case = TRUE))) {
    return("strong")
  }
  if (any(grepl("Positive-Intermediate", binding_values, ignore.case = TRUE))) {
    return("intermediate")
  }
  if (any(grepl("Positive-Low", binding_values, ignore.case = TRUE))) {
    return("weak")
  }
  if (any(grepl("Positive", binding_values, ignore.case = TRUE))) {
    return("intermediate")  # Generic "Positive"
  }

  return(NA_character_)
}

#' Get best HLA allele for a peptide (prioritizes target allele)
#'
#' @param alleles Vector of HLA alleles
#' @param target_allele Target HLA allele
#' @return Single best allele or NA
get_best_allele <- function(alleles, target_allele) {
  # Remove generic "human" and NA entries
  specific_alleles <- alleles[!is.na(alleles) &
                                tolower(alleles) != "human" &
                                alleles != ""]

  if (length(specific_alleles) == 0) {
    return(NA_character_)
  }

  # Prioritize target allele if present
  if (toupper(target_allele) %in% toupper(specific_alleles)) {
    return(target_allele)
  }

  # Otherwise return first specific allele
  return(specific_alleles[1])
}

#' Quality check for preprocessed data
#'
#' @param preprocessed_df Preprocessed data frame
#' @return Prints quality report
check_preprocessing_quality <- function(preprocessed_df) {
  cat("=== PREPROCESSING QUALITY CHECKS ===\n\n")

  # Check for required columns
  required_cols <- c("peptide_id", "evidence", "hla_allele", "binding_affinity",
                     "disease_tissue", "hla_class")
  missing_cols <- setdiff(required_cols, names(preprocessed_df))

  if (length(missing_cols) > 0) {
    cat("❌ MISSING COLUMNS:", paste(missing_cols, collapse = ", "), "\n")
  } else {
    cat("✓ All required columns present\n")
  }

  # Check data completeness
  cat("\n--- Data Completeness ---\n")
  for (col in required_cols) {
    if (col %in% names(preprocessed_df)) {
      pct_missing <- round(sum(is.na(preprocessed_df[[col]])) / nrow(preprocessed_df) * 100, 1)
      pct_present <- 100 - pct_missing
      cat(sprintf("  %s: %.1f%% present (%.1f%% missing)\n", col, pct_present, pct_missing))
    }
  }

  # Check for duplicates
  cat("\n--- Duplicate Check ---\n")
  n_duplicates <- sum(duplicated(preprocessed_df$peptide_id))
  cat("  Duplicate peptide IDs:", n_duplicates, "\n")
  if (n_duplicates > 0) {
    cat("  ⚠️  WARNING: Duplicates found! Check aggregation.\n")
  }

  # Evidence distribution
  cat("\n--- Evidence Quality ---\n")
  cat("  Evidence levels:\n")
  print(table(preprocessed_df$evidence, useNA = "ifany"))

  # Binding distribution
  cat("\n--- Binding Information ---\n")
  cat("  Binding categories:\n")
  print(table(preprocessed_df$binding_affinity, useNA = "ifany"))

  # HLA information
  cat("\n--- HLA Information ---\n")
  cat("  Peptides with specific alleles:",
      sum(!is.na(preprocessed_df$hla_allele)),
      sprintf("(%.1f%%)", sum(!is.na(preprocessed_df$hla_allele))/nrow(preprocessed_df)*100), "\n")
  cat("  Peptides matching target allele:",
      sum(preprocessed_df$has_target_allele, na.rm = TRUE), "\n")
  cat("  Class distribution:\n")
  print(table(preprocessed_df$hla_class, useNA = "ifany"))

  # Disease/tissue
  cat("\n--- Disease/Tissue Context ---\n")
  cat("  Tissue categories:\n")
  print(table(preprocessed_df$disease_tissue, useNA = "ifany"))

  cat("\n=== END QUALITY CHECKS ===\n")

  invisible(TRUE)
}

#' Summary statistics for Bayesian assessment results
#'
#' @param assessed_peptides Data frame from bayesian_peptide_assessment
#' @return Data frame with summary statistics
summarize_assessment <- function(assessed_peptides) {
  summary <- assessed_peptides %>%
    group_by(confidence_level) %>%
    summarise(
      Count = n(),
      `Mean Posterior (%)` = round(mean(posterior_prob), 1),
      `Median Posterior (%)` = round(median(posterior_prob), 1),
      `SD` = round(sd(posterior_prob), 1),
      `Min` = round(min(posterior_prob), 1),
      `Max` = round(max(posterior_prob), 1),
      .groups = "drop"
    ) %>%
    arrange(desc(`Mean Posterior (%)`))

  return(summary)
}
