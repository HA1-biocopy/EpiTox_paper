# =============================================================================
# Bayesian_peptide_assessment.R
# =============================================================================
# Bayesian Peptide Evidence Assessment System
#
# This file contains:
# 1. Main Bayesian assessment function
# 2. Visualization functions
# 3. Report generation
# 4. Complete workflow example
# =============================================================================

library(dplyr)
library(ggplot2)

# Source utility functions
source("utils_exp.R")

# =============================================================================
# MAIN BAYESIAN ASSESSMENT FUNCTION
# =============================================================================

#' Bayesian Evidence Assessment for Peptides
#'
#' Assesses peptides using Bayesian probability updating to answer:
#' "Given the quality/level of experimental evidence, how confident are we
#'  this peptide is meaningfully relevant?"
#'
#' @param peptide_data Preprocessed data frame (from preprocess_for_bayesian)
#' @param target_allele Target HLA allele
#' @param target_class Target HLA class
#' @param prior_prob Starting probability (default: 0.50 - medium prior)
#' @param likelihood_ratios Optional custom likelihood ratios
#'
#' @return Data frame with Bayesian assessment results

bayesian_peptide_assessment <- function(peptide_data,
                                        target_allele,
                                        target_class = "I",
                                        prior_prob = 0.50,
                                        likelihood_ratios = NULL) {

  cat("\n")
  cat("======================================================================\n")
  cat("BAYESIAN PEPTIDE EVIDENCE ASSESSMENT\n")
  cat("======================================================================\n\n")

  cat("Configuration:\n")
  cat("  Target allele:", target_allele, "\n")
  cat("  Target class:", target_class, "\n")
  cat("  Prior probability:", prior_prob, "(", prior_prob*100, "%)\n")
  cat("  Peptides to assess:", nrow(peptide_data), "\n\n")

  # Default likelihood ratios (empirically justified)
  if (is.null(likelihood_ratios)) {
    likelihood_ratios <- list(
      # Evidence from databases/studies
      multiple_studies = 20,    # Very strong: ≥3 independent studies
      two_studies = 10,           # Strong: 2 studies
      one_study = 5,         # Moderate: 1 study
      predicted = 0.2,          # PENALTY: No experimental evidence (< 1.0)

      # HLA allele match
      same_allele = 12,         # Very strong: exact allele match
      different_allele = 2,     # Weak: different allele
      unknown_allele = 1,       # Neutral: no allele info

      # Binding affinity (allele-specific)
      strong_binding = 8,       # Strong experimental binding
      intermediate_binding = 4, # Moderate binding
      weak_binding = 2,         # Weak but detectable
      negative_binding = 0.1,   # STRONG PENALTY: Experimentally negative
      no_binding = 1,           # Neutral: no binding data

      # Tissue/disease context
      normal_tissue = 15,       # Very strong: found in normal (HIGH RISK!)
      same_disease = 12,        # Strong: same disease context
      similar_disease = 5,      # Moderate: related disease
      unknown_tissue = 1,       # Neutral: no tissue info

      # HLA class (only when allele unknown)
      same_class = 3,           # Weak support: same class
      different_class = 0.5     # Penalty: different class
    )
  }

  cat("Likelihood Ratios:\n")
  cat("  Evidence: multiple_studies =", likelihood_ratios$multiple_studies, "\n")
  cat("  Evidence: predicted =", likelihood_ratios$predicted, "(penalty)\n")
  cat("  Allele: same =", likelihood_ratios$same_allele, "\n")
  cat("  Binding: negative =", likelihood_ratios$negative_binding, "(strong penalty)\n")
  cat("  Tissue: normal =", likelihood_ratios$normal_tissue, "\n\n")

  # Initialize output columns
  peptide_data$prior_prob <- prior_prob
  peptide_data$prior_odds <- prior_prob / (1 - prior_prob)
  peptide_data$posterior_odds <- NA
  peptide_data$posterior_prob <- NA
  peptide_data$confidence_level <- NA
  peptide_data$evidence_chain <- ""
  peptide_data$lr_product <- 1
  peptide_data$interpretation <- ""

  # Process each peptide
  cat("Processing peptides...\n")
  pb <- txtProgressBar(min = 0, max = nrow(peptide_data), style = 3)

  for (i in 1:nrow(peptide_data)) {

    odds <- peptide_data$prior_odds[i]
    lr_product <- 1
    evidence_components <- c()

    # -------------------------------------------------------------------------
    # 1. EVIDENCE LEVEL
    # -------------------------------------------------------------------------
    if (!is.na(peptide_data$evidence[i])) {
      evidence_key <- tolower(as.character(peptide_data$evidence[i]))

      if (evidence_key %in% names(likelihood_ratios)) {
        lr <- likelihood_ratios[[evidence_key]]
        odds <- odds * lr
        lr_product <- lr_product * lr
        evidence_components <- c(evidence_components,
                                 sprintf("Evidence[%s:%.1f]", evidence_key, lr))
      }
    }

    # -------------------------------------------------------------------------
    # 2. HLA ALLELE MATCH
    # -------------------------------------------------------------------------
    allele_used <- FALSE
    if (!is.na(peptide_data$hla_allele[i])) {
      if (toupper(as.character(peptide_data$hla_allele[i])) == toupper(target_allele)) {
        lr <- likelihood_ratios$same_allele
        odds <- odds * lr
        lr_product <- lr_product * lr
        evidence_components <- c(evidence_components,
                                 sprintf("Allele[same:%.1f]", lr))
      } else {
        lr <- likelihood_ratios$different_allele
        odds <- odds * lr
        lr_product <- lr_product * lr
        evidence_components <- c(evidence_components,
                                 sprintf("Allele[diff:%.1f]", lr))
      }
      allele_used <- TRUE
    }

    # -------------------------------------------------------------------------
    # 3. BINDING AFFINITY (allele-specific)
    # -------------------------------------------------------------------------
    if (!is.na(peptide_data$binding_affinity[i])) {
      binding_cat <- tolower(as.character(peptide_data$binding_affinity[i]))
      lr_key <- paste0(binding_cat, "_binding")

      if (lr_key %in% names(likelihood_ratios)) {
        lr <- likelihood_ratios[[lr_key]]
        odds <- odds * lr
        lr_product <- lr_product * lr
        evidence_components <- c(evidence_components,
                                 sprintf("Binding[%s:%.1f]", binding_cat, lr))
      }
    }

    # -------------------------------------------------------------------------
    # 4. TISSUE/DISEASE CONTEXT
    # -------------------------------------------------------------------------
    if (!is.na(peptide_data$disease_tissue[i])) {
      tissue_key <- tolower(as.character(peptide_data$disease_tissue[i]))

      if (tissue_key == "same_disease") {
        lr <- likelihood_ratios$same_disease
        tissue_label <- "same"
      } else if (tissue_key == "normal") {
        lr <- likelihood_ratios$normal_tissue
        tissue_label <- "normal"
      } else if (tissue_key == "similar") {
        lr <- likelihood_ratios$similar_disease
        tissue_label <- "similar"
      } else {
        lr <- likelihood_ratios$unknown_tissue
        tissue_label <- "unknown"
      }

      odds <- odds * lr
      lr_product <- lr_product * lr
      evidence_components <- c(evidence_components,
                               sprintf("Tissue[%s:%.1f]", tissue_label, lr))
    }

    # -------------------------------------------------------------------------
    # 5. HLA CLASS (only if allele unknown)
    # -------------------------------------------------------------------------
    if (!allele_used && !is.na(peptide_data$hla_class[i])) {
      if (toupper(as.character(peptide_data$hla_class[i])) == toupper(target_class)) {
        lr <- likelihood_ratios$same_class
        odds <- odds * lr
        lr_product <- lr_product * lr
        evidence_components <- c(evidence_components,
                                 sprintf("Class[same:%.1f]", lr))
      } else {
        lr <- likelihood_ratios$different_class
        odds <- odds * lr
        lr_product <- lr_product * lr
        evidence_components <- c(evidence_components,
                                 sprintf("Class[diff:%.1f]", lr))
      }
    }

    # -------------------------------------------------------------------------
    # FINALIZE CALCULATIONS
    # -------------------------------------------------------------------------
    peptide_data$posterior_odds[i] <- odds
    peptide_data$posterior_prob[i] <- odds / (1 + odds)
    peptide_data$lr_product[i] <- lr_product
    peptide_data$evidence_chain[i] <- paste(evidence_components, collapse = " × ")

    # Assign confidence level and interpretation
    prob <- peptide_data$posterior_prob[i]
    if (prob >= 0.80) {
      peptide_data$confidence_level[i] <- "High"
      peptide_data$interpretation[i] <- "Strong evidence - HIGH CONCERN"
    } else if (prob >= 0.50) {
      peptide_data$confidence_level[i] <- "Medium"
      peptide_data$interpretation[i] <- "Moderate evidence - MONITOR"
    } else if (prob >= 0.30) {
      peptide_data$confidence_level[i] <- "Low"
      peptide_data$interpretation[i] <- "Weak evidence - LOW PRIORITY"
    } else {
      peptide_data$confidence_level[i] <- "Very Low"
      peptide_data$interpretation[i] <- "Insufficient evidence - MINIMAL CONCERN"
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Convert probabilities to percentages
  peptide_data$posterior_prob <- round(peptide_data$posterior_prob * 100, 2)
  peptide_data$prior_prob <- round(peptide_data$prior_prob * 100, 2)

  cat("\n\n=== ASSESSMENT COMPLETE ===\n\n")

  # Summary statistics
  cat("Results Summary:\n")
  summary_table <- summarize_assessment(peptide_data)
  print(summary_table)

  cat("\n")
  cat("======================================================================\n")

  return(peptide_data)
}

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

#' Plot probability distribution by confidence level
#'
#' @param assessed_peptides Output from bayesian_peptide_assessment
#' @return ggplot object
plot_probability_distribution <- function(assessed_peptides) {
  ggplot(assessed_peptides, aes(x = posterior_prob, fill = confidence_level)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 50, linetype = "dashed", color = "black", linewidth = 1) +
    scale_fill_manual(values = c("High" = "#d73027",
                                 "Medium" = "#fee090",
                                 "Low" = "#abd9e9",
                                 "Very Low" = "#4575b4")) +
    labs(title = "Posterior Probability Distribution",
         subtitle = "Higher probability = stronger experimental evidence",
         x = "Posterior Probability (%)",
         y = "Count",
         fill = "Confidence Level") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right")
}

#' Plot prior vs posterior comparison
#'
#' @param assessed_peptides Output from bayesian_peptide_assessment
#' @return ggplot object
plot_prior_vs_posterior <- function(assessed_peptides) {
  assessed_peptides$delta <- assessed_peptides$posterior_prob - assessed_peptides$prior_prob

  ggplot(assessed_peptides, aes(x = prior_prob, y = posterior_prob)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = confidence_level, size = abs(delta)), alpha = 0.6) +
    scale_color_manual(values = c("High" = "#d73027",
                                  "Medium" = "#fee090",
                                  "Low" = "#abd9e9",
                                  "Very Low" = "#4575b4")) +
    scale_size_continuous(range = c(1, 8)) +
    labs(title = "Impact of Evidence on Belief",
         subtitle = "Distance from diagonal = strength of evidence effect",
         x = "Prior Probability (%)",
         y = "Posterior Probability (%)",
         color = "Confidence Level",
         size = "Evidence Impact") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
}

#' Plot evidence vs prediction score (if available)
#'
#' @param assessed_peptides Output from bayesian_peptide_assessment
#' @param presentation_score_col Name of column with presentation score
#' @return ggplot object or NULL
plot_evidence_vs_prediction <- function(assessed_peptides,
                                        presentation_score_col = "presentation_score") {

  if (!presentation_score_col %in% names(assessed_peptides)) {
    cat("Column", presentation_score_col, "not found. Skipping plot.\n")
    return(NULL)
  }

  ggplot(assessed_peptides, aes_string(x = presentation_score_col, y = "posterior_prob")) +
    geom_point(aes(color = confidence_level, size = lr_product), alpha = 0.6) +
    geom_hline(yintercept = 70, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 70, linetype = "dashed", color = "red") +
    annotate("text", x = 85, y = 85,
             label = "HIGH PRIORITY\n(Evidence + Prediction)",
             fontface = "bold", size = 3) +
    annotate("text", x = 85, y = 30,
             label = "Predicted\n(Needs validation)", size = 3) +
    annotate("text", x = 30, y = 85,
             label = "Evidence\n(Check prediction)", size = 3) +
    annotate("text", x = 30, y = 30,
             label = "Low priority", size = 3) +
    scale_color_manual(values = c("High" = "#d73027",
                                  "Medium" = "#fee090",
                                  "Low" = "#abd9e9",
                                  "Very Low" = "#4575b4")) +
    labs(title = "Experimental Evidence vs Computational Prediction",
         subtitle = "Quadrants guide prioritization strategy",
         x = "Presentation Score (Prediction)",
         y = "Posterior Probability (Bayesian Evidence %)",
         color = "Evidence Confidence",
         size = "Evidence Strength") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
}

# =============================================================================
# REPORT GENERATION
# =============================================================================

#' Generate complete assessment report
#'
#' @param assessed_peptides Output from bayesian_peptide_assessment
#' @param output_prefix Prefix for output files
#' @param include_plots Whether to generate plots
#' @return List with summary and high priority peptides
generate_assessment_report <- function(assessed_peptides,
                                       output_prefix = "peptide_assessment",
                                       include_plots = TRUE) {

  cat("\n=== GENERATING ASSESSMENT REPORT ===\n\n")

  # 1. Summary statistics
  summary_table <- summarize_assessment(assessed_peptides)
  openxlsx::write.xlsx(summary_table,
            paste0(output_prefix, "_summary.xlsx"),
            row.names = FALSE)
  cat("✓ Summary statistics saved\n")

  # 2. Full results
  openxlsx::write.xlsx(assessed_peptides,
            paste0(output_prefix, "_full_results.xlsx"),
            row.names = FALSE)
  cat("✓ Full results saved\n")

  # 3. High priority peptides
  high_priority <- assessed_peptides %>%
    filter(confidence_level == "High") %>%
    arrange(desc(posterior_prob))

  openxlsx::write.xlsx(high_priority,
            paste0(output_prefix, "_high_priority.xlsx"),
            row.names = FALSE)
  cat("✓ High priority peptides saved (n =", nrow(high_priority), ")\n")

  # 4. Generate plots
  if (include_plots) {
    p1 <- plot_probability_distribution(assessed_peptides)
    ggsave(paste0(output_prefix, "_distribution.png"), p1, width = 10, height = 6)
    cat("✓ Distribution plot saved\n")

    p2 <- plot_prior_vs_posterior(assessed_peptides)
    ggsave(paste0(output_prefix, "_prior_vs_posterior.png"), p2, width = 10, height = 6)
    cat("✓ Prior vs posterior plot saved\n")
  }

  cat("\n=== REPORT GENERATION COMPLETE ===\n")

  invisible(list(
    summary = summary_table,
    high_priority = high_priority
  ))
}

