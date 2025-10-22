# =============================================================================
# bayesian_peptide_assessment.R - Context-Weighted Version
# =============================================================================
# Bayesian Peptide Evidence Assessment System
#
# This version implements Option A:
# - Strong LRs for clear discrimination
# - Context-weighted evidence (allele + normal tissue priority)
# - Accept binary-like behavior for experimental evidence
# - Provide both posterior_prob AND rank_score for flexibility
# =============================================================================

library(dplyr)
library(ggplot2)

# Source utility functions
source("utils_exp.R")

# =============================================================================
# MAIN BAYESIAN ASSESSMENT FUNCTION
# =============================================================================

#' Bayesian Evidence Assessment with Context Weighting
#'
#' @param peptide_data Preprocessed data frame (from preprocess_for_bayesian)
#' @param target_allele Target HLA allele
#' @param target_class Target HLA class
#' @param prior_prob Starting probability (default: 0.50)
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
  cat("CONTEXT-WEIGHTED BAYESIAN EVIDENCE ASSESSMENT\n")
  cat("======================================================================\n\n")

  cat("Configuration:\n")
  cat("  Target allele:", target_allele, "\n")
  cat("  Target class:", target_class, "\n")
  cat("  Prior probability:", prior_prob, "(", prior_prob*100, "%)\n")
  cat("  Peptides to assess:", nrow(peptide_data), "\n\n")

  # Default likelihood ratios (Option A - Strong discrimination)
  if (is.null(likelihood_ratios)) {
    likelihood_ratios <- list(
      # Context-weighted evidence (Allele + Normal tissue priority)
      multiple_high_relevance = 12,      # 2+ studies in target context
      one_high_relevance = 8,            # 1 study in target context
      multiple_medium_relevance = 6,     # 2+ studies with target allele
      one_medium_relevance = 4,          # 1 study with target allele
      multiple_low_relevance = 3,        # 2+ studies partial match
      one_low_relevance = 2,             # 1 study partial match
      predicted = 0.3,                    # No experimental evidence

      # HLA allele match (MOST IMPORTANT per user priority)
      same_allele = 10,
      different_allele = 2,
      unknown_allele = 1,

      # Binding affinity - Two-level assessment
      # Peptide state (base)
      strong_binder = 5,
      intermediate_binder = 3,
      weak_binder = 2,
      non_binder = 1,

      # Target multipliers
      strong_on_target = 2.5,
      intermediate_on_target = 1.5,
      weak_on_target = 1.0,
      negative_on_target = 0.2,
      unknown_target = 1.0,

      # Tissue/disease context (Normal tissue = HIGH RISK per user priority)
      normal_tissue = 10,
      same_disease = 6,
      similar_disease = 3,
      unknown_tissue = 1,

      # HLA class (only when allele unknown)
      same_class = 2,
      different_class = 0.5
    )
  }

  cat("Likelihood Ratios (Context-Weighted):\n")
  cat("  Evidence: multiple_high_relevance =", likelihood_ratios$multiple_high_relevance, "\n")
  cat("  Evidence: one_high_relevance =", likelihood_ratios$one_high_relevance, "\n")
  cat("  Evidence: predicted =", likelihood_ratios$predicted, "(penalty)\n")
  cat("  Allele: same =", likelihood_ratios$same_allele, "(PRIORITY)\n")
  cat("  Binding: strong × strong_on_target = max",
      likelihood_ratios$strong_binder * likelihood_ratios$strong_on_target, "\n")
  cat("  Tissue: normal =", likelihood_ratios$normal_tissue, "(HIGH RISK)\n\n")

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
    # 1. CONTEXT-WEIGHTED EVIDENCE
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
    # 2. HLA ALLELE MATCH (MOST IMPORTANT)
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
    # 3. BINDING AFFINITY (Two-level: peptide state × target status)
    # -------------------------------------------------------------------------
    if (!is.na(peptide_data$binding_affinity[i])) {
      binding_cat <- tolower(as.character(peptide_data$binding_affinity[i]))

      # Get base LR from peptide binding state
      base_lr <- switch(binding_cat,
                        "strong" = likelihood_ratios$strong_binder,
                        "intermediate" = likelihood_ratios$intermediate_binder,
                        "weak" = likelihood_ratios$weak_binder,
                        "negative" = likelihood_ratios$non_binder,
                        1)

      # Get target multiplier (if we know target allele status)
      target_mult <- 1.0
      if (!is.na(peptide_data$has_target_allele[i]) && peptide_data$has_target_allele[i]) {
        # We have target-specific binding info
        target_mult <- switch(binding_cat,
                              "strong" = likelihood_ratios$strong_on_target,
                              "intermediate" = likelihood_ratios$intermediate_on_target,
                              "weak" = likelihood_ratios$weak_on_target,
                              "negative" = likelihood_ratios$negative_on_target,
                              1.0)
      } else {
        # No target-specific info, use neutral multiplier
        target_mult <- likelihood_ratios$unknown_target
      }

      final_binding_lr <- base_lr * target_mult
      odds <- odds * final_binding_lr
      lr_product <- lr_product * final_binding_lr
      evidence_components <- c(evidence_components,
                               sprintf("Binding[%s×%.1f=%.1f]", binding_cat, target_mult, final_binding_lr))
    }

    # -------------------------------------------------------------------------
    # 4. TISSUE/DISEASE CONTEXT (Normal tissue = HIGH RISK)
    # -------------------------------------------------------------------------
    if (!is.na(peptide_data$disease_tissue[i])) {
      tissue_key <- tolower(as.character(peptide_data$disease_tissue[i]))

      if (tissue_key == "same_disease") {
        lr <- likelihood_ratios$same_disease
        tissue_label <- "same"
      } else if (tissue_key == "normal") {
        lr <- likelihood_ratios$normal_tissue
        tissue_label <- "normal"
      } else if (tissue_key == "similar_tissue") {
        lr <- likelihood_ratios$similar_disease
        tissue_label <- "similar_tissue"
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
plot_probability_distribution <- function(assessed_peptides) {
  ggplot(assessed_peptides, aes(x = posterior_prob, fill = confidence_level)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 50, linetype = "dashed", color = "black", size = 1) +
    scale_fill_manual(values = c("High" = "#C61E19",
                                 "Medium" = "#FDE399",
                                 "Low" = "#8ECAE6",
                                 "Very Low" = "#A2C510")) +
    labs(title = "Posterior Probability Distribution (Context-Weighted)",
         subtitle = "Higher probability = stronger experimental evidence in target context",
         x = "Posterior Probability (%)",
         y = "Count",
         fill = "Confidence Level") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
}

#' Plot prior vs posterior comparison
plot_prior_vs_posterior <- function(assessed_peptides) {
  assessed_peptides$delta <- assessed_peptides$posterior_prob - assessed_peptides$prior_prob

  ggplot(assessed_peptides, aes(x = prior_prob, y = posterior_prob)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = confidence_level, size = abs(delta)), alpha = 0.6) +
    scale_color_manual(values = c("High" = "#C61E19",
                                  "Medium" = "#FDE399",
                                  "Low" = "#8ECAE6",
                                  "Very Low" = "#A2C510")) +
    scale_size_continuous(range = c(1, 8)) +
    labs(title = "Impact of Context-Weighted Evidence",
         subtitle = "Distance from diagonal = strength of evidence effect",
         x = "Prior Probability (%)",
         y = "Posterior Probability (%)",
         color = "Confidence Level",
         size = "Evidence Impact") +
    theme_minimal()
}

#' Generate complete assessment report
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

  # 4. Context relevance breakdown
  if ("context_relevance_category" %in% names(assessed_peptides)) {
    context_summary <- assessed_peptides %>%
      filter(has_experimental_data == TRUE) %>%
      group_by(context_relevance_category, confidence_level) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(desc(count))

    openxlsx::write.xlsx(context_summary,
              paste0(output_prefix, "_context_summary.xlsx"),
              row.names = FALSE)
    cat("✓ Context relevance summary saved\n")
  }

  # 5. Generate plots
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

# =============================================================================
# COMPLETE WORKFLOW EXAMPLE
# =============================================================================

# # STEP 1: Load all scripts
# source("utils_exp.R")
# source("preprocess_experimental_data.R")
# source("bayesian_peptide_assessment.R")
#
# # STEP 2: Load your data
# peptides_df <- read.csv("your_peptide_predictions.csv")
# experimental_df <- read.csv("your_iedb_data.csv")  # MUST have 'is_normal' column!
#
# # STEP 3: Preprocess (target_disease is OPTIONAL!)
# preprocessed_data <- preprocess_for_bayesian(
#   peptides_df = peptides_df,
#   experimental_df = experimental_df,
#   target_allele = "HLA-A*02:01",
#   target_disease = NULL,  # Can specify disease or leave NULL
#   target_class = "I"
# )
#
# # STEP 4: Run Bayesian assessment
# results <- bayesian_peptide_assessment(
#   peptide_data = preprocessed_data,
#   target_allele = "HLA-A*02:01",
#   target_class = "I",
#   prior_prob = 0.50
# )
#
# # STEP 5: Generate report
# report <- generate_assessment_report(
#   assessed_peptides = results,
#   output_prefix = "context_weighted_assessment",
#   include_plots = TRUE
# )
#
# # STEP 6: View top peptides by context
# top_by_context <- results %>%
#   filter(has_experimental_data == TRUE) %>%
#   arrange(desc(posterior_prob)) %>%
#   select(peptide_id, posterior_prob, confidence_level,
#          context_relevance_category, has_target_allele,
#          n_high_relevance, evidence_chain) %>%
#   head(20)
#
# print(top_by_context)
#
# # STEP 7 (OPTIONAL): Combine with similarity for smooth ranking
# results$rank_score <- results$posterior_prob * results$similarity_score  # If you have similarity_score
# results_ranked <- results %>% arrange(desc(rank_score))
