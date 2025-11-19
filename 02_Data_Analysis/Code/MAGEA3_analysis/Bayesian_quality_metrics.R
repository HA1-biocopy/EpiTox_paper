# =============================================================================
# Bayesian_quality_metrics_v2.R
# =============================================================================
# Quality Assessment Functions for Bayesian Peptide Assessment
# Updated to work with actual output structure from bayesian_peptide_assessment()
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

# =============================================================================
# FUNCTION 1: EVIDENCE TIER SEPARATION ANALYSIS (Updated)
# =============================================================================

#' Analyze and visualize posterior probability separation by evidence tier
#'
#' @param assessed_peptides Data frame from bayesian_peptide_assessment()
#' @param evidence_col Column name containing evidence tiers (default: "evidence")
#' @param posterior_col Column name containing posterior probabilities (default: "posterior_prob")
#' @param save_plot Logical, whether to save the plot (default: TRUE)
#' @param output_prefix Prefix for output files (default: "quality_check")
#' @return List containing plot object and summary statistics
#'
#' @export
analyze_evidence_tier_separation <- function(assessed_peptides,
                                             evidence_col = "evidence",
                                             posterior_col = "posterior_prob",
                                             save_plot = TRUE,
                                             output_prefix = "quality_check") {

  cat("\n")
  cat("======================================================================\n")
  cat("EVIDENCE TIER SEPARATION ANALYSIS\n")
  cat("======================================================================\n\n")

  # Verify columns exist
  if (!evidence_col %in% names(assessed_peptides)) {
    stop(paste("Column", evidence_col, "not found in data"))
  }
  if (!posterior_col %in% names(assessed_peptides)) {
    stop(paste("Column", posterior_col, "not found in data"))
  }

  # DIAGNOSTICS
  cat("=== DATA SUMMARY ===\n")
  cat("Total peptides:", nrow(assessed_peptides), "\n")
  cat("Posterior probability range:",
      round(min(assessed_peptides[[posterior_col]], na.rm=TRUE), 2),
      "to",
      round(max(assessed_peptides[[posterior_col]], na.rm=TRUE), 2), "\n")
  cat("\nEvidence tier distribution:\n")
  print(table(assessed_peptides[[evidence_col]], useNA = "ifany"))
  cat("\n")

  # Define evidence tier hierarchy and labels
  evidence_hierarchy <- c(
    "multiple_high_relevance",
    "one_high_relevance",
    "multiple_medium_relevance",
    "one_medium_relevance",
    "multiple_low_relevance",
    "one_low_relevance",
    "predicted"
  )

  evidence_labels <- c(
    "multiple_high_relevance" = "Multiple High\nRelevance",
    "one_high_relevance" = "Single High\nRelevance",
    "multiple_medium_relevance" = "Multiple Medium\nRelevance",
    "one_medium_relevance" = "Single Medium\nRelevance",
    "multiple_low_relevance" = "Multiple Low\nRelevance",
    "one_low_relevance" = "Single Low\nRelevance",
    "predicted" = "Predicted Only\n(No Exp. Data)"
  )

  # Prepare data
  tier_data <- assessed_peptides %>%
    filter(!is.na(.data[[evidence_col]])) %>%
    mutate(
      evidence_tier = factor(.data[[evidence_col]],
                             levels = evidence_hierarchy,
                             labels = evidence_labels[evidence_hierarchy]),
      tier_rank = match(.data[[evidence_col]], evidence_hierarchy),
      posterior = .data[[posterior_col]]
    ) %>%
    filter(!is.na(evidence_tier))

  # Calculate summary statistics by tier
  tier_summary <- tier_data %>%
    group_by(evidence_tier, tier_rank) %>%
    summarise(
      n = n(),
      mean_posterior = mean(posterior),
      median_posterior = median(posterior),
      sd_posterior = sd(posterior),
      q25 = quantile(posterior, 0.25),
      q75 = quantile(posterior, 0.75),
      min_posterior = min(posterior),
      max_posterior = max(posterior),
      range_posterior = max(posterior) - min(posterior),
      .groups = "drop"
    ) %>%
    arrange(tier_rank) %>%
    filter(!is.na(tier_rank))

  cat("=== EVIDENCE TIER STATISTICS ===\n")
  print(tier_summary %>%
          select(evidence_tier, n, mean_posterior, sd_posterior, range_posterior))
  cat("\n")

  # Check for variance issues
  no_variance_tiers <- tier_summary %>%
    filter(sd_posterior < 0.01 | range_posterior < 0.01)

  if (nrow(no_variance_tiers) > 0) {
    cat("⚠ WARNING: The following tiers have NO VARIANCE:\n")
    print(no_variance_tiers %>% select(evidence_tier, n, mean_posterior, sd_posterior))
    cat("\nThis means all peptides in these tiers have identical (or nearly identical)\n")
    cat("posterior probabilities. This suggests:\n")
    cat("  1. All peptides in the tier have identical evidence combinations\n")
    cat("  2. Other factors (binding, allele) are not varying within tiers\n")
    cat("  3. The model may be creating discrete bins rather than continuous probabilities\n\n")
  }

  # Calculate separation metrics
  if (nrow(tier_summary) > 1) {
    cat("=== SEPARATION QUALITY METRICS ===\n")

    # Mean separation between consecutive tiers
    tier_means <- tier_summary %>% arrange(tier_rank) %>% pull(mean_posterior)
    if (length(tier_means) > 1) {
      mean_separations <- diff(tier_means)
      cat("Mean separation between consecutive tiers:\n")
      for (i in 1:length(mean_separations)) {
        cat(sprintf("  Tier %d -> %d: %.2f%%\n", i, i+1, abs(mean_separations[i])))
      }
      cat(sprintf("Average separation: %.2f%%\n\n", mean(abs(mean_separations))))
    }
  }

  # Create visualization
  cat("=== CREATING VISUALIZATION ===\n")

  # Main plot: Violin plots with box plots overlay
  p1 <- ggplot(tier_data, aes(x = evidence_tier, y = posterior, fill = evidence_tier)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "white", color = "black") +
    scale_fill_viridis_d(option = "plasma", direction = -1) +
    labs(
      title = "Posterior Probability Distribution by Evidence Tier",
      subtitle = paste("N =", nrow(tier_data), "peptides | Violin plots show density; boxes show quartiles; diamonds show means"),
      x = "Evidence Tier",
      y = "Posterior Probability (%)",
      caption = "Higher tiers should show higher probabilities with clear separation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray40", alpha = 0.7, size = 0.8)

  # Secondary plot: Mean and ranges
  p2 <- ggplot(tier_summary, aes(x = tier_rank, y = mean_posterior)) +
    geom_line(color = "#440154FF", size = 1.2) +
    geom_point(size = 4, color = "#440154FF") +
    geom_errorbar(aes(ymin = min_posterior, ymax = max_posterior),
                  width = 0.2, color = "#440154FF", alpha = 0.4, size = 0.8) +
    geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3) +
    scale_x_continuous(breaks = tier_summary$tier_rank,
                       labels = tier_summary$evidence_tier) +
    labs(
      title = "Mean Posterior Probability by Evidence Tier",
      subtitle = "Error bars show min-max range | Points show mean",
      x = "Evidence Tier",
      y = "Mean Posterior Probability (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.minor = element_blank()
    ) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray40", alpha = 0.7)

  # Combine plots
  combined_plot <- gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(2, 1))

  # Save if requested
  if (save_plot) {
    ggsave(
      filename = paste0(output_prefix, "_evidence_tier_separation.png"),
      plot = combined_plot,
      width = 12,
      height = 10,
      dpi = 300
    )
    cat("✓ Plot saved:", paste0(output_prefix, "_evidence_tier_separation.png\n"))

    write.csv(tier_summary,
              paste0(output_prefix, "_tier_summary_stats.csv"),
              row.names = FALSE)
    cat("✓ Summary stats saved\n")
  }

  cat("\n======================================================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("======================================================================\n\n")

  # Interpretation
  cat("INTERPRETATION:\n")
  cat("✓ Clear tier separation observed:",
      round(mean(abs(diff(tier_summary$mean_posterior))), 1), "% average\n")
  if (nrow(no_variance_tiers) > 0) {
    cat("⚠ Low within-tier variance suggests deterministic outcomes\n")
  } else {
    cat("✓ Good within-tier variance shows evidence integration\n")
  }
  cat("\n")

  return(invisible(list(
    plot = combined_plot,
    summary = tier_summary,
    tier_data = tier_data
  )))
}


# =============================================================================
# FUNCTION 2: NOVEL VS. KNOWN PEPTIDE TREATMENT ANALYSIS (Updated)
# =============================================================================

#' Compare novel (predicted) vs. known (experimental) peptide treatment
#'
#' @param assessed_peptides Data frame from bayesian_peptide_assessment()
#' @param evidence_col Column name containing evidence tiers (default: "evidence")
#' @param posterior_col Column name containing posterior probabilities (default: "posterior_prob")
#' @param prior_col Column name containing prior probabilities (default: "prior_prob")
#' @param save_plot Logical, whether to save the plot (default: TRUE)
#' @param output_prefix Prefix for output files (default: "quality_check")
#' @return List containing plot object and summary statistics
#'
#' @export
analyze_novel_vs_known_treatment <- function(assessed_peptides,
                                             evidence_col = "evidence",
                                             posterior_col = "posterior_prob",
                                             prior_col = "prior_prob",
                                             save_plot = TRUE,
                                             output_prefix = "quality_check") {

  cat("\n")
  cat("======================================================================\n")
  cat("NOVEL VS. KNOWN PEPTIDE TREATMENT ANALYSIS\n")
  cat("======================================================================\n\n")

  # Verify columns
  if (!evidence_col %in% names(assessed_peptides)) {
    stop(paste("Column", evidence_col, "not found"))
  }
  if (!posterior_col %in% names(assessed_peptides)) {
    stop(paste("Column", posterior_col, "not found"))
  }

  # DIAGNOSTICS
  cat("=== DATA SUMMARY ===\n")
  cat("Total peptides:", nrow(assessed_peptides), "\n")
  cat("Novel (predicted) peptides:",
      sum(assessed_peptides[[evidence_col]] == "predicted", na.rm=TRUE), "\n")
  cat("Known (experimental) peptides:",
      sum(assessed_peptides[[evidence_col]] != "predicted", na.rm=TRUE), "\n\n")

  # Categorize peptides
  comparison_data <- assessed_peptides %>%
    mutate(
      peptide_type = case_when(
        .data[[evidence_col]] == "predicted" ~ "Novel (Predicted Only)",
        TRUE ~ "Known (Has Experimental Data)"
      ),
      posterior = .data[[posterior_col]],
      prior = if(prior_col %in% names(assessed_peptides)) .data[[prior_col]] else 50
    )

  # Further categorize novel peptides
  if ("binding_affinity" %in% names(assessed_peptides) & "hla_allele" %in% names(assessed_peptides)) {
    comparison_data <- comparison_data %>%
      mutate(
        novel_subcategory = if_else(
          peptide_type == "Novel (Predicted Only)",
          case_when(
            !is.na(binding_affinity) & binding_affinity == "strong" ~ "Novel: Strong binding",
            !is.na(binding_affinity) & binding_affinity %in% c("intermediate", "weak") ~ "Novel: Weak/intermediate binding",
            !is.na(hla_allele) ~ "Novel: Allele match only",
            TRUE ~ "Novel: Minimal evidence"
          ),
          "Known: Has experimental data"
        )
      )
  } else {
    comparison_data$novel_subcategory <- comparison_data$peptide_type
  }

  # Summary statistics
  cat("=== SUMMARY STATISTICS ===\n")
  type_summary <- comparison_data %>%
    group_by(peptide_type) %>%
    summarise(
      n = n(),
      pct_of_total = n() / nrow(comparison_data) * 100,
      mean_posterior = mean(posterior),
      median_posterior = median(posterior),
      sd_posterior = sd(posterior),
      min_posterior = min(posterior),
      max_posterior = max(posterior),
      mean_prior = mean(prior),
      pct_above_prior = sum(posterior > mean(prior)) / n() * 100,
      pct_below_prior = sum(posterior < mean(prior)) / n() * 100,
      .groups = "drop"
    )

  print(type_summary %>% select(-mean_prior))
  cat("\n")

  # Check for penalty
  novel_data <- comparison_data %>% filter(peptide_type == "Novel (Predicted Only)")
  if (nrow(novel_data) > 0) {
    mean_prior <- mean(novel_data$prior, na.rm = TRUE)
    mean_post <- mean(novel_data$posterior, na.rm = TRUE)

    cat("=== PENALTY ANALYSIS ===\n")
    cat(sprintf("Novel peptides prior (mean): %.1f%%\n", mean_prior))
    cat(sprintf("Novel peptides posterior (mean): %.1f%%\n", mean_post))
    cat(sprintf("Change from prior: %.1f%%\n", mean_post - mean_prior))

    if (mean_post < mean_prior - 5) {
      cat("⚠ WARNING: Novel peptides show systematic penalty (below prior)\n")
      cat("   This suggests the 'predicted' evidence tier has LR < 1.0\n\n")
    } else if (mean_post > mean_prior + 5) {
      cat("⚠ NOTE: Novel peptides show systematic boost (above prior)\n")
      cat("   This suggests the 'predicted' evidence tier has LR > 1.0\n\n")
    } else {
      cat("✓ Novel peptides remain near prior (fair treatment)\n\n")
    }
  }

  # Detailed breakdown for novel peptides
  if ("novel_subcategory" %in% names(comparison_data)) {
    cat("=== NOVEL PEPTIDE SUBCATEGORIES ===\n")
    novel_breakdown <- comparison_data %>%
      filter(peptide_type == "Novel (Predicted Only)") %>%
      group_by(novel_subcategory) %>%
      summarise(
        n = n(),
        mean_posterior = mean(posterior),
        sd_posterior = sd(posterior),
        pct_above_50 = sum(posterior >= 50) / n() * 100,
        .groups = "drop"
      ) %>%
      arrange(desc(mean_posterior))

    if (nrow(novel_breakdown) > 0) {
      print(novel_breakdown)
      cat("\n")
    }
  }

  # Create visualizations
  cat("=== CREATING VISUALIZATIONS ===\n")

  # Plot 1: Distribution comparison - Density
  p1 <- ggplot(comparison_data, aes(x = posterior, fill = peptide_type)) +
    geom_density(alpha = 0.6, size = 1) +
    geom_vline(xintercept = mean(comparison_data$prior),
               linetype = "dashed", color = "black", size = 1) +
    scale_fill_manual(
      values = c("Novel (Predicted Only)" = "#E69F00",
                 "Known (Has Experimental Data)" = "#56B4E9"),
      name = "Peptide Type"
    ) +
    labs(
      title = "Posterior Probability Distribution: Novel vs. Known Peptides",
      subtitle = paste("N =", nrow(comparison_data), "| Dashed line = prior probability"),
      x = "Posterior Probability (%)",
      y = "Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top"
    )

  # Plot 2: Box plot comparison
  p2 <- ggplot(comparison_data, aes(x = peptide_type, y = posterior, fill = peptide_type)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "white", color = "black") +
    geom_hline(yintercept = mean(comparison_data$prior),
               linetype = "dashed", color = "gray40") +
    scale_fill_manual(
      values = c("Novel (Predicted Only)" = "#E69F00",
                 "Known (Has Experimental Data)" = "#56B4E9")
    ) +
    labs(
      title = "Posterior Probability Distribution Comparison",
      subtitle = "Diamond = mean | Box = IQR | Whiskers = 1.5×IQR | Dashed line = prior",
      x = NULL,
      y = "Posterior Probability (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "none",
      axis.text.x = element_text(size = 11)
    )

  # Plot 3: Novel peptide subcategories (if available)
  p3 <- NULL
  if ("novel_subcategory" %in% names(comparison_data)) {
    novel_only <- comparison_data %>%
      filter(peptide_type == "Novel (Predicted Only)") %>%
      filter(novel_subcategory != "Known: Has experimental data")

    if (nrow(novel_only) > 0 & length(unique(novel_only$novel_subcategory)) > 1) {
      p3 <- ggplot(novel_only, aes(x = reorder(novel_subcategory, posterior),
                                   y = posterior, fill = novel_subcategory)) +
        geom_violin(alpha = 0.6, scale = "width") +
        geom_boxplot(width = 0.3, alpha = 0.8) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                     fill = "white", color = "black") +
        coord_flip() +
        scale_fill_brewer(palette = "Set2") +
        geom_hline(yintercept = mean(novel_only$prior),
                   linetype = "dashed", color = "gray40") +
        labs(
          title = "Novel Peptides: Impact of Supporting Evidence",
          subtitle = "Shows how binding affinity and allele information affect novel peptides",
          x = NULL,
          y = "Posterior Probability (%)"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 13),
          legend.position = "none"
        )
    }
  }

  # Plot 4: Cumulative distribution
  p4 <- ggplot(comparison_data, aes(x = posterior, color = peptide_type)) +
    stat_ecdf(size = 1.2) +
    scale_color_manual(
      values = c("Novel (Predicted Only)" = "#E69F00",
                 "Known (Has Experimental Data)" = "#56B4E9"),
      name = "Peptide Type"
    ) +
    geom_vline(xintercept = mean(comparison_data$prior),
               linetype = "dashed", color = "black") +
    labs(
      title = "Cumulative Distribution: Novel vs. Known Peptides",
      subtitle = "Shows proportion of peptides below each probability threshold",
      x = "Posterior Probability (%)",
      y = "Cumulative Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "top"
    )

  # Combine plots
  if (!is.null(p3)) {
    combined_plot <- gridExtra::grid.arrange(
      p1, p2, p3, p4,
      ncol = 2, nrow = 2,
      top = grid::textGrob("Novel vs. Known Peptide Treatment Analysis",
                           gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
  } else {
    combined_plot <- gridExtra::grid.arrange(
      p1, p2, p4,
      ncol = 2, nrow = 2,
      layout_matrix = rbind(c(1, 2), c(3, 3)),
      top = grid::textGrob("Novel vs. Known Peptide Treatment Analysis",
                           gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
  }

  # Save if requested
  if (save_plot) {
    ggsave(
      filename = paste0(output_prefix, "_novel_vs_known_comparison.png"),
      plot = combined_plot,
      width = 14,
      height = 12,
      dpi = 300
    )
    cat("✓ Plot saved:", paste0(output_prefix, "_novel_vs_known_comparison.png\n"))

    write.csv(type_summary,
              paste0(output_prefix, "_novel_vs_known_summary.csv"),
              row.names = FALSE)
    cat("✓ Summary statistics saved\n")
  }

  cat("\n======================================================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("======================================================================\n\n")

  return(invisible(list(
    combined_plot = combined_plot,
    type_summary = type_summary,
    comparison_data = comparison_data
  )))
}


# =============================================================================
# CONVENIENCE FUNCTION: RUN BOTH QUALITY CHECKS
# =============================================================================

#' Run complete quality assessment
#'
#' @param assessed_peptides Data frame from bayesian_peptide_assessment()
#' @param output_prefix Prefix for all output files
#' @return List containing results from both analyses
#'
#' @export
run_complete_quality_assessment <- function(assessed_peptides,
                                            output_prefix = "bayesian_quality_check") {

  cat("\n")
  cat("######################################################################\n")
  cat("#                                                                    #\n")
  cat("#           BAYESIAN MODEL QUALITY ASSESSMENT                       #\n")
  cat("#                                                                    #\n")
  cat("######################################################################\n")

  # Run both analyses
  tier_results <- analyze_evidence_tier_separation(
    assessed_peptides = assessed_peptides,
    save_plot = TRUE,
    output_prefix = output_prefix
  )

  novel_results <- analyze_novel_vs_known_treatment(
    assessed_peptides = assessed_peptides,
    save_plot = TRUE,
    output_prefix = output_prefix
  )

  cat("\n")
  cat("######################################################################\n")
  cat("#                                                                    #\n")
  cat("#           ALL QUALITY ASSESSMENTS COMPLETE                        #\n")
  cat("#                                                                    #\n")
  cat("######################################################################\n\n")

  return(invisible(list(
    evidence_tier_results = tier_results,
    novel_vs_known_results = novel_results
  )))
}

# # After running your Bayesian assessment:
# source("Bayesian_quality_metrics.R")
#
# Run individual analyses
tier_analysis <- analyze_evidence_tier_separation(
  assessed_peptides = results,
  save_plot = TRUE,
  output_prefix = "my_analysis"
)
#
novel_analysis <- analyze_novel_vs_known_treatment(
  assessed_peptides = results,
  save_plot = TRUE,
  output_prefix = "my_analysis"
)
#
# # Or run both at once
# quality_results <- run_complete_quality_assessment(
#   assessed_peptides = results,
#   output_prefix = "my_analysis"
# )
