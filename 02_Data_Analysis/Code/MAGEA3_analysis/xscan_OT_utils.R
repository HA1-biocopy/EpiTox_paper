# ============================================================================
# FUNCTION 1: Extract Dissimilarity Positions
# ============================================================================

extract_dissimilarity_positions <- function(target_seq, off_target_seq, target_classifications) {
  #' Extract positions where off-target differs from target
  #'
  #' @param target_seq String: target peptide sequence (e.g., "ALYNKVFLK")
  #' @param off_target_seq String: off-target peptide sequence
  #' @param target_classifications Dataframe with columns: Position, WT_AA, Category
  #' @return Dataframe with dissimilarity analysis

  # Split sequences into characters
  target_chars <- strsplit(target_seq, "")[[1]]
  off_target_chars <- strsplit(off_target_seq, "")[[1]]

  # Check if lengths match
  if (length(target_chars) != length(off_target_chars)) {
    stop("Target and off-target sequences must have the same length")
  }

  # Identify mismatches
  mismatches <- target_chars != off_target_chars

  # Create dissimilarity dataframe
  dissim_df <- data.frame(
    Position = 1:length(target_chars),
    Target_AA = target_chars,
    OffTarget_AA = off_target_chars,
    Is_Mismatch = mismatches
  )

  # Add category information from target classifications
  dissim_df <- dissim_df %>%
    left_join(target_classifications %>% select(Position, Category, Type), by = "Position")

  return(dissim_df)
}

# ============================================================================
# FUNCTION 2: Summarize Dissimilarities by Position Type
# ============================================================================

summarize_dissimilarity_by_type <- function(dissim_df) {
  #' Summarize mismatches by position category
  #'
  #' @param dissim_df Output from extract_dissimilarity_positions()
  #' @return Dataframe with counts and percentages by category

  summary <- dissim_df %>%
    group_by(Type) %>%
    summarise(
      Total_Positions = n(),
      Mismatched_Positions = sum(Is_Mismatch, na.rm = TRUE),
      Matched_Positions = sum(!Is_Mismatch, na.rm = TRUE),
      Percent_Mismatched = (Mismatched_Positions / Total_Positions) * 100,
      .groups = "drop"
    )

  return(summary)
}

# ============================================================================
# FUNCTION 3: Classify Position Types
# ============================================================================

classify_position_types <- function(target_classifications) {
  #' Classify positions into Anchor, Exposed, and Both categories
  #'
  #' @param target_classifications Dataframe with columns: Position, WT_AA, Category
  #' @return Dataframe with additional Type column

  classified <- target_classifications %>%
    mutate(
      Type = case_when(
        Category == "Classical Anchor" ~ "Anchor",
        Category == "Flexible Anchor" ~ "Both",  # High sensitivity but also permissive
        Category == "Permissive" ~ "Exposed",
        Category == "Neutral" ~ "Exposed",
        TRUE ~ "Unknown"
      )
    )

  return(classified)
}

# ============================================================================
# FUNCTION 4: Analyze Multiple Off-Targets
# ============================================================================

analyze_offtargets <- function(target_seq, target_classifications, offtargets_df) {
  #' Analyze dissimilarity patterns for multiple off-target peptides
  #'
  #' @param target_seq String: target peptide sequence
  #' @param target_classifications Dataframe with Position, WT_AA, Category
  #' @param offtargets_df Dataframe with columns: peptide, affinity
  #' @return List containing detailed and summary results

  # Add position type classification
  target_with_types <- classify_position_types(target_classifications)

  # Initialize results list
  results_list <- list()
  i = 1
  # Analyze each off-target
  for (i in 1:nrow(offtargets_df)) {
    off_target_seq <- offtargets_df$peptide[i]
    off_target_affinity <- offtargets_df$affinity[i]
    off_target_id = paste0(off_target_seq, "_", offtargets_df$uniprot[i])

    # Extract dissimilarities
    dissim <- extract_dissimilarity_positions(target_seq, off_target_seq, target_with_types)

    # Summarize by category
    summary_category <- summarize_dissimilarity_by_type(dissim)

    # Summarize by type (Anchor/Exposed/Both)
    summary_type <- dissim %>%
      group_by(Type) %>%
      summarise(
        Total_Positions = n(),
        Mismatched_Positions = sum(Is_Mismatch, na.rm = TRUE),
        Percent_Mismatched = (Mismatched_Positions / Total_Positions) * 100,
        .groups = "drop"
      )

    # Overall statistics
    total_mismatches <- sum(dissim$Is_Mismatch)
    total_positions <- nrow(dissim)
    overall_identity <- ((total_positions - total_mismatches) / total_positions) * 100

    # Store results
    results_list[[i]] <- list(
      id = off_target_id,
      peptide = off_target_seq,
      affinity = off_target_affinity,
      overall_identity = overall_identity,
      total_mismatches = total_mismatches,
      dissim_details = dissim,
      summary_by_category = summary_category,
      summary_by_type = summary_type
    )
  }

  return(results_list)
}

# ============================================================================
# FUNCTION 5: Create Summary Table for All Off-Targets
# ============================================================================

create_summary_table <- function(analysis_results) {
  #' Create a comprehensive summary table for all off-targets
  #'
  #' @param analysis_results Output from analyze_offtargets()
  #' @return Dataframe with one row per off-target

  summary_df <- data.frame()

  for (result in analysis_results) {
    # Extract anchor mismatches
    anchor_summary <- result$summary_by_type %>%
      filter(Type == "Anchor")

    if (nrow(anchor_summary) > 0) {
      anchor_mismatches <- anchor_summary$Mismatched_Positions
      anchor_percent <- anchor_summary$Percent_Mismatched
    } else {
      anchor_mismatches <- 0
      anchor_percent <- 0
    }

    # Extract exposed mismatches
    exposed_summary <- result$summary_by_type %>%
      filter(Type == "Exposed")

    if (nrow(exposed_summary) > 0) {
      exposed_mismatches <- exposed_summary$Mismatched_Positions
      exposed_percent <- exposed_summary$Percent_Mismatched
    } else {
      exposed_mismatches <- 0
      exposed_percent <- 0
    }

    # Extract both (flexible anchor) mismatches
    both_summary <- result$summary_by_type %>%
      filter(Type == "Both")

    if (nrow(both_summary) > 0) {
      both_mismatches <- both_summary$Mismatched_Positions
      both_percent <- both_summary$Percent_Mismatched
    } else {
      both_mismatches <- 0
      both_percent <- 0
    }

    #print(result)
    # Create row
    row <- data.frame(
      id = result$id,
      Overall_Identity = round(result$overall_identity, 2),
      Total_Mismatches = result$total_mismatches,
      Anchor_Mismatches = anchor_mismatches,
      Anchor_Percent_Diff = round(anchor_percent, 2),
      Exposed_Mismatches = exposed_mismatches,
      Exposed_Percent_Diff = round(exposed_percent, 2),
      Both_Mismatches = both_mismatches,
      Both_Percent_Diff = round(both_percent, 2)
    )

    summary_df <- rbind(summary_df, row)
  }

  return(summary_df)
}

# ============================================================================
# DEBUGGING FUNCTION - Check intermediate results
# ============================================================================

debug_summary_extraction <- function(analysis_results, idx = 1) {
  #' Debug helper to see what's in the results
  #'
  #' @param analysis_results Output from analyze_offtargets()
  #' @param idx Index of off-target to inspect

  cat("=== DEBUGGING RESULT", idx, "===\n\n")

  result <- analysis_results[[idx]]

  cat("Peptide:", result$peptide, "\n")
  cat("Total mismatches:", result$total_mismatches, "\n\n")

  cat("Summary by Type:\n")
  print(result$summary_by_type)

  cat("\n\nChecking filter for 'Anchor':\n")
  anchor_summary <- result$summary_by_type %>% filter(Type == "Anchor")
  cat("Rows found:", nrow(anchor_summary), "\n")
  if (nrow(anchor_summary) > 0) {
    print(anchor_summary)
  } else {
    cat("No Anchor positions found!\n")
  }

  cat("\n\nFull dissim_details:\n")
  print(result$dissim_details)

  cat("\n\nUnique Types in data:\n")
  print(unique(result$dissim_details$Type))
}



# ============================================================================
# VISUALIZATION FUNCTION
# ============================================================================

plot_mismatch_by_type <- function(summary_table) {
  #' Create grouped bar plot showing mismatches by position type
  #'
  #' @param summary_table Output from create_summary_table()
  #' @return ggplot object

  library(ggplot2)
  library(tidyr)

  # Reshape data for plotting
  plot_data <- summary_table %>%
    select(peptide, Anchor_Mismatches, Exposed_Mismatches, Both_Mismatches) %>%
    pivot_longer(cols = c(Anchor_Mismatches, Exposed_Mismatches, Both_Mismatches),
                 names_to = "Mismatch_Type",
                 values_to = "Count") %>%
    mutate(Mismatch_Type = gsub("_Mismatches", "", Mismatch_Type))

  # Create plot
  p <- ggplot(plot_data, aes(x = peptide, y = Count, fill = Mismatch_Type)) +
    geom_bar(stat = "identity", position = "dodge", color = "#2C4255", size = 0.5) +
    scale_fill_manual(values = c("Anchor" = "#A2C510",
                                 "Exposed" = "steelblue",
                                 "Both" = "purple"),
                      name = "Position Type") +
    labs(title = "Mismatches by Position Type",
         subtitle = "Comparison of off-target peptides to target",
         x = "Off-Target Peptide",
         y = "Number of Mismatches") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(face = "bold", size = 14, color = "#2C4255"),
          plot.subtitle = element_text(size = 11, color = "#666"),
          legend.position = "bottom",
          panel.grid.major.x = element_blank())

  return(p)
}

# ============================================================================
# FUNCTION 7: Alternative Stacked Bar Plot
# ============================================================================

plot_mismatch_stacked <- function(summary_table) {
  #' Create stacked bar plot showing mismatches by position type
  #'
  #' @param summary_table Output from create_summary_table()
  #' @return ggplot object

  library(ggplot2)
  library(tidyr)

  # Reshape data for plotting
  plot_data <- summary_table %>%
    select(peptide, Anchor_Mismatches, Exposed_Mismatches, Both_Mismatches) %>%
    pivot_longer(cols = c(Anchor_Mismatches, Exposed_Mismatches, Both_Mismatches),
                 names_to = "Mismatch_Type",
                 values_to = "Count") %>%
    mutate(Mismatch_Type = gsub("_Mismatches", "", Mismatch_Type),
           Mismatch_Type = factor(Mismatch_Type, levels = c("Anchor", "Both", "Exposed")))

  # Create plot
  p <- ggplot(plot_data, aes(x = peptide, y = Count, fill = Mismatch_Type)) +
    geom_bar(stat = "identity", position = "stack", color = "#2C4255", size = 0.5) +
    scale_fill_manual(values = c("Anchor" = "#A2C510",
                                 "Both" = "purple",
                                 "Exposed" = "steelblue"),
                      name = "Position Type") +
    labs(title = "Total Mismatches by Position Type",
         subtitle = "Stacked view of off-target differences from target",
         x = "Off-Target Peptide",
         y = "Number of Mismatches") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(face = "bold", size = 14, color = "#2C4255"),
          plot.subtitle = element_text(size = 11, color = "#666"),
          legend.position = "bottom",
          panel.grid.major.x = element_blank())

  return(p)
}

# ============================================================================
# FUNCTION 8: Heatmap View of Mismatches
# ============================================================================

plot_mismatch_heatmap <- function(summary_table) {
  #' Create heatmap showing mismatch patterns
  #'
  #' @param summary_table Output from create_summary_table()
  #' @return ggplot object

  library(ggplot2)
  library(tidyr)

  # Reshape data for plotting
  plot_data <- summary_table %>%
    select(peptide, Anchor_Mismatches, Exposed_Mismatches, Both_Mismatches) %>%
    pivot_longer(cols = c(Anchor_Mismatches, Exposed_Mismatches, Both_Mismatches),
                 names_to = "Mismatch_Type",
                 values_to = "Count") %>%
    mutate(Mismatch_Type = gsub("_Mismatches", "", Mismatch_Type),
           Mismatch_Type = factor(Mismatch_Type, levels = c("Anchor", "Both", "Exposed")))

  # Create plot
  p <- ggplot(plot_data, aes(x = peptide, y = Mismatch_Type, fill = Count)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = Count), color = "white", fontface = "bold", size = 5) +
    scale_fill_gradient2(low = "#F0F0F0", mid = "steelblue", high = "#A2C510",
                         midpoint = max(plot_data$Count, na.rm = TRUE) / 2,
                         name = "Mismatches") +
    labs(title = "Mismatch Pattern Heatmap",
         subtitle = "Number of mismatches by position type and peptide",
         x = "Off-Target Peptide",
         y = "Position Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(face = "bold", size = 14, color = "#2C4255"),
          plot.subtitle = element_text(size = 11, color = "#666"),
          legend.position = "right",
          panel.grid = element_blank())

  return(p)
}

# ============================================================================
# FUNCTION 9: Overall Dataset Summary Plots
# ============================================================================

plot_dataset_summary <- function(summary_table) {
  #' Create summary plots for the entire dataset
  #'
  #' @param summary_table Output from create_summary_table()
  #' @return List of ggplot objects

  library(ggplot2)
  library(gridExtra)

  # Calculate categories based on mismatch patterns
  categorized <- summary_table %>%
    mutate(
      Category = case_when(
        Anchor_Mismatches > 0 & Exposed_Mismatches > 0 ~ "Both Anchor & Exposed",
        Anchor_Mismatches > 0 & Exposed_Mismatches == 0 ~ "Anchor Only",
        Anchor_Mismatches == 0 & Exposed_Mismatches > 0 ~ "Exposed Only",
        TRUE ~ "Perfect Match"
      ),
      Risk_Level = case_when(
        Anchor_Mismatches > 0 ~ "High Risk",
        Exposed_Mismatches >= 2 ~ "Medium Risk",
        Exposed_Mismatches == 1 ~ "Low Risk",
        TRUE ~ "No Risk"
      )
    )

  # Plot 1: Distribution by mismatch category
  category_counts <- categorized %>%
    count(Category) %>%
    mutate(Percentage = (n / sum(n)) * 100)

  p1 <- ggplot(category_counts, aes(x = "", y = n, fill = Category)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 1.5) +
    coord_polar("y") +
    geom_text(aes(label = paste0(n, "\n(", round(Percentage, 1), "%)")),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 4) +
    scale_fill_manual(values = c("Both Anchor & Exposed" = "#abd9e9",
                                 "Anchor Only" = "#A2C510",
                                 "Exposed Only" = "#FBB800",
                                 "Perfect Match" = "#F0F0F0"),
                      name = "Mismatch Pattern") +
    labs(title = "Distribution of Off-Targets by Mismatch Pattern",
         subtitle = paste("Total peptides:", nrow(summary_table))) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 14, color = "#2C4255", hjust = 0.5),
          plot.subtitle = element_text(size = 11, color = "#666", hjust = 0.5),
          legend.position = "right")

  # Plot 2: Risk level distribution (bar chart)
  risk_counts <- categorized %>%
    count(Risk_Level) %>%
    mutate(Risk_Level = factor(Risk_Level,
                               levels = c("No Risk", "Low Risk", "Medium Risk", "High Risk")),
           Percentage = (n / sum(n)) * 100)

  p2 <- ggplot(risk_counts, aes(x = Risk_Level, y = n, fill = Risk_Level)) +
    geom_bar(stat = "identity", color = "#2C4255", size = 0.8) +
    geom_text(aes(label = paste0(n, " (", round(Percentage, 1), "%)")),
              vjust = -0.5, fontface = "bold", size = 4) +
    scale_fill_manual(values = c("No Risk" = "#F0F0F0",
                                 "Low Risk" = "#90EE90",
                                 "Medium Risk" = "#FFD700",
                                 "High Risk" = "#ff6b6b"),
                      name = "Risk Level") +
    labs(title = "Risk Level Distribution",
         subtitle = "Based on anchor position mismatches",
         x = "Risk Level",
         y = "Number of Peptides") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14, color = "#2C4255"),
          plot.subtitle = element_text(size = 11, color = "#666"),
          legend.position = "none",
          panel.grid.major.x = element_blank())

  # Plot 3: Distribution of total mismatches (histogram)
  p3 <- ggplot(categorized, aes(x = Total_Mismatches)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "#2C4255", alpha = 0.8) +
    geom_vline(aes(xintercept = mean(Total_Mismatches, na.rm = TRUE)),
               color = "#A2C510", linetype = "dashed", size = 1.2) +
    annotate("text", x = mean(categorized$Total_Mismatches, na.rm = TRUE),
             y = Inf, label = paste("Mean:", round(mean(categorized$Total_Mismatches, na.rm = TRUE), 1)),
             vjust = 2, color = "#A2C510", fontface = "bold") +
    labs(title = "Distribution of Total Mismatches",
         subtitle = "Frequency of mismatch counts across all off-targets",
         x = "Number of Mismatches",
         y = "Number of Peptides") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14, color = "#2C4255"),
          plot.subtitle = element_text(size = 11, color = "#666"))

  # Plot 4: Average mismatches by position type (bar chart)
  avg_mismatches <- data.frame(
    Type = c("Anchor", "Exposed", "Both"),
    Average = c(mean(summary_table$Anchor_Mismatches, na.rm = TRUE),
                mean(summary_table$Exposed_Mismatches, na.rm = TRUE),
                mean(summary_table$Both_Mismatches, na.rm = TRUE)),
    Total = c(sum(summary_table$Anchor_Mismatches, na.rm = TRUE),
              sum(summary_table$Exposed_Mismatches, na.rm = TRUE),
              sum(summary_table$Both_Mismatches, na.rm = TRUE))
  )

  p4 <- ggplot(avg_mismatches, aes(x = Type, y = Average, fill = Type)) +
    geom_bar(stat = "identity", color = "#2C4255", size = 0.8) +
    geom_text(aes(label = paste0(round(Average, 2), "\n(Total: ", Total, ")")),
              vjust = -0.5, fontface = "bold", size = 4) +
    scale_fill_manual(values = c("Anchor" = "#A2C510",
                                 "Exposed" = "steelblue",
                                 "Both" = "purple"),
                      name = "Position Type") +
    labs(title = "Average Mismatches by Position Type",
         subtitle = "Mean mismatches per peptide with total counts",
         x = "Position Type",
         y = "Average Number of Mismatches") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14, color = "#2C4255"),
          plot.subtitle = element_text(size = 11, color = "#666"),
          legend.position = "none",
          panel.grid.major.x = element_blank())

  # Return list of plots and categorized data
  return(list(
    pie_chart = p1,
    risk_distribution = p2,
    mismatch_histogram = p3,
    avg_by_type = p4,
    categorized_data = categorized
  ))
}

# ============================================================================
# FUNCTION 10: Combined Summary Dashboard
# ============================================================================

plot_summary_dashboard <- function(summary_table) {
  #' Create a 4-panel summary dashboard
  #'
  #' @param summary_table Output from create_summary_table()
  #' @return Combined grid plot

  library(gridExtra)

  # Get all plots
  summary_plots <- plot_dataset_summary(summary_table)

  # Arrange in 2x2 grid
  dashboard <- grid.arrange(
    summary_plots$pie_chart,
    summary_plots$risk_distribution,
    summary_plots$mismatch_histogram,
    summary_plots$avg_by_type,
    ncol = 2,
    nrow = 2,
    top = grid::textGrob("Off-Target Peptide Analysis Dashboard",
                         gp = grid::gpar(fontsize = 16, fontface = "bold", col = "#2C4255"))
  )

  return(dashboard)
}

# ============================================================================
# DEBUGGING FUNCTION - Check intermediate results
# ============================================================================

debug_summary_extraction <- function(analysis_results, idx = 1) {
  #' Debug helper to see what's in the results
  #'
  #' @param analysis_results Output from analyze_offtargets()
  #' @param idx Index of off-target to inspect

  cat("=== DEBUGGING RESULT", idx, "===\n\n")

  result <- analysis_results[[idx]]

  cat("Peptide:", result$peptide, "\n")
  cat("Total mismatches:", result$total_mismatches, "\n\n")

  cat("Summary by Type:\n")
  print(result$summary_by_type)

  cat("\n\nChecking filter for 'Anchor':\n")
  anchor_summary <- result$summary_by_type %>% filter(Type == "Anchor")
  cat("Rows found:", nrow(anchor_summary), "\n")
  if (nrow(anchor_summary) > 0) {
    print(anchor_summary)
  } else {
    cat("No Anchor positions found!\n")
  }

  cat("\n\nFull dissim_details:\n")
  print(result$dissim_details)

  cat("\n\nUnique Types in data:\n")
  print(unique(result$dissim_details$Type))
}
