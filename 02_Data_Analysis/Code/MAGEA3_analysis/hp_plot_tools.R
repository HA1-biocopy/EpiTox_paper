# Peptide Identification Tool Comparison: Your Tool vs CrossDome
# Install required packages if needed
# install.packages(c("ggplot2", "reshape2", "dplyr", "RColorBrewer"))

library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

# Define the comparison data
# 1 = Fully supported, 0.5 = Partially supported, 0 = Not supported
tool_comparison <- data.frame(
  Tool = c("EpiTox", "CrossDome"),
  BLOSUM_Similarity = c(1, 0),  # Update based on your tool's capabilities
  SNPs_Integration = c(1, 0),   # Update these values
  Gene_Expression = c(1, 1),
  HLA_Binding_Prediction = c(1, 1),
  Isoforms_Analysis = c(1, 0),
  Experimental_Data_Weight = c(1, 0),
  Multiple_Ranking_Scores = c(1, 0)
)

# Method 1: Heat Map (Recommended)
create_feature_heatmap <- function(data) {
  # Reshape data for ggplot
  data_long <- melt(data, id.vars = "Tool",
                    variable.name = "Feature",
                    value.name = "Support")

  # Clean up feature names for display
  data_long$Feature <- gsub("_", " ", data_long$Feature)
  data_long$Feature <- gsub("BLOSUM", "BLOSUM", data_long$Feature)
  data_long$Feature <- gsub("SNPs", "SNPs", data_long$Feature)
  data_long$Feature <- gsub("HLA", "HLA", data_long$Feature)

  # Create heat map
  ggplot(data_long, aes(x = Feature, y = Tool, fill = factor(Support))) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = ifelse(Support == 1, "✓",
                                 ifelse(Support == 0.5, "~", "✗"))),
              size = 6, fontface = "bold") +
    scale_fill_manual(values = c("0" = "#99CFE9", "0.5" = "#FBB800", "1" = "#A2C510"),
                      labels = c("Not Supported", "Partial", "Fully Supported"),
                      name = "Support Level") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.position = "bottom",
      panel.grid = element_blank()
    ) +
    labs(
      title = "Peptide Identification Tool Feature Comparison",
      subtitle = "EpiTox vs CrossDome Package",
      x = "Features",
      y = "Tools"
    )
}

# Method 2: Side-by-side Bar Chart
create_feature_barchart <- function(data) {
  data_long <- melt(data, id.vars = "Tool",
                    variable.name = "Feature",
                    value.name = "Support")

  # Clean feature names
  data_long$Feature <- gsub("_", " ", data_long$Feature)

  ggplot(data_long, aes(x = Feature, y = Support, fill = Tool)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("EpiTox" = "#A2C510", "CrossDome" = "#99CFE9")) +
    scale_y_continuous(breaks = c(0, 0.5, 1),
                       labels = c("No", "Partial", "Yes")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Feature Support Comparison",
      subtitle = "EpiTox vs CrossDome",
      x = "Features",
      y = "Support Level",
      fill = "Tool"
    )
}

# Method 3: Radar Chart (if you prefer)
create_radar_comparison <- function(data) {
  library(fmsb)

  # Prepare data for radar chart (remove Tool column, add max/min rows)
  radar_data <- data[, -1]  # Remove Tool column
  radar_data <- rbind(
    rep(1, ncol(radar_data)),    # max values
    rep(0, ncol(radar_data)),    # min values
    radar_data[1, ],             # Your tool
    radar_data[2, ]              # CrossDome
  )

  # Clean column names
  colnames(radar_data) <- gsub("_", "\n", colnames(radar_data))

  radarchart(radar_data,
             axistype = 1,
             pcol = c("#A2C510", "#99CFE9"),
             pfcol = c(rgb(0.2, 0.6, 0.9, 0.3), rgb(0.9, 0.3, 0.3, 0.3)),
             plwd = 3,
             plty = 1,
             cglcol = "grey70",
             cglty = 1,
             cglwd = 0.8,
             axislabcol = "grey50",
             vlcex = 0.8,
             caxislabels = c("No", "", "", "", "Yes"),
             title = "Peptide Tool Feature Comparison")

  legend("topright",
         legend = c("EpiTox", "CrossDome"),
         col = c("#A2C510", "#99CFE9"),
         lty = 1, lwd = 3, cex = 0.8)
}

# Generate visualizations
print("Creating heat map...")
p1 <- create_feature_heatmap(tool_comparison)
print(p1)

print("Creating bar chart...")
p2 <- create_feature_barchart(tool_comparison)
print(p2)

print("Creating radar chart...")
create_radar_comparison(tool_comparison)

# Summary statistics
cat("\n=== Feature Coverage Summary ===\n")
your_tool_coverage <- sum(tool_comparison[1, -1]) / (ncol(tool_comparison) - 1) * 100
crossdome_coverage <- sum(tool_comparison[2, -1]) / (ncol(tool_comparison) - 1) * 100

cat("EpiTox Coverage:", round(your_tool_coverage, 1), "%\n")
cat("CrossDome Coverage:", round(crossdome_coverage, 1), "%\n")

# Feature gap analysis
features <- names(tool_comparison)[-1]
your_gaps <- features[tool_comparison[1, -1] == 0]
crossdome_gaps <- features[tool_comparison[2, -1] == 0]

cat("\nFeatures missing in EpiTox:", ifelse(length(your_gaps) > 0, paste(your_gaps, collapse = ", "), "None"), "\n")
cat("Features missing in CrossDome:", ifelse(length(crossdome_gaps) > 0, paste(crossdome_gaps, collapse = ", "), "None"), "\n")
