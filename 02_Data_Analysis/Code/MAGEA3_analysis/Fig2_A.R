# Load required libraries
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggbreak)

# Assuming your dataframe is called 'df' with columns:
# - edit_distance (or mismatches): numeric
# - Wildtype: TRUE/FALSE
# - is_known_peptide: TRUE/FALSE

# ==============================================================================
# OPTION 1: Broken/Split Y-axis
# ==============================================================================

df = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  dplyr::rename(edit_distance = mismatch)
df = cutoff_4 %>%
  dplyr::rename(edit_distance = mismatch)

# ==============================================================================
# Show edit distance for all three categories
# ==============================================================================

# Panel i: Mismatch histogram - now showing three categories
df_summary <- df %>%
  mutate(
    category = case_when(
      is_known_offtarget ~ "Previously disclosed",
      Wildtype == "Yes" ~ "Wildtype (novel)",
      Wildtype == "No" ~ "SNP-derived (novel)"
    ),
    category = factor(category,
                      levels = c("Wildtype (novel)", "SNP-derived (novel)", "Previously disclosed"))
  ) %>%
  count(edit_distance, category)

# Define colors
color_palette <- c(
  "Wildtype (novel)" = "#A2C510",
  "SNP-derived (novel)" = "#2C4255",
  "Previously disclosed" = "#FBB800"
)

# Determine break points
max_wt <- max(df_summary$n[df_summary$category == "Wildtype (novel)"], na.rm = TRUE)
max_other <- max(df_summary$n[df_summary$category != "Wildtype (novel)"], na.rm = TRUE)

# Set break range (adjust these values based on your data)
break_lower <- max_other * 1.2  # Just above the tallest non-wildtype bar
break_upper <- max_wt * 0.85    # Just below the wildtype peak

p5_panel1 <- ggplot(df_summary, aes(x = edit_distance, y = n, fill = category)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_y_break(c(break_lower, break_upper), scales = 0.5) +  # This creates the break
  labs(
    x = "Sequence mismatches",
    y = "Number of peptides",
    fill = ""
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.title = element_text(size = 10, face = "bold", hjust = 0)
  ) +
  ggtitle("I")

# Panel ii: Percentage breakdown showing novel vs. previously disclosed
summary_stats <- df %>%
  group_by(Wildtype) %>%
  summarise(
    total = n(),
    known = sum(is_known_offtarget),
    novel = total - known,
    pct_known = (known / total) * 100,
    pct_novel = (novel / total) * 100,
    .groups = 'drop'
  ) %>%
  mutate(sequence_type = factor(ifelse(Wildtype == "Yes", "Wildtype", "SNP-derived"),
                                levels = c("Wildtype", "SNP-derived")))

# Reshape for stacked bar
plot_data <- summary_stats %>%
  pivot_longer(
    cols = c(known, novel),
    names_to = "status",
    values_to = "count"
  ) %>%
  mutate(
    status = factor(status,
                    levels = c("novel", "known"),
                    labels = c("Novel", "Previously disclosed")),
    percentage = (count / total) * 100
  )

p5_panel2 <- ggplot(plot_data, aes(x = sequence_type, y = percentage, fill = status)) +
  geom_col(width = 0.6) +
  scale_fill_manual(
    values = c("Novel" = "#A2C510", "Previously disclosed" = "#FBB800"),
    breaks = c("Previously disclosed", "Novel")
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 108),
    breaks = seq(0, 100, 25)
  ) +
  labs(
    x = "",
    y = "Percentage (%)",
    fill = ""
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  # Add total count labels above bars
  geom_text(
    data = summary_stats,
    aes(x = sequence_type, y = 104, label = paste0("n=", total)),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  # Add percentage and count labels on bars
  geom_text(
    data = plot_data %>% filter(percentage > 2),
    aes(
      label = ifelse(
        percentage > 7,
        sprintf("%.1f%%\n(n=%d)", percentage, count),
       sprintf("%.1f%%", percentage)
      )
    ),
    position = position_stack(vjust = 0.7),
    size = 3,
    color = ifelse(
      (plot_data %>% filter(percentage > 2))$status == "Novel",
      "grey20",
      "white"
    ),
    fontface = "bold"
  ) +
  ggtitle("II")

# Combine panels side by side
option5_plot <- p5_panel1 + p5_panel2 +
  plot_layout(widths = c(1.3, 1))

print(option5_plot)

# Save
# ggsave("figure2A_option5.pdf", option5_plot, width = 10, height = 4.5, dpi = 300)
# ggsave("figure2A_option5.png", option5_plot, width = 10, height = 4.5, dpi = 300)


# ==============================================================================
# Alternative: If the 18 are too small to see in Panel A
# ==============================================================================

# Create a version with log scale for better visibility
p5_panel1_log <- ggplot(df_summary, aes(x = edit_distance, y = n, fill = category)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = color_palette) +
  scale_y_log10(
    expand = c(0, 0),
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000")
  ) +
  labs(
    x = "Sequence mismatches",
    y = "Number of peptides (log scale)",
    fill = ""
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.title = element_text(size = 10, face = "bold", hjust = 0)
  ) +
  ggtitle("I")

# Combine with log scale version
option5_plot_log <- p5_panel1_log + p5_panel2 +
  plot_layout(widths = c(1.3, 1))

print(option5_plot_log)
#ggsave("figure2A_option5_logscale.pdf", option5_plot_log, width = 10, height = 4.5, dpi = 300)

