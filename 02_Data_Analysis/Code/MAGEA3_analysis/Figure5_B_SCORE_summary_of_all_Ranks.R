library(dplyr)
library(tidyr)
library(ggplot2)

biocopy_colors = c("#A2C510", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")

# read in annotation and score
xscan = openxlsx::read.xlsx("/Users/hoor.alhasani/Documents/Projects/MAGEA3/results/Cutoff_4/Table/Xscan_MAGEA3.xlsx") %>%
  filter(isOfftarget == "No") %>%
  distinct(.keep_all = T)

score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  # aggregate over target
  group_by(peptide) %>%
  summarise(KD = mean(KD, na.rm =T),
            across(everything(), first),  # keep first value of all columns
            .groups = 'drop'
  ) %>%
  ungroup() %>% # not filtering for n because the majority are non-binders, by removing them the number of tested change
  filter (lev_dist <5, !peptide %in% c("KVADLIHFL", "TVAELVQFV"), !peptide %in% xscan$peptide) %>% # first peptide is the QC droped one, the second is the incorrect one (TVAELVQFL)
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  arrange(lev_dist)

epitox = openxlsx::read.xlsx("../../data/full_peptide_rankings.xlsx") %>%
  mutate(peptide = gsub(".+\\_", "", id),
         uniprot = gsub("\\_.+$", "", id),
         Tested = ifelse(peptide %in% score_reference$peptide, "Yes", "No")) %>%
  filter(confidence_level != "Low") %>%
  merge(., score_reference, all.x =T, by = "peptide") %>%
  relocate(uniprot, peptide) %>%
  relocate(Tested, KD, n, lev_dist, Outcome, log2_target_FC, .after = "multibiophys_category")


ggplot(epitox %>% filter(Outcome == "Binder"),
       aes(Bi_features_score, multibiophys_similarity, col = score_disagreement)) +
  geom_point(aes(shape = Wildtype)) +
  facet_wrap(~confidence_level) +
  theme_light() +
  #geom_smooth() +
  theme(legend.position = "bottom")

total_samples <- nrow(epitox)
proportion_data <- epitox %>%
  dplyr::rename(bayesian_rank = confidence_level) %>%
  mutate(bayesian_rank = ifelse(bayesian_rank == "Medium", "Moderate", bayesian_rank),
         Outcome = ifelse(is.na(Outcome), "NA", Outcome)) %>%
  select(Wildtype, anchor_status, bayesian_rank, Bi_feature_rank,
         multibiophys_category, Tested, Outcome) %>%
  pivot_longer(cols = c(Wildtype, anchor_status, bayesian_rank,
                        Bi_feature_rank, multibiophys_category),
               names_to = "Variable",
               values_to = "Category") %>%
  group_by(Variable, Category) %>%
  summarise(
    Total_in_category = n(),
    Tested_count = sum(Tested == "Yes"),
    Binder_count = sum(Outcome == "Binder", na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Pct_of_total = (Total_in_category / total_samples) * 100,
    Pct_tested = (Tested_count / Total_in_category) * 100,
    Pct_binders_of_tested = ifelse(Tested_count > 0,
                                   (Binder_count / Tested_count) * 100, 0),
    Pct_binders_of_total = (Binder_count / Total_in_category) * 100,
    Category = factor(Category, levels = c("High", "Moderate", "Low", "Very Low", "Yes" , "No", "Intact", "anchor_mismatch", "Both Anchor & Backbone" ))
  )

# Clean up variable names
proportion_data <- proportion_data %>%
  mutate(Variable = case_when(
    Variable == "Wildtype" ~ "Wildtype",
    Variable == "anchor_status" ~ "Anchor Status",
    Variable == "bayesian_rank" ~ "Experimental Evidence",
    Variable == "Bi_feature_rank" ~ "Bi-features Rank",
    Variable == "multibiophys_category" ~ "Multi-features Rank",
    TRUE ~ Variable  # keep original if not matched
  ))
# Reshape for stacked bar chart
plot_data_stacked <- proportion_data %>%
  mutate(
    Pct_not_tested = Pct_of_total - Pct_tested,
    Pct_non_binders = Pct_tested - Pct_binders_of_total,
    Pct_binders = Pct_binders_of_total
  ) %>%
  select(Variable, Category, Pct_not_tested, Pct_non_binders, Pct_binders) %>%
  pivot_longer(cols = c(Pct_not_tested, Pct_non_binders, Pct_binders),
               names_to = "Status",
               values_to = "Percentage") %>%
  mutate(Status = factor(Status,
                         levels = c("Pct_binders", "Pct_non_binders", "Pct_not_tested"),
                         labels = c("Binders", "Non-Binders (Tested)", "Not Tested")))

# Reshape for stacked bar chart
plot_data_stacked <- proportion_data %>%
  mutate(
    # Convert Pct_tested from percentage of category to percentage of total
    Pct_tested_of_total = (Pct_tested / 100) * Pct_of_total,

    # Now all calculations are in terms of percentage of total dataset
    Pct_not_tested = Pct_of_total - Pct_tested_of_total,
    Pct_non_binders = Pct_tested_of_total - Pct_binders_of_total,
    Pct_binders = Pct_binders_of_total
  ) %>%
  select(Variable, Category, Pct_not_tested, Pct_non_binders, Pct_binders) %>%
  pivot_longer(cols = c(Pct_not_tested, Pct_non_binders, Pct_binders),
               names_to = "Status",
               values_to = "Percentage") %>%
  mutate(Status = factor(Status,
                         levels = c("Pct_binders", "Pct_non_binders", "Pct_not_tested"),
                         labels = c("Binders", "Non-Binders (Tested)", "Not Tested")))

biocopy_colors = c("#A2C510", "#99CFE9", "#FBB800", "#939597", "#2C4255", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F")

# Create the stacked bar chart
p1 <- ggplot(plot_data_stacked, aes(x = Category, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(data = plot_data_stacked %>% filter(Percentage > 2),
            aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black", fontface = "bold") +
  facet_grid(. ~ Variable, scales = "free_x", space = "free_x") +
  labs(title = "Sample Distribution: Testing Status and Outcomes",
       subtitle = paste0("Total samples: ", total_samples,
                         " | Each bar represents 100% of samples in that category"),
       x = "Category",
       y = "Percentage of Total Samples",
       fill = "Status") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "lightgray")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  scale_fill_manual(values = c("Binders" = biocopy_colors[1],
                               "Non-Binders (Tested)" = biocopy_colors[5],
                               "Not Tested" = biocopy_colors[2]))
p1

ggsave("../../data/Figure_5_B_horz.pdf", plot = p1, width = 11, height = 5)

p1 <- ggplot(plot_data_stacked, aes(x = Category, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(data = plot_data_stacked %>% filter(Percentage > 2),
            aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black", fontface = "bold") +
  labs(#title = "Sample Distribution: Testing Status and Outcomes",
       #subtitle = paste0("Total samples: ", total_samples,
      #                   " | Each bar represents 100% of samples in that category"),
       x = "",
       y = "Percentage of Total Samples",
       fill = "Status") +
  theme_bw() +
  theme(
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "lightgray")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  scale_fill_manual(values = c("Binders" = biocopy_colors[1],
                               "Non-Binders (Tested)" = biocopy_colors[5],
                               "Not Tested" = biocopy_colors[2])) +
  facet_wrap(~Variable, ncol = 1, nrow = 5, "free_x")

print(p1)

ggsave("../../data/Figure_5_B.pdf", plot = p1, width = 4.5, height = 14)

p2 <- ggplot(plot_data_stacked, aes(x = Category, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(data = proportion_data,
            aes(x = Category, y = 102, label = paste0("n=", Total_in_category), fill = NULL),
            size = 3, fontface = "bold") +
  geom_text(data = plot_data_stacked %>% filter(Percentage > 2),
            aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 2.5, color = "white", fontface = "bold") +
  facet_grid(. ~ Variable, scales = "free_x", space = "free_x") +
  labs(title = "Sample Distribution: Testing Status and Outcomes",
       subtitle = paste0("Total samples: ", total_samples,
                         " | Each bar = 100% of samples in that category"),
       x = "Category",
       y = "Percentage of Total Samples in Category",
       fill = "Status") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "lightgray")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 110)) +
  scale_fill_manual(values = c("Binders" = biocopy_colors[1],
                               "Non-Binders (Tested)" = biocopy_colors[5],
                               "Not Tested" = biocopy_colors[2]))

print(p2)


openxlsx::write.xlsx(epitox, "../../data/OT_SCORE_full_annotation.xlsx")

df = epitox
# Filter only the tested samples (since untested ones have NA in Outcome)
tested_samples <- df[df$Tested == "Yes", ]

# Create a function to generate summary statistics
create_summary <- function(data, column_name, var_name) {
  # Create contingency table
  ct <- table(data[[column_name]], data$Outcome)

  # Calculate proportions by row
  prop <- prop.table(ct, margin = 1)

  # Create result data frame
  result <- data.frame(
    Variable = var_name,
    Category = rownames(ct),
    Binder_Count = ct[, "Binder"],
    NonBinder_Count = ct[, "Non-Binder"],
    Total = rowSums(ct),
    Binder_Prop = round(prop[, "Binder"], 3),
    NonBinder_Prop = round(prop[, "Non-Binder"], 3)
  )

  return(result)
}

# Create summaries for all variables
wildtype_summary <- create_summary(tested_samples, "Wildtype", "Wildtype")
anchor_summary <- create_summary(tested_samples, "anchor_status", "Anchor Status")
bayesian_summary <- create_summary(tested_samples, "confidence_level", "Experimental Evidence")
bifeature_summary <- create_summary(tested_samples, "Bi_feature_rank", "Bi-feature Rank")
multibiophys_summary <- create_summary(tested_samples, "multibiophys_category", "Multi-biophysics")

# Combine all summaries
complete_summary <- rbind(wildtype_summary,
                          anchor_summary,
                          bayesian_summary,
                          bifeature_summary,
                          multibiophys_summary)

# Display the complete table
print(complete_summary)

# For better readability, you can also create separate tables
cat("\n=== WILDTYPE ===\n")
print(wildtype_summary)

cat("\n=== ANCHOR STATUS ===\n")
print(anchor_summary)

cat("\n=== BAYESIAN RANK ===\n")
print(bayesian_summary)

cat("\n=== BI-FEATURE RANK ===\n")
print(bifeature_summary)

cat("\n=== MULTI-BIOPHYSICS ===\n")
print(multibiophys_summary)

# Also show the overall testing statistics
cat("\n=== TESTING OVERVIEW ===\n")
testing_overview <- table(df$Tested)
print(testing_overview)

cat("\n=== OUTCOME DISTRIBUTION (Tested samples only) ===\n")
outcome_dist <- table(tested_samples$Outcome)
print(outcome_dist)

# Export to CSV
#write.csv(complete_summary, "complete_ranking_analysis.csv", row.names = FALSE)

# Optional: Create a nice formatted table with knitr
library(knitr)
kable(complete_summary,
      col.names = c("Variable", "Category", "Binders", "Non-Binders",
                    "Total", "Binder Prop.", "Non-Binder Prop."),
      caption = "Comprehensive Analysis of All Variables vs Experimental Outcome")

library(knitr)
library(kableExtra)
rownames(combined_summary) = 1:nrow(combined_summary)
kable(combined_summary,
      col.names = c("Ranking Method", "Level", "Binders", "Non-Binders",
                    "Total", "Binder %", "Non-Binder %"),
      caption = "Performance of Different Ranking Methods") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
