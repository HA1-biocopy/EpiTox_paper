library(ggplot2)
library(dplyr)
library(lessR)

peprep = read.csv("~/Documents/Projects/PEPREP/data/PEPREP_peptides_aggregation.csv")
epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  distinct(.keep_all = T) %>%
  mutate(PEPREP = ifelse(peptide %in% peprep$peptide, T, F),
         IEDB = ifelse(Peptide_IEDB == "Yes", T, F),
         category = case_when(
           PEPREP & IEDB ~ "PEPREP + IEDB",
           PEPREP ~"PEPREP",
           IEDB ~ "IEDB",
           TRUE ~ "None"),
         category = factor(category, levels = c("None", "PEPREP", "PEPREP + IEDB", "IEDB")),
         id = paste0(uniprot, "_", peptide)) %>%
  distinct(peptide, uniprot, .keep_all = T)

df = openxlsx::read.xlsx("../../data/full_peptide_rankings.xlsx") %>%
  dplyr::rename(edit_distance = mismatch) %>%
  mutate(peptide = gsub(".+\\_", "", id),
         uniprot = gsub("\\_.+", "", id),
         inPEPREP = ifelse(peptide %in% peprep$peptide, T, F),
         inIEDB = ifelse(peptide %in% epitox$peptide[epitox$Peptide_IEDB == "Yes"], T, F),
         category = case_when(
           inPEPREP & inIEDB ~ "PEPREP + IEDB",
           inPEPREP ~"PEPREP",
           inIEDB ~ "IEDB",
           TRUE ~ "None"),
         category = factor(category, levels = c("None", "PEPREP", "PEPREP + IEDB", "IEDB"))) %>%
  distinct(peptide, uniprot, .keep_all = T) %>%
  relocate (peptide, inPEPREP, inIEDB, .after = "id") %>%
  arrange (edit_distance)


df %>% select(peptide, inPEPREP) %>% distinct() %>% select (inPEPREP) %>% table
df %>% select(peptide, inIEDB) %>% distinct() %>% select (inIEDB) %>% table
df$category %>% table

# Count and calculate proportions
counts <- epitox %>%
  count(category) %>%
  mutate(prop = n / sum(n)) %>%
  mutate(label = paste0(n, " (", round(prop * 100, 1), "%)"),
         ypos = 1- (cumsum(prop) - prop / 2))

biocopy_colors = c("#A2C510", "#958BB2", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")

# Plot A
g = ggplot(counts, aes(x = 1, y = prop, fill = category)) +
  geom_bar(stat = "identity", width = 0.4, color = "white") +
  geom_text(aes(x = 1, y = ypos, label = label),
            size = 3.3, fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Proportion of peptides found in MHC databases") +
  #coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.justification = c(0.7,0),
    #legend.box.just = "left",
    #legend.key.size = unit(0.3, "cm"),      # Smaller legend boxes
    legend.box.spacing = unit(0.01, "cm") ,
    #legend.text = element_text(size = 7),    # Smaller legend text
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c("#E0E1E1", "#99CFE9", "#2C4255", "#A2C510")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave(plot = g, "../../data/Figure_3_A.pdf")

# Plot B

df = openxlsx::read.xlsx("../../data/Bayesian/MAGEA3_full_results.xlsx") %>%
  dplyr::rename(edit_distance = mismatch) %>%
  mutate(peptide = gsub(".+\\_", "", peptide_id),
         uniprot = gsub("\\_.+", "", peptide_id)) %>%
  distinct(peptide, uniprot, .keep_all = T) %>%
  relocate (uniprot, peptide, .after = "peptide_id") %>%
  arrange (edit_distance)

plot_summary = df %>%
  group_by(confidence_level) %>%
  summarise(Mu = mean(posterior_prob),
            Counts = n(),
            Prp = paste0(round(n()/nrow(df), 2)*100, "%")) %>%
  mutate(confidence_level = factor(confidence_level, levels = c("High", "Medium", "Low", "Very Low")))

g = ggplot(plot_summary, aes(x = confidence_level, y = Mu)) +
  geom_col(aes(fill = confidence_level)) +
  geom_text(aes(label = paste("n=", Counts, "\n", Prp)),
            vjust =0.5,
            size = 3) +
  labs(x = "",
       y = "Mean Posterior Probability (%)",
       title = "Mean Posterior Probability by Evidence-level") +
  theme_minimal() +
  scale_fill_manual(values = c("#00508A", "#B0DAEE", "#DDEFF8", "#F4FAFD")) +
  theme(legend.position = "none")

ggsave(plot = g, "../../data/Figure_3_C.pdf")















# ===================== based on data/full....xlsx table

PieChart(category, data = df,
         hole = 0,
         fill = c("#E0E1E1", "#99CFE9", "#A2C510", "#2C4255"),
         color="white",
         main = "",
         labels = "input",
         labels_position= "in",
         labels_size=1.5,
         labels_color = c("black", "white"),
         labels_cex = 1.5)
# Plot B
PieChart(confidence_level, data = df,
         hole = 0,
         fill = c("#A2C510", "#99CFE9", "#2C4255", "#E0E1E1"),
         color="white",
         main = "",
         labels = "input",
         labels_position= "in",
         labels_size=1.5,
         labels_color = c("black", "white"),
         labels_cex = 1.5)

# Plot C
PieChart(Bi_feature_rank, data = df,
         hole = 0,
         fill = c("#A2C510", "#99CFE9", "#2C4255", "#E0E1E1"),
         color="white",
         main = "Bi_features score",
         labels = "input",
         labels_position= "in",
         labels_size=1.5,
         labels_color = c("black", "white"),
         labels_cex = 1.5)

library(ggridges)

#
ggplot(df, aes(x = Bi_features_score, y = anchor_status, fill = anchor_status)) +
  geom_density_ridges(
   # alpha = 0.6,
    stat = "binline",
    bins = 30,
    scale = 0.9,  # Controls overlap - higher = more overlap
    rel_min_height = 0.01  # Removes very small ridges
  ) +
  theme_ridges() +
  scale_fill_manual(values = biocopy_colors) +
  geom_vline(xintercept = .6, linetype = "dashed", color = "gray30", size = 1) +  # Between High and Moderate
  geom_vline(xintercept = .31, linetype = "dashed", color = "gray30", size = 1) +  # Between Moderate and Low
  theme(legend.position = "none")





