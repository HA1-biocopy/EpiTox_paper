library(ggplot2)
library(dplyr)
library(lessR)

biocopy_colors = c("#A2C510", "#958BB2", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")

# Figure caption:
# "Groups with fewer than 5 peptides excluded from visualization"
df = openxlsx::read.xlsx("../../data/full_peptide_rankings.xlsx") %>%
  dplyr::rename(edit_distance = mismatch) %>%
  mutate(peptide = gsub(".+\\_", "", id),
         uniprot = gsub("\\_.+", "", id)) %>%
  distinct(peptide, uniprot, .keep_all = T) %>%
  relocate (peptide, .after = "id") %>%
  arrange (edit_distance)

cdf_data <- df %>%
  group_by(Bi_feature_rank, anchor_status) %>%
  arrange(Bi_features_score) %>%
  filter(n() >= 5) %>%
  mutate(
    cumulative_pct = (row_number() / n()) * 100,
    Bi_feature_rank = factor(Bi_feature_rank, levels = c("High", "Moderate", "Low"))
  ) %>%
  ungroup()


# Plot A: Bi -> CDF
g = ggplot(cdf_data, aes(x = Bi_features_score,
                         y = cumulative_pct,
                         color = anchor_status)) +
  geom_line(size = 1.5, alpha = 0.8) +
  facet_wrap(~Bi_feature_rank, ncol = 3) +
  scale_color_manual(
    values = biocopy_colors,
    name = "Anchor Status"
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray60", alpha = 0.7) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray60", alpha = 0.7) +
  labs(
    title = "Cumulative Distribution Functions (CDF)",
    #subtitle = "Bi_features_score by Risk Rank and Anchor Status",
    x = "Bi_features_score",
    y = "Cumulative Percentage (%)",
    caption = "three peptides with mismatch at the anchor positions were removed"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
g

# ======================================================
# Plot B Multi-feature

PieChart(multibiophys_category, data = df,
         hole = 0,
         fill = c("#00508A", "#F4FAFD", "#B0DAEE"),
         color="white",
         labels = "%",
         labels_position= "in",
         labels_size=1.5,
         labels_color = c("white","black","black"),
         labels_cex = 1.5,
         main = "Multi-feature score")

# ======================================================
# Plot 4 D
# Precompute counts + correlations for each facet
corr_method <- "spearman"

stats_per_facet <- df %>%
  filter(confidence_level != "Low") %>%
  group_by(confidence_level) %>%
  summarise(
    n = n(),
    corr = cor(rank_Bi_features, rank_multibiophys, method = corr_method)
  )


g = ggplot(
  df %>% filter(confidence_level != "Low"),
  aes(rank_Bi_features, rank_multibiophys)
) +
  geom_point(aes(fill = confidence_level), shape = 21, size = 2) +
  scale_fill_manual(values = c("#00508A", "#B0DAEE", "#F4FAFD", "gray70")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  theme_light() +
  theme(legend.position = "none") +
  facet_wrap(~confidence_level) +
  labs(
    x = "Rank by Bi-feature",
    y = "Rank by Multi-feature",
    title = "Faceted distribution of scoring metrics across evidence level"
  ) +
  geom_label(
    data = stats_per_facet,
    aes(
      x = Inf, y = -Inf,
      label = paste0("n = ", n, "\n\u03C1 = ", round(corr, 2))
    ),
    hjust = 1.1,
    vjust = -0.6,
    size = 4,
    label.size = 0.3,       # border thickness
    fill = "white",         # background
    alpha = 0.8
  )
g
ggsave(plot = g, filename = "../../data/Figure3_D_3Scores_comparisons.pdf", device = cairo_pdf, width = 9, height = 4.5)

# table 4
df %>% filter(known_peptide == "Known") %>%select(confidence_level) %>%
  table # gives 19 due to protein duplicates
df %>% filter(known_peptide == "Known") %>% select(confidence_level, peptide) %>%
  distinct() %>% select(confidence_level) %>% table # gives the correct 18
df %>% filter(known_peptide == "Known", confidence_level == "High") %>%
  distinct(peptide, .keep_all = T) %>% select(multibiophys_category, Bi_feature_rank) %>%
  table
df %>% filter(known_peptide == "Known", confidence_level == "Medium") %>%
  distinct(peptide, .keep_all = T) %>% select(multibiophys_category, Bi_feature_rank) %>%
  table
df %>% filter(known_peptide == "Known", confidence_level == "Very Low") %>%
  distinct(peptide, .keep_all = T) %>% select(multibiophys_category, Bi_feature_rank) %>%
  table
df %>% filter(known_peptide == "Known") %>%
  distinct(peptide, .keep_all = T) %>% select(multibiophys_category, Bi_feature_rank, confidence_level) %>%
  table

compare_scores <- df %>%
  #filter(anchor_status == "Intact") %>%  # Remove this line if you want all peptides
  group_by(confidence_level) %>%
  summarise(
    N = n(),
    Mu_bi = mean(Bi_features_score, na.rm = TRUE),
    Mu_multi = mean(multibiophys_similarity, na.rm = TRUE),
    SD_bi = sd(Bi_features_score, na.rm = TRUE),
    SD_multi = sd(multibiophys_similarity, na.rm = TRUE),
    #anchor = ifelse(anchor_status == "Intact", n(), 0)
  ) %>%
  mutate(Percent = round(N / sum(N) * 100, 2))

# View it nicely
compare_scores


# Are Bi-feature scores significantly different across evidence levels?
kruskal.test(Bi_features_score ~ confidence_level, data = df %>%
               filter(anchor_status == "Intact"))

# Pairwise comparisons
pairwise.wilcox.test(
  df %>% filter(anchor_status == "Intact") %>% pull(Bi_features_score),
  df %>% filter(anchor_status == "Intact") %>% pull(confidence_level),
  p.adjust.method = "BH"
)

# Define high-risk thresholds (adjust based on your data)
high_bi_threshold <- 0.6  # Above 60th percentile or domain knowledge
high_multi_threshold <- 0.6

priority_categories <- df %>%
  mutate(
    Bi_risk = ifelse(Bi_features_score >= high_bi_threshold, "High", "Low"),
    Multi_risk = ifelse(multibiophys_similarity >= high_multi_threshold, "High", "Low"),
    Risk_category = case_when(
      Bi_risk == "High" & Multi_risk == "High" & confidence_level == "Very Low" ~ "Priority Unknown",
      Bi_risk == "High" & Multi_risk == "High" & confidence_level %in% c("High", "Medium") ~ "Known High Risk",
      Bi_risk == "High" & Multi_risk == "Low" ~ "Sequence-driven",
      Bi_risk == "Low" & Multi_risk == "High" ~ "Physicochemical-driven",
      TRUE ~ "Lower Priority"
    )
  ) %>%
  count(Risk_category) %>%
  mutate(Percent = round(n / sum(n) * 100, 2))

priority_categories





