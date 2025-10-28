library(ggplot2)
library(ggpubr)
library(dplyr)

epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  select(ensembl_gene_id, uniprot, Gene.Names, peptide, blosum_similarity, mismatch,
         Wildtype, Peptide_HLA_Atlas, affinity, presentation_score, Rank, Ranking_score) %>%
  mutate(id = paste0(uniprot, "_", peptide)) %>%
  distinct(id, .keep_all = T)

bayesian = openxlsx::read.xlsx("/Users/hoor.alhasani/Documents/Projects/D003/Patent_Paper/paper_materials/02_Data_Analysis/data/Bayesian/MAGEA3_full_results.xlsx") %>%
  dplyr::rename(id = peptide_id) %>%
  merge(., epitox[, c("id", "Rank", "Ranking_score")], by = "id", all = T)


ggplot(bayesian, aes(confidence_level, Ranking_score)) +
  geom_boxplot(aes(fill = Rank), position = "dodge") +
  geom_jitter(aes(fill = Rank), shape = 21) +
  labs(x = "Bayesian weight probability") +
  theme_light() +
  facet_wrap(~Wildtype)

bayesian %>% group_by(confidence_level) %>% summarise(posterior_prob_mean = mean(posterior_prob),
                                                      SD = sd(posterior_prob, na.rm =T),
                                                      count = n())


score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  filter (lev_dist <5, !peptide %in% c("KVADLIHFL", "TVAELVQFV")) %>% # first peptide is the QC droped one, the second is the incorrect one (TVAELVQFL)
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  mutate(Outcome =
           case_when(is.na(Outcome) ~"No KD",
                     n<4 ~ "Bad QC",
                     TRUE ~ Outcome)) %>%
  arrange(lev_dist)


# need the xscan to annotate the data properly
xscan = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/Xscan_MAGEA3.xlsx") %>%
  filter(isOfftarget == "No")


df = merge(score_reference[, c("peptide", "KD", "Outcome", "n", "lev_dist")],
           bayesian, by = "peptide", all =T) %>%
  filter(!peptide %in% xscan$peptide) %>%
  mutate(Rank = case_when(
    is.na(Rank) ~ "Other",
    TRUE ~ Rank ),
    Rank = factor(Rank, levels = c("High", "Moderate", "Low", "Random", "Other")),
    target = ifelse(peptide == "KVAELVHFL", "Yes", "No"),
    Outcome = ifelse(is.na(Outcome), "Not tested", Outcome)
    ) %>%
  group_by(peptide) %>%
  mutate(KD = ifelse(peptide == "KVAELVHFL", mean(KD, na.rm =T), KD)) %>%
  slice(1) %>%
  ungroup()
df$Outcome %>% table

ggplot(df, aes(confidence_level, KD)) +
  geom_boxplot(aes(fill = Rank), position = "dodge") +
  geom_jitter(aes(fill = Rank), shape = 21) +
  labs(x = "Bayesian weight probability") +
  theme_light() +
  facet_wrap(~Outcome)

ggplot(df, aes(confidence_level, KD*10^5)) +
  geom_boxplot(aes(fill = Outcome), position = "dodge") +
  geom_jitter(aes(fill = Outcome, size = target), shape = 21) +
  scale_fill_manual(values = c(#"Binder" = "#C61E19",
                               "Bad QC" = "#FDE399",
                               "Ambiguous" = "#FDE399",
                               "Binder" = "#8ECAE6",
                               "Non-Binder" = "#A2C510")) +
  labs(x = "Bayesian weight probability") +
  theme_light() +
  facet_wrap(~Wildtype)

df %>% select(Outcome, confidence_level) %>% table
# Create the summary
df_binders_summary <- df %>%
  filter(Outcome == "Binder") %>%
  filter(confidence_level %in% c("High", "Medium", "Very Low")) %>%
  count(confidence_level)

# Then plot
ggplot(df_binders_summary, aes(x = confidence_level, y = n, fill = confidence_level)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5, size = 5) +
  labs(x = "Evidence Level", y = "Number of Confirmed Binders") +
  theme_light() +
  scale_fill_manual(values = c( "#FDE399", "#8ECAE6", "#A2C510")) +
  theme(legend.position = "none")


ggplot(df, aes(confidence_level, KD*10^5)) +
  geom_boxplot(aes(fill = Outcome), position = "dodge") +
  geom_jitter(aes(fill = Outcome, size = target), shape = 21) +
  scale_fill_manual(values = c(#"Binder" = "#C61E19",
    "NA" = "#FDE399",
    "Binder" = "#8ECAE6",
    "Non-Binder" = "#A2C510")) +
  ggtitle("Bi-feature ranking score vs. Experimental evidence vs KD values") +
  theme_light() +
  labs(x = "Experimental evidence", caption = "n >3") +
  facet_wrap(~Rank, scales = "free_x")


ggplot(df, aes(confidence_level, lev_dist)) +
  geom_boxplot(aes(fill = Outcome), position = "dodge") +
  geom_jitter(aes(fill = Outcome), shape = 21) +
  scale_fill_manual(values = c(#"Binder" = "#C61E19",
    "NA" = "#FDE399",
    "Binder" = "#8ECAE6",
    "Non-Binder" = "#A2C510")) +
  labs(x = "Bayesian weight probability") +
  theme_light() +
  facet_wrap(~Wildtype)


# Prepare data (exclude Not tested)
df_tested <- df %>%
  filter(Outcome %in% c("Binder", "Non-Binder", "Ambiguous", "Bad QC")) %>%
  filter(confidence_level != "Low") # only 0 tested

# Create proportional stacked bar
ggplot(df_tested, aes(x = confidence_level, fill = Outcome)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Evidence Level",
       y = "Proportion of Tested Peptides",
       title = "Experimental Binding Validation Across Evidence Categories") +
  scale_fill_manual(values = c("Binder" = "#2ECC71",
                               "Non-Binder" = "#95A5A6",
                               "Ambiguous" = "#F39C12",
                               "Bad QC" = "#E74C3C")) +
  theme_minimal()

ggplot(df_tested, aes(x = confidence_level, fill = Outcome)) +
  geom_bar(position = "dodge") +
  labs(x = "Evidence Level",
       y = "Number of Peptides",
       title = "Binding Outcomes by Evidence Category") +
  scale_fill_manual(values = c("Binder" = "#2ECC71",
                               "Non-Binder" = "#95A5A6",
                               "Ambiguous" = "#F39C12",
                               "Bad QC" = "#E74C3C")) +
  theme_minimal() +
  geom_text(stat='count', aes(label=after_stat(count)),
            position=position_dodge(width=0.9), vjust=-0.5, size=3)


# For binders only, show KD distribution
df_binders <- df %>% filter(Outcome == "Binder")

ggplot(df_binders, aes(x = confidence_level, y = KD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_y_log10() +  # KD values are usually log-scale
  labs(x = "Evidence Level",
       y = "KD (nM)",
       title = "Binding Affinity Does Not Differ by Evidence Level") +
  theme_minimal() +
  stat_compare_means()  # adds p-value from Kruskal-Wallis

# Fisher's exact test (or Chi-square if n large enough)
# Compare Binder vs Non-Binder across evidence levels

df_test <- df %>%
  filter(Outcome %in% c("Binder", "Non-Binder")) %>%
  filter(confidence_level %in% c("High", "Medium", "Very Low"))

# Create contingency table
table_test <- table(df_test$confidence_level, df_test$Outcome)

# Fisher's exact test
fisher.test(table_test, simulate.p.value = TRUE)
# or
chisq.test(table_test)


# Proportional test (prop.test)
# Comparing High vs Very Low

binders <- c(17, 15)  # High, Very Low
tested <- c(93, 224)  # High, Very Low

prop.test(binders, tested)

# For all three groups
binders_all <- c(17, 3, 15)  # High, Medium, Very Low
tested_all <- c(93, 41, 224)

prop.test(binders_all, tested_all)


# Check for issues:
table(df$confidence_level, df$Rank)

# Common causes of NA:
# 1. Expected frequencies < 5
# 2. Missing values

# Try Fisher's exact instead:
fisher.test(table(df$confidence_level, df$Rank),
            simulate.p.value = TRUE)

# Or remove Low category (only 1 peptide):
df_clean <- df %>% filter(confidence_level != "Low")
chisq.test(table(df_clean$confidence_level, df_clean$Rank))
