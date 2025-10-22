epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  select(ensembl_gene_id, uniprot, Gene.Names, peptide, blosum_similarity, mismatch,
         Wildtype, Peptide_HLA_Atlas, affinity, presentation_score, Rank, Ranking_score) %>%
  distinct(.keep_all = T) %>%
  mutate(id = paste0(uniprot, "_", peptide)) %>%
  dplyr::rename()%>%
  #merge(., experimental_df, by = "id", all.x =T) %>%
  distinct(peptide, uniprot, .keep_all = T)

bayesian = openxlsx::read.xlsx("/Users/hoor.alhasani/Documents/Projects/D003/Patent_Paper/paper_materials/02_Data_Analysis/data/Bayesian/MAGEA3_full_results.xlsx") %>%
  #mutate(peptide = gsub(".+\\_", "", peptide_id)) %>%
  dplyr::rename(id = peptide_id) %>%
  select(id, evidence, disease_tissue, evidence_chain, confidence_level, interpretation, lr_product) %>%
  merge(., epitox, by = "id", all.x = T)

score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  filter (!is.na(KD), n>3) %>%
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  arrange(lev_dist)


# need the xscan to annotate the data properly

df = merge(score_reference[, c("peptide", "KD", "Outcome", "n", "lev_dist")],
           bayesian, by = "peptide", all =T) %>%
  mutate(Rank = ifelse(is.na(Rank), "Other", Rank),
    Rank = factor(Rank, levels = c("High", "Moderate", "Low", "Random", "Other")))


ggplot(df, aes(confidence_level, KD*10^5)) +
  geom_boxplot(aes(fill = Outcome), position = "dodge") +
  geom_jitter(aes(fill = Outcome), shape = 21) +
  scale_fill_manual(values = c(#"Binder" = "#C61E19",
                               "NA" = "#FDE399",
                               "Binder" = "#8ECAE6",
                               "Non-Binder" = "#A2C510")) +
  labs(x = "Bayesian weight probability") +
  theme_light() +
  facet_wrap(~Wildtype)

ggplot(df, aes(confidence_level, KD*10^5)) +
  geom_boxplot(aes(fill = Outcome), position = "dodge") +
  geom_jitter(aes(fill = Outcome), shape = 21) +
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
