library(dplyr)

SCORE = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  arrange(lev_dist)

cutoff_5 = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_5/Tables/HPA_genes_nTPM.xlsx") %>%
  filter (peptide %in% SCORE$peptide) %>%
  distinct(chr, uniprot, Gene.Names, peptide, mismatch, Wildtype, .keep_all = T) %>%
  merge(., SCORE, by = "peptide", all.y =T) %>%
  arrange(mismatch) %>%
  mutate(Pep_n = NULL, rep_pep = NULL, labguru_id = NULL,
         HLA_Atlas_allele = NULL, Protein_HLA_Atlas = NULL,
         IEDB_allele = NULL, Protein_IEDB = NULL, binder = NULL, presented = NULL,
         Peptide_HLA_Atlas = ifelse(Peptide_HLA_Atlas == "Yes", 1, 0),
         Peptide_IEDB = ifelse(Peptide_IEDB == "Yes", 1, 0),
         Wildtype = ifelse(Wildtype == "Yes", 1, 0),
         Wildtype = ifelse(Wildtype == "Yes", 1, 0),
         Wildtype = ifelse(Wildtype == "Yes", 1, 0),
         Outcome = case_when(
           Outcome == "Binder" ~ 1,
           Outcome == "Non-Binder" ~ 0,
           Outcome == "Ambiguous" ~ 0.5,
           TRUE ~ NA
           )
         )

pairs(cutoff_5[,c(10:18, 59:ncol(cutoff_5))], upper.panel = NULL, col = "steelblue1")
ggplot(cutoff_5, aes(KD, blosum_similarity)) +
  geom_point(aes(col = as.factor(Outcome)), alpha = 0.7) +
  #geom_jitter() +
  theme_light()

pairs(cutoff_5[,c(19:40)], upper.panel = NULL, col = "steelblue1")
pairs(cutoff_5[,c(41:58)], upper.panel = NULL, col = "steelblue1")




