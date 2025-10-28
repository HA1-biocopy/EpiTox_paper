
epitox = openxlsx::read.xlsx("../../data/full_peptide_rankings.xlsx") %>%
  mutate(peptide = gsub(".+\\_", "", id),
         uniprot = gsub("\\_.+$", "", id),
         id = NULL) %>%
  relocate(uniprot, peptide)

score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  # aggregate over target
  filter (lev_dist <5, !peptide %in% c("KVADLIHFL", "TVAELVQFV")) %>% # first peptide is the QC droped one, the second is the incorrect one (TVAELVQFL)
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  merge(., epitox, by = "peptide") %>%
  # mutate(Outcome =
  #          case_when(is.na(Outcome) ~"Tested - No KD",
  #                    n<4 ~ "Bad QC",
  #                    TRUE ~ Outcome)) %>%
  arrange(lev_dist)

epitox = epitox %>%
  merge(., score_reference[score_reference$n >3, c("labguru_id", "peptide", "KD", "n", "lev_dist", "Outcome", "log2_target_FC", "n_slides_QC_passed")], by = "peptide", all.x = T) %>%
  mutate(Tested = ifelse(!is.na(labguru_id), "Yes", "No")) %>%
  relocate(Tested, Outcome, .after = "peptide")

epitox %>%
  arrange(rank_multibiophys) %>%
  filter(is.na(Tested))
  select(confidence_level) %>%
  table()


epitox %>%
  arrange(rank_multibiophys) %>%
  filter(rank_Bi_features <500, !peptide %in% score_reference$peptide) %>%
  head()

ggplot(score_reference %>% filter(Outcome == "Binder"),
       aes(Bi_features_score, multibiophys_distance, col = score_disagreement)) +
  geom_point(aes(shape = Wildtype)) +
  facet_wrap(~confidence_level) +
  theme_light() +
  #geom_smooth() +
  theme(legend.position = "bottom")
