library(ggplot2)
library(dplyr)

cutoff_5 = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_5/Tables/PrediTopes_MAGEA3.xlsx") %>%
  distinct(chr, uniprot, Gene.Names, peptide, mismatch, Wildtype, .keep_all = T) %>%
  arrange(mismatch)

# read in patents

sheets <- openxlsx::getSheetNames("~/Documents/Projects/MAGEA3/patent.xlsx")
patents_list = lapply(sheets, function(sheet){
  patent = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/patent.xlsx", sheet = sheet) %>%
    mutate(mismatch = stringdist::stringdist("KVAELVHFL", peptide, method = "hamming"),
           source = sheet) %>%
    distinct(gene, peptide, source, .keep_all = T) %>%
    filter(!is.na(peptide))
})

names(patents_list) = sheets

patents_df = do.call("rbind", patents_list) %>%
  distinct(gene, peptide, source, .keep_all = T) %>%
  mutate(EpiTox = ifelse(peptide %in% cutoff_5$peptide, "Yes", "No")) %>%
  filter(!is.na(peptide)) %>%
  select(-source, -Note) %>%
  distinct (.keep_all = T) %>%
  arrange(mismatch)
rownames(patents_df) = 1:nrow(patents_df)

patents_df = merge(patents_df, cutoff_5[, c("peptide", "Wildtype")], by = "peptide", all.x =T) %>%
  arrange(mismatch)

x = patents_df$peptide %>% unique()

cutoff_5[cutoff_5$peptide %in% x,]

venn_list = list(EpiTox = patents_df$peptide[patents_df$EpiTox == "Yes"],
                 Patent = patents_df$peptide)

ggvenn::ggvenn(venn_list, show_percentage = F,
               fill_color = c("#A2C510", "#FBB800", "#99CFE9", "#2C4255"), stroke_color = "white",
               set_name_size = 3)
write.csv (patents_df, "../../data/known_patents_OTs.csv", row.names = F)

#==============================

cutoff_5 = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_5/Tables/HPA_genes_nTPM.xlsx") %>%
  distinct(chr, uniprot, Gene.Names, peptide, mismatch, Wildtype, .keep_all = T) %>%
  arrange(mismatch)

score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  arrange(lev_dist) %>%
  select(labguru_id, peptide, KD, n, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  arrange(KD)

df = merge(score_reference, cutoff_5[, c("uniprot", "Gene.Names", "peptide", "Wildtype",
                        "processing_score", "affinity", "blosum_similarity",
                        "Rank", "Ranking_score")], by = "peptide", all.x = T) %>%
  distinct(.keep_all = T)


ggplot(df, aes(affinity, KD, col = Wildtype)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_light()

ggplot(df, aes(blosum_similarity, KD, col = Wildtype)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_light()




