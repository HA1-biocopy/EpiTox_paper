library(ggplot2)
library(dplyr)


epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  select(ensembl_gene_id, uniprot, Gene.Names, peptide, blosum_similarity, mismatch,
         Wildtype, Peptide_HLA_Atlas, affinity, presentation_score, Rank, Ranking_score) %>%
  mutate(id = paste0(uniprot, "_", peptide)) %>%
  distinct(id, .keep_all = T)


dups_peptides = epitox$peptide[duplicated(epitox$peptide)]
dups_df = epitox %>%
  filter(peptide %in% dups_peptides, Peptide_IEDB == "Yes") %>%
  mutate(id = paste0(uniprot, "_", peptide)) %>%
  relocate(id)


bayesian = openxlsx::read.xlsx("/Users/hoor.alhasani/Documents/Projects/D003/Patent_Paper/paper_materials/02_Data_Analysis/data/Bayesian/MAGEA3_full_results.xlsx") %>%
  dplyr::rename(id = peptide_id) %>%
  filter(id %in% dups_df$id)

openxlsx::write.xlsx(bayesian, "../../data/peptide_identity_bayesian.xlsx")

ggplot(bayesian, aes(peptide, lr_product, col = confidence_level)) +
  geom_point() +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()


bayesian %>%
  arrange(peptide) %>%
  select(id, confidence_level)

bayesian %>%
  arrange(peptide) %>%
  select(id, confidence_level, Gene.Names) %>%
  distinct(.keep_all = T)

# P50281_AVHELGHAL           Medium                              MMP14
# Q9Y5R2_AVHELGHAL         Very Low                      MMP24 MT5MMP
# P51512_AVHELGHAL         Very Low                MMP16 C8orf57 MMPX2
# P51511_AVHELGHAL         Very Low                              MMP15
#
# Q15173_GVAELLEIL         Very Low                            PPP2R5B
# Q15172_GVAELLEIL           Medium                            PPP2R5A
# Q16537_GVAELLEIL         Very Low                            PPP2R5E
#
# Q9H361_KVDEAVAVL         Very Low                PABPC3 PABP3 PABPL3
# Q13310_KVDEAVAVL             High                  PABPC4 APP1 PABP4
# P11940_KVDEAVAVL             High      PABPC1 PAB1 PABP PABP1 PABPC2
#
# Q8NFP9_SVVELVMLL           Medium          NBEA BCL8B KIAA1544 LYST2
# P50851_SVVELVMLL         Very Low                 LRBA BGL CDC4L LBA
#

