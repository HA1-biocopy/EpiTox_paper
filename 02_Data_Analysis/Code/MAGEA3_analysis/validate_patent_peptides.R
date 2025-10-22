library(ggplot2)
library(dplyr)

# load score
# load epitox
# load crossdome results
# load patent-peptides

epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  dplyr::rename(edit_distance = mismatch) %>%
  distinct(uniprot, peptide, Wildtype, .keep_all = T)

# mark the peptides according to being detected or not by crossdome
# extract gene expression and compare:
# gxp,
# hla-info
# KD value
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
  mutate(peptide = ifelse(peptide == "TVAELVQFV", "TVAELVQFL", peptide),
         EpiTox = ifelse(peptide %in% epitox$peptide, "Yes", "No")
         ) %>%
  filter(!is.na(peptide),
          mismatch != 9) %>%
  select(-source, -Note) %>%
  distinct (peptide, .keep_all = T) %>%
  arrange(mismatch)
rownames(patents_df) = 1:nrow(patents_df)

patents_df = merge(patents_df, epitox[, c("peptide", "Wildtype")], by = "peptide", all.x =T) %>%
  arrange(mismatch)%>%
  distinct (peptide, .keep_all = T)

biocopy_colors = c("#A2C510", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")


# just for the plots
#results = merge(results, cutoff_4[, c("id", "presentation_score")], by.x = "peptide_id", by.y = "id", all.x = T)
PieChart(EpiTox, data = patents_df, hole = 0,
         fill = c("white", "#958BB2"),
         color="white",
         main = paste0("Disclosed peptides: ", nrow(patents_df)),
         labels = "input",
         labels_position= "in",
         labels_size=1.5,
         labels_color = c("black", "white"),
         labels_cex = 1.5)

ggplot(epitox, aes(edit_distance, fill = Wildtype)) +
  geom_bar(position = "dodge") +
  theme_light() +
  scale_fill_manual(values = c("Yes" = "#A2C510", "No" = "#2C4255")) +
  scale_y_log10() +
  coord_flip() +
  theme(legend.position = "bottom")

