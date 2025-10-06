
#devtools::install_github("antuneslab/crossdome", build_vignettes = TRUE)
library(crossdome)
library(ggplot2)
library(dplyr)

EpiTox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_5/Tables/HPA_genes_nTPM.xlsx") %>%
  distinct(chr, uniprot, Gene.Names, peptide, mismatch, Wildtype, .keep_all = T) %>%
  arrange(mismatch)



database <- cross_background(off_targets = 'KVAELVHFL', allele = "HLA-A*02:01")
result <- cross_compose(query = 'KVAELVHFL', background = database)
crossdome_rslts = result@result

rs <- cross_expression_matrix(result, pvalue_threshold = 0.005)
cross_expression_plot(object = rs)

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
  mutate(CrossDome = ifelse(peptide %in% crossdome_rslts$subject, "Yes", "No"),
         EpiTox = ifelse(peptide %in% EpiTox$peptide, "Yes", "No")) %>%
  filter(!is.na(peptide)) %>%
  select(-source, -Note) %>%
  distinct (.keep_all = T) %>%
  arrange(mismatch)
rownames(patents_df) = 1:nrow(patents_df)

patents_df = merge(patents_df, crossdome_rslts[, c("subject", "n_mismatch", "relatedness_score")], by.x = "peptide", by.y = "subject", all.x =T) %>%
  arrange(mismatch)


# EpiTox vs. CrossDome

df = crossdome_rslts %>%
  select(subject, n_mismatch) %>%
  dplyr::rename(mismatch = n_mismatch,
                peptide = subject) %>%
  filter(mismatch < 6) %>%
  mutate(Tool = "CrossDome") %>%
  distinct(.keep_all = T)


df = EpiTox %>%
  filter(Rank == "High") %>%
  select(peptide, mismatch) %>%
  distinct(.keep_all = T) %>%
  mutate(Tool = "EpiTox") %>%
  rbind(., df)


ggplot(df, aes(mismatch, fill = Tool)) +
  geom_bar(position=position_dodge(), col = "#2C4255") +
  scale_fill_manual(values = c("#FBB800", "#A2C510", "#99CFE9", "#2C4255")) +
  theme_light() +
  scale_y_log10() +
  theme(legend.position = "bottom") +
  labs(y = "# peptides (log10)")



# venn for the patent
venn_list = list(EpiTox = unique(patents_df$peptide[patents_df$EpiTox == "Yes"]),
                 CrossDome = unique(patents_df$peptide[patents_df$CrossDome == "Yes"]),
                 Patent = unique(patents_df$peptide))

ggvenn::ggvenn(venn_list, show_percentage = F,
               fill_color = c("#A2C510", "#FBB800", "#99CFE9", "#2C4255"), stroke_color = "white",
               set_name_size = 3)
