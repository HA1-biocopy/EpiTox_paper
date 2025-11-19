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
  mutate(Gene_name = gsub("[[:punct:]]", "", gene),
         gene = NULL,
         Note = NULL,
         peptide = ifelse(peptide == "TVAELVQFV", "TVAELVQFL", peptide)) %>%
  filter(!is.na(peptide), !peptide %in% c("TVAELVQFV", "FLWGPRALV")) %>%
  distinct (.keep_all = T) %>%
  group_by(peptide) %>%
  summarise(across(everything(), ~ paste(unique(.x), collapse = ", "))) %>%
  relocate(Gene_name) %>%
  arrange(mismatch, Gene_name)
rownames(patents_df) = 1:nrow(patents_df)

openxlsx::write.xlsx(patents_df, "../../data/patents_supplemantry_table3.xlsx")
