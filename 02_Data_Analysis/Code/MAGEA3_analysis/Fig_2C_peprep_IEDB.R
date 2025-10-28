library(dplyr)
library(ggplot2)

biocopy_colors = c("#A2C510", "#958BB2", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")

#############################################################################
# panel one: peptides proportion in peprep and iedb
#############################################################################


peprep = read.csv("~/Documents/Projects/PEPREP/data/PEPREP_peptides_aggregation.csv")

epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  distinct(.keep_all = T) %>%
  mutate(PEPREP = ifelse(peptide %in% peprep$peptide, T, F),
         IEDB = ifelse(Peptide_IEDB == "Yes", T, F),
         category = case_when(
           PEPREP & IEDB ~ "PEPREP + IEDB",
           PEPREP ~"PEPREP",
           IEDB ~ "IEDB",
           TRUE ~ "None"),
         category = factor(category, levels = c("None", "PEPREP", "PEPREP + IEDB", "IEDB")),
         id = paste0(uniprot, "_", peptide)) %>%
  distinct(peptide, uniprot, .keep_all = T)


# Count and calculate proportions
counts <- epitox %>%
  count(category) %>%
  mutate(prop = n / sum(n)) %>%
  #arrange(match(category, c("None", "IEDB", "PEPREP + IEDB", "PEPREP"))) %>%
  mutate(label = paste0(n, " (", round(prop * 100, 1), "%)"),
         ypos = 1- (cumsum(prop) - prop / 2))


# Plot
ggplot(counts, aes(x = 1, y = prop, fill = category)) +
  geom_bar(stat = "identity", width = 0.4, color = "white") +
  geom_text(aes(x = 1, y = ypos, label = label),
            size = 3.3, fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Proportion of peptides found in MHC databases") +
  #coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.justification = c(0.7,0),
    #legend.box.just = "left",
    #legend.key.size = unit(0.3, "cm"),      # Smaller legend boxes
    legend.box.spacing = unit(0.01, "cm") ,
    #legend.text = element_text(size = 7),    # Smaller legend text
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = biocopy_colors) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


ggplot(counts, aes(x = 1, y = prop, fill = category)) +
  geom_bar(stat = "identity", width = 0.4, color = "white") +
  geom_text(aes(x = 1, y = ypos, label = label),
            size = 3.3, fontface = "bold") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = NULL) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    #legend.justification = c(0.7,0),
    #legend.box.just = "left",
    #legend.key.size = unit(0.3, "cm"),      # Smaller legend boxes
    legend.box.spacing = unit(0.01, "cm") ,
    #legend.text = element_text(size = 7),    # Smaller legend text
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = biocopy_colors)

epitox %>% filter( PEPREP == T & IEDB == F) %>% nrow

#############################################################################
# panel two: peptides distribution accross tissues including normal
#############################################################################

files = list.files(path = "~/Documents/Projects/MAGEA3/results/Cutoff_4/PrediTopes/DB_annotation/",
                   pattern = "*sequence.csv", full.names = T)
experimental_df = lapply(files, function(f){
  if (file.size(f) == 0) return(NULL)
  first_line <- readLines(f, n = 1, warn = FALSE)
  if (length(first_line) == 0 || all(trimws(first_line) == "")) return(NULL)

  x = read.csv(f) %>%
    select(parent_source_antigen_name, curated_source_antigen.accession, linear_sequence, reference_id,
           qualitative_measure, mhc_class, mhc_allele_name,
           assay_names, disease_names) %>%
    dplyr::rename(
      peptide = linear_sequence,
      study_ref = reference_id,
      disease = disease_names,
      Experimental_method = assay_names,
      Binding = qualitative_measure,
      hla_class = mhc_class,
      hla_allele = mhc_allele_name)
  return(x)

}) %>%
  do.call("rbind",.) %>%
  mutate(uniprot = gsub("\\.\\d+", "", curated_source_antigen.accession),
         disease = gsub("\\[\\'(.+)\\'\\]", "\\1", disease),
         curated_source_antigen.accession = NULL,
         id = paste0(uniprot, "_", peptide),
         uniprot = NULL,
         is_normal = case_when(
           disease == "healthy" ~ TRUE,
           TRUE ~ FALSE
         )) %>%
  relocate(id)







#############################################################################
