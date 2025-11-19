library(ggplot2)
library(dplyr)
library(ggpubr)  # for stat_compare_means

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
         uniprot = NULL, peptide = NULL,
         is_normal = case_when(
           disease == "healthy" ~ TRUE,
           TRUE ~ FALSE
         )) %>%
  relocate(id)

results <- read.xlsx("../../data/Bayesian/MAGEA3_full_results.xlsx")

# Calculate counts
counts <- results %>%
  group_by(confidence_level) %>%
  summarise(n = n(), .groups = 'drop')

# Define the correct order for confidence levels
level_order <- c("High", "Medium", "Low", "Very Low")

# Create labels with proper ordering
df <- results %>%
  left_join(counts, by = "confidence_level") %>%
  mutate(confidence_level = factor(confidence_level, levels = level_order),
         confidence_label = paste0(confidence_level, "\n(n=", n, ")"),
         confidence_label = factor(confidence_label,
                                   levels = paste0(level_order, "\n(n=",
                                                   counts$n[match(level_order, counts$confidence_level)], ")")))


# Create the plot with statistical tests
ggplot(df, aes(x = confidence_label, y = affinity)) +
  geom_boxplot(aes(fill = confidence_label)) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(results$affinity) * 1.15,
                     label.x = 1.5) +  # Overall test
  stat_compare_means(comparisons = list(c(1, 2), c(2,4), c(1,4)),
                     method = "wilcox.test",
                     label = "p.signif") +  # Pairwise comparisons
  theme_light() +
  labs(x = "Evidence level",
       y = "log10(Affinity)",
       title = "Target HLA-allele predicted affinity vs. Bayesian classification",
       caption = "Significance was assessed with wilcoxon test") +
  scale_fill_manual(values = c("#00508A", "#B0DAEE", "#F4FAFD", "#F4FAFD")) +
  theme(legend.position = "none") +
  scale_y_log10()

ggplot(experimental_df %>%
         filter(), aes(disease, fill = Binding)) +
  geom_bar(stat = "count", position = "dodge") +
  coord_flip()
