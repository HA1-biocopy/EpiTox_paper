# Bayesian Peptide Evidence Assessment
# ============================================
#
# This script assesses off-target peptides using Bayesian probability updating
# Goal: "How critical is the off-target based on experimental evidence?"
#
# Author: Hoor. Al-Hasani
# Date: 2025
# ----------------
library(lessR)
library(dplyr)
library(tidyr)

# ============================================
# # STEP 1: Load data
# loading IEDB annotation
# ---------------------------
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
         uniprot = NULL, peptide = NULL) %>%
  relocate(id)


# loading EpiTox's peptide list
# --------------------------------
cutoff_4 = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  select(ensembl_gene_id, uniprot, Gene.Names, peptide, blosum_similarity, mismatch, Wildtype, Peptide_HLA_Atlas) %>%
  distinct(.keep_all = T) %>%
  mutate(id = paste0(uniprot, "_", peptide)) %>%
  #merge(., experimental_df, by = "id", all.x =T) %>%
  distinct(.keep_all = T)
# ============================================


# ============================================
# # STEP 2. Load all scripts
source("utils_exp.R")
source("Bayesian_preprocess_experimental_data.R")
source("Bayesian_peptide_assessment.R")

# # STEP 2. Preprocess
preprocessed <- preprocess_for_bayesian(cutoff_4, experimental_df, target_allele = "HLA-A*02:01",
                                        target_disease = "bladder", target_class = "I")

#   target_disease = "melanoma",
#   target_class = "I"
# ============================================
# # STEP 3. Assess
biocopy_colors = c("#A2C510", "#99CFE9", "#FBB800", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")

results <- bayesian_peptide_assessment(preprocessed, target_allele, target_class)

PieChart(confidence_level, data = results, hole = 0,
         fill = biocopy_colors,
         color="white",
         main = paste0("Total peptides: ", nrow(results)),
         labels_cex = 0.6)

# ============================================
# # STEP 4. Report
DIR = "/Users/hoor.alhasani/Documents/Projects/D003/Patent_Paper/paper_materials/02_Data_Analysis/data/Bayesian/MAGEA3"
generate_assessment_report(results,
                           output_prefix = DIR)
# ============================================
#


# =============================================================================
# COMPLETE WORKFLOW EXAMPLE
# =============================================================================

# # STEP 1: Load all scripts
# source("utils_exp.R")
# source("preprocess_experimental_data.R")
# source("bayesian_peptide_assessment.R")
#
# # STEP 2: Load your data
# peptides_df <- read.csv("your_peptide_predictions.csv")
# experimental_df <- read.csv("your_iedb_data.csv")
#
# # STEP 3: Preprocess
# preprocessed_data <- preprocess_for_bayesian(
#   peptides_df = peptides_df,
#   experimental_df = experimental_df,
#   target_allele = "HLA-A*02:01",
#   target_disease = "melanoma",
#   target_class = "I"
# )
#
# # STEP 4: Run Bayesian assessment
# results <- bayesian_peptide_assessment(
#   peptide_data = preprocessed_data,
#   target_allele = "HLA-A*02:01",
#   target_class = "I",
#   prior_prob = 0.50
# )
#
# # STEP 5: Generate report
# report <- generate_assessment_report(
#   assessed_peptides = results,
#   output_prefix = "my_peptide_assessment",
#   include_plots = TRUE
# )
#
# # STEP 6: View top peptides
# top_peptides <- results %>%
#   arrange(desc(posterior_prob)) %>%
#   select(peptide_id, posterior_prob, confidence_level,
#          interpretation, evidence_chain, has_target_allele, n_studies) %>%
#   head(20)
#
# print(top_peptides)
