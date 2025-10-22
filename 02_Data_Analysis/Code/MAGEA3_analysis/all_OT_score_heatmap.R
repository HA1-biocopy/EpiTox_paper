library(dplyr)
library(ggplot2)

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

epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  distinct(uniprot, peptide, mismatch, .keep_all = T) %>%
  relocate(Rank, .after = "peptide") %>%
  mutate(id = paste0(uniprot, "_", peptide),
         is_healthy = ifelse(id %in% experimental_df$id[experimental_df$is_normal], "Yes", "No")) %>%
  arrange(mismatch, Rank)




epitox$Peptide_IEDB %>% table
epitox$is_healthy %>% table
epitox[, c("Peptide_HLA_Atlas", "is_healthy")] %>% table
epitox[, c("Peptide_HLA_Atlas", "is_healthy", "IEDB_allele")] %>% table

score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  filter (peptide %in% epitox$peptide) %>% #, Outcome == "Binder") %>%
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  merge(., epitox[,c("uniprot", "peptide", "Wildtype", "blosum_similarity", "Rank")], by = "peptide", all.x = T) %>%
  arrange(lev_dist)


expr_mat <- as.matrix(epitox[, 27:ncol(epitox)])
rownames(expr_mat) <- paste(epitox$uniprot, epitox$peptide, sep = "_")
row_annot <- data.frame(Wildtype = epitox$Wildtype,
                        Rank =epitox$Rank,
                        SCORE = epitox$peptide %in% score_reference$peptide,
                        AB_binder = epitox$peptide %in% score_reference$peptide[score_reference$Outcome == "Binder"]) %>%
  mutate(SCORE = ifelse(SCORE, "Yes", "No"),
         AB_binder = ifelse(AB_binder, "Yes", "No"))

rownames(row_annot) <- paste(epitox$uniprot, epitox$peptide, sep = "_")


# Compared to target
target_gxp = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  filter(peptide == "KVAELVHFL") %>%
  distinct(peptide, .keep_all = T)
target_gxp = target_gxp[, 27:ncol(target_gxp)]
rownames(target_gxp) = "P43362_KVAELVHFL"

target_gxp_df = t(target_gxp)  %>% as.data.frame()
target_gxp_df$Tissue = rownames(target_gxp_df)
rownames(target_gxp_df) = 1:nrow(target_gxp_df)

# if target_gxp is a one-row data frame
df <- sweep(log2(expr_mat+1), 2, log2(as.numeric(target_gxp[1, ])+1), "-")
mat = as.matrix(df)

ann_colors <- list(
  Wildtype = c(
    Yes = "#A2C510",    # green
    No = "#99CFE9"  # red
  ),
  SCORE = c(
    Yes = "#2C4255",    # green
    No = "white"  # red
  ),
  AB_binder = c(
    Yes = "#2C4255",    # green
    No = "white"  # red
  ),
  Rank= c(
    High = "#FBB800",
    Low = "white",
    Moderate = "#FEF1CC",
    Random = "white"
  )
)

# group_order <- order(row_annot$Rank)
# mat_ordered <- mat[group_order, ]
# row_annot_ordered <- row_annot[group_order, ]


# gaps <- cumsum(table(row_annot_ordered$Rank))
# pheatmap(mat,
#          cluster_rows = FALSE,    # disable clustering if you want groups preserved
#          cluster_cols = TRUE,
#          annotation_row = row_annot_ordered,
#          gaps_row = gaps,         # visually separate groups
#          annotation_colors = ann_colors,
#          color = colorRampPalette(c("white", "red"))(100),
#          breaks = seq(0, 10, length.out = 101),
#          legend_breaks = c(0, 1, 2, 3, 5, 10),
#          legend_labels = c("same", "2x", "4x", "8x", "32x", "1024x"),
#          legend_title = "log2 ratio",
#          na_col = "grey")


library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(c(min(mat), 0, max(mat)),
                      c("#6B7B88", "white", "#A2C510"))

cht = Heatmap(t(mat),
        name = "log2(OT) - log2(target)",
        cluster_rows = F,
        cluster_columns = TRUE,
        col = col_fun,
        heatmap_legend_param = list(
          direction = "horizontal",   # ← makes the legend bar horizontal
          legend_width = unit(4, "cm")), # optional, adjust width
        top_annotation = HeatmapAnnotation(df = row_annot, col = ann_colors),
        column_title_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 1),
        column_split = row_annot$Rank, # ← separates by annotation
        row_names_gp = gpar(fontsize = 7),
)

draw(cht,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legends = TRUE)










