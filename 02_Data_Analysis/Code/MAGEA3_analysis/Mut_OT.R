

cutoff_4 = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  filter(Wildtype == "No") %>%
  distinct(chr, uniprot, Gene.Names, peptide, mismatch, Wildtype, .keep_all = T) %>%
  arrange(mismatch)

epitox = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_5/Tables/PrediTopes_MAGEA3.xlsx") %>%
  filter(Wildtype == "No") %>%
  distinct(chr, uniprot, Gene.Names, peptide, mismatch, Wildtype, .keep_all = T) %>%
  merge(., cutoff_4[,c("peptide", "Rank")], by = "peptide", all.x = T) %>%
  relocate(Rank, .after = "peptide") %>%
  arrange(mismatch, Rank)

score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  filter (peptide %in% epitox$peptide, Outcome == "Binder") %>%
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  merge(., cutoff_4[,c("peptide", "blosum_similarity", "Rank")], by = "peptide", all.x = T) %>%
  arrange(lev_dist)


# Quick calculation to prioritize:
#For each off-target, estimate: (Expression level) / KD
# This gives you a relative "sink" effect - which off-target will actually compete for your antibody in patients.

ggplot(cutoff_4 %>% filter(affinity <600)) +
  geom_point(aes(blosum_similarity, affinity, col = Rank)) + theme_light()

ggplot(cutoff_4) +
  geom_point(aes(blosum_similarity, affinity, col = Rank)) +
  theme_light()

mut_ot = cutoff_4 %>%
  filter(peptide %in% score_reference$peptide)
mut_ot_df = mut_ot[, 27:ncol(mut_ot)]
rownames(mut_ot_df) = paste0(mut_ot$uniprot, "_", mut_ot$peptide)

library(pheatmap)

# Convert to matrix
mat <- as.matrix(mut_ot_df)

# Heatmap
pheatmap(mat,
         scale = "row",          # scale each gene across tissues
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE)

# Compared to target
target_gxp = cutoff_5 %>%
  filter(peptide == "KVAELVHFL") %>%
  distinct(peptide, .keep_all = T)
target_gxp = target_gxp[, 27:ncol(target_gxp)]
rownames(target_gxp) = "P43362_KVAELVHFL"

target_gxp_df = t(target_gxp)  %>% as.data.frame()
target_gxp_df$Tissue = rownames(target_gxp_df)
rownames(target_gxp_df) = 1:nrow(target_gxp_df)

# if target_gxp is a one-row data frame
df <- sweep(log2(mut_ot_df+1), 2, log2(as.numeric(target_gxp[1, ])+1), "-")
mat = as.matrix(df)

# Heatmap
pheatmap(df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c( "white", "red"))(100),
         breaks = seq(0, 10, length.out = 101),
         #legend_breaks = c(0, 1, 2, 3, 4, 5),
         legend_labels = c("same", "2x", "4x", "8x"),
         legend_title = "log2 ratio",
         na_col = "grey")

pheatmap(df,
         cluster_rows = T,
         cluster_cols = T,
         color = colorRampPalette(c("white", "red"))(100),
         breaks = seq(0, 10, length.out = 101),
         legend_breaks = c(0, 1, 2, 3, 5, 10),
         legend_labels = c("same", "2x", "4x", "8x", "32x", "1024x"),
         legend_title = "log2 ratio",
         na_col = "grey")

##########################################################
peptides = c("FQAELVHPA", "IVAYLVHYV") #unique(score_reference$peptide)

df = openxlsx::read.xlsx("~/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  filter(uniprot %in% epitox$uniprot[epitox$peptide %in% peptides]) %>%
  distinct(chr, uniprot, Gene.Names, peptide, mismatch, Wildtype, .keep_all = T)

score_reference = openxlsx::read.xlsx("/Volumes/lab/03_HighSCORE/P048/P048_01_Results Triplicates/TTP-230_new/PPB-68_kinetic_table.xlsx") %>%
  dplyr::rename(peptide = peptide_sequence) %>%
  filter (peptide %in% df$peptide) %>% #, Outcome == "Binder") %>%
  select(labguru_id, peptide, KD, n, lev_dist, Outcome, log2_target_FC, n_slides_QC_passed) %>%
  merge(., df[,c("uniprot", "peptide", "Wildtype", "blosum_similarity", "Rank")], by = "peptide", all.x = T) %>%
  arrange(lev_dist)

expr_mat <- as.matrix(df[, 27:ncol(df)])
rownames(expr_mat) <- paste(df$uniprot, df$peptide, sep = "_")
row_annot <- data.frame(Wildtype = df$Wildtype, Rank =df$Rank,
                        SCORE = df$peptide %in% score_reference$peptide)
row_annot$SCORE = ifelse(row_annot$SCORE, "Yes", "No")
rownames(row_annot) <- paste(df$uniprot, df$peptide, sep = "_")


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
    Yes = "#A2C510",    # green
    No = "white"  # red
  ),
  Rank= c(
    High = "#FBB800",
    Low = "white",
    Moderate = "white"
  )
)

pheatmap(
  mat,
  annotation_row = row_annot,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_colors = ann_colors,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("white", "firebrick3"))(100),
  breaks = seq(0, 5, length.out = 101),
  fontsize_row = 6
)
pheatmap(mat,
         cluster_rows = T,
         cluster_cols = T,
         annotation_row = row_annot,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("white", "red"))(100),
         breaks = seq(0, 10, length.out = 101),
         legend_breaks = c(0, 1, 2, 3, 5, 10),
         legend_labels = c("same", "2x", "4x", "8x", "32x", "1024x"),
         legend_title = "log2 ratio",
         na_col = "grey")




