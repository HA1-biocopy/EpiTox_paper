library(dplyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(ggseqlogo)


# prepare data and annotation
df = openxlsx::read.xlsx("/Users/hoor.alhasani/Documents/Projects/MAGEA3/results/Cutoff_4/Table/HPA_genes_nTPM.xlsx") %>%
  mutate(id = paste0(uniprot, "_", peptide)) %>%
  relocate(id)

#x = openxlsx::read.xlsx("../../data/OT_SCORE_full_annotation.xlsx")
COLS = colnames(df)[c(1:3,7,28:ncol(df))]
epitox = openxlsx::read.xlsx("../../data/OT_SCORE_full_annotation.xlsx") %>%
  select(id, Wildtype, mismatch, blosum_similarity, anchor_status, known_peptide,
         confidence_level, Bi_feature_rank, multibiophys_category,
         Tested, Outcome, log2_target_FC) %>%
  dplyr::rename(Multi_feature_rank = multibiophys_category,
                Experimental_Evidence = confidence_level) %>%
  merge(., df[,COLS], by = "id", all.x =T) %>%
  mutate(Experimental_Evidence = ifelse(Experimental_Evidence == "Medium", "Moderate", Experimental_Evidence)) %>%
  filter(Outcome == "Binder") %>%
  arrange(anchor_status, Wildtype, Bi_feature_rank, Multi_feature_rank,
          Experimental_Evidence, log2_target_FC)

target = epitox %>%
  filter (id == "P43362_KVAELVHFL")

epitox = epitox %>%
  filter(id != "P43362_KVAELVHFL")

# ============================================================================================================
biocopy_colors = c("#A2C510", "#FBB800", "#99CFE9", "#939597", "#C61E19", "#438D99", "#958BB2", "#6B7B88",
                   "#338232", "#F08000", "#3373A1", "#64686A", "#D14B47", "#98D0BC", "#4F3D7F", "#2C4255")
# ============================================================================================================
#source("complexHM_utlis.R")

top_annotation <- function(gxp = NA,  to_heatmap = NA, scale = F){

  if(scale){
    gxp_colSums = colSums(scale(t(to_heatmap)), na.rm = T)
  }else{
    gxp_colSums = colSums(t(to_heatmap), na.rm = T)
  }
  top_ha = HeatmapAnnotation(
    anchor_status = gxp$anchor_status,
    WT = gxp$Wildtype,
    Known_OT = gxp$known_peptide,
    #Mismatch = anno_barplot(gxp$mismatch, gp = gpar(col = "white", fill = "#8ECAE6"), height = unit(1, "cm")),
    #BLOSUM62 = anno_barplot(gxp$blosum_similarity,  gp = gpar( fill = "#EAD7A4"), height = unit(1, "cm")),
    log2_target_KD_FC = anno_barplot(gxp$log2_target_FC,  gp = gpar( fill = "#EAD7A4"), height = unit(1, "cm")),
    Bi_feature_rank = gxp$Bi_feature_rank,
    Multi_feature_rank = gxp$Multi_feature_rank,
    Experimental_Evidence = gxp$Experimental_Evidence,
    gxp = anno_barplot(
      gxp_colSums,
      bar_width = 1, gp = gpar(fill = "#438D99", lwd = 0.2, col ="black"),
      height = unit(1, "cm")
    ),
    col = list(
      anchor_status = c("anchor_mismatch" = "#338232", "Intact" = "#4F3D7F",
                        "Both Anchor & Backbone" = "#FBB800"),
      WT = c("No" = "#2C4255", "Yes" = "#EAECEE"),
      Known_OT = c("Known" = "#2C4255", "unkwon" = "#EAECEE"),
      Bi_feature_rank = c("High" = "#00508A", "Moderate" = "#99B9D0", "Low" = "#EAECEE"),
      Multi_feature_rank = c("High" = "#00508A", "Moderate" = "#99B9D0", "Low" = "#EAECEE"),
      Experimental_Evidence = c("High" = "#00508A", "Moderate" = "#99B9D0", "Very Low" = "#EAECEE")
    ),
    gap = unit(c(1,0,2,2,0,0,2), "mm"),
    border = c(Bi_feature_rank =T, Multi_feature_rank =T, Experimental_Evidence =T)
  )
  return (top_ha)
}






expr_mat = epitox[, 16:ncol(epitox)]
df <- sweep(log2(expr_mat+1), 2, log2(as.numeric(target[2, 16:ncol(epitox)])+1), "-") # took the isoform with highest gxp
mat = as.matrix(df)
rownames(mat) = epitox$id
# Define your custom breakpoints for better granularity above zero
col_prot_fun = colorRamp2(
  c(min(mat, na.rm = TRUE), 0, 1, 2, 4, 6, 8, max(mat, na.rm = TRUE)),
  c("white", "#CCDCE8", "#F6F9E7", "#ECF3CF", "#DAE89F", "#C7DC70", "#B5D140", "#A2C510")
)
prot_top_ha = top_annotation(gxp = epitox, to_heatmap = mat)

hpa_ch_plot = mat %>%
  t() %>%
  Heatmap(name = "Gene expression (nTPM)",
          col = col_prot_fun,
          rect_gp = gpar(col = "white", lwd = 0.2),
          row_gap = unit(2, "mm"), border = TRUE,
          #row_split = tissues_prot$Organ,
          row_title_rot = 75,
          row_title_gp = gpar(fontsize = 9),
          row_names_gp = gpar(fontsize = 9),
          column_split = epitox$known_peptide,
          cluster_columns = F,
          column_names_gp = gpar(fontsize = 6),
          show_column_dend = F,
          column_names_rot = 45,
          na_col = "gray",
          top_annotation = prot_top_ha,
          #left_annotation = row_ha,
          heatmap_legend_param = list(
            direction = "horizontal",
            legend_width = unit(15, "cm"),
            at = c(min(mat, na.rm = TRUE), 0, 1, 2, 4, 6, 8, max(mat, na.rm = TRUE)),
            labels = c(round(min(mat, na.rm = TRUE)), "0", "1", "2", "4", "6", "8",
                       round(max(mat, na.rm = TRUE))),
            title = "log2 FC"
          ),
  )

# generate heatmap
pdf("../../data/Figure_4_C_Binders_GXP_nTPM.pdf", width = 14, height = 10)
draw(hpa_ch_plot, column_title = "Confirmed binders according to highSCORE",
     merge_legend = TRUE, heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()



# Figure D

peptides = gsub(".+\\_", "", epitox$id)
ggseqlogo(peptides)

non_binder = x$peptide[which(x$Outcome == "Non-Binder" & x$Tested == "Yes")]
ggseqlogo(non_binder, ncol =1, method = "prob") +
  scale_fill_manual(values = c("Polar" = "#FBB800", "Neutral" = "#99CFE9", "Acidic" = "#C61E19", "Basic" = "#A2C510", "Hydrophobic" = "#438D99"))




