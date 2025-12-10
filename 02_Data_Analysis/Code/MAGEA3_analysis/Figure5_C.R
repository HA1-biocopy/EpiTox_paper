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
  mutate(Experimental_Evidence = ifelse(Experimental_Evidence == "Medium", "Moderate", Experimental_Evidence),
         genes = sub("\\s.*", "", Gene.Names)) %>%
  filter(Outcome == "Binder") %>%
  arrange(anchor_status, Wildtype, Bi_feature_rank,
          Multi_feature_rank, Experimental_Evidence, log2_target_FC)

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
    #Known_MAGEA3 = gxp$known_peptide,
    anchor_status = gxp$anchor_status,
    WT = gxp$Wildtype,
    #Mismatch = anno_barplot(gxp$mismatch, gp = gpar(col = "white", fill = "#8ECAE6"), height = unit(1, "cm")),
    #BLOSUM62 = anno_barplot(gxp$blosum_similarity,  gp = gpar( fill = "#EAD7A4"), height = unit(1, "cm")),
    log2_target_KD_FC = anno_barplot(gxp$log2_target_FC,  gp = gpar( fill = "#6696B9"), height = unit(1.3, "cm")),
    Total = anno_barplot(
      gxp_colSums,
      bar_width = 1, gp = gpar(fill = "#338232", lwd = 0.2, col ="black"),
      height = unit(1.3, "cm")
    ),
    col = list(
      anchor_status = c("anchor_mismatch" = "#F08000", "Intact" = "#438D99",
                        "Both Anchor & Backbone" = "#EAD7A4"),
      WT = c("No" = "#2C4255", "Yes" = "#F0F0F0")#,
      #Known_MAGEA3 = c("Known" = "#958BB2", "unkwon" = "#DCD8E5")
    ),
    gap = unit(c(1,1,2,2,0,0,2), "mm"),
    border = c(Bi_feature_rank =T, Multi_feature_rank =T, Experimental_Evidence =T)
  )
  return (top_ha)
}

bottom_annotation <- function(gxp = NA,  to_heatmap = NA, scale = F){
  top_ha = HeatmapAnnotation(
    Bi_feature_rank = gxp$Bi_feature_rank,
    Multi_feature_rank = gxp$Multi_feature_rank,
    Experimental_Evidence = gxp$Experimental_Evidence,
    col = list(
      Bi_feature_rank = c("High" = "#00508A", "Moderate" = "#B0DAEE", "Low" = "#F4FAFD"),
      Multi_feature_rank = c("High" = "#00508A", "Moderate" = "#B0DAEE", "Low" = "#F4FAFD"),
      Experimental_Evidence = c("High" = "#00508A", "Moderate" = "#B0DAEE", "Very Low" = "#F4FAFD")
    ),
    #gap = unit(c(1,0,2,2,0,0,2), "mm"),
    border = c(Bi_feature_rank =T, Multi_feature_rank =T, Experimental_Evidence =T)
  )
  return (top_ha)
}


tissues_prot = data.frame(tissue = colnames(epitox[16:(ncol(epitox)-1)])) %>%
  mutate(category = case_when(
    tissue %in% c("adipose.tissue", "skin") ~ "Connective / Integumentary",
    tissue %in% c("adrenal.gland", "thyroid.gland", "parathyroid.gland", "pancreas") ~ "Endocrine",
    tissue %in% c("heart.muscle", "skeletal.muscle", "smooth.muscle") ~ "Muscular",
    tissue %in% c("cerebral.cortex", "choroid.plexus") ~ "Nervous",
    tissue %in% c("lung") ~ "Respiratory",
    tissue %in% c("liver", "gallbladder") ~ "Hepatic / Biliary",
    tissue %in% c("kidney", "urinary.bladder") ~ "Urinary",
    tissue %in% c("spleen", "thymus", "lymph.node", "tonsil","appendix", "bone.marrow") ~ "Immune", #"Immune / Lymphatic",
    tissue %in% c("colon", "duodenum", "esophagus", "rectum", "salivary.gland", "small.intestine", "stomach", "tongue") ~ "Digestive",
    tissue %in% c("breast", "cervix", "endometrium", "epididymis", "fallopian.tube", "ovary", "placenta", "prostate", "seminal.vesicle", "testis") ~ "Reproductive",
    #tissue %in% c("appendix", "bone.marrow") ~ "Immune / Hematopoietic",
    TRUE ~ "Other"
  ))



expr_mat = epitox[, 16:(ncol(epitox)-1)]
# log2( (expr_my_gene + 1) / (expr_target + 1) ) — i.e. log₂ fold change with a +1 pseudocount
df <- sweep(log2(expr_mat+1), 2, log2(as.numeric(target[2, 16:(ncol(target)-1)])+1), "-") # took the isoform with highest gxp
mat = as.matrix(df)
rownames(mat) = paste0(epitox$genes, "_", gsub(".+\\_", "", epitox$id))
peptides = gsub(".+\\_", "", epitox$id)
# Define your custom breakpoints for better granularity above zero
col_prot_fun = colorRamp2(
  c(min(mat, na.rm = TRUE), 0, 1, 2, 4, 6, max(mat, na.rm = TRUE)),
  c("#64686A", "#F6F9E7", "#D6E6D6", "#ADCDAD", "#85B484", "#5C9B5B", "#338232")
)
prot_top_ha = top_annotation(gxp = epitox, to_heatmap = mat)
prot_bottom_ha = bottom_annotation(gxp = epitox, to_heatmap = mat)

hpa_ch_plot = mat %>%
  t() %>%
  Heatmap(name = "Gene expression (nTPM)",
          col = col_prot_fun,
          rect_gp = gpar(col = "black", lwd = 0.2),
          row_gap = unit(2, "mm"), border = TRUE,
          row_split = tissues_prot$category,
          row_title_rot = 0,
          row_title_gp = gpar(fontsize = 9),
          row_names_gp = gpar(fontsize = 10),
          column_split = epitox$known_peptide,
          cluster_columns = F,
          column_names_gp = gpar(fontsize = 9),
          show_column_dend = F,
          show_row_dend = F,
          column_names_rot = 45,
          na_col = "gray",
          top_annotation = prot_top_ha,
          bottom_annotation = prot_bottom_ha,
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
draw(hpa_ch_plot, column_title = "Confirmed binders according to highSCORE",
           merge_legend = TRUE, heatmap_legend_side = "bottom",
           annotation_legend_side = "bottom")

pdf("../../data/Figure_6_C_Binders_GXP_nTPM_tissueSplit.pdf", width = 14, height = 12)
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




