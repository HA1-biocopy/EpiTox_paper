top_annotation <- function(gxp = NA,  to_heatmap = NA, scale = F){

  if(scale){
    gxp_colSums = colSums(scale(t(to_heatmap)), na.rm = T)
  }else{
    gxp_colSums = colSums(t(to_heatmap), na.rm = T)
  }
  top_ha = HeatmapAnnotation(
                             anchor_status = gxp$anchor_status,
                             #WT = gxp$Wildtype,
                             Mismatch = anno_barplot(gxp$mismatch, gp = gpar(col = "white", fill = "#FBB800"), height = unit(1, "cm")),
                             BLOSUM62 = anno_barplot(gxp$blosum_similarity,  gp = gpar( fill = "#EAD7A4"), height = unit(2, "cm")),
                             log2_target_KD_FC = anno_barplot(gxp$log2_target_FC,  gp = gpar( fill = "#EAD7A4"), height = unit(1, "cm")),
                             Known_OT = gxp$known_peptide,
                             Bi_feature_rank = gxp$Bi_feature_rank,
                             Multi_feature_rank = gxp$Multi_feature_rank,
                             Experimental_Evidence = gxp$Experimental_Evidence,
                             gxp = anno_barplot(
                               gxp_colSums,
                               bar_width = 1, gp = gpar(fill = "#438D99", lwd = 0.2),
                               height = unit(1.5, "cm")
                             ),
                             # col = list(anchor_status = c("High" = "#C61E19", "Medium" = "#D14B47", "Low" = "#E8A5A3", "Random" = "#F4D2D1"),
                             #            WT = c("Yes" = "#2C4255", "No" = "#EAECEE"),
                             #            Known_OT = c("Yes" = "#FBB800", "No" = "#FEF1CC"),
                             #            Bi_feature_rank = c("Yes" = "#2C4255", "No" = "#EAECEE"),
                             #            Multi_feature_rank = c("Yes" = "#FBB800", "No" = "#FEF1CC"),
                             #            Experimental_Evidence = c("Yes" = "#2C4255", "No" = "#8ECAE6")
                             # )
                             #border = c(KD = TRUE, Mismatch = FALSE, IEDB = T, IEDB_allele = T, HLA_Atlas = T, HLA_Atlas_allele =T, WT =F, Binder =T, Presented =T),
                             #gap = unit(c(1,2,2,2,2,0,2,0,2,0,2), "mm"),
                             # annotation_legend_param = list(
                             #   Binder = list(at = c("High", "Medium", "Low", "Random"), nrow=1, direction = "horizontal", title_position = "lefttop"),
                             #   Presented = list(at = c("Very high", "High", "Random"), nrow=1, direction = "horizontal", title_position = "lefttop"),
                             #   HLA_Atlas = list(title = "Peptide in DBs", nrow=1, direction = "horizontal", title_position = "lefttop"),
                             #   IEDB_allele = list(title = allele, nrow=1, direction = "horizontal", title_position = "lefttop")
                             # ),
                             # show_legend = c("IEDB" =FALSE, "HLA_Atlas_allele" = F),
                             # annotation_label = c("Wildtype","#Mismatch", "#Peptides", "#Replicates",
                             #                      "Affinity (molar)", "Affinity (cat)", "Presented (cat)",
                             #                      "in IEDB", allele, "in HLA-Atlas", allele,
                             #                      "GXP-colSums")
  )
  return (top_ha)
}
