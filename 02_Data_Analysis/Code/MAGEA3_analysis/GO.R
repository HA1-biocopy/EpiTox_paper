# Install required packages (run once)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "AnnotationDbi"))
install.packages(c("ggplot2", "dplyr"))

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse, etc.
library(enrichplot)
library(ggplot2)
library(dplyr)

# Your protein/gene list (replace with your data)
# Can be gene symbols, Entrez IDs, Ensembl IDs, or UniProt IDs
gene_list <- peptide_data$uniprot[peptide_data$Rank == "High"] %>% unique()

# Convert to Entrez IDs if needed (from gene symbols)
gene_entrez <- bitr(gene_list,
                    fromType = "UNIPROT",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# GO Enrichment Analysis
go_enrichment <- enrichGO(gene = gene_entrez$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          ont = "ALL",  # "BP", "MF", "CC", or "ALL"
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)

# View results
head(go_enrichment, 10)
go_enrichment_df = as.data.frame(go_enrichment)

# Check which specific proteins are in the HDAC binding term
hdac_genes <- go_enrichment@result %>%
  filter(grepl("histone deacetylase", Description, ignore.case=TRUE)) %>%
  pull(geneID)

print(hdac_genes)  # See which of your 3 proteins are involved


# Export results
write.csv(as.data.frame(go_enrichment), "GO_enrichment_results.csv")

# KEGG Pathway Analysis
kegg_enrichment <- enrichKEGG(gene = gene_entrez$ENTREZID,
                              organism = 'hsa')

# Convert KEGG IDs to gene symbols
kegg_enrichment <- setReadable(kegg_enrichment, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_enrichment_df = as.data.frame(kegg_enrichment)

write.csv(as.data.frame(kegg_enrichment), "KEGG_enrichment_results.csv")

# Visualizations
# 1. Dot plot
dotplot(go_enrichment, showCategory=15)
ggsave("GO_dotplot.pdf", width=10, height=8)

# 2. Bar plot
barplot(go_enrichment, showCategory=15)
ggsave("GO_barplot.pdf", width=10, height=8)

# 3. Gene-Concept Network
cnetplot(go_enrichment, categorySize="pvalue",
         foldChange=NULL, showCategory = 7)
ggsave("GO_cnetplot.pdf", width=12, height=10)

# 4. Enrichment Map (shows relationships between terms)
emapplot(pairwise_termsim(go_enrichment), showCategory=20)
ggsave("GO_emapplot.pdf", width=12, height=10)
