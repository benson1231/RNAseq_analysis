# load library for package ----------------------------------------------
library(magrittr)
library(clusterProfiler) 
library(enrichplot)
library(ggplot2)      # we use ggplot2 to add x axis labels (ex: ridgeplot)
library(pathview)
library(cowplot)
library(DOSE)

source("RNAseq_function.R")


# load data ---------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1/metal_anno"


# use function run GSEA and KEGG ------------------------------------------
gse <- gsea_run(data_path = data_path,
                file="ip_Y_V_S_CO_0_deg.xlsx", all_gene = T
                )
keg <- kegg_run(data_path=data_path,
                file = "ip_Y_V_S_CO_0_deg.xlsx", all_gene = T
                )



# run GSEA -------------------------------------------------------------------
df <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_CO_0_deg.xlsx"))
organism <- "org.Hs.eg.db"
# $csv file's colume namm of log2 fold change
original_gene_list <- df$M
# $csv file's colume namm of ENSEMBL ID 
names(original_gene_list) <- df$ENSEMBL
# omit any NA values and sort the list in decreasing order
gsea_gene_list <- na.omit(original_gene_list) %>% 
  sort(., decreasing = TRUE)
#run GSEA
gse <- gseGO(geneList = gsea_gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

# gsea plot --------------------------------------------------------------------
# dot plot(gse)
dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(.~.sign)
ggsave(paste0(file_name,"gse_dot",'.png'), height = 12, width = 10, units = "in")

# Enrichment map(gse)
Enrichment_map_gsea <- pairwise_termsim(gse)
emapplot(Enrichment_map_gsea)
ggsave(paste0(file_name,"gse_Enrichment_w",'.png'), height = 12, width = 12, units = "in")

# UpSet Plot
upsetplot(gse)

# Category Netplot(gse)
# categorySize can be either 'pvalue' or 'geneNum'
cat <- DOSE::setReadable(gse, "org.Hs.eg.db", keyType = "ENSEMBL")  # mapping geneID to gene Symbol
cnetplot(cat, categorySize = "pvalue", foldChange = gsea_gene_list, showCategory = 3)
ggsave(paste0("gse_Category",".png"), height = 8, width = 10, units = "in")

# Ridgeplot(gse)
ridgeplot(gse, showCategory = 15) + labs(x = "enrichment distribution")
ggsave(paste0(file_name,"gse_Ridgeplot_w",'.png'), height = 15, width = 10, units = "in")

# GSEA plot(gse)
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
ggsave(paste0(file_name,"gse_GSEA_w",'.png'), height = 8, width = 10, units = "in")
# multiple GSEA plot
for(i in 1:5){
  gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = 1)
  ggsave(paste0(file_name,"gse_GSEA_w", i, ".png"), height = 8, width = 10, units = "in")
}

# plot GSEA result in gseaplot2()
enrichplot::gseaplot2(gsea_d_bap, geneSetID = 1, title = gsea_d_bap$Description[1])
enrichplot::gseaplot2(gse, geneSetID = 1:3, subplots = 1:3)
# overlap GSEA plot
enrichplot::gseaplot2(gse, geneSetID = 1:3, pvalue_table = TRUE,
                      color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

# gene rank
gsearank(gse, 1, title = gse[1, "Description"])

# PubMed trend of enriched terms
terms <- gse$Description[1:3]
enrichplot::pmcplot(terms, 2012:2022, proportion = FALSE)
ggsave(paste0(file_name,".gse_Pubmed",'.png'), height = 8, width = 10, units = "in")




# 若需轉換ID種類 --------------------------------------------------------
library(dplyr)
new_df <- clusterProfiler::bitr(df$ENSEMBL, 
                                fromType = "ENSEMBL", 
                                toType = c("ENTREZID","SYMBOL"), 
                                OrgDb = "org.Hs.eg.db") %>% 
  right_join(.,df, by = "ENSEMBL")
gene.df <- bitr(df$ENSEMBL, fromType = "ENSEMBL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
colnames(gene.df) <- c("ENSEMBL","ENTREZID")
merged_df <- dplyr::left_join(x = df, y = gene.df)
# run kegg --------------------------------------------------------------
# $csv file's colume namm of log2 fold change
kegg_gene_list <- merged_df$M
# Name vector with ENTREZ ids
names(kegg_gene_list) <- merged_df$ENTREZID
# omit any NA values and sort the list in decreasing order
kegg_gene_list<- na.omit(kegg_gene_list) %>% 
  sort(., decreasing = TRUE)
# Run KEGG
# change KEGG Organism Code from https://www.genome.jp/kegg/catalog/org_list.html 
kegg_organism <-  "hsa"  # human is "hsa"
keg <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
# kegg plot --------------------------------------------------------------------
# pairwise_termsim
kmat <- enrichplot::pairwise_termsim(keg)

# Heatmap-like functional classification
heat1 <- enrichplot::heatplot(kmat, showCategory=5)
heat2 <- enrichplot::heatplot(keg, foldChange=kegg_gene_list, showCategory=5)
heat2
cowplot::plot_grid(heat1, heat2, ncol=1, labels=LETTERS[1:2])
ggsave(paste0(file_name,".kegg_func_Heatmap",'.png'), height = 8, width = 10, units = "in")

# Tree plot
tree1 <- treeplot(kmat)
tree1
tree2 <- treeplot(kmat, cluster.params = list(method = "average"))
aplot::plot_list(tree1, tree2, tag_levels='A')
ggsave(paste0(file_name,".kegg_tree_plot",'.png'), height = 8, width = 12, units = "in")

# UpSet Plot
upsetplot(keg)

# Gene-Concept Network
## convert gene ID to Symbol
net <- DOSE::setReadable(keg, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(net, categorySize="pvalue", color.params = list(foldChange = kegg_gene_list))
# circular form
cnetplot(net, color.params = list(foldChange = kegg_gene_list, edge = TRUE), circular = TRUE) 
ggsave(paste0(file_name,".kegg_network",'.png'), height = 12, width = 12, units = "in")

# Bar Plot???
bar <- names(kegg_gene_list)[abs(kegg_gene_list) > 2] %>% 
        DOSE::enrichDGN() %>% 
        barplot(., showCategory=20)
bar
# dot plot
dotplot(keg, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ggsave(paste0(file_name,".kegg_dot",'.png'), height = 12, width = 10, units = "in")

# Encrichment map
emapplot(kmat)
ggsave(paste0(file_name,".kegg_Encrichment",'.png'), height = 10, width = 10, units = "in")

# Ridgeplot
ridgeplot(keg, showCategory = 20) + labs(x = "enrichment distribution")
ggsave(paste0(file_name,".kegg_Ridgeplot",'.png'), height = 20, width = 12, units = "in")

# GSEA Plot
#Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(keg, by = "all", title = keg$Description[1], geneSetID = 1)
ggsave(paste0(file_name,".kegg_GSEA",'.png'), height = 8, width = 10, units = "in")
# multiple GSEA plot
for(i in 1:5){
  gseaplot(keg, by = "all", title = gse$Description[i], geneSetID = 1)
  ggsave(paste0(file_name,".keg_GSEA", i, ".png"), height = 8, width = 10, units = "in")
}

# plot GSEA result in gseaplot2()
enrichplot::gseaplot2(keg, geneSetID = 1, title = keg$Description[1])
enrichplot::gseaplot2(keg, geneSetID = 1:3, subplots = 1:3)
# multiple GSEA plot
for(i in 1:3){
  plot <- gseaplot2(keg, geneSetID = i, subplots = 1:3, title = gse$Description[i]) %>% print()
  ggsave(paste0(file_name,".keg_GSEA", i, ".png"), height = 8, width = 10, units = "in")
}
# overlap GSEA plot
enrichplot::gseaplot2(keg, geneSetID = 1:3, pvalue_table = TRUE,
                      color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "line")


# The gsearank() plot the ranked list of genes belong to the specific gene set.
rank <- lapply(1:3, function(i) {
  anno <- keg[i, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  gsearank(keg, i, keg[i, 2]) + xlab(NULL) +ylab(NULL) +
    annotate("text", 10000, keg[i, "enrichmentScore"] * .75, label = lab, hjust=0, vjust=0)
})
plot_grid(plotlist = rank, ncol=1)

# KEGG(for single pathway) 
setwd("./plot")

kegg_pathway_id <- "hsa04115"
pathview(gene.data = gsea_gene_list, pathway.id = kegg_pathway_id, 
         species = kegg_organism, gene.idtype = gene.idtype.list[3], kegg.native = T)
# KEGG(for multiple pathways)
setwd("./plot")
keg_ddr <- c("hsa03030","hsa03410","hsa03420","hsa03430",
             "hsa03440","hsa03450","hsa03460")
keg_cyp <- c("hsa00980","hsa00982","hsa05204")
keg_egfr <- c("hsa01521")
for (kegg_pathway_id in keg_cyp) {
   pathview(gene.data = gsea_gene_list, pathway.id = kegg_pathway_id, 
            species = kegg_organism, gene.idtype = gene.idtype.list[3], kegg.native = T)
}







