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

# use function to run GSEA ----------------------------------------------
gse <- gsea_run("ip_Y_V_S_AS_0_deg.xlsx", all_gene = T)
gse <- gsea_run("ip_Y_V_S_CO_0_deg.xlsx", all_gene = T)
gse <- gsea_run("ip_Y_V_S_LCD_0_deg.xlsx", all_gene = T)
gse<- gsea_run("ip_Y_V_S_HCD_0_deg.xlsx", all_gene = T)
gse <- gsea_run("ip_Y_V_S_BAP_0_deg.xlsx", all_gene = T)
gse <- gsea_run("ip_Y_V_S_AS_BAP_0_deg.xlsx", all_gene = T)
gse <- gsea_run("ip_Y_V_S_CO_BAP_0_deg.xlsx", all_gene = T)
gse <- gsea_run("ip_Y_V_S_LCD_BAP_0_deg.xlsx", all_gene = T)
gse <- gsea_run("ip_Y_V_S_HCD_BAP_0_deg.xlsx", all_gene = T)

# use function to run KEGG ------------------------------------------------
keg <- kegg_run("ip_Y_V_S_AS_0_deg.xlsx", all_gene = T)
keg <- kegg_run("ip_Y_V_S_CO_0_deg.xlsx", all_gene = T)
keg <- kegg_run("ip_Y_V_S_LCD_0_deg.xlsx", all_gene = T)
keg <- kegg_run("ip_Y_V_S_HCD_0_deg.xlsx", all_gene = T)
keg <- kegg_run("ip_Y_V_S_BAP_0_deg.xlsx", all_gene = T)
keg <- kegg_run("ip_Y_V_S_AS_BAP_0_deg.xlsx", all_gene = T)
keg <- kegg_run("ip_Y_V_S_CO_BAP_0_deg.xlsx", all_gene = T)
keg <- kegg_run("ip_Y_V_S_LCD_BAP_0_deg.xlsx", all_gene = T)
keg <- kegg_run("ip_Y_V_S_HCD_BAP_0_deg.xlsx", all_gene = T)

# GSEA plots --------------------------------------------------------------
# dot plot(gse)
dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(.~.sign)
# Enrichment map(gse)
Enrichment_map_gsea <- pairwise_termsim(gse)
emapplot(Enrichment_map_gsea)
# UpSet Plot
upsetplot(gse)
# Category Netplot(gse)
# categorySize can be either 'pvalue' or 'geneNum'
cat <- DOSE::setReadable(gse, "org.Hs.eg.db", keyType = "ENSEMBL")  # mapping geneID to gene Symbol
cnetplot(cat, categorySize = "pvalue", showCategory = 3,
         color.params = list(foldChange = gse@geneList))
# filtered data(ORA)
cnetplot(cat, categorySize = "pvalue", showCategory = 3,
         color.params = list(foldChange = gse@geneList), circular=T)
# Ridgeplot(gse)
ridgeplot(gse, showCategory = 15) + labs(x = "enrichment distribution")
# GSEA plot(gse)
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
# plot GSEA result in gseaplot2()
enrichplot::gseaplot2(gse, geneSetID = 1, title = gse$Description[1])
enrichplot::gseaplot2(gse, geneSetID = 1:3, subplots = 1:3)
# overlap GSEA plot
enrichplot::gseaplot2(gse, geneSetID = 1:3, pvalue_table = TRUE,
                      color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
# gene rank
gsearank(gse, 1, title = gse[1, "Description"])


# KEGG plots --------------------------------------------------------------
# pairwise_termsim
kmat <- enrichplot::pairwise_termsim(keg)

# Heatmap-like functional classification
enrichplot::heatplot(kmat, showCategory=5, foldChange = keg@geneList)
# Tree plot
treeplot(kmat, cluster.params = list(method = "average"))
# UpSet Plot
upsetplot(keg)
# Gene-Concept Network
## convert gene ID to Symbol
net <- DOSE::setReadable(keg, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(net, categorySize="pvalue", showCategory = 3,
         color.params = list(foldChange = keg@geneList))
# circular form
cnetplot(net, color.params = list(foldChange = keg@geneList, edge = TRUE), 
         showCategory = 3, circular = TRUE) 
# dot plot
dotplot(keg, showCategory = 10, title = "Enriched Pathways" , split=".sign") + 
  facet_grid(.~.sign)
# Encrichment map
emapplot(kmat)
# Ridgeplot
ridgeplot(keg, showCategory = 20) + labs(x = "enrichment distribution")
# GSEA Plot
#Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(keg, by = "all", title = keg$Description[1], geneSetID = 1)
# plot GSEA result in gseaplot2()
enrichplot::gseaplot2(keg, geneSetID = 1, title = keg$Description[1])
enrichplot::gseaplot2(keg, geneSetID = 1:3, subplots = 1:3)
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


# KEGG pathway analysis ---------------------------------------------------
setwd("./plot")

# KEGG(for single pathway) 
kegg_pathway_id <- "hsa04115"
pathview(gene.data = gse@geneList, pathway.id = kegg_pathway_id, 
         species = "hsa", gene.idtype = gene.idtype.list[3], kegg.native = T)

# KEGG(for multiple pathways)
keg_ddr <- c("hsa03030","hsa03410","hsa03420","hsa03430",
             "hsa03440","hsa03450","hsa03460")
keg_cyp <- c("hsa00980","hsa00982","hsa05204")
keg_egfr <- c("hsa01521")
for (kegg_pathway_id in keg_ddr) {
  pathview(gene.data = gse@geneList, pathway.id = kegg_pathway_id, 
           species = "hsa", gene.idtype = gene.idtype.list[3], kegg.native = T)
}







