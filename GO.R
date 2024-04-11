# library -----------------------------------------------------------------
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
source("RNAseq_function.R")

# create_enrich_list ------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/V"
list <- create_enrich_list("ip_Y_V_S_HCD_0_deg.xlsx",type = "ENTREZID")
gene <- names(list)[abs(list) > 1]

# GO classification -------------------------------------------------------
ugo <- groupGO(gene     = gene,
               OrgDb    = "org.Hs.eg.db",
               ont      = "MF",
               level    = 3,
               readable = TRUE)

head(ugo)
categorys <- ugo@result %>% arrange(desc(Count)) %>% pull(Description)
barplot(ugo,showCategory = categorys[1:10],x = "GeneRatio")

# GO over-representation analysis -----------------------------------------
ego <- enrichGO(gene          = gene,
                universe      = names(list),
                OrgDb         = "org.Hs.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
barplot(ego, x = "GeneRatio")
dotplot(ego)
goplot(ego)

# GO Gene Set Enrichment Analysis -----------------------------------------
gse <- gseGO(geneList     = list,  # desc order list
             OrgDb        = "org.Hs.eg.db",
             ont          = "BP",
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE)
head(gse)
dotplot(gse)
goplot(gse)

# Reactome enrichment analysis --------------------------------------------
# Reactome pathway over-representation analysis
de_pathway <- enrichPathway(gene=gene, pvalueCutoff = 0.05, readable=TRUE)
head(de_pathway)
dotplot(de_pathway)

# Reactome pathway gene set enrichment analysis
gse_pathway <- gsePathway(list, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)
head(gse_pathway)
dotplot(gse_pathway)

# Pathway Visualization
# change colors
red_white_blue_colors <- c("blue", "white", "red")
rwb_palette <- colorRampPalette(red_white_blue_colors)
n_colors <- 10
colors <- rwb_palette(n_colors)
# change work space
setwd("./plot")

viewPathway(gse_pathway@result$Description[1], 
            readable = TRUE, 
            foldChange = list) + 
  ggtitle(gse_pathway@result$Description[1])+
  scale_color_gradientn(name = "fold change", colors=colors, na.value = "#E5C494")

for(i in 1:3){
  viewPathway(gse_pathway@result$Description[i], 
              readable = TRUE, 
              foldChange = list) + 
    ggtitle(gse_pathway@result$Description[i]) +
    scale_color_gradientn(name = "fold change", colors=colors, na.value = "#E5C494")
  
  cat("-> saving image ",i,"\n")
  ggsave(paste0(gse_pathway@result$Description[i],'.png'), 
         height = 12, width = 12, units = "in")
}

viewPathway("Xenobiotics", readable = TRUE, foldChange = list) +
  scale_color_gradientn(name = "fold change", colors=colors, na.value = "#E5C494")

ggsave("Xenobiotics.png", 
       height = 12, width = 12, units = "in")
