# 0.load library for package ----------------------------------------------
library(tidyverse)
library(clusterProfiler) 
library(enrichplot)
library(ReactomePA)
library(decoupleR)
library(pathview)
library(cowplot)
library(DOSE)
library(ComplexHeatmap)
library(ggVennDiagram)
library(colorRamp2)
library(KEGGREST)
library(openxlsx)
library(ggrepel)
library(edgeR)

source("RNAseq_function.R")

# 1.setting -----------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/Vitro"
name_df <- read.xlsx("/Users/benson/Documents/project/RNA-seq1-3/data/name.xlsx")
# color setting
red_white_blue_colors <- c("blue", "white", "red")
rwb_palette <- colorRampPalette(red_white_blue_colors)
n_colors <- 10
colors <- rwb_palette(n_colors)

# 2.input raw_count data ----------------------------------------------------
### load annotation data
gene_df <- "/Users/benson/Documents/project/RNA-seq1-3/data/anno_gene.RDS" %>% 
  readRDS() %>%
  select(ENSEMBL,SYMBOL)

### load raw reads counts 
mycount_77 <- "/Users/benson/Documents/project/RNA-seq1-3/data/myTMM.RDS" %>% 
  readRDS() %>% as.data.frame()
mycount_df <- mycount_77 %>% 
  select(all_of(name_df$group_name[27:39]))
abbr_count <- mycount_77 %>% 
  select(all_of(name_df$group_name[27:39])) %>% 
  setNames(name_df$abbreviate[27:39])
# count for decoupleR
gene_count <- mycount_df %>% 
  rownames_to_column("ENSEMBL") %>% 
  left_join(gene_df,"ENSEMBL") %>% 
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL")

# 3.MD plot -----------------------------------------------------------------
plot_MD("ip_Y_V_S_HCD_BAP_0_deg.xlsx")

# 4.MDS -------------------------------------------------------------------
### euclidean distance
euclidean_dist <- dist(t(abbr_count), method = "euclidean")
distance_matrix <- as.matrix(euclidean_dist)

# euclidean distance heatmap
ComplexHeatmap::Heatmap(distance_matrix, name = "euclidean_dist")

### Create a DGEList object
y <- DGEList(count) 
# set factor
sampleinfo$sample <- factor(sampleinfo$sample)
sampleinfo$treatment <- factor(sampleinfo$treatment)
sampleinfo$culture <- factor(sampleinfo$culture)
sampleinfo$cell <- factor(sampleinfo$cell)
levels(sampleinfo$sample)
levels(sampleinfo$treatment)
levels(sampleinfo$culture)
levels(sampleinfo$cell)
# color
col.sample <- c("black")[sampleinfo$sample]
col.treatment <- c('#FFE66F','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000',
                   '#E0E0E0','#FF9797','#ADADAD','#005AB5','#000079','#0080FF',
                   '#0000E3')[sampleinfo$treatment]
col.culture <- c("#0072E3","darkblue")[sampleinfo$culture]
col.cell <- c('#FF9224','#6F00D2')[sampleinfo$cell]
# sample
plotMDS(y,col=col.sample ,xlab = "PC1",ylab = "PC2")
title("sample")
### treatment(carcinogen)
plotMDS(y,col=col.treatment,xlab = "PC1",ylab = "PC2")
legend("topleft",fill=c('#FFE66F','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000',
                        '#E0E0E0','#FF9797','#ADADAD','#005AB5','#000079','#0080FF',
                        '#0000E3'),
       legend=levels(sampleinfo$treatment))
title("Carcinogen Treatment")
### culture(long/short-term)
plotMDS(y,col=col.culture,xlab = "PC1",ylab = "PC2")
legend("topleft",fill=c("darkblue","#0072E3"),legend=levels(sampleinfo$culture))
title("Protocol(Short-term/Long-term)")
### cell type
plotMDS(y,col=col.cell,xlab = "PC1",ylab = "PC2")
legend("topleft",fill=c("orange","#6F00D2"),legend=levels(sampleinfo$cell))
title("Cell type")

# 5.ORA ---------------------------------------------------------------------
### enrichKEGG(KEGG pathway)
keg_result <- run_keg_path("ip_Y_V_S_CO_BAP_0_deg.xlsx", dir = "up",log_crit = 1)
# bar plot
barplot(keg_result, showCategory=10, font.size = 9, x = "GeneRatio", label_format = 40,
        title = "")   # change
# bar plot with qscore
mutate(keg_result, qscore = -log(p.adjust, base=10)) %>% 
  barplot(., showCategory=10, font.size = 9, x = "qscore", label_format = 40,
          title = "")   # change

### enrichPathway(Reactome pathway)
react_result <- run_reactome("ip_Y_V_S_HCD_BAP_0_deg.xlsx",dir = "up")
# bar plot
barplot(react_result, showCategory=10, font.size = 9, x = "GeneRatio", label_format = 40,
        title = "") 
# bar plot with qscore
mutate(react_result, qscore = -log(p.adjust, base=10)) %>% 
  barplot(., showCategory=10, font.size = 9, x = "qscore", label_format = 40,
          title = "")   # change

# 6. ORA several groups with heatmap -----------------------------------------
file_list <- c("ip_Y_V_S_CO_0_deg.xlsx", "ip_Y_V_S_BAP_0_deg.xlsx", "ip_Y_V_S_CO_BAP_0_deg.xlsx")
group_names <- c("CO", "BAP", "CO_BAP")
### KEGG
plot_heatmap(file_list, group_names, analysis = "kegg", dir="up", title = "")
### Reactome
plot_heatmap(file_list, group_names, analysis = "reactome", dir="up", title = "")


# 7.GSEA analysis -----------------------------------------------------------
### KEGG
keg <- kegg_run("ip_Y_V_S_HCD_BAP_0_deg.xlsx")  # change
### dot plot
dotplot(keg, showCategory = 10, label_format=50, 
        title = "Enriched Pathways", split=".sign") + 
  facet_grid(.~.sign)
### Gene-Concept Network
net <- DOSE::setReadable(keg, 'org.Hs.eg.db', 'ENTREZID') # convert gene ID to Symbol
# View(net@result %>% rownames_to_column("pathway_ID"))
cnetplot(net, categorySize="pvalue", 
         showCategory = net@result$Description[c(1)],  # change
         color.params = list(foldChange = net@geneList[abs(net@geneList)>1]  # change
         )) +
  scale_color_gradientn(name = "logFC", colors=colors, 
                        na.value = "#E5C494", limits = c(-2, 2))
### GSEA plot
ID <- 39
enrichplot::gseaplot2(keg, geneSetID = ID, title = keg$Description[ID])
### Tree plot
kmat <- enrichplot::pairwise_termsim(keg)
treeplot(kmat, cluster.params = list(method = "average"))

# 8.get DEG list ------------------------------------------------------------
file_name <- "ip_Y_V_S_HCD_BAP_0_deg.xlsx"  # change
group <- "CD"
# up_all
DEG <- get_deg(file_name, log_crit = c(1,-1),dir = "up",type = "SYMBOL")
draw_from_list(list = DEG, groups = group, id = "SYMBOL", show_row_names = F)
# up_100
top <- get_deg(file_name, log_crit = c(1,-1),dir = "up",type = "SYMBOL", top = 100)
draw_from_list(list = top, groups = group, id = "SYMBOL")
# down_all
DEG <- get_deg(file_name, log_crit = c(1,-1),dir = "down",type = "SYMBOL")
draw_from_list(list = DEG, groups = group, id = "SYMBOL", show_row_names = F)
# down_100
top <- get_deg(file_name, log_crit = c(1,-1),dir = "down",type = "SYMBOL", top = 100)
draw_from_list(list = top, groups = group, id = "SYMBOL")

# 9.heatmap -----------------------------------------------------------------
# single
list <- c("CHRNA4", "CHRNB4", "VDR", "PGR", "ALK", "CYP1A1")
draw_from_list(list =list , groups = "ALL", id = "SYMBOL")

# complex
list1 <- c("EGFR","ALK","ROS1","BRAF","MET","RET","Her2","KRAS","TP53","PTEN",
          "ERBB2", "HRAS", "NRAS","STK11","NTRK1", "NTRK2","NTRK3")
list2 <- c("MT1A", "MT1B", "MT1E", "MT1F", "MT1G", "MT1H", "MT1M", "MT1X","MT2A",
           "TP53","NFKB1","NQO1","GCLC","MCM2","ALK","NPM1","YAP1","JUN","EGFR")
p1 <- draw_from_list(list = list1,
                     groups = "ALL", show_row_names = T,
                     id = "SYMBOL",label_num = F,anno = T,title = "")
p2 <- draw_from_list(list = list2,
                     groups = "ALL", show_row_names = T,
                     id = "SYMBOL",label_num = F,anno = F,title = "")
p1 %v% p2

# 10.cluster analysis ------------------------------------------------------
file_name <- "ip_Y_V_S_CO_BAP_0_deg.xlsx"
group <- "only_HCD"
k_value <- 5
draw_heatmap(file_name, groups = group, log_crit = 1, km = k_value)
# run cluster analysis
cluster_result <- draw_heatmap(file_name, groups = group, log_crit = 1, 
                               row_km = T, km=k_value, return_cluster = T)
# get cluster gene SYMBOL
cluster1 <- names(cluster_result [cluster_result == 1])
cluster2 <- names(cluster_result [cluster_result == 2])
cluster3 <- names(cluster_result [cluster_result == 3])
cluster4 <- names(cluster_result [cluster_result == 4])
cluster5 <- names(cluster_result [cluster_result == 5])
cluster6 <- names(cluster_result [cluster_result == 6])
# heatmap with control/DMSO
p1 <- draw_from_list(cluster1, groups = group,title = "C1")
p2 <- draw_from_list(cluster2, groups = group,anno = F,title = "C2")
p3 <- draw_from_list(cluster3, groups = group,anno = F,title = "C3")
p4 <- draw_from_list(cluster4, groups = group,anno = F,title = "C4")
p5 <- draw_from_list(cluster5, groups = group,anno = F,title = "C5")

p1 %v% p2 %v% p3 %v% p4 %v% p5
# gene SYMBOL output
writeLines(unlist(cluster1), "cluster_genes.txt")
draw_from_list(list = cluster4, groups = group, id = "SYMBOL",show_row_names = T,
               title = "C4")

# 11.venn diagram ------------------------------------------------------------
CO_up <- get_deg("ip_Y_V_S_CO_0_deg.xlsx", log_crit = c(1,-1),dir = "up")
BAP_up <- get_deg("ip_Y_V_S_BAP_0_deg.xlsx", log_crit = c(1,-1),dir = "up")
CO_BAP_up <- get_deg("ip_Y_V_S_CO_BAP_0_deg.xlsx", log_crit = c(1,-1),dir = "up")
# venn diagram list
venn_list <- list(CO = CO_up,
                  BAP = BAP_up,
                  CO_BAP =CO_BAP_up)
# venn diagram
ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#FF2D2D")
ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#6A6AFF")
# upset
ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0, force_upset = T)
# get venn diagram result(gene list)
venn_result <- process_region_data(Venn(venn_list))
venn_result
draw_from_list(list = venn_result$item[[7]], groups = "CO")

### use divenn2.0
get_divenn("ip_Y_V_S_CO_0_deg.xlsx")
get_divenn("ip_Y_V_S_CO_BAP_0_deg.xlsx")
# https://divenn.tch.harvard.edu/v2/


# 12.pathway and TFs analysis using decoupleR --------------------------------
net <- "/Users/benson/Documents/project/RNA-seq1-3/data/net.RDS" %>% readRDS()
net_TF <- "/Users/benson/Documents/project/RNA-seq1-3/data/net_TF.RDS" %>% readRDS()

### heatmap
# pathway heatmap
run_path_heatmap(gene_count)
# TFs heatmap
run_TF_heatmap(gene_count)

### pathway for specific group
file_name <- "ip_Y_V_S_CO_BAP_0_deg.xlsx"   # change
# pathway barplot
run_pathway(file_name, title = abb(file_name))
# pathway specific genes MD plot
path <- plot_pathway(file_name, "p53", title = abb(file_name))
# top50 heatmap
path_top50 <- path %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(ID)
draw_from_list(list = path_top50, groups = "CD",
               id = "SYMBOL",label_num = F,anno = T,title = "")

### TFs for specific group
file_name <- "ip_Y_V_S_CO_BAP_0_deg.xlsx"   # change
# TFs barplot
run_TF(file_name, title = abb(file_name))
# TFs specific genes MD plot
TF <- plot_TF(file_name, "HIF1A", title = abb(file_name))
# ggsave("TF.jpeg")
# top50 heatmap
TF_top50 <- TF %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(ID)
draw_from_list(list = TF_top50,
               groups = "CO",
               id = "SYMBOL",label_num = F,anno = T,title = "")
## for loop
# for (i in name_df$file_name[29:39]) {
#   plot_TF(i, "HIF1A", title = abb(i))
#   ggsave(paste(abb(i),"HIF1A.jpeg"))
#   dev.off()
# }









