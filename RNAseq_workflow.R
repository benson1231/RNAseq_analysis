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
library(KEGGREST)

source("RNAseq_function.R")

# 0.setting -----------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/ST"
output_dir <- "/Users/benson/Documents/project/RNA-seq1-3/output"
name_df <- openxlsx::read.xlsx("/Users/benson/Documents/project/RNA-seq1-3/data/name.xlsx")
# color setting
red_white_blue_colors <- c("blue", "white", "red")
rwb_palette <- colorRampPalette(red_white_blue_colors)
n_colors <- 10
colors <- rwb_palette(n_colors)

# 1.input raw_count data ----------------------------------------------------
### load annotation data
gene_df <- "/Users/benson/Documents/project/RNA-seq1-3/data/anno_gene.RDS" %>% 
  readRDS() %>%
  select(ENSEMBL,SYMBOL)

### load raw reads counts
mylogcount_df <- "/Users/benson/Documents/project/RNA-seq1-3/data/nor_logcounts.RDS" %>% 
  readRDS() %>% as.data.frame() %>% select(name_df$group_name) %>% setNames(name_df$abbreviate)
mycount_df <- "/Users/benson/Documents/project/RNA-seq1-3/data/nor_counts.RDS" %>% 
  readRDS() %>% as.data.frame()
abbr_count <- mycount_df %>% select(name_df$group_name) %>% setNames(name_df$abbreviate)
# count for decoupleR
gene_count <- mycount_df %>% 
  rownames_to_column("ENSEMBL") %>% 
  left_join(gene_df,"ENSEMBL") %>% 
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL") %>% 
  select(name_df$group_name) %>% setNames(name_df$abbreviate)

# 2.bar plot --------------------------------------------------------------
draw_bar("ANGPT2")
draw_bar("ANKRD1")

# 3.MD plot -----------------------------------------------------------------
# only DEG
plot_MD("ip_D_V_S_HCD_0_deg.xlsx", only_DE = T)
# with non-significant
plot_MD("ip_Y_V_S_HCD_BAP_0_deg.xlsx")

# 4.MDS -------------------------------------------------------------------
### euclidean distance
euclidean_dist <- dist(t(mylogcount_df), method = "euclidean")
distance_matrix <- as.matrix(euclidean_dist) 
# euclidean distance heatmap
ComplexHeatmap::Heatmap(distance_matrix, name = "euclidean_dist")

### Create a DGEList object
y <- DGEList(abbr_count) 
abbr_sampleinfo <- data.frame(row = names(abbr_count),
                         sample = names(abbr_count), 
                         treatment = str_sub(names(abbr_count), start=5),
                         culture = str_sub(names(abbr_count), start=3,end=3),
                         cell    = str_sub(names(abbr_count), start=1,end=1)) %>%
  column_to_rownames('row') 
head(abbr_sampleinfo)
# set factor
abbr_sampleinfo$sample <- factor(abbr_sampleinfo$sample)
abbr_sampleinfo$treatment <- factor(abbr_sampleinfo$treatment)
abbr_sampleinfo$culture <- factor(abbr_sampleinfo$culture)
abbr_sampleinfo$cell <- factor(abbr_sampleinfo$cell)
levels(abbr_sampleinfo$sample)
levels(abbr_sampleinfo$treatment)
levels(abbr_sampleinfo$culture)
levels(abbr_sampleinfo$cell)
# color
col.sample <- c("black")[abbr_sampleinfo$sample]
col.treatment <- c('#FFE66F','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000',
                   '#E0E0E0','#FF9797','#ADADAD','#005AB5','#000079','#0080FF',
                   '#0000E3')[abbr_sampleinfo$treatment]
col.culture <- c("#0072E3","darkblue")[abbr_sampleinfo$culture]
col.cell <- c('#00DB00','#FF88C2','#0080FF', '#FF8000')[abbr_sampleinfo$cell]
# # sample
# plotMDS(y,col=col.sample ,xlab = "Dimension 1",ylab = "Dimension 2")
# title("sample")
### treatment(carcinogen)
plotMDS(y,col=col.treatment,xlab = "Dimension 1",ylab = "Dimension 2")
legend("topleft",fill=c('#FFE66F','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000',
                        '#E0E0E0','#FF9797','#ADADAD','#005AB5','#000079','#0080FF',
                        '#0000E3'),
       legend=levels(abbr_sampleinfo$treatment))
title("Carcinogen Treatment")
# ### culture(long/short-term)
# plotMDS(y,col=col.culture,xlab = "Dimension 1",ylab = "Dimension 2")
# legend("topleft",fill=c("darkblue","#0072E3"),legend=levels(abbr_sampleinfo$culture))
# title("Protocol(Short-term/Long-term)")
### cell type
plotMDS(y,col=col.cell,xlab = "Dimension 1",ylab = "Dimension 2")
legend("topleft",fill=c('#00DB00','#FF88C2','#0080FF', '#FF8000'),legend=levels(abbr_sampleinfo$cell))
title("Cell type")

# 5.ORA ---------------------------------------------------------------------
### enrichKEGG(KEGG pathway)
keg_result <- run_keg_path("ip_Y_V_S_CO_BAP_0_deg.xlsx", dir = "up",log_crit = 1)
view(keg_result %>% as.data.frame() %>% rownames_to_column("new_ID"))
# bar plot
barplot(keg_result, showCategory=10, font.size = 9, x = "GeneRatio", label_format = 40,
        title = "")   # change
# bar plot with qscore
mutate(keg_result, qscore = -log(p.adjust, base=10)) %>% 
  barplot(., showCategory=10, font.size = 9, x = "qscore", label_format = 40,
          title = "HCD")   # change
# heatmap
# list <- force_list(keg_result@result$geneID[4])
# draw_from_list(list =list , groups = "CO", id = "SYMBOL")
list <- get_kegg_list(path_ID = keg_result@result$ID[6]) 
draw_from_list(list =list , groups = "CD", id = "SYMBOL",title = "Mineral absorption")
kegg_ID_order <- 3  # change
group_order <- 6  # change
plotMD_kegg(name_df$file_name[group_order],keg_result@result$ID[kegg_ID_order],
            title = name_df$abbreviate[group_order],with_line = F,only_DE = F)
# bar plot with qscore
mutate(keg_result, qscore = -log(p.adjust, base=10)) %>% 
  barplot(., showCategory=10, font.size = 9, x = "qscore", label_format = 40,
          title = "CO")   # change

### enrichPathway(Reactome pathway)
react_result <- run_reactome("ip_D_V_S_HCD_0_deg.xlsx",dir = "up")
# view(react_result %>% as.data.frame() %>% rownames_to_column("new_ID"))
# heatmap
list <- force_list(react_result@result$geneID[6])
draw_from_list(list =list , groups = "ALL", id = "SYMBOL")
# bar plot
barplot(react_result, showCategory=10, font.size = 9, x = "GeneRatio", label_format = 40,
        title = "") 
# bar plot with qscore
mutate(react_result, qscore = -log(p.adjust, base=10)) %>% 
  barplot(., showCategory=10, font.size = 9, x = "qscore", label_format = 40,
          title = "")   # change

# 6.ORA several groups with heatmap -----------------------------------------
file_list <- c("ip_W_V_S_HCD_0_deg.xlsx", "ip_D_V_S_HCD_0_deg.xlsx", 
               "ip_Y_V_S_HCD_0_deg.xlsx","ip_L_V_S_HCD_0_deg.xlsx")
group_names <- c("WT", "Del19", "YAP","L858R")
### KEGG
p1 <- plot_heatmap(file_list, group_names, analysis = "kegg", dir="down", col_title = "HCD",row_title = "KEGG")
### Reactome
p2 <- plot_heatmap(file_list, group_names, analysis = "reactome", dir="down", row_title = "Reactome")
p1%v%p2

# 7.GSEA analysis -----------------------------------------------------------
### KEGG
keg <- kegg_run("ip_Y_V_S_HCD_0_deg.xlsx")  # change
### dot plot
dotplot(keg, showCategory = 10, label_format=50, 
        title = "Enriched Pathways", split=".sign") + 
  facet_grid(.~.sign)
### Gene-Concept Network
keg_gene <- DOSE::setReadable(keg, 'org.Hs.eg.db', 'ENTREZID') # convert gene ID to Symbol
View(keg_gene@result %>% rownames_to_column("pathway_ID"))
cnetplot(keg_gene, categorySize="pvalue", 
         showCategory = keg_gene@result$Description[c(16)],  # change
         color.params = list(foldChange = keg_gene@geneList[abs(keg_gene@geneList)>1]  # change
         )) +
  scale_color_gradientn(name = "logFC", colors=colors, 
                        na.value = "#E5C494", limits = c(-2, 2))
### GSEA plot
ID <- 15
enrichplot::gseaplot2(keg, geneSetID = ID, title = keg$Description[ID])
### Tree plot
kmat <- enrichplot::pairwise_termsim(keg)
treeplot(kmat, cluster.params = list(method = "average"))
# heatmap
list <- force_list(keg_gene@result$core_enrichment[15])
draw_from_list(list =list , groups = "ALL", id = "SYMBOL")

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
list <- c("CHRNA4", "CHRNB4", "VDR", "PGR", "ALKBH4", "CYP1A1", "CYP1B1","HIF1A","YAP1")
list <- c("GADD45B","FHIT","CDKN1A","GADD45G","BAD","HRAS","BAX","MAP2K2","CDK4","AKT1","ALKBH4")
list <- c("TP53","YAP1","LATS1","LATS2","MTF1","MT1F")
draw_from_list(list = list , groups = "CON", id = "SYMBOL")

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
group <- "only_CO"
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
writeLines(unlist(cluster1), file.path(output_dir,"cluster_genes.txt"))
draw_from_list(list = cluster4, groups = group, id = "SYMBOL",show_row_names = T,
               title = "C4")

# 11.venn diagram ------------------------------------------------------------
W_up <- get_deg("ip_W_V_S_AS_0_deg.xlsx", log_crit = c(1,-1),dir = "down")
L_up <- get_deg("ip_L_V_S_AS_0_deg.xlsx", log_crit = c(1,-1),dir = "down")
D_up <- get_deg("ip_D_V_S_AS_0_deg.xlsx", log_crit = c(1,-1),dir = "down")
Y_up <- get_deg("ip_Y_V_S_AS_0_deg.xlsx", log_crit = c(1,-1),dir = "down")
# venn diagram list
venn_list <- list(W = W_up,
                  L = L_up,
                  D = D_up,
                  Y = Y_up)
# venn diagram
ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#FF2D2D")+ggtitle("AS")
ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#6A6AFF")+ggtitle("AS")
# upset
ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0, force_upset = T)
# get venn diagram result(gene list)
venn_result <- process_region_data(Venn(venn_list))
venn_result
draw_from_list(list = venn_result$item[[15]], groups = "AS")

### use divenn2.0
get_divenn("ip_Y_V_S_CO_0_deg.xlsx",output_name = "CO.xlsx")
get_divenn("ip_Y_V_S_CO_BAP_0_deg.xlsx",output_name = "CO_BAP.xlsx")
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
file_name <- "ip_D_V_S_HCD_0_deg.xlsx"   # change
# pathway barplot
run_pathway(file_name, title = abb(file_name))
# pathway specific genes MD plot
path <- plot_pathway(file_name, "Hypoxia", title = abb(file_name))
# top50 heatmap
path_top50 <- path %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(ID)
draw_from_list(list = path_top50, groups = "CD",
               id = "SYMBOL",label_num = F,anno = T,title = "")

### TFs for specific group
file_name <- "ip_W_V_S_AS_BAP_0_deg.xlsx"   # change
# TFs barplot
run_TF(file_name, title = abb(file_name))
# TFs specific genes MD plot
TF <- plot_TF(file_name, "SP1", title = abb(file_name))
# ggsave("TF.jpeg")
# top50 heatmap
TF_top50 <- TF %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(ID)
draw_from_list(list = TF_top50,
               groups = "CO",
               id = "SYMBOL",label_num = F,anno = T,title = "")
# # for loop
# for (i in name_df$file_name[29:39]) {
#   plot_TF(i, "HIF1A", title = abb(i))
#   ggsave(paste(abb(i),file.path(output_dir,"HIF1A.jpeg")))
#   dev.off()
# }









