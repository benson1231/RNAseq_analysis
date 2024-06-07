### 0.load library for package ----------------------------------------------
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
library(openxlsx)
library(ggrepel)
library(edgeR)
library(TCGAbiolinks)
# load my function
source("RNAseq_function.R")
# pre-load data
wd <- "/Users/benson/Documents/project/RNA-seq1-3/"
setwd(wd)
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/ST"
output_dir <- "/Users/benson/Documents/project/RNA-seq1-3/output"
name_df <- openxlsx::read.xlsx("/Users/benson/Documents/project/RNA-seq1-3/data/name.xlsx")
# color setting
red_white_blue_colors <- c("blue", "white", "red")
rwb_palette <- grDevices::colorRampPalette(red_white_blue_colors)
n_colors <- 10
colors <- rwb_palette(n_colors)
clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
          'Del19' = '#0066FF', 'YAP' = '#CA8EFF')

### 1-1.input raw_count data ----------------------------------------------------
# load annotation data
gene_df <- "/Users/benson/Documents/project/RNA-seq1-3/data/anno_gene.RDS" %>% 
  readRDS() %>% dplyr::select(ENSEMBL,SYMBOL)
# load log raw reads counts(logCPM)
mylogCPM <- "/Users/benson/Documents/project/RNA-seq1-3/data/nor_logcounts.RDS" %>% 
  readRDS() %>% as.data.frame() %>% dplyr::select(name_df$group_name) %>% 
  setNames(name_df$abbreviate) 
mylogCPM_SYMBOL <- "/Users/benson/Documents/project/RNA-seq1-3/data/nor_logcounts.RDS" %>% 
  readRDS() %>% as.data.frame() %>% dplyr::select(name_df$group_name) %>% 
  setNames(name_df$abbreviate) %>% rownames_to_column("ENSEMBL") %>% 
  left_join(gene_df,"ENSEMBL") %>% 
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL")
# load raw reads counts(CPM)
mycount_df <- "/Users/benson/Documents/project/RNA-seq1-3/data/nor_counts.RDS" %>% 
  readRDS() %>% as.data.frame()
# clear column name(abbreviate)
abbr_count <- mycount_df %>% 
  dplyr::select(name_df$group_name) %>% setNames(name_df$abbreviate)
# row name from ENSEMBL to SYMBOL for decoupleR
gene_count <- abbr_count %>% 
  rownames_to_column("ENSEMBL") %>% 
  left_join(gene_df,"ENSEMBL") %>% 
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL")
# TCGA patients' mutation rate
muta_rate <- "/Users/benson/Documents/project/RNA-seq1-3/data/muta_rate.RDS" %>% 
  readRDS()

### 2-1.heatmap ----------------------------------------------------------
gene_list <- c("ROX1","YAP1","LATS1","LATS2","TEAD")
draw_from_list(gene_list, groups = "CO")
draw_from_list(gene_list, groups = "CD", col_cluster = T, row_cluster = T)
# heatmap for lists
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


### 2-2.bar plot  -----------------------------------------------------------
# bar plot of gene expression(CPM)
draw_bar("RRAD", group = "ALL")

### 2-3.MD plot --------------------------------------------------------
# log fold changes(differences, D) versus average log values(means, M)
plot_MD("ip_L_V_S_AS_BAP_0_deg.xlsx", only_DE = T, plot_type = 1)
# with non-significant genes(it may take while)
plot_MD("ip_L_V_S_HCD_BAP_0_deg.xlsx")

### 3-1.euclidean distance(use logCPM) ----------------------------
euclidean_dist <- dist(t(mylogCPM), method = "euclidean")
distance_matrix <- as.matrix(euclidean_dist) 
# euclidean distance heatmap
ComplexHeatmap::Heatmap(distance_matrix, name = "euclidean_dist")

### 3-2.highly variable genes ------------------------------------
# preparing sample annotation
abbr_sampleinfo <- data.frame(row = names(mylogCPM),
                              sample = names(mylogCPM), 
                              treatment = str_sub(names(mylogCPM), start=3),
                              cell = str_sub(names(mylogCPM), start=1,end=1)) %>%
  column_to_rownames('row') 
head(abbr_sampleinfo)
# calculate variance
var_genes <- apply(mylogCPM, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- mylogCPM[select_var,] %>% as.matrix()
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
# plot heatmap of 500 high variable genes
draw_from_list(rownames(highly_variable_lcpm),groups = "ALL",id = "ENSEMBL",
               show_row_names = F, title = "Top 500 most variable genes across samples",
               col_cluster = F,filter = F,row_km = 5)

### highly variable genes in gene SYMBOL
var_genes <- apply(gene_count, 1, var) %>% sort(,decreasing = T)
head(var_genes, 50)

### 3-3.MDS ----------------------------------
# Create a DGEList object
DEG_obj <- edgeR::DGEList(abbr_count) 
# set factor
abbr_sampleinfo$sample <- factor(abbr_sampleinfo$sample)
abbr_sampleinfo$treatment <- factor(abbr_sampleinfo$treatment)
abbr_sampleinfo$cell <- factor(abbr_sampleinfo$cell)
levels(abbr_sampleinfo$sample)
levels(abbr_sampleinfo$treatment)
levels(abbr_sampleinfo$cell)
# color
col.sample <- c("black")[abbr_sampleinfo$sample]
col.treatment <- c('#FFE66F','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000',
                   '#E0E0E0','#FF9797','#ADADAD','#005AB5','#000079','#0080FF',
                   '#0000E3')[abbr_sampleinfo$treatment]
col.cell <- c('#00DB00','#FF88C2','#0080FF', '#FF8000')[abbr_sampleinfo$cell]
# treatment(carcinogen)
plotMDS(DEG_obj,col=col.treatment,xlab = "Dimension 1",ylab = "Dimension 2")
title("Carcinogen Treatment")
legend("topleft",fill=c('#FFE66F','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000',
                        '#E0E0E0','#FF9797','#ADADAD','#005AB5','#000079','#0080FF',
                        '#0000E3'),
       legend=levels(abbr_sampleinfo$treatment))
# cell type
plotMDS(DEG_obj,col=col.cell,xlab = "Dimension 1",ylab = "Dimension 2")
title("Cell type")
legend("topleft",fill=c('#00DB00','#FF88C2','#0080FF', '#FF8000'),legend=levels(abbr_sampleinfo$cell))

### 3-4.3D PCA -------------------------------------------------------------
library(plotly)
# PCA
count_mat <- t(abbr_count)
prin_comp <- prcomp(count_mat, rank. = 3)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components <- components %>% rownames_to_column("group")
# calculate PC cover ratio
explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
pc1 <- round(100 * explained_variance_ratio[1],1)
pc2 <- round(100 * explained_variance_ratio[2],1)
pc3 <- round(100 * explained_variance_ratio[3],1)
tit <- paste0('Explained Variance PC1 to PC3 ',pc1,"%/",pc2,"%/",pc3,"%")
# set colors
groups_color <- rep(c('#FFFF37','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000','#E0E0E0',
                      '#FF9797','#ADADAD','#005AB5','#000079','#0080FF','#0000E3'),4)
clone_color <- rep(c('#F4A7B7','#FBE251','#A5DEE4','#FF1493'),each=13)
# 3D PCA plot
fig <- plotly::plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~group, 
               colors = clone_color) %>%
  plotly::add_markers(size = 12)  
fig %>% layout(title = tit, scene = list(bgcolor = "#e5ecf6"))

# another
plot3D::scatter3D(components$PC1, components$PC2, components$PC3,
                  pch = 16,xlab = pc1,ylab = pc2,zlab=pc3)

### 4-1.DEG number in each group --------------------------------------------
num <- c(3:13, 16:26, 29:39, 42:52)
# get DEG number in each groups
DEG_num <- data.frame(file = character(), up_regulation = numeric(),
                     down_regulation = numeric(), stringsAsFactors = FALSE)
for (i in name_df$file_name[num]) {
  # 讀取 Excel 文件
  data <- readxl::read_xlsx(file.path(data_path, i))
  # 過濾 M column > 1 的 row (up-regulation)
  filtered_data_greater_than_1 <- data %>% filter(M > 1) %>% 
    filter(geneBiotype=="protein_coding") %>% 
    dplyr::select(SYMBOL, M) %>% 
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) 
  up_regulation <- nrow(filtered_data_greater_than_1)
  # 過濾 M column < -1 的 row (down-regulation)
  filtered_data_less_than_minus_1 <- data %>% filter(M < -1) %>% 
    filter(geneBiotype=="protein_coding") %>% 
    dplyr::select(SYMBOL, M) %>% 
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) 
  down_regulation <- nrow(filtered_data_less_than_minus_1)
  # 將结果添加到 DEG_num
  DEG_num <- rbind(DEG_num, data.frame(file = abb(i,type = "file"), 
                                       up_regulation = up_regulation,
                                       down_regulation = down_regulation,
                                       stringsAsFactors = FALSE))
}
# 查看结果
head(DEG_num)
DEG_num$group <- name_df$abbreviate[num]
# 將數據轉成長格式
DEG_num_long <- DEG_num %>%
  mutate(down_regulation = -down_regulation) %>%
  tidyr::pivot_longer(cols = c(up_regulation, down_regulation),
                      names_to = "count_type", values_to = "count")
head(DEG_num_long)
DEG_num_long$group <- factor(DEG_num_long$group, levels = DEG_num$group)
DEG_num_long$count_type <- factor(DEG_num_long$count_type, 
                                 levels = c("up_regulation","down_regulation"))
# 雙向bar chart
ggplot(DEG_num_long, aes(x = group, y = count, fill = count_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = abs(count)), vjust = ifelse(DEG_num_long$count > 0, -0.3, 1.3), position = position_dodge(width = 0.5)) +
  scale_y_continuous(labels = abs) +  # 使 y 轴标签显示为正值
  labs(title = "DEG", x = "Short-term Groups", y = "Genes", fill = "Type") +
  scale_fill_manual(values = c("up_regulation" = "red", "down_regulation" = "darkblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### 4-2.get DEGs and plot in each groups ---------------------------------------
top <- 100
logcrit <- 1
# 初始化
DEG_df <- data.frame(num = 1:top)
for (i in name_df$file_name[num]) {
  # 讀取 Excel 文件
  data <- readxl::read_xlsx(file.path(data_path, i))
  # 过滤 M 列大于 1 的行
  filtered_data_greater_than_1 <- data %>% filter(M > logcrit) %>% 
    filter(geneBiotype=="protein_coding") %>%  
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>%
    arrange(desc((M))) %>% head(top) %>% pull(SYMBOL)
  # 过滤 M 列小于 -1 的行
  filtered_data_less_than_minus_1 <- data %>% filter(M < -logcrit) %>% 
    filter(geneBiotype=="protein_coding") %>% 
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>%
    arrange(desc((M))) %>% head(top) %>% pull(SYMBOL)
  # 補足未滿
  if (length(filtered_data_greater_than_1) < top) {
    missing_rows <- rep(NA, top - length(filtered_data_greater_than_1))
    filtered_data_greater_than_1 <- c(filtered_data_greater_than_1, missing_rows)
  }
  if (length(filtered_data_less_than_minus_1) < top) {
    missing_rows <- rep(NA, top - length(filtered_data_less_than_minus_1))
    filtered_data_less_than_minus_1 <- c(filtered_data_less_than_minus_1, missing_rows)
  }
  # 創建新的列名
  up_colname <- paste0(abb(i,"file"), "_up")
  down_colname <- paste0(abb(i,"file"), "_down")
  # 將排序后的 M 列添加到 DEG_df
  DEG_df[[up_colname]] <- filtered_data_greater_than_1
  DEG_df[[down_colname]] <- filtered_data_less_than_minus_1
}
head(DEG_df)
DEG <- DEG_df[,-1]
# openxlsx::write.xlsx(DEG,"data/DEG.xlsx")
# DEG <- readxl::read_excel("data/DEG.xlsx") %>% as.data.frame()
# 對數據框中的每一列計數
de_count_all <- table(unlist(DEG)) %>% sort(decreasing = T)
head(de_count_all,50)

### get DEG in each clone
w_up <- DEG[, grep("^W.*_up$", colnames(DEG), value = TRUE)] %>% unlist() %>% table()%>% as.data.frame()
l_up <- DEG[, grep("^L.*_up$", colnames(DEG), value = TRUE)] %>% unlist() %>% table()%>% as.data.frame()
d_up <- DEG[, grep("^D.*_up$", colnames(DEG), value = TRUE)] %>% unlist() %>% table()%>% as.data.frame()
y_up <- DEG[, grep("^Y.*_up$", colnames(DEG), value = TRUE)] %>% unlist() %>% table()%>% as.data.frame()
# count DEG and max_different in 4 clone
deg_count_df <- w_up %>% 
  full_join(l_up, by=".") %>% 
  full_join(d_up, by=".") %>% 
  full_join(y_up, by=".") %>% 
  setNames(c("gene","WT","L858R","Del19","YAP")) %>% 
  replace_na(list(WT = 0, L858R = 0, Del19 = 0, YAP = 0)) %>% 
  mutate(max_diff= pmax(WT, L858R, Del19, YAP) - pmin(WT, L858R, Del19, YAP)) %>%
  mutate(total_count=WT+L858R+Del19+YAP) %>% 
  filter(total_count > 20) %>% 
  arrange(desc(max_diff))
head(deg_count_df)
# 将数据从宽格式转换为长格式
deg_count_df$gene <- factor(deg_count_df$gene, levels = deg_count_df$gene)
df_long <- deg_count_df %>% .[1:10,] %>% 
  pivot_longer(cols = c("WT", "L858R", "Del19", "YAP"), 
               names_to = "category", 
               values_to = "value")
# 绘制条形图
ggplot(df_long, aes(x = gene, y = value, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "", x = "gene", y = "count number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = clone)
draw_bar(deg_count_df$gene[1])

### use variance to get highly variable genes in each clone
clone_var_genes <- deg_count_df %>% 
  column_to_rownames("gene") %>% .[,1:4] %>% 
  apply(., 1, var) %>% 
  sort(,decreasing = T) %>% 
  as.data.frame() %>% 
  setNames("variance") %>% rownames_to_column("gene")
deg_count_df <- deg_count_df %>% left_join(., clone_var_genes, by = "gene")
head(deg_count_df, 10)
draw_bar(clone_var_genes$gene[4])

# # DEG count df
# deg_count <- data.frame(gene=enrich_de)
# deg_count$gene.Var1 <- factor(deg_count$gene.Var1, levels = deg_count$gene.Var1)
# # 绘制 count 的柱状图
# ggplot(deg_count[1:30,], aes(x = gene.Var1, y = gene.Freq)) +
#   geom_bar(stat = "identity", fill = "skyblue") +
#   theme_minimal() +
#   labs(title = "count of DEG in all groups", x = "gene Symbol", y = "Count") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # bar plot
# draw_bar(names(enrich_de)[10],group = "ALL")

### 4-3.get all DEG list ---------------------------------------------
all_deg_list <- c()
for(i in name_df$file_name[num]){
  # 讀取 Excel 文件
  data <- readxl::read_xlsx(file.path(data_path, i))
  # 过滤 M 列大于 1 的行
  deg <- data %>% filter(abs(M) > 1) %>% 
    filter(geneBiotype=="protein_coding") %>%  
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>%
    arrange(desc(M)) %>% pull(SYMBOL)
  all_deg_list <- c(all_deg_list,deg)
}
all_deg_list <- unique(all_deg_list)
# saveRDS(all_deg_list,"data/all_deg_list.RDS")
all_deg_list <- file.path(wd,"data/all_deg_list.RDS") %>% readRDS
length(all_deg_list)

### 5.venn diagram ----------------------------------------------------------
# plot 4 clone venn diagram
AZA_venn <- multi_venn(name_df$file_name[c(3,16,29,42)],dir = "up", title = "AZA")
DAC_venn <- multi_venn(name_df$file_name[c(4,17,30,43)],dir = "up", title = "DAC")
AS_venn <- multi_venn(name_df$file_name[c(5,18,31,44)],dir = "up", title = "AS")
CO_venn <- multi_venn(name_df$file_name[c(6,19,32,45)],dir = "up", title = "CO")
LCD_venn <- multi_venn(name_df$file_name[c(7,20,33,46)],dir = "up", title = "LCD")
HCD_venn <- multi_venn(name_df$file_name[c(8,21,34,47)],dir = "up", title = "HCD")
BAP_venn <- multi_venn(name_df$file_name[c(9,22,35,48)],dir = "up", title = "BAP")
AS_BAP_venn <- multi_venn(name_df$file_name[c(10,23,36,49)],dir = "up", title = "AS_BAP")
CO_BAP_venn <- multi_venn(name_df$file_name[c(11,24,37,50)],dir = "up", title = "CO_BAP")
LCD_BAP_venn <- multi_venn(name_df$file_name[c(12,25,38,51)],dir = "up", title = "LCD_BAP")
HCD_BAP_venn <- multi_venn(name_df$file_name[c(13,26,39,52)],dir = "up", title = "HCD_BAP")
# upset plot
HCD_BAP_venn <- multi_venn(name_df$file_name[c(13,26,39,52)],dir = "up", upset = T)
# get venn list gene
print(HCD_venn)
venn_list <- AS_venn$item[[15]]
draw_from_list(venn_list, groups = "AS")
plot_MD(name_df$file_name[5], plot_type = 1,with_line = F)

# ### use divenn2.0
# get_divenn("ip_Y_V_S_CO_0_deg.xlsx",output_name = "CO.xlsx")
# get_divenn("ip_Y_V_S_CO_BAP_0_deg.xlsx",output_name = "CO_BAP.xlsx")
# # https://divenn.tch.harvard.edu/v2/

### 6.get DEG list ------------------------------------------------------------
file_name <- name_df$file_name[6]  # change
group <- "CO"
# up_all
DEG <- get_deg(file_name, log_crit = c(1,-1),dir = "up",type = "SYMBOL")
draw_from_list(list = DEG, groups = group, id = "SYMBOL", show_row_names = F)
# up_100
top <- get_deg(file_name, log_crit = c(1,-1),dir = "up",type = "SYMBOL", top = 50)
draw_from_list(list = top, groups = group, id = "SYMBOL")
# down_all
DEG <- get_deg(file_name, log_crit = c(1,-1),dir = "down",type = "SYMBOL")
draw_from_list(list = DEG, groups = group, id = "SYMBOL", show_row_names = F)
# down_100
top <- get_deg(file_name, log_crit = c(1,-1),dir = "down",type = "SYMBOL", top = 50)
draw_from_list(list = top, groups = group, id = "SYMBOL")

### 7.ORA several groups with heatmap -----------------------------------------
file_list <- c("ip_W_V_S_HCD_BAP_0_deg.xlsx", "ip_L_V_S_HCD_BAP_0_deg.xlsx", 
               "ip_D_V_S_HCD_BAP_0_deg.xlsx","ip_Y_V_S_HCD_BAP_0_deg.xlsx")
group_names <- c("WT", "L858R", "Del19","YAP")
### KEGG
title <- "HCD_BAP"
dir <- "up"
p1 <- plot_heatmap(file_list, group_names, analysis = "kegg", dir=dir, col_title = title,
                   row_title = "KEGG", top = 20)
### Reactome
p2 <- plot_heatmap(file_list, group_names, analysis = "reactome", dir=dir, 
                   row_title = "Reactome", top = 20)
p1 %v% p2

### 8.ORA pathway in single group ----------------------------------------------
### enrichKEGG(KEGG pathway)
keg_result <- run_keg_path(name_df$file_name[10], dir = "up",log_crit = 1)
view(keg_result %>% as.data.frame() %>% rownames_to_column("new_ID"))
# bar plot
barplot(keg_result, showCategory=10, font.size = 9, x = "GeneRatio", label_format = 40,
        title = "")   # change
# bar plot with qscore
mutate(keg_result, qscore = -log(p.adjust, base=10)) %>% 
  barplot(., showCategory=10, font.size = 9, x = "qscore", label_format = 40,
          title = "")   # change
# heatmap
# list <- force_list(keg_result@result$geneID[4])
# draw_from_list(list =list , groups = "CO", id = "SYMBOL")
list <- get_kegg_list(path_ID = keg_result@result$ID[6]) 
draw_from_list(list =list , groups = "CD", id = "SYMBOL",title = "")
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

### 9.GSEA analysis -----------------------------------------------------------
### KEGG
keg <- kegg_run("ip_D_V_S_HCD_0_deg.xlsx")  # change
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

### 10.cluster analysis ------------------------------------------------------
file_name <- "ip_Y_V_S_HCD_BAP_0_deg.xlsx"
group <- "only_CD"
k_value <- 5
draw_heatmap(groups = group, log_crit = 1, km = k_value,list = list)
# run cluster analysis
cluster_result <- draw_heatmap(file_name,groups = group, log_crit = 1, 
                               row_km = T, km=k_value, return_cluster = T)

# draw heatmap from 4 clone DEGs
list1 <- get_deg(name_df$file_name[13])
list2 <- get_deg(name_df$file_name[26])
list3 <- get_deg(name_df$file_name[39])
list4 <- get_deg(name_df$file_name[52])
list <- c(list1,list2,list3,list4) %>% unique()
length(list)
cluster_result <- draw_heatmap(groups = group, log_crit = 1, list = list,
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

### 11.pathway and TFs analysis using decoupleR --------------------------------
net <- "/Users/benson/Documents/project/RNA-seq1-3/data/net.RDS" %>% readRDS()
net_TF <- "/Users/benson/Documents/project/RNA-seq1-3/data/net_TF.RDS" %>% readRDS()

### pathway 
run_path_heatmap(gene_count)
# pathway for specific group
file_name <- "ip_L_V_S_AS_BAP_0_deg.xlsx"   # change
# pathway barplot
run_pathway(file_name, title = abb(file_name))
# pathway specific genes MD plot
path <- plot_pathway(file_name, "Hypoxia", title = abb(file_name))
# top50 heatmap
path_top50 <- path %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(ID)
draw_from_list(list = path_top50, groups = "CD",
               id = "SYMBOL",label_num = F,anno = T,title = "")
plot_MD(file_name, list = path_top50, plot_type = 2)

### TFs 
run_TF_heatmap(gene_count)
# TFs for specific group
file_name <- "ip_Y_V_S_HCD_BAP_0_deg.xlsx"   # change
# TFs barplot
run_TF(file_name, title = abb(file_name, type = "file"))
tf_candidate <- run_TF(file_name, title = abb(file_name, type = "file"), plot = F)
# TFs specific genes MD plot
tf <- "NFKB"
TF <- plot_TF(file_name, TF = tf, title = "regulons of")
# ggsave("TF.jpeg")
# top50 heatmap
TF_top50 <- TF %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(ID)
draw_from_list(list = TF_top50,
               groups = "CD",
               id = "SYMBOL",label_num = F,anno = T,title = "")

### target gene's TFs
target_gene <- "ANKRD1"
target_gene_TF <- net_TF %>% filter(target==target_gene) %>% pull(source)
draw_from_list(c(target_gene_TF,target_gene),groups = "CD")

### 12.RTN for TFs GSEA ---------------------------------------------------------
library(RTN)
# column annotation and interesting TFs
colAnnotation <- openxlsx::read.xlsx("/Users/benson/Documents/project/RNA-seq1-3/data/sampleinfo.xlsx",rowNames = T)
tf_candidate <- run_TF(file_name, title = abb(file_name),plot = F) %>% 
  pull(source) %>% head(6)
print(tf_candidate)
### TNI object
rtni <- RTN::tni.constructor(expData = as.matrix(abbr_count), 
                             regulatoryElements = tf_candidate, 
                             colAnnotation = colAnnotation, 
                             rowAnnotation = gene_df)

rtni <- tni.permutation(rtni, nPermutations = 100)  # Please set nPermutations >= 1000
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni)
tni.regulon.summary(rtni)
tni.regulon.summary(rtni, regulatoryElements = "JUN")
regulons <- tni.get(rtni, what = "regulons.and.mode", idkey = "SYMBOL")
head(regulons$JUN)

### TNA object and GSEA
FC_list <- get_list(file_name,type = "ENSEMBL")
head(FC_list)
DEG_list <- get_deg(file_name,type = "ENSEMBL")
head(DEG_list)
# create TNA object
rtna <- tni2tna.preprocess(object = rtni, 
                           phenotype = FC_list, 
                           hits = DEG_list, 
                           phenoIDs = gene_df)
# Run the MRA method
rtna <- tna.mra(rtna)
# Get MRA results;
mra <- tna.get(rtna, what="mra", ntop = -1)  # 'ntop = -1' will return all results, regardless of a threshold
head(mra)
### Run the GSEA method
rtna <- tna.gsea1(rtna, nPermutations=100)  # Please set nPermutations >= 1000
# Get GSEA results
gsea1 <- tna.get(rtna, what="gsea1", ntop = -1)
head(gsea1)
# Plot GSEA results
tna.plot.gsea1(rtna, labPheno="abs(log2 fold changes)", ntop = -1, 
               filepath = "./output")
### Run the GSEA-2T method
rtna <- tna.gsea2(rtna, nPermutations = 100)  # Please set nPermutations >= 1000
# Get GSEA-2T results
gsea2 <- tna.get(rtna, what = "gsea2", ntop = -1)
head(gsea2$differential)
# Plot GSEA-2T results
tna.plot.gsea2(rtna, labPheno="log2 fold changes", tfs="TP53",filepath = "./output")
for (i in tf_candidate) {
  tna.plot.gsea2(rtna, labPheno="log2 fold changes", tfs=i,filepath = "./output")
}

### regulon activity heatmap 
# Compute regulon activity for individual samples
rtni1st <- tni.gsea2(rtni, regulatoryElements = tfs)
metabric_regact <- tni.get(rtni1st, what = "regulonActivity")
# Get sample attributes from the 'rtni1st' dataset
metabric_annot <- tni.get(rtni1st, "colAnnotation")
# Get treat and clone attributes for pheatmap
treat <- c("AZA","DAC","AS","CO","LCD","HCD","BAP","AS_BAP","CO_BAP",
           "LCD_BAP","HCD_BAP")
clone <- c("W","L","D","Y")
t_annot <- metabric_annot[,treat]
c_annot <- metabric_annot[,clone]
# Plot regulon activity profiles
pheatmap(t(metabric_regact$dif), 
         main="Short-term",
         annotation_col = t_annot, 
         show_colnames = T, annotation_legend = FALSE, 
         clustering_method = "ward.D2", fontsize_row = 6,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")
pheatmap(t(metabric_regact$dif), 
         main="Short-term",
         annotation_col = c_annot, 
         show_colnames = T, annotation_legend = FALSE, 
         clustering_method = "ward.D2", fontsize_row = 6,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")

### 13.SPP ---------------------------------------------------------------------
### download data from SPP website
# http://www.signalingpathways.org/ominer/query.jsf
spp_TF <- readxl::read_excel("/Users/benson/Downloads/SRX_results_file.xlsx") %>% 
  tidyr::separate_rows(`IPAGS|BSM|Other AGS`, sep = "\\|") %>% 
  dplyr::distinct(`IPAGS|BSM|Other AGS`) %>%
  pull(`IPAGS|BSM|Other AGS`)
head(spp_TF)
# get DEG from excel
deg <- get_deg(name_df$file_name[10],dir = "all",log_crit = 1)
head(deg)
# select SPP TFs in DEG(intersect)
target_gene <- "ANKRD1"
de_spp_TF <- intersect(spp_TF, deg)
list <- c(target_gene, de_spp_TF)
list <- factor(list,levels = list)
draw_from_list(list,groups = "CD", title = paste("TFs of",target_gene),row_cluster = FALSE)
plot_MD(name_df$file_name[8], list = tfs,only_DE = F,title = paste("TFs of",target_gene),plot_type = 2)

### CollecTRI TF
target_gene <- "AREG"
target_gene_tf <- net_TF %>% filter(target==target_gene)
deg <- get_deg(file_name)
target_gene_de_tf  <- intersect(net_TF$source, deg)
draw_from_list(target_gene_de_tf,"CO")
plot_TF(file_name,"JUN",plot = T)
plot_TF(file_name,"TP53",plot = T,plot_type = 2)
plot_TF(file_name,"FOS",plot = T)
plot_MD(file_name,list="TP53",only_DE = F,plot_type = 2)

### load CollecTRI network TFs
de_TF <- intersect(net_TF$source, deg)
length(de_TF)
draw_from_list(de_all_TF,groups = "CD")
# all TFs on CollecTRI network
net_TF_list_all <- net_TF %>% distinct(source) %>% unlist()
draw_from_list(net_TF_list_all,groups = "ALL",show_row_names = F,col_cluster = T,row_km = 4)

### 14.get KEGG gene list in DEG ----------------------------------------------
# 獲取所有hsa（人類）的pathway信息
pathways_list <- KEGGREST::keggList("pathway", "hsa")
kegg_id_df <- data.frame(num=c(1:length(pathways_list)),
                         hsa = names(pathways_list), 
                         pathway = pathways_list)
# get specific pathway for our study
pathway_df <- kegg_id_df[c(90,91,92,108:114,133,164,170,268,270,291),]
pathway_id <- pathway_df$hsa[2]

# get KEGG gene list
kegg_list <- get_kegg_list(pathway_id)
pathway_name <- get_kegg_list(pathway_id, return_list = F)
list <- intersect(kegg_list, all_deg_list)
# heatmap
draw_from_list(list, groups = "ALL",title = pathway_name)
# MD plot
plot_MD(name_df$file_name[8], title = pathway_name, list = list, plot_type = 2)
# bar plot of gene expression
draw_bar("GSTM5")


### 15.TCGA regulons -----------------------------------------------------------
# ### get mutation rate in TCGA-LUAD patients 
# mutation_TCGA <- maf_object@data %>% dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>% 
#   mutate("patient" = str_sub(.$Tumor_Sample_Barcode,0,12))
# head(mutation_TCGA)
# # get number of patients 
# patient_num <- mutation_TCGA %>%
#   dplyr::distinct(patient) %>% 
#   pull(patient) %>% 
#   length()
# # count mutation rate in different patients
# muta_rate <- mutation_TCGA %>% 
#   distinct(patient,Hugo_Symbol) %>% 
#   group_by(Hugo_Symbol) %>% 
#   summarise(count = n()) %>% 
#   arrange(desc(count)) %>% 
#   mutate(mutation_ratio=count/patient_num) %>% 
#   setNames(c("ID","count","mutation_ratio"))
# muta_rate$ID <- factor(muta_rate$ID, levels = muta_rate$ID)
# head(muta_rate)
# # saveRDS(muta_rate,"muta_rate.RDS")

### load mutation rate table
muta_rate <- "/Users/benson/Documents/project/RNA-seq1-3/data/muta_rate.RDS" %>% 
  readRDS()
### filter genes we interesting 
# get TFs' regulons
tf <- "JUN"   # TFs in our study
plot_TF(name_df$file_name[13], tf, title = abb(file_name), plot_type = 1)
tf_regu_list1 <- tf_regu_df %>% arrange(desc(abs(logFC))) %>% pull(ID)
# top50 heatmap
draw_from_list(list = tf_regu_list1[1:50], groups = "CD",
               id = "SYMBOL",label_num = F,anno = T,title = "")

### get total TF's regulons
total_regulons <- net_TF %>% filter(source == tf) %>% pull(target)
length(total_regulons)
# get specific TF regulons' mutation rate
target_tfs_regulon <- muta_rate %>% filter(ID %in% total_regulons)
# 绘制 mutation_ratio 的柱状图
ggplot(target_tfs_regulon[1:10,], aes(x = ID, y = `mutation_ratio(%)`)) +
  geom_bar(stat = "identity", fill = "salmon") +
  theme_minimal() +
  labs(title = paste("regulons of",tf), x = "Gene Symbol", y = "Mutation Rate (%)")
# 绘制 count 的柱状图
ggplot(target_tfs_regulon[1:10,], aes(x = ID, y = `mutation_ratio(%)`)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = paste("regulons of",tf), x = "gene Symbol", y = "Count")


### 16.TCGA-LUAD RNA-seq and survival ---------------------------------------
### RNA-seq data 
# load TCGA-LUAD RNA-seq count
tcga_count <- "/Users/benson/Documents/project/RNA-seq1-3/data/tcga_count.RDS" %>% 
  readRDS()
# load TCGA-LUAD patient clinical data
clin.LUAD <- "/Users/benson/Documents/project/RNA-seq1-3/data/clin_LUAD.RDS" %>% 
  readRDS()
# 抓出亞洲女性與非亞洲女性id
asian_id <- clin.LUAD %>% 
  as.data.frame()%>% 
  filter(race=="asian" & gender=="female")%>% 
  pull(submitter_id)
non_asian_id <- clin.LUAD %>% 
  as.data.frame()%>% 
  filter(race!="asian" & gender=="female") %>% 
  pull(submitter_id)
### draw box plot
draw_TCGA_boxplot("RAB7B")

### survival
# load TCGA-LUAD patient data
maf_object <- "/Users/benson/Documents/project/RNA-seq1-3/data/maf_object.RDS" %>%
  readRDS()
# load TCGA-LUAD patient clinical data
clin.LUAD <- "/Users/benson/Documents/project/RNA-seq1-3/data/clin_LUAD.RDS" %>% 
  readRDS()
# TCGA-RNAseq data
tcga_count <- "/Users/benson/Documents/project/RNA-seq1-3/data/tcga_count.RDS" %>% 
  readRDS()

# survival curve of mutation (in 'output' file)
draw_muta_survival("SERPINE1")
# survival curve of gene expression (in 'output' file)
draw_TCGA_survival("CPA4", population = "ALL")

### 17.2020 cellpress -------------------------------------------------------
cell_info <- readxl::read_excel("/Users/benson/Documents/project/RNA-seq1-3/data/2020cell_info.xlsx") %>% 
  as.data.frame() %>% 
  mutate(Stage=ifelse(Stage %in% c("IA","IB"), "stage I",Stage)) %>% 
  mutate(Stage=ifelse(Stage %in% c("IIA","IIB"), "stage II",Stage)) %>% 
  mutate(Stage=ifelse(Stage %in% c("IIIA","IIIB"), "stage III",Stage)) %>%
  mutate(Stage=ifelse(Stage %in% c("IV"), "stage IV",Stage)) %>% 
  rename(Smoking_Status=`Smoking Status`)
  
cell_count <- readxl::read_excel("/Users/benson/Documents/project/RNA-seq1-3/data/2020cell.xlsx") %>% 
  as.data.frame() %>% group_by(gene) %>% summarize(across(where(is.numeric), sum)) %>% 
  column_to_rownames("gene")
cell_count$Median <- NULL

# draw boxplot by 'Gender','Smoking Status','Stage' or 'EGFR_Status'
gene <- names(enrich_de[1])
gene <- "KRT80"
print(gene)
draw_cell_boxplot(gene = gene, by="EGFR_Status")
draw_cell_boxplot(gene = gene, by="Stage")
draw_cell_boxplot(gene = gene, by="Gender")
draw_cell_boxplot(gene = gene, by="Smoking_Status")

# signature
sig <- get_deg(name_df$file_name[13]) %>% head(10)
title <- ""
draw_cell_boxplot(sig = sig, by="EGFR_Status", title = title)
draw_cell_boxplot(sig = sig, by="Stage", title = title)
draw_cell_boxplot(sig = sig, by="Gender", title = title)
draw_cell_boxplot(sig = sig, by="Smoking_Status", title = title)
# ggsave("output/1.png")

### 18.Spearman's correlation for signature ----------------------------------
### Spearman's correlation
draw_Spearman_bar(sig_num = 5, group = "AS",top = 30)
draw_Spearman_bar(sig_num = 52, top = 10, count_df = "cell")

### signature
draw_boxplot(sig_num = 10, group = "ALL", top = 20)

# 19.hallmark in cancer ---------------------------------------------------
hallmark_df <- read.csv("/Users/benson/Desktop/hallmark/W_AS_hallmark", sep = "\t")

