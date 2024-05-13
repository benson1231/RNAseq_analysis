### We load the required packages -------------------------------------------
library(decoupleR)
library(OmnipathR)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(ggrepel)
library(colorRamp2)

source("RNAseq_function.R")

### count data input --------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/Vitro"
# annotation df
gene_df <- "/Users/benson/Documents/project/RNA-seq1-3/anno_gene.RDS" %>% 
  readRDS() %>% select(ENSEMBL,SYMBOL,ENTREZID)
anno_df <- "/Users/benson/Documents/project/RNA-seq1-3/anno_gene.RDS" %>% 
  readRDS()
# groups
group_names <- c("ip_L_V_L_CON", "ip_L_V_L_DMS", "ip_L_V_L_AZA", "ip_L_V_L_DAC", "ip_L_V_L_AS",
                 "ip_L_V_L_CO", "ip_L_V_L_LCD", "ip_L_V_L_HCD", "ip_L_V_L_BAP", "ip_L_V_L_AS_BAP",
                 "ip_L_V_L_CO_BAP", "ip_L_V_L_LCD_BAP", "ip_L_V_L_HCD_BAP", "ip_Y_V_S_CON",
                 "ip_Y_V_S_DMS", "ip_Y_V_S_AZA", "ip_Y_V_S_DAC", "ip_Y_V_S_AS", "ip_Y_V_S_CO",
                 "ip_Y_V_S_LCD", "ip_Y_V_S_HCD", "ip_Y_V_S_BAP", "ip_Y_V_S_AS_BAP", "ip_Y_V_S_CO_BAP",
                 "ip_Y_V_S_LCD_BAP", "ip_Y_V_S_HCD_BAP")
# count df
mycount_77 <- "/Users/benson/Documents/project/RNA-seq1-3/myTMM.RDS" %>% 
  readRDS() %>% as.data.frame()
mycount_df <- mycount_77 %>% 
  select(.,all_of(group_names))
count <- mycount_df %>% 
  rownames_to_column("ENSEMBL") %>% 
  left_join(gene_df,"ENSEMBL") %>% 
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL")
# design
design <- colnames(mycount_df) %>% as.data.frame() %>% setNames("design")
design

### pathway heatmap ---------------------------------------------------------
### load pathway data
# net <- get_progeny(organism = 'human', top = 500)
# saveRDS(net,"net.RDS")
net <- "/Users/benson/Documents/project/RNA-seq1-3/net.RDS" %>% readRDS()
net
sample_acts <- run_mlm(mat=count, net=net, .source='source', .target='target',
                       .mor='weight', minsize = 5)
sample_acts
# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Scale per feature
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

### heatmap
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 



### two groups comparing pathway --------------------------------------------
### read file
file_name <- "ip_Y_V_S_HCD_BAP_0_deg.xlsx"   # change
deg_FC <- get_df(file_name ,de=T,log_crit = c(1,-1)) %>% 
  group_by(gene) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "gene")
deg_df <- get_df(file_name ,all = T) %>% group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL") %>% 
  rownames_to_column("ID") %>% 
  .[,c(1,5,6,7,8)] %>% 
  setNames(c("ID","treat","control","M","D"))%>% 
  mutate(Average_expression =(treat+control)/2) %>% 
  .[,c(1,4,6)] %>% 
  setNames(c("ID","logFC","Average_expression"))

### run MLM analysis
contrast_acts <- run_mlm(mat=deg_FC, net=net, .source='source', .target='target',
                         .mor='weight', minsize = 5)
contrast_acts

### pathway enrich plot
ggplot(contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") + ggtitle("Y-ST-Co_BAP")   # change title

### Specific pathway related genes
pathway <- "EGFR"   # change
logFC_criteria <- 1   # change

### filter pathway specific genes
path <- net %>%
  filter(source == pathway) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target') %>% 
  left_join(deg_df,"ID") %>% 
  mutate(color = if_else(logFC > logFC_criteria, '1', color)) %>%
  mutate(color = if_else(logFC < logFC_criteria & logFC > negative(logFC_criteria), '2', color)) %>%   
  mutate(color = if_else(logFC < negative(logFC_criteria), '3', color)) 
### MD plot
ggplot(path, aes(x = log(Average_expression), y = logFC, color = color)) + geom_point() +
  scale_colour_manual(values = c("red","grey","royalblue3")) +
  geom_label_repel(aes(label = ID)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(pathway)
### top50 heatmap
path_top50 <- path %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(ID)
draw_from_list(list = path_top50,
               groups = "ALL",
               id = "SYMBOL",label_num = F,anno = T,title = "")


### TFs heatmap -------------------------------------------------------------
### load pathway data
# net_TF <- get_collectri(organism='human', split_complexes=FALSE)
# saveRDS(net_TF,"net_TF.RDS")
net_TF <- "/Users/benson/Documents/project/RNA-seq1-3/net_TF.RDS" %>% readRDS()
net_TF
sample_TF <- run_ulm(mat=count, net=net_TF, .source='source', .target='target',
                       .mor='mor', minsize = 5)
sample_TF

# plot
n_tfs <- 25
# Transform to wide matrix
sample_TF_mat <- sample_TF %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- sample_TF %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
sample_TF_mat <- sample_TF_mat[,tfs]

# Scale per sample
sample_TF_mat <- scale(sample_TF_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

### heatmap
pheatmap(sample_TF_mat, border_color = NA, color=my_color, breaks = my_breaks) 


### two groups comparing TFs ------------------------------------------------
file_name <- "ip_Y_V_S_HCD_0_deg.xlsx"
de_FC <- get_df(file_name ,de=T) %>% group_by(gene) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "gene")
de <- get_df(file_name ,all = T) %>% group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL") %>% 
  rownames_to_column("ID") %>% 
  .[,c(1,5,6,7,8)] %>% 
  setNames(c("ID","treat","control","M","D"))%>% 
  mutate(Average_expression =(treat+control)/2) %>% 
  .[,c(1,4,6)] %>% setNames(c("ID","logFC","Average_expression"))
### Run ULM
contrast_TF <- run_ulm(mat=de_FC[, "logFC", drop=FALSE], net=net_TF, .source='source', .target='target',
                         .mor='mor', minsize = 5)
contrast_TF

# Filter top TFs in both signs
f_contrast_TF <- contrast_TF %>%
  mutate(rnk = NA)
msk <- f_contrast_TF$score > 0
f_contrast_TF[msk, 'rnk'] <- rank(-f_contrast_TF[msk, 'score'])
f_contrast_TF[!msk, 'rnk'] <- rank(-abs(f_contrast_TF[!msk, 'score']))
tfs <- f_contrast_TF %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_TF <- f_contrast_TF %>%
  filter(source %in% tfs)

### TF enrich plot
ggplot(f_contrast_TF, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("TFs") + ggtitle("")   # change

### Specific TF related genes -----------------------------------------------
tf <- 'HIF1A'   # change
logFC_criteria <- 1    # change

# filter TF related genes
df <- net_TF %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target') %>% 
  left_join(de,"ID") %>% 
  mutate(color = if_else(logFC > logFC_criteria, '1', color)) %>%
  mutate(color = if_else(logFC < logFC_criteria & logFC > negative(logFC_criteria), '2', color)) %>%   
  mutate(color = if_else(logFC < negative(logFC_criteria), '3', color)) 

### MD plot
ggplot(df, aes(x = log(Average_expression), y = logFC, color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","grey","royalblue3")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(tf)
### top50 heatmap
TF_top50 <- df %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(ID)
draw_from_list(list = TF_top50,
               groups = "ALL",
               id = "SYMBOL",label_num = F,anno = T,title = "")
