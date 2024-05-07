## We load the required packages
library(decoupleR)
library(OmnipathR)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/pw_bk.html

inputs_dir <- system.file("extdata", package = "decoupleR")
data <- readRDS(file.path(inputs_dir, "bk_data.rds"))

# input -------------------------------------------------------------------
count <- mycount_df %>% rownames_to_column("ENSEMBL") %>% left_join(gene_df,"ENSEMBL") %>% 
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL")

design <- colnames(mycount_df) %>% as.data.frame() %>% setNames("design")
design


net <- get_progeny(organism = 'human', top = 500)
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
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 


# new ---------------------------------------------------------------------
# example
deg <- get_df("ip_Y_V_S_CO_0_deg.xlsx",de=T) %>% group_by(gene) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "gene")

contrast_acts <- run_mlm(mat=deg, net=net, .source='source', .target='target',
                         .mor='weight', minsize = 5)
contrast_acts
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
  xlab("Pathways")

pathway <- 'p53'

df <- net %>%
  filter(source == pathway) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')
inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[c('logFC')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(weight > 0 & logFC > 0, '1', color)) %>%
  mutate(color = if_else(weight > 0 & logFC < 0, '2', color)) %>%
  mutate(color = if_else(weight < 0 & logFC > 0, '2', color)) %>%
  mutate(color = if_else(weight < 0 & logFC < 0, '1', color))

ggplot(df, aes(x = weight, y = logFC, color = color)) + geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(pathway)

# TF ----------------------------------------------------------------------
de <- get_df("ip_Y_V_S_CO_0_deg.xlsx",all = T) %>% group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames(., var = "SYMBOL") %>% .[,6:7]
net_TF <- get_collectri(organism='human', split_complexes=FALSE)
net_TF
sample_TF <- run_ulm(mat=count, net=net_TF, .source='source', .target='target',
                       .mor='mor', minsize = 5)
sample_TF

### plot
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

# Plot
pheatmap(sample_TF_mat, border_color = NA, color=my_color, breaks = my_breaks) 

### logFC
# Run ulm
contrast_TF <- run_ulm(mat=deg[, "logFC", drop=FALSE], net=net_TF, .source='source', .target='target',
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

# Plot
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
  xlab("TFs")

# 
tf <- 'HIF1A'

df <- net_TF %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(de),rownames(df)))
df <- df[inter, ]
df[,'logfc'] <- de[inter, ]$M
df[,'D'] <- de[inter, ]$D
df <- df %>%
  mutate(color = if_else(mor > 0 & logfc > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & logfc < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & logfc > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & logfc < 0, '1', color))

ggplot(df, aes(x = logfc, y = log(D), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(tf)
