### library -----------------------------------------------------------------
library(clusterProfiler)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)

source("RNAseq_function.R")

### data path ---------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/V"


# function ----------------------------------------------------------------
get_path <- function(file_name, group_name){
  genelist <- get_list(file_name)
  gene <- names(genelist)[genelist > 1]
  kk <- enrichKEGG(gene         = gene,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  df <- kk@result[,c("Description","pvalue")] %>% 
    filter(pvalue < 0.01) %>% 
    setNames(c("Description",group_name))
  return(df)
}

df1 <- get_path("ip_Y_V_S_CO_0_deg.xlsx", group_name = "CO")
df2 <- get_path("ip_Y_V_S_BAP_0_deg.xlsx", group_name = "BAP")
df3 <- get_path("ip_Y_V_S_CO_BAP_0_deg.xlsx", group_name = "CO_BAP")
df <- df1 %>% full_join(df2,"Description") %>% full_join(df3,"Description") %>% 
  column_to_rownames("Description") %>% as.matrix() 
mat <- log10(df) %>% negative()

### heatmap
col_fun = colorRamp2(c(0, 5), c("white", "red"))
ComplexHeatmap::Heatmap(mat,na_col = "grey",cluster_rows = F,cluster_columns = F,
                        col = col_fun, name = "-log10(P)" ,
                        row_names_gp = gpar(fontsize = 7), 
                        column_names_gp = gpar(fontsize = 10), 
                        column_names_rot = 45)


# details -----------------------------------------------------------------
### up1
genelist <- get_list("ip_Y_V_S_BAP_0_deg.xlsx")
gene <- names(genelist)[genelist > 1]
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
df1 <- kk@result[,c(4,8)] %>% filter(p.adjust<0.05) %>% setNames(c("Description","BAP"))
### up2
genelist <- get_list("ip_Y_V_S_CO_0_deg.xlsx")
gene <- names(genelist)[genelist > 1]
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
df2 <- kk@result[,c(4,8)] %>% filter(p.adjust<0.05) %>% setNames(c("Description","HCD"))
### up3
genelist <- get_list("ip_Y_V_S_CO_BAP_0_deg.xlsx")
gene <- names(genelist)[genelist > 1]
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
df3 <- kk@result[,c(4,8)] %>% filter(p.adjust<0.05) %>% setNames(c("Description","HCD_BAP"))

### merge
df <- df1 %>% full_join(df2,"Description") %>% full_join(df3,"Description") %>% 
  column_to_rownames("Description") %>% as.matrix()
mat <- log10(df) %>% negative()

### heatmap
col_fun = colorRamp2(c(0, 5), c("white", "red"))
ComplexHeatmap::Heatmap(mat,na_col = "grey",cluster_rows = F,cluster_columns = F,
                        col = col_fun, name = "-log10(P)" ,
                        row_names_gp = gpar(fontsize = 7), 
                        column_names_gp = gpar(fontsize = 10), 
                        column_names_rot = 45)
