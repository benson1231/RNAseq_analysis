### library -----------------------------------------------------------------
library(clusterProfiler)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)

source("RNAseq_function.R")

### pre-load ---------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/V"
col_fun = colorRamp2(c(0, 5), c("white", "red"))

### function ----------------------------------------------------------------
get_path <- function(file_name, group_name, dir="up", log_crit=1, return_keg_result = F){
  if (!(dir %in% c("up","down"))) {
    stop("Invalid type. Allowed values are 'up' or 'down'.")
  }
  genelist <- get_list(file_name)
  if(dir=="up"){
    gene <- names(genelist)[genelist >= log_crit]
  } else{
    gene <- names(genelist)[genelist <= negative(log_crit)]
  }
  # run kegg ORA
  kk <- enrichKEGG(gene         = gene,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  # ENTREZID to SYMBOL
  net <- DOSE::setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
  if(return_keg_result==T){
    return(net@result)
  }
  # get pvalue
  df <- net@result[,c("Description","pvalue")] %>% 
    filter(pvalue < 0.01) %>% 
    setNames(c("Description",group_name))
  
  cat(c(" -> '",group_name,dir,"'regulation","-> log2FC criteria:",log_crit,"\n"))
  return(df)
}
### get kegg result data
keg_df <- get_path("ip_Y_V_S_HCD_0_deg.xlsx", group_name = "HCD",return_keg_result = T)

### run kegg and get pvalue
df1 <- get_path("ip_Y_V_S_CO_0_deg.xlsx", group_name = "CO")
df2 <- get_path("ip_Y_V_S_BAP_0_deg.xlsx", group_name = "BAP")
df3 <- get_path("ip_Y_V_S_CO_BAP_0_deg.xlsx", group_name = "CO_BAP")
combined_df <- df1 %>% 
  full_join(df2,"Description") %>% 
  full_join(df3,"Description") %>% 
  column_to_rownames("Description") %>% as.matrix() 
# -log(P)
mat <- log10(combined_df) %>% negative()
# heatmap
ComplexHeatmap::Heatmap(mat,na_col = "grey",cluster_rows = F,cluster_columns = F,
                        col = col_fun, name = "-log10(P)" ,
                        row_names_gp = gpar(fontsize = 7), 
                        column_names_gp = gpar(fontsize = 8), 
                        column_names_rot = 45,column_title = "")


# details -----------------------------------------------------------------
### up1
genelist <- get_list("ip_Y_V_S_BAP_0_deg.xlsx")
gene <- names(genelist)[genelist > 1]
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
net <- DOSE::setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
df1 <- net@result[,c(4,8)] %>% filter(p.adjust<0.05) %>% setNames(c("Description","BAP"))
net <- DOSE::setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
### up2
genelist <- get_list("ip_Y_V_S_CO_0_deg.xlsx")
gene <- names(genelist)[genelist > 1]
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
net <- DOSE::setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
df2 <- net@result[,c(4,8)] %>% filter(p.adjust<0.05) %>% setNames(c("Description","HCD"))
### up3
genelist <- get_list("ip_Y_V_S_CO_BAP_0_deg.xlsx")
gene <- names(genelist)[genelist > 1]
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
df3 <- net@result[,c(4,8)] %>% filter(p.adjust<0.05) %>% setNames(c("Description","HCD_BAP"))

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
                        column_names_rot = 45,column_title = "")
