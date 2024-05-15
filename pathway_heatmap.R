### library -----------------------------------------------------------------
library(clusterProfiler)
library(pathfindR)
library(enrichplot)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)

source("RNAseq_function.R")

### pre-load ---------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/V"
col_fun = colorRamp2(c(0, 5), c("white", "red"))

### run KEGG pathway by clusterProfiler's enrichKEGG  -----------------------
keg_result <- run_keg_path("ip_Y_V_S_HCD_BAP_0_deg.xlsx")
### bar plot
barplot(keg_result, showCategory=20, font.size = 9, x = "GeneRatio", label_format = 40,
        title = "")   # change

### compare groups and plot heatmap
file_list <- c("ip_Y_V_S_CO_0_deg.xlsx", "ip_Y_V_S_BAP_0_deg.xlsx", "ip_Y_V_S_CO_BAP_0_deg.xlsx")
group_names <- c("CO", "BAP", "CO_BAP")
plot_heatmap(file_list, group_names, analysis = "kegg",title = "")

### run reactome by enrichPathway -----------------------------------------
### Reactome pathway over-representation analysis
CO <- run_reactome("ip_Y_V_S_CO_0_deg.xlsx",dir = "up")
### bar plot
barplot(CO, showCategory=20, font.size = 9, x = "GeneRatio", label_format = 40,
        title = "HCD") 

### compare groups and plot heatmap
file_list <- c("ip_Y_V_S_CO_0_deg.xlsx", "ip_Y_V_S_BAP_0_deg.xlsx", "ip_Y_V_S_CO_BAP_0_deg.xlsx")
group_names <- c("CO", "BAP", "CO_BAP")
plot_heatmap(file_list, group_names, analysis = "reactome",title = "")

### run pathfind ------------------------------------------------------------
run_pathfind <- function(file,
                         log_crit = c(1,-1),
                         dir="all",  # all/up/down
                         db="KEGG" # “KEGG”, “Reactome”, “BioCarta”, “GO-All”, “GO-BP”, “GO-CC”, “GO-MF”, “cell_markers”
){
  if (!exists("run_pathfindR")) {
    stop("run_pathfindR function is not available.")
  }
  df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
    filter(.,!(SYMBOL%in% c("havana","ensembl_havana","havana_tagene")),
           geneBiotype=="protein_coding") %>% 
    select(SYMBOL, M) %>% 
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), mean)) %>% 
    mutate(adj.P.Val=0.05) %>% 
    setNames(c("Gene.symbol","logFC","adj.P.Val")) %>% 
    as.data.frame()
  
  if(dir=="all"){
    df <- df %>% 
      filter(logFC >log_crit[1]| logFC < log_crit[2]) 
    cat(c(" -> abs(log2FC) larger than", log_crit[1],"->",nrow(df),"genes.\n"))
  } else if(dir=="up"){
    df <- df %>% 
      filter(logFC > log_crit[1])
    cat(c(" -> log2FC larger than", log_crit[1],"->",nrow(df),"genes.\n"))
  } else if(dir=="down"){
    df <- df %>% 
      filter(logFC < log_crit[2])
    cat(c(" -> log2FC smaller than", log_crit[2],"->",nrow(df),"genes.\n"))
  }
  # run pathfind
  output_df <- pathfindR::run_pathfindR(df, p_val_threshold = 0.05, gene_sets = db)
  cat(c(" ->",db, "analysis completed"))
  return(output_df)
}

CO_df <- run_pathfind("ip_Y_V_S_CO_0_deg.xlsx",dir = "up")
BAP_df <- run_pathfind("ip_Y_V_S_BAP_0_deg.xlsx",dir = "up")
CO_BAP_df <- run_pathfind("ip_Y_V_S_CO_BAP_0_deg.xlsx",dir = "up")

df1 <- CO_df[,c(2,6)] %>% head(10) %>% setNames(c("Description","CO"))
df2 <- BAP_df[,c(2,6)] %>% head(10) %>% setNames(c("Description","BAP"))
df3 <- CO_BAP_df[,c(2,6)] %>% head(10) %>% setNames(c("Description","CO_BAP"))

combined_df <- df1 %>% 
  full_join(df2,"Description") %>% 
  full_join(df3,"Description") %>% 
  column_to_rownames("Description") %>% as.matrix()
mat <- log10(combined_df) %>% negative()

### heatmap
col_fun = colorRamp2(c(0, 5), c("white", "red"))
ComplexHeatmap::Heatmap(mat,na_col = "grey",cluster_rows = F,cluster_columns = F,
                        col = col_fun, name = "-log10(P)" ,
                        row_names_gp = gpar(fontsize = 8), 
                        column_names_gp = gpar(fontsize = 10), 
                        column_names_rot = 45,column_title = "")


