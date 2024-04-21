# library -----------------------------------------------------------------
library(pathfindR)
library(tidyverse)

# get_pathfind function ---------------------------------------------------
get_pathfind <- function(file,
                         log_crit=1
                         ){
  df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
    select(SYMBOL, M) %>% 
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), mean)) %>% 
    filter(abs(M)>1) %>% 
    mutate(adj.P.Val=0.05) %>% 
    setNames(c("Gene.symbol","logFC","adj.P.Val")) %>% 
    as.data.frame()
  return(df)
}

# input and filter data ---------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/V"

AS_df <- get_pathfind("ip_Y_V_S_AS_0_deg.xlsx")
CO_df <- get_pathfind("ip_Y_V_S_CO_0_deg.xlsx")
LCD_df <- get_pathfind("ip_Y_V_S_LCD_0_deg.xlsx")
HCD_df <- get_pathfind("ip_Y_V_S_HCD_0_deg.xlsx")
BAP_df <- get_pathfind("ip_Y_V_S_BAP_0_deg.xlsx")
AS_BAP_df <- get_pathfind("ip_Y_V_S_AS_BAP_0_deg.xlsx")
CO_BAP_df <- get_pathfind("ip_Y_V_S_CO_BAP_0_deg.xlsx")
LCD_BAP_df <- get_pathfind("ip_Y_V_S_LCD_BAP_0_deg.xlsx")
HCD_BAP_df <- get_pathfind("ip_Y_V_S_HCD_BAP_0_deg.xlsx")

# run pathfindR -----------------------------------------------------------
# The available gene sets in pathfindR are “KEGG”, “Reactome”, “BioCarta”, 
# “GO-All”, “GO-BP”, “GO-CC” and “GO-MF”
# pin_name_path argument can be one of “Biogrid”, “STRING”, “GeneMania”, 
# “IntAct”, “KEGG”, “mmu_STRING”
output_df1 <- run_pathfindR(CO_df, p_val_threshold = 0.05, gene_sets = "KEGG")
output_df2 <- run_pathfindR(BAP_df, p_val_threshold = 0.05, gene_sets = "KEGG")
output_df3 <- run_pathfindR(CO_BAP_df, p_val_threshold = 0.05, gene_sets = "KEGG")

### enrichment dot plot
enrichment_chart(
  result_df = output_df1,
  top_terms = 15
)
# change color parameters
vertex_cols <- c(`Common term` = "#FCCA46", `A-only term` = "#00E3E3", 
                 `B-only term` = "#019858", `Up gene` = "#FF2D2D", 
                 `Down gene` = "#6A6AFF", `Conflicting gene` = "gray")

### combine two data sets
combined_df <- combine_pathfindR_results(
  result_A = output_df1,
  result_B = output_df3,
  plot_common = F
)
# default plot
combined_results_graph(
  combined_df, 
  use_description = T
) + scale_colour_manual(values = vertex_cols, name = NULL)
# select specific Term_Description
combined_results_graph(
  combined_df, 
  use_description = T,
  selected_terms = combined_df$Term_Description[85]
) + scale_colour_manual(values = vertex_cols, name = NULL)

### subset B_only
b_only <- combined_df %>% subset(status=="B only")
combined_results_graph(
  combined_df, 
  use_description = T
) + scale_colour_manual(values = vertex_cols, name = NULL)
# select specific Term_Description
combined_results_graph(
  combined_df, 
  use_description = F, 
  selected_terms = combined_df$ID[2]
) + scale_colour_manual(values = vertex_cols, name = NULL)

# cluster analysis --------------------------------------------------------
# The available gene sets in pathfindR are “KEGG”, “Reactome”, “BioCarta”, 
# “GO-All”, “GO-BP”, “GO-CC” and “GO-MF”
# pin_name_path argument can be one of “Biogrid”, “STRING”, “GeneMania”, 
# “IntAct”, “KEGG”, “mmu_STRING”
output_df <- run_pathfindR(HCD_df, p_val_threshold = 0.05, gene_sets = "GO-All")
### Hierarchical Clustering
output_clustered <- cluster_enriched_terms(output_df,
                                           plot_dend = FALSE, 
                                           plot_clusters_graph = FALSE
                                           )
# plotting only selected clusters for better visualization
selected_clusters <- subset(output_clustered, Cluster %in% 1:5)
enrichment_chart(selected_clusters, plot_by_cluster = TRUE)
