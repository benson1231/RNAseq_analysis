# library -----------------------------------------------------------------
library(pathfindR)
library(tidyverse)

# input and filter data
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/V"
co <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_CO_0_deg.xlsx")) %>% 
  select(SYMBOL, M) %>% 
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), mean)) %>% 
  filter(abs(M)>1) %>% 
  mutate(adj.P.Val=0.05) %>% 
  setNames(c("Gene.symbol","logFC","adj.P.Val")) %>% 
  as.data.frame()
head(co)

cobap <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_CO_BAP_0_deg.xlsx")) %>% 
  select(SYMBOL, M) %>% 
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), mean)) %>% 
  filter(abs(M)>1) %>% 
  mutate(adj.P.Val=0.05) %>% 
  setNames(c("Gene.symbol","logFC","adj.P.Val")) %>% 
  as.data.frame()
head(cobap)

# run pathfindR
# The available gene sets in pathfindR are “KEGG”, “Reactome”, “BioCarta”, 
# “GO-All”, “GO-BP”, “GO-CC” and “GO-MF”
output_df1 <- run_pathfindR(co, p_val_threshold = 0.05, gene_sets = "KEGG")
output_df2 <- run_pathfindR(cobap, p_val_threshold = 0.05, gene_sets = "KEGG")
# change color parameters
vertex_cols <- c(`Common term` = "#FCCA46", `A-only term` = "#00E3E3", 
                 `B-only term` = "#019858", `Up gene` = "#FF2D2D", 
                 `Down gene` = "#6A6AFF", `Conflicting gene` = "gray")

# combine two data sets
combined_df <- combine_pathfindR_results(
  result_A = output_df1,
  result_B = output_df2,
  plot_common = F
) 
head(combined_df,2)
combined_results_graph(combined_df, use_description = T,
                       selected_terms = combined_df$Term_Description[124]) + 
  scale_colour_manual(values = vertex_cols, name = NULL)
# subset B_only
b_only <- combined_df %>% subset(status=="B only")
combined_results_graph(combined_df, use_description = F, 
                       selected_terms = b_only$ID[1]) +
  scale_colour_manual(values = vertex_cols, name = NULL)
# specific pathway
combined_results_graph(combined_df, use_description = F, 
                       selected_terms = combined_df$ID[2]) +
  scale_colour_manual(values = vertex_cols, name = NULL)

# out put data
combined_results_graph(
  combined_df, 
  use_description = T,
  selected_terms = combined_df$Term_Description[1:4]
) + scale_colour_manual(values = vertex_cols, name = NULL) +
  ggtitle("")
