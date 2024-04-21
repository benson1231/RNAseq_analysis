# library -----------------------------------------------------------------
library(ggVennDiagram)
library(COmplexHeatmap)
library(tidyverse)
library(openxlsx)

source("RNAseq_function.R")

# data pathway ------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/V"

# up ----------------------------------------------------------------------
BAP_up <- get_deg(file = "ip_Y_V_S_BAP_0_deg.xlsx",
                  log_crit = c(1,-1),dir = "up")
CO_up <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
                   log_crit = c(1,-1),dir = "up")
CO_BAP_up <- get_deg(file = "ip_Y_V_S_CO_BAP_0_deg.xlsx",
                       log_crit = c(1,-1),dir = "up")
CO_gene <- list(CO = CO_up,
                BAP = BAP_up,
                CO_BAP =CO_BAP_up
                )
 
ggVennDiagram(CO_gene,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#FF2D2D")

CO_list <- process_region_data(Venn(CO_gene))
# draw heatmap
draw_from_list(list = CO_list$item[[3]],groups = "CO",id = "ENSEMBL")

# 將venn diagram結果輸出成excel
venn_to_excel(CO_list, name = "CO_up")

# down ----------------------------------------------------------------------
BAP_down <- get_deg(file = "ip_Y_V_S_BAP_0_deg.xlsx",
                    log_crit = c(1,-1),dir = "down")
CO_down <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
                   log_crit = c(1,-1),dir = "down")
CO_BAP_down <- get_deg(file = "ip_Y_V_S_CO_BAP_0_deg.xlsx",
                       log_crit = c(1,-1),dir = "down")
CO_gene <- list(CO = CO_down,
                BAP = BAP_down,
                CO_BAP =CO_BAP_down)

ggVennDiagram(CO_gene,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#6A6AFF")

CO_list <- process_region_data(Venn(CO_gene))
# draw heatmap
draw_from_list(list = CO_list$item[[3]],groups = "CD",id = "ENSEMBL")
# 將venn diagram結果輸出成excel
venn_to_excel(CO_list, name = "CO_down")


# # metal -------------------------------------------------------------------
# # up
# CO_up <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
#                   log_crit = c(1,-1),dir = "up")
# CO_up <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
#                   log_crit = c(1,-1),dir = "up")
# CO_up <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
#                       log_crit = c(1,-1),dir = "up")
# CO_up <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
#                   log_crit = c(1,-1),dir = "up")
# BAP_up <- get_deg(file = "ip_Y_V_S_BAP_0_deg.xlsx",
#                   log_crit = c(1,-1),dir = "up")
# gene <- list(CO = CO_up,
#              CO = CO_up,
#              CO =CO_up,
#              CO = CO_up,
#              BAP = BAP_up
# )
# 
# ggVennDiagram(gene,label_percent_digit = 1,label_alpha = 0, force_upset = T)
# # down
# CO_down <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
#                   log_crit = c(1,-1),dir = "down")
# CO_down <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
#                   log_crit = c(1,-1),dir = "down")
# CO_down <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
#                       log_crit = c(1,-1),dir = "down")
# CO_down <- get_deg(file = "ip_Y_V_S_CO_0_deg.xlsx",
#                   log_crit = c(1,-1),dir = "down")
# BAP_down <- get_deg(file = "ip_Y_V_S_BAP_0_deg.xlsx",
#                   log_crit = c(1,-1),dir = "down")
# gene <- list(CO = CO_down,
#              CO = CO_down,
#              CO =CO_down,
#              CO = CO_down,
#              BAP = BAP_down
# )
# 
# ggVennDiagram(gene,label_percent_digit = 1,label_alpha = 0, force_upset = T)
# 
# CO_list <- process_region_data(Venn(gene))
# draw_from_list(list = CO_list$item[[7]],groups = "CD",id = "ENSEMBL")
