# up ----------------------------------------------------------------------
BAP_up <- get_deg(data_path = data_path,file = "ip_Y_V_S_BAP_0_deg.xlsx",
                    log_crit = c(1,-1),dir = "up")
LCD_up <- get_deg(data_path = data_path,file = "ip_Y_V_S_LCD_0_deg.xlsx",
                   log_crit = c(1,-1),dir = "up")
LCD_BAP_up <- get_deg(data_path = data_path,file = "ip_Y_V_S_LCD_BAP_0_deg.xlsx",
                       log_crit = c(1,-1),dir = "up")
LCD_gene <- list(LCD = LCD_up,
                BAP = BAP_up,
                LCD_BAP =LCD_BAP_up
                )
 
ggVennDiagram(LCD_gene,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#FF2D2D")

LCD_list <- process_region_data(Venn(LCD_gene))
draw_from_list(list = LCD_list$item[[7]],groups = "CD",id = "ENSEMBL")

# down ----------------------------------------------------------------------
BAP_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_BAP_0_deg.xlsx",
                    log_crit = c(1,-1),dir = "down")
LCD_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_LCD_0_deg.xlsx",
                   log_crit = c(1,-1),dir = "down")
LCD_BAP_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_LCD_BAP_0_deg.xlsx",
                       log_crit = c(1,-1),dir = "down")
LCD_gene <- list(LCD = LCD_down,
                BAP = BAP_down,
                LCD_BAP =LCD_BAP_down)

ggVennDiagram(LCD_gene,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#6A6AFF")

LCD_list <- process_region_data(Venn(LCD_gene))
draw_from_list(list = LCD_list$item[[7]],groups = "CD",id = "ENSEMBL")

# metal -------------------------------------------------------------------
AS_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_AS_0_deg.xlsx",
                  log_crit = c(1,-1),dir = "down")
CO_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_CO_0_deg.xlsx",
                  log_crit = c(1,-1),dir = "down")
LCD_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_LCD_0_deg.xlsx",
                      log_crit = c(1,-1),dir = "down")
HCD_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_HCD_0_deg.xlsx",
                  log_crit = c(1,-1),dir = "down")
BAP_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_BAP_0_deg.xlsx",
                  log_crit = c(1,-1),dir = "down")
gene <- list(AS = AS_down,
             CO = CO_down,
             LCD =LCD_down,
             HCD = HCD_down,
             BAP = BAP_down
)

ggVennDiagram(gene,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#6A6AFF")

LCD_list <- process_region_data(Venn(gene))
draw_from_list(list = LCD_list$item[[7]],grodowns = "CD",id = "ENSEMBL")
