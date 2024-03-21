# library -----------------------------------------------------------------
library(magrittr)
library(tidyverse)
library(stringr)
library(ComplexHeatmap)

# load raw reads counts ---------------------------------------------------------------
raw_counts_df <- read.csv("/Users/benson/Documents/project/RNA-seq1/mycounts_f.txt")

# df<- raw_counts_df[,grepl("ip_L_V_L|ip_Y_V_S", colnames(mycount_df))]
mycount_df <- raw_counts_df %>% select(.,ip_L_V_L_CON,ip_L_V_L_DMS,ip_L_V_L_AZA,ip_L_V_L_DAC,
                            ip_L_V_L_AS,ip_L_V_L_CO,ip_L_V_L_LCD,ip_L_V_L_HCD,
                            ip_L_V_L_BAP, ip_L_V_L_AS_BAP, ip_L_V_L_CO_BAP,
                            ip_L_V_L_LCD_BAP,ip_L_V_L_HCD_BAP,
                            ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
                            ip_Y_V_S_AS,ip_Y_V_S_CO,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
                            ip_Y_V_S_BAP, ip_Y_V_S_AS_BAP, ip_Y_V_S_CO_BAP,
                            ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP) 
rownames(mycount_df) <- raw_counts_df$X

df_log <- mycount_df %>% na.omit() %>% t() %>% scale(scale = T) %>% t() %>% na.omit()
table(is.na(df_log))

# heatmap with annotation ---------------------------------------------------------
col <- colnames(mycount_df)

# agent name
agent <- factor (
  str_replace_all(col, c("ip_L_V_L_|ip_Y_V_S_"='')),
  levels=c('CON','DMS',"AZA","DAC",'AS',"CO","LCD","HCD","BAP",
           "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))

# clone name
clone <- factor( 
  str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                       "D" = "DEL19", "Y" = "YAP")),
  levels=c('WT','L858R',"DEL19","YAP"))

# draw hp
ha <- HeatmapAnnotation(agent = agent, clone = clone,
                       col = list(agent=c('CON'='#E0E0E0', 'DMS'='#ADADAD', 'AZA'='#FFD2D2', 'DAC'='#FF9797',
                                          'AS'='#FFFF37','CO'='#FF5151','LCD'='#0080FF',
                                          'HCD'='#005AB5', 'BAP'='#00DB00', 'AS_BAP'='#FFDC35', 'CO_BAP'= '#EA0000',
                                          'LCD_BAP'= '#0000E3', "HCD_BAP"="#000079"),
                                  clone=c('WT'="#FF2D2D",'L858R'="#FF9224",
                                          "DEL19"= "#66B3FF","YAP"="#2828FF")))

ComplexHeatmap::Heatmap(df_log, top_annotation = ha, cluster_columns = F, 
                        show_row_names = F, show_column_names = F,
                        name = "Z-score"
                        )

# DEG heatmap -----------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1/raw/Y_L_V"
co <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_CO_0_deg.xlsx"))
co <- co %>% filter(abs(M)>4)
list <- co$ENSEMBL
gene_df <- co[,c("ENSEMBL","SYMBOL")] 
mycount_df$ENSEMBL <- rownames(mycount_df)
new_df <- mycount_df %>% left_join(.,gene_df, by="ENSEMBL")
mat <- mycount_df %>% filter(.,rownames(mycount_df) %in% list) %>% 
  left_join(.,gene_df, by="ENSEMBL") %>% group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% 
  column_to_rownames(., var = "SYMBOL") %>% 
  select(ip_L_V_L_CON,ip_L_V_L_DMS,ip_L_V_L_CO,ip_L_V_L_BAP,ip_L_V_L_CO_BAP)

mat_scale <- mat %>% t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()


col <- colnames(mat1)
# agent name
agent <- factor (
  str_replace_all(col, c("ip_L_V_L_|ip_Y_V_S_"='')),
  levels=c('CON','DMS',"AZA","DAC",'AS',"CO","LCD","HCD","BAP",
           "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))

# clone name
clone <- factor( 
  
  str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                       "D" = "DEL19", "Y" = "YAP")),
  
  levels=c('WT','L858R',"DEL19","YAP"))

# draw hp
ha <- HeatmapAnnotation(agent = agent, clone = clone,
                        col = list(agent=c('CON'='#E0E0E0', 'DMS'='#ADADAD', 'AZA'='#FFD2D2', 'DAC'='#FF9797',
                                           'AS'='#FFFF37','CO'='#FF5151','LCD'='#0080FF',
                                           'HCD'='#005AB5', 'BAP'='#00DB00', 'AS_BAP'='#FFDC35', 'CO_BAP'= '#EA0000',
                                           'LCD_BAP'= '#0000E3', "HCD_BAP"="#000079"),
                                   clone=c('WT'="#FF2D2D",'L858R'="#FF9224",
                                           "DEL19"= "#66B3FF","YAP"="#2828FF")))

ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                        show_row_names = T, show_column_names = F,
                        name = "Z-score"
)


# draw heatmap function ---------------------------------------------------
draw_heatmap <- function(data_path=data_path,
                         file=file,
                         log_crit=3,
                         groups=groups
                         ){
  data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
  data <- data %>% filter(abs(M)>log_crit)  # filter log2FC criteria
  list <- data$ENSEMBL  # 抓出差異ensembl id
  gene_df <- data[,c("ENSEMBL","SYMBOL")] # 
  mycount_df$ENSEMBL <- rownames(mycount_df)
  new_df <- mycount_df %>% left_join(.,gene_df, by="ENSEMBL")
  mat <- mycount_df %>% filter(.,rownames(mycount_df) %in% list) %>% 
    left_join(.,gene_df, by="ENSEMBL") %>% group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>% 
    column_to_rownames(., var = "SYMBOL")

  if(groups == "AS"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AS,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP)
  } else if(groups == "CO"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_CO,ip_Y_V_S_BAP,ip_Y_V_S_CO_BAP)
  } else if(groups == "CD"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
                               ip_Y_V_S_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else if(groups == "BAP"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP
                               ,ip_Y_V_S_CO_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else{
    cat(c("<- please type in groups", "\n"))
    break
  }
  cat(c("<- filtering data", "\n"))
  
  mat_scale <- data_mat %>% t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()
  cat(c("<- annotation", "\n"))
  
  col <- colnames(mat_scale)
  
  # agent name
  agent <- factor (
    str_replace_all(col, c("ip_L_V_L_|ip_Y_V_S_"='')),
    levels=c('CON','DMS',"AZA","DAC",'AS',"CO","LCD","HCD","BAP",
             "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
  
  # clone name
  clone <- factor( 
    str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                         "D" = "DEL19", "Y" = "YAP")),
    levels=c('WT','L858R',"DEL19","YAP"))
  
  # draw hp
  ha <- HeatmapAnnotation(agent = agent, clone = clone,
                          col = list(agent=c('CON'='#E0E0E0', 'DMS'='#ADADAD', 'AZA'='#FFD2D2', 'DAC'='#FF9797',
                                             'AS'='#FFFF37','CO'='#FF5151','LCD'='#0080FF',
                                             'HCD'='#005AB5', 'BAP'='#00DB00', 'AS_BAP'='#FFDC35', 'CO_BAP'= '#EA0000',
                                             'LCD_BAP'= '#0000E3', "HCD_BAP"="#000079"),
                                     clone=c('WT'="#FF2D2D",'L858R'="#FF9224",
                                             "DEL19"= "#66B3FF","YAP"="#2828FF")))
  cat(c("<- drawing heatmap", "\n"))
  
  ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                          show_row_names = T, show_column_names = F,
                          name = "Z-score"
  )
  
}
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1/raw/Y_L_V"
draw_heatmap(data_path = data_path ,file = "ip_Y_V_S_BAP_0_deg.xlsx",groups = "AS",log_crit = 4)
