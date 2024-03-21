library(magrittr)
library(tidyverse)
library(stringr)
library(ComplexHeatmap)
df <- read.csv("/Users/benson/Documents/project/RNA-seq1/mycount_df.csv")

# my_df <- df[,grepl("ip_L_V_L", colnames(mycount_df))]
mycount_df <- df %>% select(.,ip_L_V_L_CON,ip_L_V_L_DMS,ip_L_V_L_AZA,ip_L_V_L_DAC,
                            ip_L_V_L_AS,ip_L_V_L_CO,ip_L_V_L_LCD,ip_L_V_L_HCD,
                            ip_L_V_L_BAP, ip_L_V_L_AS_BAP, ip_L_V_L_CO_BAP,
                            ip_L_V_L_LCD_BAP,ip_L_V_L_HCD_BAP,
                            ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
                            ip_Y_V_S_AS,ip_Y_V_S_CO,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
                            ip_Y_V_S_BAP, ip_Y_V_S_AS_BAP, ip_Y_V_S_CO_BAP,
                            ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP) 

df_log <- mycount_df %>% na.omit() %>% t() %>% scale(scale = T) %>% t() %>% na.omit()
table(is.na(df_log))

# with annotation ---------------------------------------------------------
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

#draw hp
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
