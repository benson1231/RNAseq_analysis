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

y_log <- mycount_df %>% na.omit() %>% t() %>% scale(scale = T) %>% t() %>% na.omit()
table(is.na(y_log))


plot <- ComplexHeatmap::Heatmap(y_log, show_row_names = F, na_col = "black",
                                cluster_columns = F)
plot



# with annotation ---------------------------------------------------------
col <- colnames(mycount_df)

# agent name
agent = factor (
  str_replace_all(col, c("ip_L_V_L_|ip_Y_V_S_"='')),
  levels=c('CON','DMS',"AZA","DAC",'AS',"CO","LCD","HCD","BAP",
           "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
# clone name
clone = factor( 
  
  str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                       "D" = "DEL19", "Y" = "YAP")),
  
  levels=c('WT','L858R',"DEL19","YAP"))

#draw hp
ha = HeatmapAnnotation(agent = agent, clone = clone,
                       col = list(agent=c('Con'='#c44dff', 'DMS'='#ff0066', 'AZA'='#ff99cc', 'DAC'='#66ffff',
                                          'AS'='#b3ff66','CO'='#000000','LCD'='#999999',
                                          'HCD'='#008ae6', 'BAP'='#ffff00', 'AS_BAP'='#663300', 'CO_BAP'= '#008000',
                                          'LCD_BAP'= 'pink', "HCD_BAP"="orange"),
                                  clone=c('WT'="red",'L858R'="orange",
                                          "DEL19"= "yellow","YAP"="green")))

ComplexHeatmap::Heatmap(y_log, top_annotation = ha, cluster_columns = F, 
                        show_row_names = F, show_column_names = F
                        )
