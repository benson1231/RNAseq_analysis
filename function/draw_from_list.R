### function
draw_from_list <- function(list,
                           groups=groups,
                           id="ENSEMBL"  # ENSEMBL/SYMBOL
                           ){
  if(id=="ENSEMBL"){
    list <- list  
    mycount_df$ENSEMBL <- rownames(mycount_df)
    new_df <- mycount_df %>% left_join(.,gene_df, by="ENSEMBL")
    mat <- mycount_df %>% filter(.,rownames(mycount_df) %in% list) %>% 
      left_join(.,gene_df, by="ENSEMBL") %>% group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>% 
      column_to_rownames(., var = "SYMBOL")
  }else if(id=="SYMBOL"){
    list <- list  
    mycount_df$ENSEMBL <- rownames(mycount_df)
    new_df <- mycount_df %>% 
      left_join(.,gene_df, by="ENSEMBL") %>% 
      group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>% 
      na.omit() %>% 
      column_to_rownames(.,var="SYMBOL") 
    mat <- new_df %>% 
      filter(.,rownames(new_df) %in% list)
  }else{
    cat(c("<- error, check gene id", "\n"))
    break
  }
  
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

### run function
list <- c("CYP1A1","MTF1","MTG1")
draw_from_list(list = list,groups = "AS",id = "SYMBOL")
