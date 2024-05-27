# function for RNA-seq data 
# def
CON <- c(1,2, 14,15, 27,28, 40,41)
AS <- c(1,2,9,5,10, 14,15,22,18,23, 27,28,35,31,36, 40,41,48,44,49)
only_AS <- c(9,5,10, 22,18,23, 35,31,36, 48,44,49)
CO <- c(1,2,9,6,11, 14,15,22,19,24, 27,28,35,32,37, 40,41,48,45,50)
only_CO <- c(9,6,11, 22,19,24, 35,32,37, 48,45,50)
LCD <- c(1,2,9,7,12, 14,15,22,20,25, 27,28,35,33,38, 40,41,48,46,51)
only_LCD <- c(9,7,12, 22,20,25, 35,33,38, 48,46,51)
HCD <- c(1,2,9,8,13, 14,15,22,21,26, 27,28,35,34,39, 40,41,48,47,52)
only_HCD <- c(9,8,13, 22,21,26, 35,34,39, 48,47,52)
CD <- c(1,2,9,7,8,12,13, 14,15,22,20,21,25,26, 27,28,35,33,34,38,39, 40,41,48,46,47,51,52)
only_CD <- c(9,7,8,12,13, 22,20,21,25,26, 35,33,3438,39, 48,46,47,51,52)
BAP <- c(1,2,9,10,11,12,13, 14,15,22,23,24,25,26, 27,28,35,36,37,38,39, 40,41,48,49,50,51,52)
# abb ---------------------------------------------------------------------
abb <- function(name,type="group"){
  if(!(type%in%c("file","group"))){
    stop("'type' must be 'file' or 'group'")
  }
  if(type=="file"){
    df <- name_df %>% filter(file_name==name) %>% pull(abbreviate)
  } else{
    df <- name_df %>% filter(group_name==name) %>% pull(abbreviate)
  }
  return(df)
}

# force_list --------------------------------------------------------------
force_list <- function(x){
  input_string <- x
  output_list <- unlist(strsplit(input_string, split = "/"))
  return(output_list)
}

# # filter_gene -------------------------------------------------------------
# filter_gene <- function(file_name){
#   df <- readxl::read_xlsx(file.path(data_path, file_name))
#   sum <- df[,c(1,10:11)] %>% column_to_rownames("ENSEMBL") %>% 
#     setNames(c("group1","group2")) %>% 
#     mutate("avg" = (group1 + group2)/2) %>% 
#     mutate("g1"=group1 > 1) %>% 
#     mutate("g2"=group2 > 1) %>% 
#     filter(avg > 2)
#   df <- df %>% filter(ENSEMBL %in% rownames(sum)) 
#   return(df)
# }

# draw_heatmap ------------------------------------------------------------
draw_heatmap <- function(file=file,
                         log_crit=3,
                         groups="ALL",
                         show_row_names = FALSE,
                         log_scale=FALSE,
                         row_km=T,
                         km=0,
                         return_cluster=F,
                         title=""
                         ){
  if(!(groups%in%c("AS","CO","LCD","HCD","CD","BAP","only_AS","only_CO",
                  "only_LCD","only_HCD","only_CD","CON","ALL"))){
    stop('"groups" must be "AS","CO","LCD","HCD","CD","BAP","only_AS","only_CO",
         "only_LCD","only_HCD","only_CD","CON" or "ALL".')
  }
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
  data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
  data <- data %>% filter(abs(M)>log_crit,
                          !(SYMBOL%in% c("havana","ensembl_havana","havana_tagene")))  
  cat(c(" -> log2FC criteria is", log_crit,"\n"))
  
  list <- data$ENSEMBL  # 抓出差異ensembl id
  gene_df <- data[,c("ENSEMBL","SYMBOL")] # 
  mycount_df$ENSEMBL <- rownames(mycount_df)
  new_df <- mycount_df %>% left_join(.,gene_df, by="ENSEMBL")
  mat <- mycount_df %>% filter(.,rownames(mycount_df) %in% list) %>% 
    left_join(.,gene_df, by="ENSEMBL") %>% group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>% 
    column_to_rownames(., var = "SYMBOL")
  
  # select groups
  if(groups == "AS"){
    data_mat <- mat %>% select(all_of(name_df$group_name[AS]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#00DB00',
                          'AS' = '#FFFF37', 'AS_BAP' = '#FFDC35'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "CO"){
    data_mat <- mat %>% select(all_of(name_df$group_name[CO]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#00DB00',
                          'CO' = '#FF5151', 'CO_BAP' = '#EA0000'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "CD"){
    data_mat <- mat %>% select(all_of(name_df$group_name[CD]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#00DB00',
                          'LCD' = '#0080FF', 'HCD' = '#005AB5', 
                          'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "BAP"){
    data_mat <- mat %>% select(all_of(name_df$group_name[BAP]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#00DB00',
                          'AS_BAP' = '#FFDC35', 'CO_BAP' = '#EA0000',
                          'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "ALL"){
    data_mat <- mat %>% select(all_of(name_df$group_name))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'AZA' = '#FFD2D2', 'DAC' = '#FF9797',
                          'BAP' = '#00DB00','AS' = '#FFFF37', 'CO' = '#FF5151',
                          'LCD' = '#0080FF', 'HCD' = '#005AB5', 
                          'AS_BAP' = '#FFDC35', 'CO_BAP' = '#EA0000',
                          'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_CO"){
    data_mat <- mat %>% select(all_of(name_df$group_name[only_CO]))
    ann <- list(agent = c('BAP' = '#00DB00','CO' = '#FF5151', 'CO_BAP' = '#EA0000'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_AS"){
    data_mat <- mat %>% select(all_of(name_df$group_name[only_AS]))
    ann <- list(agent = c('BAP' = '#00DB00', 'AS' = '#FFFF37', 'AS_BAP' = '#FFDC35'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_LCD"){
    data_mat <- mat %>% select(all_of(name_df$group_name[only_LCD]))
    ann <- list(agent = c('BAP' = '#00DB00', 'LCD' = '#0080FF',  'LCD_BAP' = '#0000E3'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_HCD"){
    data_mat <- mat %>% select(name_df$group_name[only_HCD])
    ann <- list(agent = c('BAP' = '#00DB00','HCD' = '#005AB5', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_CD"){
    data_mat <- mat %>% select(all_of(name_df$group_name[only_CD]))
    ann <- list(agent = c('BAP' = '#00DB00','LCD' = '#0080FF', 'HCD' = '#005AB5',
                          'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "CON"){
    data_mat <- mat %>% select(all_of(name_df$group_name[CON]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else{
    return(NULL)
  }
  
  if(log_scale==TRUE){
    mat_scale <- data_mat %>% log2() %>% 
      t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()
  }else{
    mat_scale <- data_mat %>% t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()
  }
  
  ### heat-map argument
  col <- colnames(mat_scale)
  
  # agent name
  agent <- factor (
    str_replace_all(col, c("ip_W_V_S_|ip_L_V_S_|ip_D_V_S_|ip_Y_V_S_"='')),
    levels=c('CON','DMS',"AZA","DAC","BAP",'AS',"CO","LCD","HCD",
             "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
  
  # clone name
  clone <- factor( 
    str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                         "D" = "Del19", "Y" = "YAP")),
    levels=c('WT','L858R',"Del19","YAP"))
  
  # draw hp
  ha <- HeatmapAnnotation(agent = agent, clone = clone,
                          col = ann)

  heat <- ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                                  show_row_names = show_row_names, 
                                  show_column_names = F, row_title = title,
                                  name = "Z-score",row_km = km)
  cat(c( " ->", nrow(mat_scale), "genes"))
  
  if(return_cluster==T){
    group = kmeans((mat_scale), centers = km)$cluster
    return(group)
  } else{
    heat
  }
}

# draw_from_list ----------------------------------------------------------
draw_from_list <- function(list,
                           groups="ALL",
                           id="SYMBOL",
                           cluster=TRUE,
                           anno=TRUE,
                           title="",
                           show_row_names = TRUE, # ENSEMBL/SYMBOL
                           label_num=F,
                           log_scale=FALSE
                           ){
  if(!(groups%in%c("AS","CO","LCD","HCD","CD","BAP","only_AS","only_CO",
                  "only_LCD","only_HCD","only_CD","CON","ALL"))){
    stop('"groups" must be "AS","CO","LCD","HCD","CD","BAP","only_AS","only_CO",
         "only_LCD","only_HCD","only_CD","CON" or "ALL".')
  }
  cat(c(" -> input list with",length(list),"genes","\n"))
  if(id=="ENSEMBL"){
    list <- list  
    mat <- mycount_df %>%
      mutate(.,ENSEMBL=rownames(mycount_df)) %>% 
      left_join(.,gene_df, by="ENSEMBL") %>% 
      filter(.,rownames(mycount_df) %in% list,
             !(SYMBOL%in% c("havana","ensembl_havana","havana_tagene"))) %>% 
      group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>% na.omit() %>%  
      column_to_rownames(., var = "SYMBOL")
  }else if(id=="SYMBOL"){
    list <- list  
    mat <- mycount_df %>% 
      mutate(.,ENSEMBL=rownames(mycount_df)) %>% 
      left_join(.,gene_df, by="ENSEMBL") %>% 
      group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>% 
      na.omit() %>% 
      filter(.,SYMBOL %in% list,
             !(SYMBOL %in% c("havana","ensembl_havana","havana_tagene"))) %>% 
      column_to_rownames(.,var="SYMBOL") 
  }else{
    cat(c("<- error: check gene id", "\n"))
    return(NULL)
  }
  cat(c(" -> get",length(rownames(mat)),"genes matrix. Scaling and plotting.","\n"))
  
  # select groups
  if(groups == "AS"){
    data_mat <- mat %>% select(all_of(name_df$group_name[AS]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#00DB00',
                          'AS' = '#FFFF37', 'AS_BAP' = '#FFDC35'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "CO"){
    data_mat <- mat %>% select(all_of(name_df$group_name[CO]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#00DB00',
                          'CO' = '#FF5151', 'CO_BAP' = '#EA0000'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "CD"){
    data_mat <- mat %>% select(all_of(name_df$group_name[CD]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#00DB00',
                          'LCD' = '#0080FF', 'HCD' = '#005AB5', 
                          'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "BAP"){
    data_mat <- mat %>% select(all_of(name_df$group_name[BAP]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#00DB00',
                          'AS_BAP' = '#FFDC35', 'CO_BAP' = '#EA0000',
                          'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "ALL"){
    data_mat <- mat %>% select(all_of(name_df$group_name))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'AZA' = '#FFD2D2', 'DAC' = '#FF9797',
                          'BAP' = '#00DB00','AS' = '#FFFF37', 'CO' = '#FF5151',
                          'LCD' = '#0080FF', 'HCD' = '#005AB5', 
                          'AS_BAP' = '#FFDC35', 'CO_BAP' = '#EA0000',
                          'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_CO"){
    data_mat <- mat %>% select(all_of(name_df$group_name[only_CO]))
    ann <- list(agent = c('BAP' = '#00DB00','CO' = '#FF5151', 'CO_BAP' = '#EA0000'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_AS"){
    data_mat <- mat %>% select(all_of(name_df$group_name[only_AS]))
    ann <- list(agent = c('BAP' = '#00DB00', 'AS' = '#FFFF37', 'AS_BAP' = '#FFDC35'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_LCD"){
    data_mat <- mat %>% select(all_of(name_df$group_name[only_LCD]))
    ann <- list(agent = c('BAP' = '#00DB00', 'LCD' = '#0080FF',  'LCD_BAP' = '#0000E3'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_HCD"){
    data_mat <- mat %>% select(name_df$group_name[only_HCD])
    ann <- list(agent = c('BAP' = '#00DB00','HCD' = '#005AB5', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "only_CD"){
    data_mat <- mat %>% select(all_of(name_df$group_name[only_CD]))
    ann <- list(agent = c('BAP' = '#00DB00','LCD' = '#0080FF', 'HCD' = '#005AB5',
                          'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else if(groups == "CON"){
    data_mat <- mat %>% select(all_of(name_df$group_name[CON]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD'),
                clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
                          'Del19' = '#FF8800', 'YAP' = '#FF44AA'))
  } else{
    return(NULL)
  }
  
  # scaling
  if(log_scale==TRUE){
    mat_scale <- data_mat %>% log2() %>%  t() %>% scale() %>% t() %>% as.matrix() %>% na.omit()
  } else{
    mat_scale <- data_mat  %>%  t() %>% scale() %>% t() %>% as.matrix() %>% na.omit()
  }
  
  # color
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  # heat-map argument
  col <- colnames(mat_scale)
  
  # agent name
  agent <- factor (
    str_replace_all(col, c("ip_W_V_S_|ip_L_V_S_|ip_D_V_S_|ip_Y_V_S_"='')),
    levels=c('CON','DMS',"AZA","DAC","BAP",'AS',"CO","LCD","HCD",
             "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
  
  # clone name
  clone <- factor( 
    str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                         "D" = "Del19", "Y" = "YAP")),
    levels=c('WT','L858R',"Del19","YAP"))
  
  # draw hp
  ha <- HeatmapAnnotation(agent = agent, clone = clone,
                          col = ann)
  
  if(anno==T){
    if(label_num==T){
      ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                              show_row_names = show_row_names,
                              show_column_names = F,cluster_rows = cluster,
                              col = col_fun,
                              name = "Z-score", row_title = title,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.1f", mat_scale[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
    }else{
      ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                              show_row_names = show_row_names,col = col_fun,
                              show_column_names = F,cluster_rows = cluster,
                              name = "Z-score", row_title = title
                              )
    }
  } else{
    if(label_num==T){
      ComplexHeatmap::Heatmap(mat_scale, cluster_columns = F, 
                              show_row_names = show_row_names,col = col_fun,
                              show_column_names = F,cluster_rows = cluster,
                              name = "Z-score", row_title = title,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.1f", mat_scale[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
    }else{
      ComplexHeatmap::Heatmap(mat_scale, cluster_columns = F, 
                              show_row_names = show_row_names,col = col_fun,
                              show_column_names = F,cluster_rows = cluster,
                              name = "Z-score", row_title = title,
                              )
    }
  }
}

# get_deg --------------------------------------------------------
get_deg <- function(file=file,
                    log_crit = c(1,-1),
                    dir="all",  # all/up/down
                    type="SYMBOL",
                    top="all"
                    ){
  
  # 檢查 type 是否有效
  if (!(type %in% c("ENTREZID", "ENSEMBL", "SYMBOL"))) {
    stop("Invalid type. Allowed values are 'ENTREZID','ENSEMBL'or 'SYMBOL'.")
  }
  
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
  raw_data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
  
  # 根據不同的類型選擇相應的欄位名稱
  if(type == "ENTREZID") {
    id_col <- "ENTREZID"
  } else if(type == "ENSEMBL") {
    id_col <- "ENSEMBL"
  } else{
    id_col <- "SYMBOL"
  }
  
  # 處理資料
  data <- raw_data %>%
    as.data.frame() %>% 
    filter(.,!(SYMBOL%in% c("havana","ensembl_havana","havana_tagene")),
           geneBiotype=="protein_coding") %>% 
    select(all_of(id_col),M) %>%
    group_by(across(all_of(id_col))) %>%
    summarize(across(where(is.numeric), mean)) %>%
    arrange(desc(M))
  
  if(top=="all"){
    if(dir=="all"){
      data <- data %>% filter(M >log_crit[1]| M < log_crit[2])
      cat(c(" -> abs(log2FC) larger than", log_crit[1],"->",length(data$M),"genes.\n"))
    }else if(dir=="up"){
      data <- data %>% filter(M > log_crit[1])
      cat(c(" -> log2FC larger than", log_crit[1],"->",length(data$M),"genes.\n"))
    }else if(dir=="down"){
      data <- data %>% filter(M < log_crit[2])
      cat(c(" -> log2FC smaller than", log_crit[2],"->",length(data$M),"genes.\n"))
    }else{
      cat(c("<- error, check direction", "\n"))
      return(NULL)
    }
  } else{
    if(dir=="all"){
      data <- data %>% filter(M >log_crit[1]| M < log_crit[2]) %>% 
        arrange(desc(abs(M))) %>% .[1:top,]
      cat(c(" -> abs(log2FC) larger than",log_crit[1],"->",length(data$M),"genes.\n"))
    }else if(dir=="up"){
      data <- data %>% filter(M > log_crit[1]) %>% 
        arrange(desc(abs(M))) %>% .[1:top,]
      cat(c(" -> log2FC larger than", log_crit[1],"->",length(data$M),"genes.\n"))
    }else if(dir=="down"){
      data <- data %>% filter(M < log_crit[2]) %>% 
        arrange(desc(abs(M))) %>% .[1:top,]
      cat(c(" -> log2FC smaller than", log_crit[2],"->",length(data$M),"genes.\n"))
    }else{
      cat(c("<- error, check direction", "\n"))
      return(NULL)
    }
  }

  # 抓出差異gene id
  if(type=="ENTREZID"){
    list <- data$ENTREZID %>% na.omit()
  }else if(type=="ENSEMBL"){
    list <- data$ENSEMBL
  }else{
    list <- data$SYMBOL
  }
  
  return(list)
}

# get_list ------------------------------------------------------
get_list <- function(file_name, 
                     type = "ENTREZID"
                     ) {
  # 檢查 type 是否有效
  if (!(type %in% c("ENTREZID", "ENSEMBL"))) {
    stop("Invalid type. Allowed values are 'ENTREZID' or 'ENSEMBL'.")
  }
  # 讀取資料
  raw_df <- readxl::read_xlsx(file.path(data_path, file_name)) %>% as.data.frame()
  
  # 根據不同的類型選擇相應的欄位名稱
  if (type == "ENTREZID") {
    id_col <- "ENTREZID"
  } else {
    id_col <- "ENSEMBL"
  }
  
  # 處理資料
  df <- raw_df %>% as.data.frame() %>% 
    select(all_of(id_col),M) %>%
    group_by(across(all_of(id_col))) %>%
    summarize(across(where(is.numeric), mean)) %>%
    arrange(desc(M))  %>%
    dplyr::rename(.,logFC=M) # 確保將 'M' 列重命名為 'logFC'
  
  # 創建列表
  enrich_list <- df$logFC 
  names(enrich_list) <- df[[id_col]]
  enrich_list <-  na.omit(enrich_list)
  
  return(enrich_list)
}

# gsea_run ----------------------------------------------------------------
gsea_run <- function(file,
                     method="GSEA",
                     list,
                     list_id="ensembl",   # ensembl/symbol
                     gsea_term="ALL"  # "BP","MF","CC"
                     ){
  # 檢查 type 是否有效
  if (!(list_id %in% c("ensembl", "symbol"))) {
    stop("Invalid type. Allowed values are 'ensembl' or 'symbol'")
  }
  if (!(method %in% c("GSEA", "ORA"))) {
    stop("Invalid type. Allowed methods are 'GSEA' or 'ORA'")
  }
  
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
  if(method=="GSEA"){
    df <- readxl::read_xlsx(file.path(data_path, file))
    cat(c(" >- input all gene and run GSEA","\n"))
  } else {
    if(list_id=="ensembl"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        filter(ENSEMBL %in% list)
    } else {
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        filter(SYMBOL %in% list)
    } 
    cat(c(" -> input selected gene and run ORA","\n"))
  } 
  
  # set organism for human
  organism <- "org.Hs.eg.db"
  cat(c(" -> organism type is", organism,"\n"))
  
  # $csv file's colume namm of log2 fold change
  original_gene_list <- df$M
  # $csv file's colume namm of ENSEMBL ID 
  names(original_gene_list) <- df$ENSEMBL
  # omit any NA values and sort the list in decreasing order
  gsea_gene_list <- na.omit(original_gene_list) %>% 
    sort(., decreasing = TRUE)
  cat(c(" -> length of the gene list is",length(gsea_gene_list),"\n"))
  #run GSEA
  cat(c(" -> running GSEA in", gsea_term ,"GSEA-term\n"))
  gse <- gseGO(geneList = gsea_gene_list, 
               ont = gsea_term, 
               keyType = "ENSEMBL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
  return(gse)
}


# kegg_run ----------------------------------------------------------------
kegg_run <- function(file,
                     method="GSEA",
                     list,
                     list_id="ensembl"   # ensembl/symbol
                     ){
  # 檢查 type 是否有效
  if (!(list_id %in% c("ensembl", "symbol"))) {
    stop("Invalid type. Allowed values are 'ensembl' or 'symbol'")
  }
  if (!(method %in% c("GSEA", "ORA"))) {
    stop("Invalid type. Allowed methods are 'GSEA' or 'ORA'")
  }
  
  # select gene in list
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
  if(method=="GSEA"){
    df <- readxl::read_xlsx(file.path(data_path, file))
    cat(c(" >- input all gene and run KEGG","\n"))
  } else {
    if(list_id=="ensembl"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        filter(ENSEMBL %in% list)
    } else {
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        filter(SYMBOL %in% list)
    } 
    cat(c(" >- input selected gene and run ORA","\n"))
  } 
  
  # change KEGG Organism Code from https://www.genome.jp/kegg/catalog/org_list.html 
  kegg_organism <-  "hsa"  # human is "hsa"
  cat(c(" -> organism type is", kegg_organism,"\n"))
  
  # selcet log2FC value
  kegg_gene_list <- df$M
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- df$ENTREZID
  # omit any NA values and sort the list in decreasing order
  kegg_gene_list<- na.omit(kegg_gene_list) %>% 
    sort(., decreasing = TRUE)
  cat(c(" -> length of the gene list is",length(kegg_gene_list),"\n"))
  # Run KEGG
  cat(c(" -> running KEGG", "\n"))
  kegg <- gseKEGG(geneList     = kegg_gene_list,
                  organism     = kegg_organism,
                  nPerm        = 10000,
                  minGSSize    = 3,
                  maxGSSize    = 800,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "none",
                  keyType       = "ncbi-geneid")
  return(kegg)
}

# get_df ------------------------------------------------------------------
get_df <- function(file,
                   de=F,
                   list=c(),
                   dir="all",
                   log_crit=c(1,-1),
                   with_D=FALSE,
                   all=FALSE
                   ){
  if(all==TRUE){
    df <- readxl::read_xlsx(file.path(data_path, file)) %>% as.data.frame()
    return(df)
  }
  
  if(de==TRUE){
    if(dir=="all"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(abs(M) > log_crit[1]) %>% 
        select(M, SYMBOL) %>% 
        setNames(c("logFC","gene"))
    }else if(dir=="up"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(M > log_crit[1]) %>% 
        select(M, SYMBOL) %>% 
        setNames(c("logFC","gene"))
    }else if(dir=="down"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(M < log_crit[2]) %>% 
        select(M, SYMBOL) %>% 
        setNames(c("logFC","gene"))
    }else{
      cat(c("<- error, check direction", "\n"))
      return(NULL)
    }
    
    return(df)
    
  } else {
    
    if(with_D==T){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(ENSEMBL %in% list) %>% 
        select(M,D, SYMBOL) %>% 
        setNames(c("logFC","D","gene"))
    } else{
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(ENSEMBL %in% list) %>% 
        select(M, SYMBOL) %>% 
        setNames(c("logFC","gene"))
    }
    
    return(df)
  }
}

# venn_to_excel -----------------------------------------------------------
# 將venn diagram結果輸出成excel
venn_to_excel <- function(venn_list, name) {
  # list中資料依序加入data frame
  max_length <- max(sapply(venn_list$item, length))
  data_df <- data.frame(matrix(NA, ncol = length(venn_list$item), nrow = max_length))
  for (i in 1:length(venn_list$item)) {
    sublist <- venn_list$item[[i]]
    data_df[, i] <- c(sublist, rep(NA, max_length - length(sublist)))
    colnames(data_df) <- venn_list$name
  }
  # 設置欄寬
  wb <- openxlsx::createWorkbook()
  addWorksheet(wb, "Sheet1")
  writeData(wb, "Sheet1", data_df)
  for (i in 1:ncol(data_df)) {
    setColWidths(wb, sheet = 1, cols = i, widths = "auto")
  }
  saveWorkbook(wb, paste0(name, ".xlsx"))
  cat(paste("-> Output completed. File name is ", paste0(name, ".xlsx")))
}

# get_kegg_list -----------------------------------------------------------
get_kegg_list <- function(path_ID
                          ){
  # 获取指定通路的基因列表
  pathway_df <- KEGGREST::keggGet(path_ID)[[1]]
  pathway_genes <- pathway_df$GENE
  cat(c(" -> pathway name:", pathway_df$NAME,"\n"))
  # 提取基因名称
  gene_names <- sapply(strsplit(pathway_genes, ";"), `[`, 1)
  # 删除奇数索引的元素
  gene_list_even <- gene_names[seq_along(gene_names) %% 2 == 0]
  cat(c(" ->",length(gene_list_even),"genes involved.\n"))
  return(gene_list_even)
}

# get_anno ----------------------------------------------------------------
get_anno <- function(list,
                     file,
                     id_type){
  # 檢查 type 是否有效
  if (!(id_type %in% c("ENSEMBL","SYMBOL"))) {
    stop("Invalid type. Allowed values are 'ENSEMBL','SYMBOL'.")
  }
  # 讀檔獲取df
  if(id_type=="ENSEMBL"){
    df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
      as.data.frame() %>% 
      filter(ENSEMBL %in% list)
  } else if(id_type=="SYMBOL"){
    df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
      as.data.frame() %>% 
      filter(SYMBOL %in% list)
  }
}


# ensembl_to_symbol -------------------------------------------------------
ensembl_to_symbol <- function(list){
  df <- gene_df %>% filter(ENSEMBL %in% list) %>% pull(SYMBOL)
}
# negative ----------------------------------------------------------------
negative <- function(x) {
  return(-x)
}

# run_keg_path -----------------------------------------------------------
run_keg_path <- function(file_name, dir="up", log_crit=1){
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
  message("Running KEGG analysis")
  kk <- enrichKEGG(gene         = gene,
                   organism     = 'hsa',
                   pvalueCutoff = 0.1)
  # ENTREZID to SYMBOL
  net <- DOSE::setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
  return(net)
}

# run_reactome ------------------------------------------------------------
run_reactome <- function(file_name, dir="up",log_crit=1,
                         return_df = F){
  if (!(dir %in% c("up","down"))) {
    stop("Invalid type. Allowed values are 'up' or 'down'.")
  }
  genelist <- get_list(file_name)
  if(dir=="up"){
    gene <- names(genelist)[genelist >= log_crit]
  } else{
    gene <- names(genelist)[genelist <= negative(log_crit)]
  }
  # Reactome pathway over-representation analysis
  message("Running reactome analysis")
  de_pathway <- enrichPathway(gene=gene, pvalueCutoff = 0.1, readable=TRUE)
  return(de_pathway)
}
# get_p -------------------------------------------------------------------
get_p <- function(x, group_name, top=10, p_value_cutoff = 0.01){
  pvalue <- x@result %>% .[,c("Description","pvalue")] %>% 
    filter(pvalue < p_value_cutoff) %>% 
    setNames(c("Description",group_name)) %>% head(top)
  return(pvalue)
}

# plot_heatmap ------------------------------------------------------------
plot_heatmap <- function(file_list, group_names, analysis, dir="up",col_title="", row_title="") {
  # 檢查參數是否匹配
  if (length(file_list) != length(group_names)) {
    stop("file_list and group_names must have the same length")
  }
  if (!(analysis %in% c("kegg", "reactome"))) {
    stop("Invalid type. Allowed values are 'kegg' or 'reactome'.")
  }
  # 讀取並處理每個文件
  if(analysis=="kegg"){
    dfs <- lapply(file_list, function(file) run_keg_path(file, log_crit = 1, dir = dir))
  } else{
    dfs <- lapply(file_list, function(file) run_reactome(file, log_crit = 1, dir = dir))
  }
  processed_dfs <- mapply(get_p, dfs, group_name = group_names, SIMPLIFY = FALSE)
  
  # 合併所有數據框
  combined_df <- Reduce(function(x, y) full_join(x, y, by = "Description"), processed_dfs) %>% 
    column_to_rownames("Description") %>% 
    as.matrix()
  
  # 計算 -log10(P) 值
  mat <- log10(combined_df) %>% -.
  
  # 繪製熱圖
  if(dir=="up"){
    col_fun = colorRamp2(c(1, 5), c("white", "red"))
  } else{
    col_fun = colorRamp2(c(1, 5), c("white", "blue"))
  }
  ComplexHeatmap::Heatmap(mat, 
                          na_col = "grey", 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE,
                          col = col_fun, 
                          name = "-log10(P)",
                          row_names_gp = gpar(fontsize = 7), 
                          column_names_gp = gpar(fontsize = 8),
                          column_title_gp = gpar(fontsize = 10),
                          column_names_rot = 45, 
                          column_title = col_title, row_title = row_title,
                          heatmap_legend_param = list(   # 控制color bar在0~5之間
                            at = c(1, 2, 3, 4, 5),  
                            labels = c("1", "2", "3", "4", "5")  ))
}

# run_path_heatmap --------------------------------------------------------
run_path_heatmap <- function(count){
  count <- count %>% as.matrix()
  message("running MLM analysis")
  sample_acts <- run_mlm(mat=count, net=net, .source='source', .target='target',
                         .mor='weight', minsize = 5)
  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()
  
  # Scale per feature
  sample_acts_mat <- scale(sample_acts_mat)
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))
  ### heatmap
  p <- pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)
  return(p)
}

# run_TF_heatmap ----------------------------------------------------------
run_TF_heatmap <- function(count){
  count <- count %>% as.matrix()
  message("running ULM analysis")
  sample_TF <- run_ulm(mat=count, net=net_TF, .source='source', .target='target',
                       .mor='mor', minsize = 5)
  # plot
  n_tfs <- 25
  # Transform to wide matrix
  sample_TF_mat <- sample_TF %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()
  
  # Get top tfs with more variable means across clusters
  tfs <- sample_TF %>%
    group_by(source) %>%
    summarise(std = sd(score)) %>%
    arrange(-abs(std)) %>%
    head(n_tfs) %>%
    pull(source)
  sample_TF_mat <- sample_TF_mat[,tfs]
  
  # Scale per sample
  sample_TF_mat <- scale(sample_TF_mat)
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))
  
  ### heatmap
  pheatmap(sample_TF_mat, border_color = NA, color=my_color, breaks = my_breaks) 
}

# run_pathway -------------------------------------------------------------
run_pathway <- function(file_name,title=""){
  deg_FC <- get_df(file_name ,de=T,log_crit = c(1,-1)) %>% 
    group_by(gene) %>%
    summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
    column_to_rownames(., var = "gene")
  
  ### run MLM analysis
  message("running MLM analysis")
  contrast_acts <- run_mlm(mat=deg_FC, net=net, .source='source', .target='target',
                           .mor='weight', minsize = 5)
  
  ### pathway enrich plot
  ggplot(contrast_acts, aes(x = reorder(source, score), y = score)) + 
    geom_bar(aes(fill = score), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
                         mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
          axis.text.y = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    xlab("Pathways") + ggtitle(title)   # change title
}

# plot_pathway  ----------------------------------------------------------
plot_pathway <- function(file_name, pathway, title="",logFC_criteria = 1){
  deg_df <- get_df(file_name ,all = T) %>% group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
    column_to_rownames(., var = "SYMBOL") %>% 
    rownames_to_column("ID") %>% 
    .[,c(1,5,6,7,8)] %>% 
    setNames(c("ID","treat","control","M","D"))%>% 
    mutate(Average_expression =(treat+control)/2) %>% 
    .[,c(1,4,6)] %>% 
    setNames(c("ID","logFC","Average_expression"))
  ### filter pathway specific genes
  path <- net %>%
    filter(source == pathway) %>%
    arrange(target) %>%
    mutate(ID = target, color = "3") %>%
    column_to_rownames('target') %>% 
    left_join(deg_df,"ID") %>% 
    mutate(color = if_else(logFC > logFC_criteria, '1', color)) %>%
    mutate(color = if_else(logFC < logFC_criteria & logFC > negative(logFC_criteria), '2', color)) %>%   
    mutate(color = if_else(logFC < negative(logFC_criteria), '3', color)) 
  ### MD plot
  p <- ggplot(path, aes(x = log(Average_expression), y = logFC, color = color)) + geom_point() +
    scale_colour_manual(values = c("red","grey","royalblue3")) +
    geom_label_repel(aes(label = ID)) + 
    theme_minimal() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ggtitle(paste(title,pathway)) 
  print(p)
  return(path)
}

# run_TF ------------------------------------------------------------------
run_TF <- function(file_name, title="", n_tfs=25){
  ### two groups comparing TFs
  de_FC <- get_df(file_name ,de=T) %>% group_by(gene) %>%
    summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
    column_to_rownames(., var = "gene")
  
  ### Run ULM
  contrast_TF <- run_ulm(mat=de_FC[, "logFC", drop=FALSE], net=net_TF, .source='source', .target='target',
                         .mor='mor', minsize = 5)
  
  # Filter top TFs in both signs
  f_contrast_TF <- contrast_TF %>%
    mutate(rnk = NA)
  msk <- f_contrast_TF$score > 0
  f_contrast_TF[msk, 'rnk'] <- rank(-f_contrast_TF[msk, 'score'])
  f_contrast_TF[!msk, 'rnk'] <- rank(-abs(f_contrast_TF[!msk, 'score']))
  tfs <- f_contrast_TF %>%
    arrange(rnk) %>%
    head(n_tfs) %>%
    pull(source)
  f_contrast_TF <- f_contrast_TF %>%
    filter(source %in% tfs)
  
  ### TF enrich plot
  ggplot(f_contrast_TF, aes(x = reorder(source, score), y = score)) + 
    geom_bar(aes(fill = score), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
                         mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
          axis.text.y = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    xlab("TFs") + ggtitle(title) 
}

# plot_TF -----------------------------------------------------------------
plot_TF <- function(file_name, TF, title="", logFC_criteria = 1){
  de <- get_df(file_name ,all = T) %>% group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
    column_to_rownames(., var = "SYMBOL") %>% 
    rownames_to_column("ID") %>% 
    .[,c(1,5,6,7,8)] %>% 
    setNames(c("ID","treat","control","M","D"))%>% 
    mutate(Average_expression =(treat+control)/2) %>% 
    .[,c(1,4,6)] %>% setNames(c("ID","logFC","Average_expression"))
  ### Specific TF related genes 
  tf <- TF 
  logFC_criteria <- logFC_criteria  
  
  # filter TF related genes
  df <- net_TF %>%
    filter(source == tf) %>%
    arrange(target) %>%
    mutate(ID = target, color = "3") %>%
    column_to_rownames('target') %>% 
    left_join(de,"ID") %>% 
    mutate(color = if_else(logFC > logFC_criteria, '1', color)) %>%
    mutate(color = if_else(logFC < logFC_criteria & logFC > negative(logFC_criteria), '2', color)) %>%   
    mutate(color = if_else(logFC < negative(logFC_criteria), '3', color)) 
  
  ### MD plot
  p <- ggplot(df, aes(x = log(Average_expression), y = logFC, color = color, size=abs(mor))) +
    geom_point() +
    scale_colour_manual(values = c("red","grey","royalblue3")) +
    geom_label_repel(aes(label = ID, size=1)) + 
    theme_minimal() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ggtitle(paste(title,tf))
  print(p)
  return(df)
}





# get_divenn --------------------------------------------------------------
get_divenn <- function(file_name, top=50,output_name="up_down.xlsx"){
  df <- read.xlsx(file.path(data_path, file_name)) %>% 
    filter(., geneBiotype==" protein_coding") %>% 
    arrange(desc(abs(M))) %>% 
    head(top) %>% 
    mutate(direaction = case_when(
      M >= 1 ~ 1,
      M < 1 & M > -1 ~ 3,
      M <= -1 ~ 2 )) %>% 
    filter(., direaction != 3) %>% 
    select(SYMBOL, direaction)
  output_dir <- output_dir
  write.xlsx(df, file.path(output_dir,output_name))
  cat(c(" -> ",output_dir,"\n -> ",output_name))
}


# plot_MD -----------------------------------------------------------------
plot_MD <- function(file_name, title="", logFC_criteria = 1, only_DE=F){
  df <- get_df(file_name ,all = T) %>% group_by(SYMBOL) %>%
    filter(.,!(SYMBOL%in% c("havana","ensembl_havana","havana_tagene")),
           geneBiotype=="protein_coding") %>% 
    summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
    column_to_rownames(., var = "SYMBOL") %>% 
    rownames_to_column("ID") %>% 
    .[,c(1,5,6,7,8)] %>% 
    setNames(c("ID","treat","control","M","D"))%>% 
    mutate(Average_expression =(treat+control)/2) %>% 
    .[,c(1,4,6)] %>% setNames(c("ID","logFC","Average_expression")) %>% 
    mutate(ID = ID, color = "3") %>% 
    mutate(color = if_else(logFC > logFC_criteria, '1', color)) %>%
    mutate(color = if_else(logFC < logFC_criteria & logFC > negative(logFC_criteria), '3', color)) %>%   
    mutate(color = if_else(logFC < negative(logFC_criteria), '2', color))
  
  up <- df %>% filter(color==1) %>% nrow()
  down <- df %>% filter(color==2) %>% nrow()
  non <- df %>% filter(color==3) %>% nrow()
  cat(c(" -> up:",up," -> down:",down, " -> non:",non,"\n"))
  if(only_DE==T){
    df <- df %>% filter(color == 1 | color == 2)
  }
  
  ### MD plot
  message("plotting")
  ggplot(df, aes(x = log(Average_expression), y = logFC, color = color)) +
    geom_point() +
    scale_colour_manual(values = c("red","royalblue3","grey")) +
    geom_label_repel(aes(label = ID, size=1)) + 
    theme_minimal() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ggtitle(title) +
    geom_hline(yintercept = c(-logFC_criteria, logFC_criteria), 
             color = "dodgerblue", linetype = 'solid', linewidth = 0.7)
}

# plotMD_kegg -------------------------------------------------------------
plotMD_kegg <- function(file_name, ID, title="",logFC_criteria = 1, only_DE=T,
                        with_line=T){
  # 获取指定通路的基因列表
  pathway_df <- KEGGREST::keggGet(ID)[[1]]
  pathway_genes <- pathway_df$GENE
  cat(c(" -> pathway name:", pathway_df$NAME,"\n"))
  # 提取基因名称
  gene_names <- sapply(strsplit(pathway_genes, ";"), `[`, 1)
  # 删除奇数索引的元素
  gene_list_even <- gene_names[seq_along(gene_names) %% 2 == 0]
  cat(c(" ->",length(gene_list_even),"genes involved.\n"))
  
  # filter data
  cat(c(" -> ",file_name,"/n"))
  df <- get_df(file_name ,all = T) %>% group_by(SYMBOL) %>%
    filter(.,!(SYMBOL%in% c("havana","ensembl_havana","havana_tagene")),
           geneBiotype=="protein_coding") %>% 
    filter(.,SYMBOL %in% gene_list_even)%>%
    summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
    column_to_rownames(., var = "SYMBOL") %>% 
    rownames_to_column("ID") %>% 
    .[,c(1,5,6,7,8)] %>% 
    setNames(c("ID","treat","control","M","D"))%>% 
    mutate(Average_expression =(treat+control)/2) %>% 
    .[,c(1,4,6)] %>% setNames(c("ID","log2FC","Average_expression")) %>% 
    mutate(ID = ID, color = "3") %>% 
    mutate(color = if_else(log2FC > logFC_criteria, '1', color)) %>%
    mutate(color = if_else(log2FC < logFC_criteria & log2FC > negative(logFC_criteria), '3', color)) %>%   
    mutate(color = if_else(log2FC < negative(logFC_criteria), '2', color))
  
  up <- df %>% filter(color==1) %>% nrow()
  down <- df %>% filter(color==2) %>% nrow()
  non <- df %>% filter(color==3) %>% nrow()
  cat(c(" -> up:",up," -> down:",down, " -> non:",non,"\n"))
  if(only_DE==T){
    df <- df %>% filter(color==1|color==2)
  }
  if(up==0){
    
  }
  
  ### MD plot
  message("drawing MD plot")
  p <- ggplot(df, aes(x = log(Average_expression), y = log2FC, color = color)) +
    geom_point() +
    scale_colour_manual(values = c("1" = "red", "2" = "royalblue3", "3" = "grey")) +
    geom_label_repel(aes(label = ID, size=1)) + 
    theme_minimal() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ggtitle(paste(title, "  ",pathway_df$NAME))
  
  if(with_line==T){
    p + geom_hline(yintercept = c(-logFC_criteria, logFC_criteria), 
                   color = "dodgerblue", linetype = 'solid', linewidth = 0.7)
  } else {
    p
  }
}

# draw_bar ----------------------------------------------------------------
draw_bar <- function(gene){
  mat <- abbr_count %>%
    rownames_to_column("ENSEMBL") %>% 
    left_join(.,gene_df, by="ENSEMBL") %>% 
    filter(.,SYMBOL == gene) %>% 
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>%  
    column_to_rownames(., var = "SYMBOL") 
  mat_t <- t(mat) %>%  as.data.frame() %>% rownames_to_column("ID") %>% 
    setNames(c("ID","CPM"))
  mat_t$ID <- factor(mat_t$ID, levels = mat_t$ID)
  # 绘制条形图
  clone = c('WT' = '#0080FF', 'L858R' = '#00FF00',
            'Del19' = '#FF69B4', 'YAP' = '#FF8000')
  mat_t$clone <- rep(c("WT", "L858R", "Del19", "YAP"), each = 13)
  # 绘制条形图
  ggplot(mat_t, aes(x = ID, y = CPM, fill = clone)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = clone) +  # 使用自定义颜色
    labs(title = paste0(gene," (RNA-seq)"),
         x = "Groups",
         y = "Expression Level(CPM)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}
