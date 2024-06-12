# function for RNA-seq data ---------------------------------------------------
# pre-def
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
only_CD <- c(9,7,8,12,13, 22,20,21,25,26, 35,33,34,38,39, 48,46,47,51,52)
BAP <- c(1,2,9,10,11,12,13, 14,15,22,23,24,25,26, 27,28,35,36,37,38,39, 40,41,48,49,50,51,52)
# abb ---------------------------------------------------------------------
abb <- function(name,type="group"){
  if(!(type%in%c("file","group"))){
    stop("'type' must be 'file' or 'group'")
  }
  if(type=="file"){
    df <- name_df %>% filter(file_name==name) %>% pull(abbreviate)
  } else{
    df <- name_df %>% filter(abbreviate==name) %>% pull(abbreviate)
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
draw_heatmap <- function(file=NULL,
                         log_crit=3,
                         groups="ALL",
                         show_row_names = FALSE,
                         log_scale=T,
                         row_km=T,
                         km=0,
                         return_cluster=F,
                         title="",
                         list=NULL
                         ){
  if(!(groups%in%c("AS","CO","LCD","HCD","CD","BAP","only_AS","only_CO",
                  "only_LCD","only_HCD","only_CD","CON","ALL"))){
    stop('"groups" must be "AS","CO","LCD","HCD","CD","BAP","only_AS","only_CO",
         "only_LCD","only_HCD","only_CD","CON" or "ALL".')
  }
  
  if(log_scale==T){
    CPM_df <- mylogCPM
  } else {
    CPM_df <- myCPM
  }

  mat <- CPM_df %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% 
    left_join(.,gene_df, by="ENSEMBL") %>% 
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>% 
    column_to_rownames(., var = "SYMBOL") %>% na.omit()
  
  if(!is.null(file)){
    cat(c(" -> load data from",file.path(data_path, file),"\n"))
    data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
    data <- data %>% filter(abs(M)>log_crit,
                            !(SYMBOL%in% c("havana","ensembl_havana","havana_tagene")))  
    cat(c(" -> log2FC criteria is", log_crit,"\n"))
    deg <- data$SYMBOL 
    mat <- mat %>% filter(rownames(.)%in% deg)  # 抓出差異ensembl id
    if(!is.null(list)){
      mat <- mat %>% filter(rownames(.)%in% list)
      cat(" -> filter file DEGs in 'list'.\n")
    } 
  } else{
      if(!is.null(list)){
        mat <- mat %>% filter(rownames(.)%in% list)
        cat(" -> filter genes in 'list'.\n")
      } else {
        stop(" error: both 'file' and 'list' missing.")
      }
  }
  
  # select groups
  if(groups == "AS"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[AS]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#33FF33',
                          'AS' = '#FFFF37', 'AS_BAP' = '#FFDC35'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "CO"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[CO]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#33FF33',
                          'CO' = '#FF5151', 'CO_BAP' = '#CC0000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "CD"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[CD]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#33FF33',
                          'LCD' = '#1E90FF', 'HCD' = '#0000CD', 
                          'LCD_BAP' = '#0088A8', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "BAP"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[BAP]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#33FF33',
                          'AS_BAP' = '#FFDC35', 'CO_BAP' = '#CC0000',
                          'LCD_BAP' = '#0088A8', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "ALL"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'AZA' = '#FFD2D2', 'DAC' = '#FFA488',
                          'BAP' = '#33FF33','AS' = '#FFFF37', 'CO' = '#FF5151',
                          'LCD' = '#1E90FF', 'HCD' = '#0000CD', 
                          'AS_BAP' = '#FFDC35', 'CO_BAP' = '#CC0000',
                          'LCD_BAP' = '#0088A8', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_CO"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[only_CO]))
    ann <- list(agent = c('BAP' = '#33FF33','CO' = '#FF5151', 'CO_BAP' = '#CC0000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_AS"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[only_AS]))
    ann <- list(agent = c('BAP' = '#33FF33', 'AS' = '#FFFF37', 'AS_BAP' = '#FFDC35'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_LCD"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[only_LCD]))
    ann <- list(agent = c('BAP' = '#33FF33', 'LCD' = '#1E90FF',  'LCD_BAP' = '#0088A8'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_HCD"){
    data_mat <- mat %>% dplyr::select(name_df$abbreviate[only_HCD])
    ann <- list(agent = c('BAP' = '#33FF33','HCD' = '#0000CD', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_CD"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[only_CD]))
    ann <- list(agent = c('BAP' = '#33FF33','LCD' = '#1E90FF', 'HCD' = '#0000CD',
                          'LCD_BAP' = '#0088A8', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "CON"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[CON]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else{
    return(NULL)
  }
  
  mat_scale <- data_mat %>% t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()
  
  ### heat-map argument
  col <- colnames(mat_scale)
  
  # agent name
  agent <- factor (
    str_replace_all(col, c("^W_|^L_|^D_|^Y_"='')),
    levels=c('CON','DMS',"AZA","DAC","BAP",'AS',"CO","LCD","HCD",
             "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
  
  # clone name
  clone <- factor( 
    str_replace_all(str_sub(col, 1,1), c("W" = "WT", "L" = "L858R", 
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
                           row_cluster=TRUE,
                           anno=TRUE,
                           title="",
                           show_row_names = TRUE, # ENSEMBL/SYMBOL
                           label_num=F,
                           log_scale=T,
                           col_cluster=FALSE,
                           row_km = 0,
                           filter=T
                           ){
  if(!(groups%in%c("AS","CO","LCD","HCD","CD","BAP","only_AS","only_CO",
                  "only_LCD","only_HCD","only_CD","CON","ALL"))){
    stop('"groups" must be "AS","CO","LCD","HCD","CD","BAP","only_AS","only_CO",
         "only_LCD","only_HCD","only_CD","CON" or "ALL".')
  }
  if(log_scale==T){
    CPM_df <- mylogCPM %>% as.data.frame()
  } else {
    CPM_df <- myCPM %>% as.data.frame()
  }
  cat(c(" -> input list with",length(list),"genes","\n"))
  if(id=="ENSEMBL"){
    list <- list  
    if(filter==F){
      mat <- CPM_df %>%
        mutate(.,ENSEMBL=rownames(CPM_df)) %>% 
        filter(.,rownames(CPM_df) %in% list)
    } else{
      mat <- CPM_df %>%
        mutate(.,ENSEMBL=rownames(CPM_df)) %>% 
        left_join(.,gene_df, by="ENSEMBL") %>% 
        filter(.,rownames(CPM_df) %in% list,
               !(SYMBOL%in% c("havana","ensembl_havana","havana_tagene"))) %>% 
        group_by(SYMBOL) %>%
        summarize(across(where(is.numeric), sum)) %>% na.omit() %>%  
        column_to_rownames(., var = "SYMBOL")
    }
  }else if(id=="SYMBOL"){
    list <- list  
    mat <- CPM_df %>% 
      mutate(.,ENSEMBL=rownames(CPM_df)) %>% 
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
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[AS]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#33FF33',
                          'AS' = '#FFFF37', 'AS_BAP' = '#FFDC35'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "CO"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[CO]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#33FF33',
                          'CO' = '#FF5151', 'CO_BAP' = '#CC0000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "CD"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[CD]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#33FF33',
                          'LCD' = '#1E90FF', 'HCD' = '#0000CD', 
                          'LCD_BAP' = '#0088A8', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "BAP"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[BAP]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'BAP' = '#33FF33',
                          'AS_BAP' = '#FFDC35', 'CO_BAP' = '#CC0000',
                          'LCD_BAP' = '#0088A8', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "ALL"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'AZA' = '#FFD2D2', 'DAC' = '#FFA488',
                          'BAP' = '#33FF33','AS' = '#FFFF37', 'CO' = '#FF5151',
                          'LCD' = '#1E90FF', 'HCD' = '#0000CD', 
                          'AS_BAP' = '#FFDC35', 'CO_BAP' = '#CC0000',
                          'LCD_BAP' = '#0088A8', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_CO"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[only_CO]))
    ann <- list(agent = c('BAP' = '#33FF33','CO' = '#FF5151', 'CO_BAP' = '#CC0000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_AS"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[only_AS]))
    ann <- list(agent = c('BAP' = '#33FF33', 'AS' = '#FFFF37', 'AS_BAP' = '#FFDC35'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_LCD"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[only_LCD]))
    ann <- list(agent = c('BAP' = '#33FF33', 'LCD' = '#1E90FF',  'LCD_BAP' = '#0088A8'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_HCD"){
    data_mat <- mat %>% dplyr::select(name_df$abbreviate[only_HCD])
    ann <- list(agent = c('BAP' = '#33FF33','HCD' = '#0000CD', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "only_CD"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[only_CD]))
    ann <- list(agent = c('BAP' = '#33FF33','LCD' = '#1E90FF', 'HCD' = '#0000CD',
                          'LCD_BAP' = '#0088A8', 'HCD_BAP' = '#000000'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else if(groups == "CON"){
    data_mat <- mat %>% dplyr::select(all_of(name_df$abbreviate[CON]))
    ann <- list(agent = c('CON' = '#E0E0E0', 'DMS' = '#ADADAD'),
                clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
                          'Del19' = '#0066FF', 'YAP' = '#CA8EFF'))
  } else{
    return(NULL)
  }
  
  # scaling
  mat_scale <- data_mat  %>%  t() %>% scale() %>% t() %>% as.matrix() %>% na.omit()
  
  # color
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  # heat-map argument
  col <- colnames(mat_scale)
  
  # agent name
  agent <- factor(
    str_replace_all(col, c("^W_|^L_|^D_|^Y_"='')),
    levels=c('CON','DMS',"AZA","DAC","BAP",'AS',"CO","LCD","HCD",
             "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
  
  # clone name
  clone <- factor( 
    str_replace_all(str_sub(col, 1,1), c("W" = "WT", "L" = "L858R", 
                                         "D" = "Del19", "Y" = "YAP")),
    levels=c('WT','L858R',"Del19","YAP"))
  
  # draw hp
  ha <- HeatmapAnnotation(agent = agent, clone = clone,
                          col = ann)
  
  if(anno==T){
    if(label_num==T){
      ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, 
                              cluster_columns = col_cluster, 
                              show_row_names = show_row_names,
                              show_column_names = F,cluster_rows = row_cluster,
                              col = col_fun, row_km = row_km ,
                              name = "Z-score", row_title = title,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.1f", mat_scale[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
    }else{
      ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, 
                              cluster_columns = col_cluster, 
                              show_row_names = show_row_names,col = col_fun,
                              show_column_names = F,cluster_rows = row_cluster,
                              name = "Z-score", row_title = title,row_km = row_km
                              )
    }
  } else{
    if(label_num==T){
      ComplexHeatmap::Heatmap(mat_scale, cluster_columns = col_cluster, 
                              show_row_names = show_row_names,col = col_fun,
                              show_column_names = F,cluster_rows = row_cluster,
                              name = "Z-score", row_title = title,row_km = row_km ,
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.1f", mat_scale[i, j]), x, y, gp = gpar(fontsize = 10))
                              })
    }else{
      ComplexHeatmap::Heatmap(mat_scale, cluster_columns = col_cluster, 
                              show_row_names = show_row_names,col = col_fun,
                              show_column_names = F,cluster_rows = row_cluster,
                              name = "Z-score", row_title = title,row_km = row_km 
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
  cat(c(" -> load data from", file,"\n"))
  raw_data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔

  # 處理資料
  data <- raw_data %>%
    as.data.frame() %>%
    filter(geneBiotype == "protein_coding") %>%
    dplyr::select(!!sym(type), M) %>%
    group_by(!!sym(type)) %>%
    summarize(across(where(is.numeric), sum))
  
  if(top=="all"){
    if(dir=="all"){
      data <- data %>% filter(M >log_crit[1]| M < log_crit[2])
      cat(c(" -> abs(log2FC) larger than", log_crit[1],"->",nrow(data),"genes.\n"))
    }else if(dir=="up"){
      data <- data %>% filter(M > log_crit[1])
      cat(c(" -> log2FC larger than", log_crit[1],"->",nrow(data),"genes.\n"))
    }else if(dir=="down"){
      data <- data %>% filter(M < log_crit[2])
      cat(c(" -> log2FC smaller than", log_crit[2],"->",nrow(data),"genes.\n"))
    }else{
      cat(c("<- error, check direction", "\n"))
      return(NULL)
    }
  } else{
    if(dir=="all"){
      data <- data %>% filter(M >log_crit[1]| M < log_crit[2]) %>% 
        arrange(desc(abs(M))) %>% .[1:top,]
      cat(c(" -> abs(log2FC) larger than",log_crit[1],"->",nrow(data),"genes.\n"))
    }else if(dir=="up"){
      data <- data %>% filter(M > log_crit[1]) %>% 
        arrange(desc(abs(M))) %>% .[1:top,]
      cat(c(" -> log2FC larger than", log_crit[1],"->",nrow(data),"genes.\n"))
    }else if(dir=="down"){
      data <- data %>% filter(M < log_crit[2]) %>% 
        arrange(desc(abs(M))) %>% .[1:top,]
      cat(c(" -> log2FC smaller than", log_crit[2],"->",nrow(data),"genes.\n"))
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
  if (!(type %in% c("ENTREZID", "ENSEMBL","SYMBOL"))) {
    stop("Invalid type. Allowed values are 'ENTREZID', 'ENSEMBL' or 'SYMBOL'.")
  }
  # 讀取資料
  raw_df <- readxl::read_xlsx(file.path(data_path, file_name)) %>% as.data.frame()
  
  # 根據不同的類型選擇相應的欄位名稱
  if (type == "ENTREZID") {
    id_col <- "ENTREZID"
  } else if (type == "ENSEMBL"){
    id_col <- "ENSEMBL"
  } else {
    id_col <- "SYMBOL"
  }
  
  # 處理資料
  df <- raw_df %>% as.data.frame() %>% 
    filter(geneBiotype=="protein_coding") %>%
    dplyr::select(all_of(id_col),M) %>%
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
    cat(c(" -> raw data from ",file,"\n"))
    return(df)
  }
  
  if(de==TRUE){
    if(dir=="all"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(geneBiotype=="protein_coding") %>%
        filter(abs(M) > log_crit[1]) %>% 
        dplyr::select(M, SYMBOL) %>% 
        setNames(c("logFC","gene"))
    }else if(dir=="up"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(geneBiotype=="protein_coding") %>%
        filter(M > log_crit[1]) %>% 
        dplyr::select(M, SYMBOL) %>% 
        setNames(c("logFC","gene"))
    }else if(dir=="down"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(geneBiotype=="protein_coding") %>%
        filter(M < log_crit[2]) %>% 
        dplyr::select(M, SYMBOL) %>% 
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
        filter(geneBiotype=="protein_coding") %>%
        filter(ENSEMBL %in% list) %>% 
        dplyr::select(M,D, SYMBOL) %>% 
        setNames(c("logFC","D","gene"))
    } else{
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        as.data.frame() %>% 
        filter(geneBiotype=="protein_coding") %>%
        filter(ENSEMBL %in% list) %>% 
        dplyr::select(M, SYMBOL) %>% 
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
get_kegg_list <- function(path_ID, return_list=T
                          ){
  # 获取指定通路的基因列表
  pathway_df <- KEGGREST::keggGet(path_ID)[[1]]
  pathway_genes <- pathway_df$GENE
  # 提取基因名称
  gene_names <- sapply(strsplit(pathway_genes, ";"), `[`, 1)
  # 删除奇数索引的元素
  gene_list_even <- gene_names[seq_along(gene_names) %% 2 == 0]
  cat(c(" ->",length(gene_list_even),"genes involved.\n"))
  if(return_list==T){
    cat(c(" -> pathway name:", pathway_df$NAME,"\n"))
    return(gene_list_even)
  } else {
    return(pathway_df$NAME)
  }
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
get_p <- function(x, abbreviate, top=0, p_value_cutoff = 0.01){
  pvalue <- x@result %>% .[,c("Description","pvalue")] %>% 
    filter(pvalue < p_value_cutoff) %>% 
    setNames(c("Description",abbreviate))
  if(top != 0){
    pvalue <- pvalue %>% head(top) 
  }
  return(pvalue)
}

# plot_heatmap ------------------------------------------------------------
plot_heatmap <- function(file_list, abbreviates, analysis, dir="up",
                         col_title="", row_title="",top=20) {
  # 檢查參數是否匹配
  if (length(file_list) != length(abbreviates)) {
    stop("file_list and abbreviates must have the same length")
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
  processed_dfs <- mapply(get_p, dfs, abbreviate = abbreviates, SIMPLIFY = FALSE)
  
  # 合併所有數據框
  combined_df <- Reduce(function(x, y) full_join(x, y, by = "Description"), processed_dfs) %>% 
    column_to_rownames("Description") %>% as.matrix()
  # 計算na值
  non_NA_count <- rowSums(!is.na(combined_df)) %>% sort(.,decreasing = T)
  
  # 計算 -log10(P) 值
  mat <- log10(combined_df) %>% -.
  mat <- mat[names(non_NA_count), ] %>% head(top)
  rownames(mat) <- factor(rownames(mat), levels = names(non_NA_count))
  
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
run_TF <- function(file_name, title="", n_tfs=25, plot=T){
  ### two groups comparing TFs
  cat(c(" -> read file",file_name,"\n"))
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
    filter(source %in% tfs) %>% 
    arrange(desc(score))
  if(plot==F){
    return(f_contrast_TF)
  }
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
plot_TF <- function(file_name, TF, title="", logFC_criteria = 1, only_DE=F,
                    with_line=T, plot_type=1, plot=T){
  if(!(plot_type%in%c(1,2))){
    stop(" -> 'plot_type' must be '1' or '2'. 1: MD plot. 2: mutation rate M plot.")
  }
  cat(" -> TF:",TF,"\n")
  tf <- net_TF %>% setNames(c("ID","regulon","mor"))
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
    mutate(ID = target, color = "non", dir="4") %>%
    column_to_rownames('target') %>% 
    left_join(de,"ID") %>% 
    mutate(color = if_else(logFC > logFC_criteria, 'up', color)) %>%
    mutate(color = if_else(logFC < logFC_criteria & logFC > negative(logFC_criteria), 'non', color)) %>%   
    mutate(color = if_else(logFC < negative(logFC_criteria), 'down', color)) %>% 
    mutate(dir = ifelse(mor > 0 & color == 'up', "1", dir)) %>% 
    mutate(dir = ifelse(mor < 0 & color == 'down', "2", dir)) %>% 
    mutate(dir = ifelse(mor > 0 & color == 'down', "3", dir)) %>%  
    mutate(dir = ifelse(mor < 0 & color == 'up', "3", dir))
  if(only_DE==T){
    df <- df %>% filter(color != 'non')
    cat(" -> only plot DEG\n")
  }
  ### MD plot
  if(plot_type==1){
    p <- ggplot(df, aes(x = log(Average_expression), y = logFC, color = dir, size=abs(mor))) +
      geom_point() +
      scale_colour_manual(values = c("1" = "red", "2" = "royalblue3",
                                     "3"="#CA8EFF","4"="gray")) +
      geom_label_repel(aes(label = ID, size=1)) + 
      theme_minimal() +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0, linetype = 'dotted') +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      ggtitle(paste(title, "Regulons of",tf))
  } else{
    df <- df %>% left_join(muta_rate,by = "ID")
    p <- ggplot(df, aes(x = `mutation_ratio(%)`, y = logFC, color = dir)) +
      geom_point() +
      scale_colour_manual(values = c("1" = "red", "2" = "royalblue3",
                                     "3"="#CA8EFF","4"="gray")) +
      geom_label_repel(aes(label = ID, size=1)) + 
      theme_minimal() +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0, linetype = 'dotted') +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      ggtitle(paste(title, "Regulons of",tf))
  }
  if(with_line==T){
    p <- p + geom_hline(yintercept = c(-logFC_criteria, logFC_criteria), 
                        color = "dodgerblue", linetype = 'solid', linewidth = 0.7)
  }
  if(plot==T){
    print(p)
  } else{
    return(df)
  }
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
    dplyr::select(SYMBOL, direaction)
  output_dir <- output_dir
  write.xlsx(df, file.path(output_dir,output_name))
  cat(c(" -> ",output_dir,"\n -> ",output_name))
}


# plot_MD -----------------------------------------------------------------
plot_MD <- function(file_name, title="", logFC_criteria = 1, only_DE=T,
                    list=c(),with_line=T,plot_type=1){
  if(!(plot_type%in%c(1,2))){
    stop(" -> 'plot_type' must be '1' or '2'. 1: MD plot. 2: logFC and TCGA-mutation-rate plot.")
  }
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
  # if list was given, filter genes in list 
  if(length(list)!=0){
    df <- df %>% filter(ID %in% list)
    cat(c(" -> filter ",length(list)," genes in input list.\n"))
  }
  # calculate number of up/down regulation
  up <- df %>% filter(color==1) %>% nrow()
  down <- df %>% filter(color==2) %>% nrow()
  non <- df %>% filter(color==3) %>% nrow()
  cat(c(" -> up:",up," -> down:",down, " -> non:",non,"\n"))
  if(only_DE==T){
    df <- df %>% filter(color == 1 | color == 2)
    cat(" -> only plot DEG\n")
  }
  ### MD plot
  if(plot_type==1){
    cat(" -> plotting MD plot\n")
    p <- ggplot(df, aes(x = log(Average_expression), y = logFC, color = color)) +
      geom_point() +
      scale_colour_manual(values = c("1" = "red", "2" = "royalblue3", "3" = "grey")) +
      geom_label_repel(aes(label = ID, size=1)) + 
      theme_minimal() +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0, linetype = 'dotted') +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      ggtitle(title) 
  } else{
    df <- df %>% left_join(muta_rate,by = "ID")
    cat(" -> plotting logFC and TCGA-mutation-rate plot\n")
    p <- ggplot(df, aes(x = `mutation_ratio(%)`, y = logFC, color = color))  +
      geom_point() +
      scale_colour_manual(values = c("1" = "red", "2" = "royalblue3", "3" = "grey")) +
      geom_label_repel(aes(label = ID, size=1)) + 
      theme_minimal() +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0, linetype = 'dotted') +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      ggtitle(title) 
  }
  
  if(with_line==F){
    p
  } else{
    p + geom_hline(yintercept = c(-logFC_criteria, logFC_criteria), 
                  color = "dodgerblue", linetype = 'solid', linewidth = 0.7)
  }
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
    filter(geneBiotype=="protein_coding") %>%
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
    cat(" -> only plot DEG")
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
draw_bar <- function(gene, group="ALL"){
  if(!(group %in% c("AS","CO","CD","ALL","CON"))){
    stop(" -> group must be 'AS','CO','CD','ALL' or 'CON'.")
  }
  mat <- myCPM %>%
    rownames_to_column("ENSEMBL") %>% 
    left_join(.,gene_df, by="ENSEMBL") %>% 
    filter(.,SYMBOL == gene) %>% 
    group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>%  
    column_to_rownames(., var = "SYMBOL") 
  if(group=="AS"){
    mat <- mat %>% dplyr::select(name_df$abbreviate[AS])
  } else if(group=="CO"){
    mat <- mat %>% dplyr::select(name_df$abbreviate[CO])
  } else if(group=="CD"){ 
    mat <- mat %>% dplyr::select(name_df$abbreviate[CD])
  } else if(group=="CON"){
    mat <- mat %>% dplyr::select(name_df$abbreviate[CON])
  } else {
    mat
  }
  mat_t <- t(mat) %>%  as.data.frame() %>% rownames_to_column("ID") %>% 
    setNames(c("ID","CPM"))
  mat_t$ID <- factor(mat_t$ID, levels = mat_t$ID)
  # 绘制条形图
  clone = c('WT' = '#AAAAAA', 'L858R' = '#FF4500',
            'Del19' = '#0066FF', 'YAP' = '#CA8EFF')
  mat_t$clone <- rep(c("WT", "L858R", "Del19", "YAP"), each = nrow(mat_t)/4)
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


# draw_survival -----------------------------------------------------------
draw_muta_survival <- function(gene){
  mutation_samples <- subset(maf_object@data, Hugo_Symbol == gene)$Tumor_Sample_Barcode
  # 提取了樣本ID的前12個字符，用於後續臨床樣本配對
  mutation_samples_short <- substr(mutation_samples, 1, 12)
  length(mutation_samples_short)
  head(mutation_samples_short)
  # 新增紀錄突變狀態的欄位
  clin.LUAD$mutation_status <- ifelse(clin.LUAD$submitter_id %in% mutation_samples_short, "Mutated", "Wildtype")
  # plotting
  TCGAbiolinks::TCGAanalyze_survival(
    data = clin.LUAD,
    clusterCol = "mutation_status",
    main = "TCGA Set\n LUAD",
    height = 10,
    width=10,
    legend = gene, 
    filename = paste0("output/",gene,"_survival.pdf")
  )
}

# draw_TCGA_survival ------------------------------------------------------
draw_TCGA_survival <- function(gene, population = ALL){
  if(population=="ALL"){
    tcga_count <- tcga_count
  } else if(population=="female"){
    list <- clin.LUAD %>% filter(gender=="female") %>% pull(submitter_id)
    tcga_count <- tcga_count %>% .[,(substr(colnames(tcga_count), 1, 12) %in% list)]
  } else if(population=="asian"){
    list <- clin.LUAD %>% filter(race=="asian") %>% pull(submitter_id)
    tcga_count <- tcga_count %>% .[,(substr(colnames(.), 1, 12) %in% list)]
  }
  exp <- tcga_count[rownames(tcga_count) == gene, ] %>% t() %>% 
    as.data.frame() %>% setNames("gene_expression") %>% 
    mutate(ID=substr(rownames(.), 1, 12)) %>% 
    left_join(clin.LUAD, by = c("ID" = "submitter_id")) %>% 
    mutate(gene_exp_level = ifelse(gene_expression >= median(.[["gene_expression"]]), "high","low"))
  
  TCGAanalyze_survival(
    data = exp,
    clusterCol = "gene_exp_level",
    main = paste("TCGA",population),
    height = 10,
    width=10,
    legend = gene,
    filename = paste("output/",gene,"survival.pdf") 
  )
}

# multi_venn --------------------------------------------------------------
multi_venn <- function(file_list, dir="up", upset=FALSE, title=""){
  # get 4 clone
  W <- get_deg(file_list[1], log_crit = c(1,-1), dir = dir,type = "SYMBOL")
  L <- get_deg(file_list[2], log_crit = c(1,-1), dir = dir,type = "SYMBOL")
  D <- get_deg(file_list[3], log_crit = c(1,-1), dir = dir,type = "SYMBOL")
  Y <- get_deg(file_list[4], log_crit = c(1,-1), dir = dir,type = "SYMBOL")
  # venn diagram list
  venn_list <- list(W = W, L = L, D = D, Y = Y)
  # venn diagram
  if(dir=="up"){
    if(upset==T){
      p <- ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0, force_upset = T)
    } else {
      p <- ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0) +
        scale_fill_gradient(low="white",high = "#FF2D2D")+ ggtitle(title)
    }
  } else if(dir=="down"){
    if(upset==T){
      p <- ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0, force_upset = T)
    } else {
      p <- ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0) +
        scale_fill_gradient(low="white",high = "#6A6AFF")+ ggtitle(title)
    }
  } else {
    if(upset==T){
      p <- ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0, force_upset = T)
    } else {
      p <- ggVennDiagram(venn_list,label_percent_digit = 1,label_alpha = 0) +
        scale_fill_gradient(low="white",high = "orange")+ ggtitle(title)
    }
  }
  print(p)
  # get venn list
  venn_result <- process_region_data(Venn(venn_list))
  return(venn_result)
}

# draw_TCGA_boxplot ------------------------------------------------------
draw_TCGA_boxplot <- function(gene){
  mat <- tcga_count %>% filter(rownames(.) == gene) %>% 
    t() %>% as.data.frame() %>% rownames_to_column("sample") %>% 
    setNames(c("sample","count")) %>% mutate(ID=substr(sample, 1, 12)) %>% 
    mutate(group="") %>% 
    mutate(group=ifelse(ID %in% asian_id,"asian",group)) %>% 
    mutate(group=ifelse(ID %in% non_asian_id,"non_asian",group))
  # 提取非亚洲人和亚洲人的样本计数
  non_asian_counts <- mat$count[mat$group == "non_asian"]
  asian_counts <- mat$count[mat$group == "asian"]
  # 执行 t 检验
  t_test_result <- t.test(non_asian_counts, asian_counts)
  # 绘制基因表达量的箱线图
  p <- ggplot(mat, aes(x = group, y = log(count))) +
    geom_boxplot()+
    labs(title = paste0(gene," (TCGA RNA expression)"),
         x = "Groups",
         y = "Expression Level (log.count)") +
    theme_minimal()
  # 在箱线图上添加 t 检验结果
  p + geom_text(aes(x = 1.25, y = 7, 
                    label = paste("p-value:", format(t_test_result$p.value, digits = 3))),
                hjust = 0, vjust = 0, color = "black")
}

# draw_cell_boxplot --------------------------------------------------------
draw_cell_boxplot <- function(gene=NULL, by=NULL, sig=NULL, top=50, 
                              title=NULL){
  if(is.null(sig)){
    if(!(by %in% c("Gender","Smoking_Status","Stage","EGFR_Status"))){
      stop("error: 'by' must be 'Gender','Smoking_Status','Stage','EGFR_Status'")
    }
    df <- cell_count %>% dplyr::filter(rownames(.)==gene) %>% t() %>%  
      as.data.frame() %>% 
      setNames("count") %>% rownames_to_column("ID") %>% left_join(cell_info,"ID")
    # 根據by參數設置因子
    df[[by]] <- as.factor(df[[by]])
    # 執行單因素ANOVA
    formula <- as.formula(paste("count ~", `by`))
    anova_result <- aov(formula, data = df)
    
    # 如果有顯著差異，進行Tukey HSD檢定
    tukey_result <- TukeyHSD(anova_result)
    
    # 提取與第一個組別比較的p值
    first_group <- levels(df[[by]])[1]
    comparisons <- sapply(levels(df[[by]])[-1], function(group) {
      ifelse(first_group < group, paste0(group, "-", first_group), paste0(group, "-", first_group))
    })
    p_values <- sapply(comparisons, function(comp) {
      tukey_result[[by]][comp, "p adj"]
    })
    
    max_counts <- df %>%
      group_by(!!sym(by)) %>%
      summarise(max_count = max(count))
    # 創建包含 p 值的數據框
    p_labels <- data.frame(
      Group = levels(df[[by]])[-1],
      p_value = round(p_values, 4)
    ) %>%
      rename(!!by := Group) %>% 
      left_join(max_counts, by = by) %>%
      mutate(y_position = max_count + 0.5)
    
    ggplot(df, aes(x = !!sym(by), y = count)) +
      geom_boxplot() +
      labs(title = paste0(gene, " (RNA expression in Chen, Yi-Ju et al.)"),
           x = "Groups",
           y = "log2FC") +
      theme_minimal() + 
      geom_hline(yintercept = 0, color = "dodgerblue", 
                 linetype = 'solid', linewidth = 0.7) +
      geom_hline(yintercept = 1, color = "red", 
                 linetype = 'solid', linewidth = 0.7) +
      geom_text(data = p_labels, aes(x = !!sym(by), y = y_position, 
                                     label = paste0("p = ", p_value)),
                vjust = -0.5, color = "black")
    
  } else{
    cat(c(" -> get",length(sig), "gene list for signature\n"))
    # 过滤数据，提取所需基因
    cell_info <- cell_info %>% as.data.frame()
    filtered_counts <- cell_count %>% dplyr::filter(rownames(.) %in% sig)
    long_data <- filtered_counts %>% as.data.frame() %>% mutate(gene=rownames(.)) %>% 
      pivot_longer(cols = -gene, names_to = "ID", values_to = "value") %>% 
      left_join(cell_info,"ID") 
    
    sorted_data <- long_data %>%
      mutate(Gender = factor(Gender, levels = c("Male", "Female"))) %>% 
      mutate(Smoking_Status = factor(Smoking_Status, levels = c("Nonsmoke", "Ex-smoker","Current_Smoker"))) %>% 
      mutate(Stage = factor(Stage, levels = c("stage I", "stage II","stage III","stage IV"))) %>%
      mutate(EGFR_Status = factor(EGFR_Status, levels = c("WT", "L858R","exon19del","L858R.exon19del","others"))) %>% 
      arrange(!!sym(by))
    sorted_data$ID <- factor(sorted_data$ID, levels = unique(sorted_data$ID))
    
    if(by=="Gender"){
      col <- c('Male' = 'white', 'Female' = 'grey40')
    } else if(by=="Smoking_Status"){
      col <- c("Nonsmoke"='white', "Ex-smoker"='grey',"Current_Smoker"='black')
    } else if(by=="Stage"){
      col <- c('stage I' = 'white','stage II' = 'grey85','stage III' = 'grey60','stage IV' = 'black')
    } else if(by=='EGFR_Status'){
      col <- c('WT' = 'grey', 'L858R' = '#FF4500',
               'exon19del' = '#0066FF','L858R.exon19del'='yellow','others' = 'black')
    }
    
    # Create the boxplot plot
    ggplot(sorted_data, aes(x = ID, y = value, fill = !!sym(by))) +
      geom_boxplot() +
      labs(title = paste(title ,"signature in Chen, Yi-Ju et al."),
           x = "patients",
           y = "Expression Value (log2FC)") +
      theme_minimal() +  
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_continuous(limits = c(-4, 4)) +
      geom_hline(yintercept = 0, 
                 color = "black", linetype = 'solid', linewidth = 0.7)+
      geom_hline(yintercept = 1, 
                 color = "red", linetype = 'solid', linewidth = 0.7) + 
      scale_fill_manual(values = col)
  }
}


# draw_Spearman_bar -------------------------------------------------------
draw_Spearman_bar <- function(sig_num=NULL,group="ALL",top=50,count_df="my"){
  if(!(count_df %in% c("my","cell"))){
    stop(" error: count must be 'my' or 'cell'\n")
  }
  # 函数：计算组别之间的Spearman相关系数
  compute_spearman <- function(group1, group2) {
    cor(group1, group2, method = "spearman")
  }
  
  if(count_df == "my"){
    sig <- get_deg(name_df$file_name[sig_num], type = "ENSEMBL") %>% head(top)
    cat(c(" -> get",name_df$abbreviate[sig_num], "top",top,"signature\n"))
    filtered_counts <- abbr_count %>% filter(rownames(.) %in% sig)
    group1 <- filtered_counts[,name_df$abbreviate[sig_num]]
  } else{
    my_df <- readxl::read_excel(file.path(data_path,name_df$file_name[sig_num])) %>% 
      dplyr::select(SYMBOL,M) %>% setNames(c("SYMBOL",name_df$abbreviate[sig_num]))
    sig <- get_deg(name_df$file_name[sig_num], type = "SYMBOL") %>% head(top)
    cat(c(" -> get",name_df$abbreviate[sig_num], "top",top,"signature\n"))
    filtered_counts <- cell_count %>% rownames_to_column("SYMBOL") %>% 
    left_join(.,my_df,by="SYMBOL") %>%filter(SYMBOL %in% sig) %>% 
      column_to_rownames("SYMBOL")
    group1 <- filtered_counts[,name_df$abbreviate[sig_num]]
    filtered_counts[,name_df$abbreviate[sig_num]] <- NULL
  }
  
  # 定义组别名称列表
  abbreviates <- colnames(filtered_counts)
  
  # 初始化存储相关系数的矩阵
  spearman_corr_matrix <- matrix(nrow = length(abbreviates), ncol = 1)
  rownames(spearman_corr_matrix) <- abbreviates
  colnames(spearman_corr_matrix) <- "signature"
  
  # 计算组别之间的Spearman相关系数
  for (i in 1:length(abbreviates)) {
    group2 <- filtered_counts[, abbreviates[i]]
    spearman_corr_matrix[i, 1] <- compute_spearman(group1, group2)
  }
  
  if(count_df=="my"){
    data <- spearman_corr_matrix %>% as.data.frame() %>% setNames("Spearman") %>% 
      rownames_to_column("group") %>% mutate(clone=substr(group,1,1))
    data$clone <- factor(data$clone, levels = c("W","L","D","Y"))
    data$group <- factor(data$group, levels = name_df$abbreviate)
    clone = c('W' = '#AAAAAA', 'L' = '#FF4500',
              'D' = '#0066FF', 'Y' = '#CA8EFF')
    
    if(group=="AS"){
      data <- data %>% dplyr::filter(group %in% name_df$abbreviate[c(5,18,31,44)])
    } else if(group=="CO"){
      data <- data %>% dplyr::filter(group %in% name_df$abbreviate[c(6,19,32,45)])
    } else if(group=="CD"){ 
      data <- data %>% dplyr::filter(group %in% name_df$abbreviate[c(7,8,20,21,33,34,46,47)])
    } else if(group=="ALL"){
      data
    }
    
    ggplot(data, aes(x = group, y = Spearman, fill = clone)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = clone) +  # 使用自定义颜色
      labs(title = paste(name_df$abbreviate[sig_num],"signature (RNA-seq)"),
           x = "Groups",
           y = "Spearman's correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else{
    data <- spearman_corr_matrix %>% as.data.frame() %>% setNames("Spearman") %>% 
      rownames_to_column("ID") %>% left_join(cell_info,by = "ID") %>% 
      arrange(desc(Spearman))
    
    w <- data %>% filter(EGFR_Status=="WT") %>% pull(ID)
    l <- data %>% filter(EGFR_Status=="L858R") %>% pull(ID)
    d <- data %>% filter(EGFR_Status=="exon19del") %>% pull(ID)
    dl <- data %>% filter(EGFR_Status=="L858R.exon19del") %>% pull(ID)
    other <- data %>% filter(EGFR_Status=="others") %>% pull(ID)
    
    EGFR_Status <- c('WT' = 'grey', 'L858R' = '#FF4500',
                     'exon19del' = '#0066FF','L858R.exon19del'='yellow','others' = 'black')
    data$ID <- factor(data$ID,levels = c(w,l,d,dl,other))
    
    ggplot(data, aes(x = ID, y = Spearman, fill = EGFR_Status)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = EGFR_Status) +  # 使用自定义颜色
      labs(title = paste(name_df$abbreviate[sig_num],"signature (RNA-seq)"),
           x = "Groups",
           y = "Spearman's correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
}


# draw_boxplot -------------------------------------------------------------
draw_boxplot <- function(sig_num,group="ALL", top=50){
  sig <- get_deg(name_df$file_name[sig_num], type = "ENSEMBL",dir = "up") %>% 
    head(top)
  cat(c(" -> get",name_df$abbreviate[sig_num], "signature\n"))
  
  # 过滤数据，提取所需基因
  filtered_counts <- abbr_count %>% filter(rownames(.) %in% sig)
  
  if(group=="AS"){
    mat <- filtered_counts %>% dplyr::select(name_df$abbreviate[AS])
    long_data <- mat %>% as.data.frame() %>% mutate(gene=rownames(.)) %>% 
      pivot_longer(cols = -gene, names_to = "condition", values_to = "value")
    long_data$condition <- factor(long_data$condition,levels = name_df$abbreviate[AS])
  } else if(group=="CO"){
    mat <- filtered_counts %>% dplyr::select(name_df$abbreviate[CO])
    long_data <- mat %>% as.data.frame() %>% mutate(gene=rownames(.)) %>% 
      pivot_longer(cols = -gene, names_to = "condition", values_to = "value")
    long_data$condition <- factor(long_data$condition,levels = name_df$abbreviate[CO])
  } else if(group=="CD"){ 
    mat <- filtered_counts %>% dplyr::select(name_df$abbreviate[CD])
    long_data <- mat %>% as.data.frame() %>% mutate(gene=rownames(.)) %>% 
      pivot_longer(cols = -gene, names_to = "condition", values_to = "value")
    long_data$condition <- factor(long_data$condition,levels = name_df$abbreviate[CD])
  } else {
    mat <- filtered_counts
    long_data <- mat %>% as.data.frame() %>% mutate(gene=rownames(.)) %>% 
      pivot_longer(cols = -gene, names_to = "condition", values_to = "value")
    long_data$condition <- factor(long_data$condition,levels = name_df$abbreviate)
  }
  
  long_data <- long_data %>% mutate(clone=substr(condition,1,1))
  long_data$clone <- factor(long_data$clone, levels = c("W","L","D","Y"))
  clone = c('W' = '#AAAAAA', 'L' = '#FF4500',
            'D' = '#0066FF', 'Y' = '#CA8EFF')
  
  # Create the violin plot
  ggplot(long_data, aes(x = condition, y = log(value), fill = clone)) +
    geom_boxplot() +
    labs(title = paste(name_df$abbreviate[sig_num],"signature"),
         x = "Groups",
         y = "Gene expression") +
    theme_minimal() +
    scale_fill_manual(values = clone) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
