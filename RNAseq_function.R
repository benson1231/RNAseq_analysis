# function for visualization of RNA-seq data 
# draw_heatmap ------------------------------------------------------------
draw_heatmap <- function(file=file,
                         log_crit=3,
                         groups=groups,
                         show_row_names = TRUE
                         ){
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
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
  } else if(groups == "ALL"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
                               ip_Y_V_S_AS,ip_Y_V_S_CO,
                               ip_Y_V_S_LCD,ip_Y_V_S_HCD,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP,
                               ip_Y_V_S_CO_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else{
    cat(c("<- please type in groups", "\n"))
    return(NULL)
  }
  
  cat(c(" -> filtering data", "\n"))
  
  mat_scale <- data_mat %>% t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()
  cat(c(" -> annotation", "\n"))
  
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
  cat(c(" -> drawing heatmap", "\n"))
  
  ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                          show_row_names = show_row_names, 
                          show_column_names = F,
                          name = "Z-score"
  )
  
}

# draw_from_list ----------------------------------------------------------
draw_from_list <- function(list,
                           groups=groups,
                           id="ENSEMBL",
                           cluster=TRUE,
                           anno=TRUE,
                           title="",
                           show_row_names = TRUE# ENSEMBL/SYMBOL
){
  if(id=="ENSEMBL"){
    list <- list  
    mycount_df$ENSEMBL <- rownames(mycount_df)
    new_df <- mycount_df %>% left_join(.,gene_df, by="ENSEMBL")
    mat <- mycount_df %>% filter(.,rownames(mycount_df) %in% list) %>% 
      left_join(.,gene_df, by="ENSEMBL") %>% group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
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
    return(NULL)
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
  } else if(groups == "ALL"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
                               ip_Y_V_S_AS,ip_Y_V_S_CO,
                               ip_Y_V_S_LCD,ip_Y_V_S_HCD,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP,
                               ip_Y_V_S_CO_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else{
    cat(c("<- please type in groups", "\n"))
    reture(NULL)
  }
  cat(c(" -> filtering data", "\n"))
  
  mat_scale <- data_mat %>% t() %>% scale() %>% t() %>% as.matrix() %>% na.omit()
  cat(c(" -> annotation", "\n"))
  
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
                                             'AS'='#FFFF37','CO'='#FF7575','LCD'='#0080FF',
                                             'HCD'='#005AB5', 'BAP'='#00DB00', 'AS_BAP'='#FFD306', 'CO_BAP'= '#FF0000',
                                             'LCD_BAP'= '#0000E3', "HCD_BAP"="#000079"),
                                     clone=c('WT'="#FF2D2D",'L858R'="#FF9224",
                                             "DEL19"= "#66B3FF","YAP"="#2828FF")))
  cat(c(" -> drawing heatmap", "\n"))
  
  if(anno==T){
    ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                            show_row_names = show_row_names,
                            show_column_names = F,cluster_rows = cluster,
                            name = "Z-score", row_title = title
    )
  }else{
    ComplexHeatmap::Heatmap(mat_scale, cluster_columns = F, 
                            show_row_names = show_row_names, 
                            show_column_names = F,
                            cluster_rows = cluster,
                            name = "Z-score", row_title = title
    )
  }
}

# get_deg --------------------------------------------------------
get_deg <- function(file=file,
                    log_crit = c(1,-1),
                    dir="all"  # all/up/down
                    
){
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
  data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
  if(dir=="all"){
    data <- data %>% filter(M >log_crit[1]| M < log_crit[2])
  }else if(dir=="up"){
    data <- data %>% filter(M > log_crit[1])
  }else if(dir=="down"){
    data <- data %>% filter(M < log_crit[2])
  }else{
    cat(c("<- error, check direction", "\n"))
    return(NULL)
  }
  list <- data$ENSEMBL 
  return(list)# 抓出差異ensembl id
}

# gsea_run ----------------------------------------------------------------
gsea_run <- function(file,
                     all_gene=FALSE, 
                     list,
                     list_id="ensembl"   # ensembl/symbol
                     ){
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
  if(all_gene==TRUE){
    df <- readxl::read_xlsx(file.path(data_path, file))
    cat(c(" >- input all gene and run GSEA","\n"))
  } else {
    if(list_id=="ensembl"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        filter(ENSEMBL %in% list)
    } else if(list_id=="symbol"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        filter(SYMBOL %in% list)
    } else{
      cat(c(" <- error: check list_id","\n"))
      return(NULL)
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
  cat(c(" -> running GSEA", "\n"))
  gse <- gseGO(geneList = gsea_gene_list, 
               ont ="ALL", 
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
                     all_gene=FALSE,
                     list,
                     list_id="ensembl"   # ensembl/symbol
                     ){
  # select gene in list
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
  if(all_gene==TRUE){
    df <- readxl::read_xlsx(file.path(data_path, file))
    cat(c(" >- input all gene and run KEGG","\n"))
  } else {
    if(list_id=="ensembl"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        filter(ENSEMBL %in% list)
    } else if(list_id=="symbol"){
      df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
        filter(SYMBOL %in% list)
    } else{
      cat(c(" <- error: check list_id","\n"))
      return(NULL)
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

