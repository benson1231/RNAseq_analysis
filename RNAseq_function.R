# function for visualization of RNA-seq data 
# draw_heatmap ------------------------------------------------------------
draw_heatmap <- function(file=file,
                         log_crit=3,
                         groups="ALL",
                         show_row_names = TRUE
                         ){
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
  cat(c(" -> filtering data", "\n"))
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
    cat(c("<- check groups", "\n"))
    return(NULL)
  }
  
  mat_scale <- data_mat %>% t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()
  
  ### heat-map argument
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
  
  ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                          show_row_names = show_row_names, 
                          show_column_names = F,
                          name = "Z-score"
  )
  
}

# draw_from_list ----------------------------------------------------------
draw_from_list <- function(list,
                           groups="ALL",
                           id="ENSEMBL",
                           cluster=TRUE,
                           anno=TRUE,
                           title="",
                           show_row_names = TRUE# ENSEMBL/SYMBOL
                           ){
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
  
  # scaling
  mat_scale <- data_mat %>% t() %>% scale() %>% t() %>% as.matrix() %>% na.omit()
  
  # heat-map argument
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
                    dir="all",  # all/up/down
                    type="ENSEMBL",
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
  df <- raw_df %>%
    select(all_of(id_col),M) %>%
    group_by(across(all_of(id_col))) %>%
    summarize(across(where(is.numeric), mean)) %>%
    arrange(desc(M)) %>% 
    setNames(c(all_of(id_col),"logFC"))
  
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
                   log_crit=c(1,-1)
                   ){
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
    df <- readxl::read_xlsx(file.path(data_path, file)) %>% 
      as.data.frame() %>% 
      filter(ENSEMBL %in% list) %>% 
      select(M, SYMBOL) %>% 
      setNames(c("logFC","gene"))
    
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
get_kegg_list <- function(path_name
                          ){
  # 獲取指定通路的基因列表
  pathway_df <- KEGGREST::keggGet(path_name)[[1]]
  pathway_genes <- pathway_df$GENE
  cat(c(" -> pathway name:", pathway_df$NAME,"\n"))
  # 提取基因名
  gene_names <- sapply(strsplit(pathway_genes, ";"), `[`, 1)
  # 刪除奇數索引之元素
  gene_list_even <- gene_names[seq_along(gene_names) %% 2 == 0]
  cat(c(" ->",length(gene_list_even),"genes involved.\n"))
  return(gene_list_even)
}

