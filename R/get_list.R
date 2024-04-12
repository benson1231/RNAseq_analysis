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
    select(all_of(id_col), M) %>%
    group_by(across(all_of(id_col))) %>%
    summarize(across(where(is.numeric), mean)) %>%
    arrange(desc(M)) 
  
  # 創建列表
  enrich_list <- df$M 
  names(enrich_list) <- df[[id_col]]
  enrich_list <-  na.omit(enrich_list)
  
  return(enrich_list)
}
