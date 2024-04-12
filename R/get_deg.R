get_deg <- function(file=file,
                    log_crit = c(1,-1),
                    dir="all",  # all/up/down
                    type="ENSEMBL"
                    ){
  
  # 檢查 type 是否有效
  if (!(type %in% c("ENTREZID", "ENSEMBL"))) {
    stop("Invalid type. Allowed values are 'ENTREZID' or 'ENSEMBL'.")
  }
  
  cat(c(" -> load data from",file.path(data_path, file),"\n"))
  data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
  
  # filter log2FC criteria
  if(dir=="all"){
    data <- data %>% filter(M >log_crit[1]| M < log_crit[2])
    cat(c(" -> abs(log2FC) larger than", log_crit[1],"\n"))
  }else if(dir=="up"){
    data <- data %>% filter(M > log_crit[1])
    cat(c(" -> log2FC larger than", log_crit[1],"\n"))
  }else if(dir=="down"){
    data <- data %>% filter(M < log_crit[2])
    cat(c(" -> log2FC smaller than", log_crit[2],"\n"))
  }else{
    cat(c("<- error, check direction", "\n"))
    return(NULL)
  }
  
  # 抓出差異gene id
  if(type=="ENTREZID"){
    list <- data$ENTREZID %>% na.omit()
  }else{
    list <- data$ENSEMBL
  }
  
  return(list)
}
