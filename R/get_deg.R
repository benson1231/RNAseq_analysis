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
