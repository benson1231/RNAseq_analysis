### function
get_deg <- function(data_path=data_path,
                    file=file,
                    log_crit = c(1,-1),
                    dir="all"  # all/up/down
                    ){
  data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
  if(dir=="all"){
    data <- data %>% filter(abs(M)>log_crit[1])
  }else if(dir=="up"){
    data <- data %>% filter(M > log_crit[1])
  }else if(dir=="down"){
    data <- data %>% filter(abs(M)< log_crit[2])
  }else{
    cat(c("<- error, check direction", "\n"))
  }
  data <- data %>% filter(abs(M)>log_crit)  # filter log2FC criteria
  list <- data$ENSEMBL 
  return(list)# 抓出差異ensembl id
}

### run function
BAP_all <- get_deg(data_path = data_path,file = "ip_Y_V_S_BAP_0_deg.xlsx",
                   log_crit = c(1,-1),dir = "up")
