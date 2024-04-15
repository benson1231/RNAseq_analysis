get_df <- function(file,
                   de=F,
                   list,
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
