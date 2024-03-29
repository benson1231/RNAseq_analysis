gsea_run <- function(data_path=data_path,
                     file,
                     all_gene=FALSE, 
                     list,
                     list_id="ensembl"   # ensembl/symbol
                     ){
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
    cat(c(" >- input selected gene and run ORA","\n"))
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
  cat(c(" -> running KEGG", "\n"))
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
