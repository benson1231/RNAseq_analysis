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
