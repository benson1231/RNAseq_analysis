# library -----------------------------------------------------------------
library(SummarizedExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)

# # get TCGA-LUAD RNA-seq data
# query_rna <- GDCquery(
#   project = "TCGA-LUAD",
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "STAR - Counts"
# )
# # 下载RNA-Seq数据
# GDCdownload(query_rna)
# rna_data <- GDCprepare(query_rna)
# saveRDS(rna_data,"rna_data.RDS")

# load TCGA-LUAD RNA-seq data ----------------------------------------
rna_data <- "/Users/benson/Documents/project/RNA-seq1-3/data/rna_data.RDS" %>% 
  readRDS()
# load TCGA-LUAD patient clinical data
clin.LUAD  <- "/Users/benson/Documents/project/RNA-seq1-3/data/clin_LUAD.RDS" %>% 
  readRDS()

# 假设种族信息在 "race_list" 列
# 过滤出东亚人群的数据
asian_id <- clin.LUAD  %>% 
  as.data.frame()%>% 
  filter(race=="asian" & gender=="female")%>% 
  pull(submitter_id)
non_asian_id <- clin.LUAD  %>% 
  as.data.frame()%>% 
  filter(race!="asian" & gender=="female") %>% 
  pull(submitter_id)

# 从RNA-Seq数据中提取东亚人群的样本
tcga_rna <- rna_data[, substr(colnames(rna_data), 1, 12) %in% c(asian_id, non_asian_id)]
head(tcga_rna)

# 提取基因表达量
gene_expression <- SummarizedExperiment::assay(tcga_rna)
head(gene_expression)

# get annotation
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = substr(rownames(gene_expression), 1, 15), 
                                      columns = "SYMBOL", keytype = "ENSEMBL")

tcga_count <- gene_expression %>% as.data.frame() %>% mutate(ENSEMBL=substr(rownames(.), 1, 15)) %>% 
  left_join(.,gene_symbols,by = "ENSEMBL",relationship = "many-to-many") %>% 
  group_by(SYMBOL) %>% summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
  column_to_rownames("SYMBOL")
# saveRDS(tcga_count,"data/tcga_count.RDS")

# run TCGA-LUAD
tcga_count <- "/Users/benson/Documents/project/RNA-seq1-3/data/tcga_count.RDS" %>% 
  readRDS()
# load TCGA-LUAD patient clinical data
clin.LUAD  <- "/Users/benson/Documents/project/RNA-seq1-3/data/clin_LUAD.RDS" %>% 
  readRDS()
# 过滤出东亚人群的数据
asian_id <- clin.LUAD  %>% 
  as.data.frame()%>% 
  filter(race=="asian" & gender=="female")%>% 
  pull(submitter_id)
length(asian_id)
non_asian_id <- clin.LUAD  %>% 
  as.data.frame()%>% 
  filter(race!="asian" & gender=="female") %>% 
  pull(submitter_id)
length(non_asian_id)

# 假设你选择了第一个数字作为键类型
draw_TCGA_boxplot <- function(gene){
  mat <- tcga_count %>% filter(rownames(.)==gene) %>% 
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
    labs(title = paste0(gene," (TCGA RNA-seq)"),
         x = "Groups",
         y = "Expression Level (log(count))") +
    theme_minimal()
  # 在箱线图上添加 t 检验结果
  p + geom_text(aes(x = 1.25, y = 7, 
                    label = paste("p-value:", format(t_test_result$p.value, digits = 3))),
                hjust = 0, vjust = 0, color = "black")
  
}
draw_TCGA_boxplot("SERPINE1")

