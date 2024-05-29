### library -----------------------------------------------------------------
library(RTN)
library(TCGAbiolinks)
library(maftools)
library(survival)
library(DT)

# http://www.signalingpathways.org/ominer/query.jsf
# https://www.bioconductor.org/packages/release/bioc/vignettes/RTN/inst/doc/RTN.html

### OS analysis -------------------------------------------------------------
# # 查詢肺腺癌("TCGA-LUAD")
# query <- GDCquery(
#   project = "TCGA-LUAD",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
# )
# # 只需下載一次
# GDCdownload(query)
# # get data
# maf <- GDCprepare(query)
# maf_object <- read.maf(maf = maf)
# saveRDS(maf_object,"maf_object.RDS")
maf_object <- "/Users/benson/Documents/project/RNA-seq1-3/data/maf_object.RDS" %>%
  readRDS()
# # get clinical data
# clin.LUAD <- GDCquery_clinic("TCGA-LUAD", "clinical")
# saveRDS(clin.LUAD,"clin_LUAD.RDS")
clin.LUAD <- "/Users/benson/Documents/project/RNA-seq1-3/data/clin_LUAD.RDS" %>% 
  readRDS()
# set gene you interest
gene_of_interest <- "EGFR"
# get mutation sample ID
mutation_samples <- subset(maf_object@data, Hugo_Symbol == gene_of_interest)$Tumor_Sample_Barcode
# 提取了樣本ID的前12個字符，用於後續臨床樣本配對
mutation_samples_short <- substr(mutation_samples, 1, 12)
length(mutation_samples_short)
head(mutation_samples_short)
# 新增紀錄突變狀態的欄位
clin.LUAD$mutation_status <- ifelse(clin.LUAD$submitter_id %in% mutation_samples_short, "Mutated", "Wildtype")
# plotting
TCGAanalyze_survival(
  data = clin.LUAD,
  clusterCol = "mutation_status",
  main = "TCGA Set\n LUAD",
  height = 10,
  width=10,
  legend = gene_of_interest, 
  filename = paste0("output/",gene_of_interest,"_survival.pdf")
)

### TNI object -----------------------------------------------------------
# Input 1: 'expData', a named gene expression matrix (genes on rows, samples on cols); 
# Input 2: 'regulatoryElements', a vector listing genes regarded as TFs
# Input 3: 'rowAnnotation', an optional data frame with gene annotation
# Input 4: 'colAnnotation', an optional data frame with sample annotation
tfs <- c("JUN","FOS","SP1","HSF2","HIF1A","AP1","ESR1")
abbr_mat <- abbr_count %>% as.matrix()
colAnnotation <- openxlsx::read.xlsx("/Users/benson/Documents/project/RNA-seq1-3/data/sampleinfo.xlsx",rowNames = T)
rtni <- tni.constructor(expData = abbr_mat, 
                        regulatoryElements = tfs, 
                        colAnnotation = colAnnotation, 
                        rowAnnotation = gene_df)
# p.s. alternatively, 'expData' can be a 'SummarizedExperiment' object
# Please set nPermutations >= 1000
rtni <- tni.permutation(rtni, nPermutations = 100)
rtni <- tni.bootstrap(rtni)
rtni <- tni.dpi.filter(rtni)
tni.regulon.summary(rtni)
tni.regulon.summary(rtni, regulatoryElements = "JUN")
regulons <- tni.get(rtni, what = "regulons.and.mode", idkey = "SYMBOL")
head(regulons$JUN)

# TNA object and GSEA -----------------------------------------------------
# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype
FC_list <- get_list(name_df$file_name[13],type = "ENSEMBL")
head(FC_list)
DEG_list <- get_deg(name_df$file_name[13],type = "ENSEMBL")
head(DEG_list)
rtna <- tni2tna.preprocess(object = rtni, 
                           phenotype = FC_list, 
                           hits = DEG_list, 
                           phenoIDs = gene_df)
# Run the MRA method
rtna <- tna.mra(rtna)
# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra <- tna.get(rtna, what="mra", ntop = -1)
head(mra)
# Run the GSEA method
# Please set nPermutations >= 1000
rtna <- tna.gsea1(rtna, nPermutations=100)
# Get GSEA results
gsea1 <- tna.get(rtna, what="gsea1", ntop = -1)
head(gsea1)
# Plot GSEA results
tna.plot.gsea1(rtna, labPheno="abs(log2 fold changes)", ntop = -1, filepath = "./output")
# Run the GSEA-2T method
# Please set nPermutations >= 1000
rtna <- tna.gsea2(rtna, nPermutations = 100)
# Get GSEA-2T results
gsea2 <- tna.get(rtna, what = "gsea2", ntop = -1)
head(gsea2$differential)
# Plot GSEA-2T results
tna.plot.gsea2(rtna, labPheno="log2 fold changes", tfs="JUN",filepath = "./output")


### regulon heatmap ---------------------------------------------------------------------
# A list of transcription factors of interest (here 36 risk-associated TFs)
tfs <- c("JUN","FOS","SP1","HSF2","HIF1A","AP1","ESR1")
# Compute regulon activity for individual samples
rtni1st <- tni.gsea2(rtni, regulatoryElements = tfs)
metabric_regact <- tni.get(rtni1st, what = "regulonActivity")
# Get sample attributes from the 'rtni1st' dataset
metabric_annot <- tni.get(rtni1st, "colAnnotation")
# Get ER+/- and PAM50 attributes for pheatmap
treat <- c("AZA","DAC","AS","CO","LCD","HCD","BAP","AS_BAP","CO_BAP",
           "LCD_BAP","HCD_BAP")
clone <- c("W","L","D","Y")
t_annot <- metabric_annot[,treat]
c_annot <- metabric_annot[,clone]
# Plot regulon activity profiles
pheatmap(t(metabric_regact$dif), 
         main="Short-term",
         annotation_col = t_annot, 
         show_colnames = T, annotation_legend = FALSE, 
         clustering_method = "ward.D2", fontsize_row = 6,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")
pheatmap(t(metabric_regact$dif), 
         main="Short-term",
         annotation_col = c_annot, 
         show_colnames = T, annotation_legend = FALSE, 
         clustering_method = "ward.D2", fontsize_row = 6,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")


### SPP ---------------------------------------------------------------------
### download data from SPP
spp_TF <- readxl::read_excel("/Users/benson/Downloads/SRX_results_file.xlsx") %>% 
  tidyr::separate_rows(`IPAGS|BSM|Other AGS`, sep = "\\|") %>% 
  dplyr::distinct(`IPAGS|BSM|Other AGS`) %>%
  pull(`IPAGS|BSM|Other AGS`)
head(spp_TF)
deg <- get_deg(name_df$file_name[8],dir = "all",log_crit = 1)
head(deg)
# select SPP TFs in DEG
de_spp_TF <- intersect(spp_TF, deg)
draw_from_list(de_spp_TF,groups = "CD")

### load CollecTRI network TFs
de_all_TF <- intersect(net_TF$source, deg)
length(de_all_TF)
draw_from_list(de_all_TF,groups = "CD")

### TCGA regulons -----------------------------------------------------------
mutation_TCGA <- maf_object@data %>% select(Hugo_Symbol, Tumor_Sample_Barcode) %>% 
  mutate("patient" = str_sub(.$Tumor_Sample_Barcode,0,12))
head(mutation_TCGA)
# get number of patients 
patient_num <- mutation_TCGA %>%
  dplyr::distinct(patient) %>% 
  pull(patient) %>% 
  length()
# count mutation rate in different patients
muta_rate <- mutation_TCGA %>% 
  distinct(patient,Hugo_Symbol) %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% 
  mutate(mutation_ratio=count/patient_num) %>% 
  setNames(c("ID","count","mutation_ratio"))
muta_rate$ID <- factor(muta_rate$ID, levels = muta_rate$ID)
head(muta_rate)
# saveRDS(muta_rate,"muta_rate.RDS")
muta_rate <- "/Users/benson/Documents/project/RNA-seq1-3/data/muta_rate.RDS" %>% 
  readRDS()
# filter genes we interesting 
# all of TFs
all_tf <- muta_rate %>% filter(ID %in% net_TF$source)
# TFs in our study
tf <- "JUN"
tf_df <- plot_TF(name_df$file_name[13], tf, title = abb(file_name),plot = F)
# ggsave("TF.jpeg")
# top50 heatmap
regulons_list <- tf_df %>% arrange(desc(abs(logFC))) %>% head(100) %>% pull(ID)
draw_from_list(list = regulons_list, groups = "CD",
               id = "SYMBOL",label_num = F,anno = T,title = "")
regulons <- muta_rate %>% filter(ID %in% regulons_list)

# 绘制 mutation_ratio 的柱状图
ggplot(regulons[1:10,], aes(x = ID, y = `mutation_ratio(%)`)) +
  geom_bar(stat = "identity", fill = "salmon") +
  theme_minimal() +
  labs(title = paste("regulons of",tf), x = "Hugo Symbol", y = "Mutation Ratio")
# 绘制 count 的柱状图
ggplot(regulons[1:10,], aes(x = ID, y = `mutation_ratio(%)`)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = paste("regulons of",tf), x = "Hugo Symbol", y = "Count")

# plot logFC_TCGA mutation plot
tf_df <- plot_TF(file_name, "JUN",plot = F)
regulons_list <- tf_df %>% arrange(desc(abs(logFC))) %>% pull(ID)
length(regulons_list)
plot_mutaRate(name_df$file_name[8],title = tf,list=regulons_list,with_line = F)

