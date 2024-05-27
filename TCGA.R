# library -----------------------------------------------------------------
library(TCGAbiolinks)
library(maftools)
library(survival)
library(DT)

# query details -------------------------------------------------------------------
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
# GDCdownload(query)
maf <- GDCprepare(query)

query_clin <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML"
)
# GDCdownload(query_clin)
clinical <- GDCprepare_clinic(query_clin, clinical.info = "patient")

# 将突变数据转化为maftools对象
maf_object <- read.maf(maf = maf)
datatable(getSampleSummary(maf_object),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
oncoplot(maf_object, top = 10, removeNonMutated = TRUE)
titv = titv(maf_object, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)

# 从maftools对象中提取样本突变信息
mutation_data <- getSampleSummary(maf_object)

# 选择感兴趣的基因，例如TP53
gene_of_interest <- "KRAS"
mutation_samples <- subset(maf_object@data, Hugo_Symbol == gene_of_interest)$Tumor_Sample_Barcode
# 提取样本ID的前12个字符
mutation_samples_short <- substr(mutation_samples, 1, 12)
length(mutation_samples_short)
head(mutation_samples_short)

# 提取样本ID
clinical$sample <- substr(clinical$bcr_patient_barcode, 1, 12)

# 标记临床数据中的突变状态
clinical$mutation_status <- ifelse(clinical$sample %in% mutation_samples_short, "Mutated", "Wildtype")
# 检查突变状态列的水平
table(clinical$mutation_status)

# 创建生存对象
clinical$OS.time <- as.numeric(clinical$days_to_death)
clinical$OS.event <- ifelse(clinical$vital_status == "Dead", 1, 0)

# 创建生存分析对象
surv_object <- Surv(time = clinical$OS.time, event = clinical$OS.event)

# 确保mutation_status列为因子类型并有两个水平
clinical$mutation_status <- as.factor(clinical$mutation_status)
if (length(levels(clinical$mutation_status)) > 1) {
  # 创建生存分析模型
  fit <- survfit(surv_object ~ mutation_status, data = clinical)
  
  # 绘制生存曲线
  plot(fit, col = c("blue", "red"), lty = 1:2, xlab = "Days", ylab = "Overall Survival Probability")
  legend("topright", legend = c("Wildtype", "Mutated"), col = c("blue", "red"), lty = 1:2)
  title(paste0("Survival Curve by Mutation Status of ",gene_of_interest))
  # 使用Cox回归模型评估突变状态对生存率的影响
  cox_fit <- coxph(surv_object ~ mutation_status, data = clinical)
  summary(cox_fit)
} else {
  cat("mutation_status 列只有一个水平，无法进行Cox回归分析。\n")
}


# new ---------------------------------------------------------------------
# 查詢肺腺癌("TCGA-LUAD")
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
# 只需下載一次
GDCdownload(query)
# get data
maf <- GDCprepare(query)
maf_object <- read.maf(maf = maf)
# set gene you interest
gene_of_interest <- "TP53"
# get mutation sample ID
mutation_samples <- subset(maf_object@data, Hugo_Symbol == gene_of_interest)$Tumor_Sample_Barcode
# 提取了樣本ID的前12個字符，用於後續臨床樣本配對
mutation_samples_short <- substr(mutation_samples, 1, 12)
length(mutation_samples_short)
head(mutation_samples_short)
# get clinical data
clin.LUAD <- GDCquery_clinic("TCGA-LUAD", "clinical")
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

