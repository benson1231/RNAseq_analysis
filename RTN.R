library(RTN)

# Input 1: 'expData', a named gene expression matrix (genes on rows, samples on cols); 
# Input 2: 'regulatoryElements', a vector listing genes regarded as TFs
# Input 3: 'rowAnnotation', an optional data frame with gene annotation
# Input 4: 'colAnnotation', an optional data frame with sample annotation
tfs <- c("JUN","FOS","SP1","HSF2","HIF1A","AP1","ESR1")
abbr_mat <- abbr_count %>% as.matrix()
colAnnotation <- openxlsx::read.xlsx("/Users/benson/Documents/project/RNA-seq1-3/data/sampleinfo.xlsx",rowNames = T)
rtni <- tni.constructor(expData = abbr_mat, 
                        regulatoryElements = tfs, colAnnotation = colAnnotation, 
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


# new ---------------------------------------------------------------------
library(RTN)
library(Fletcher2013b)
library(pheatmap)
# Load 'rtni1st' data object, which includes regulons and expression profiles
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



