# library -----------------------------------------------------------------
library(NOISeq)
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(tidyverse)

# load raw count ----------------------------------------------------------
raw_count <- read.csv("/Users/benson/Documents/raw_data/RNA-seq1-3/df_total.txt") %>% 
  column_to_rownames("X")
length <- raw_count[,1:2]
raw_count_df <- raw_count %>% .[,-c(1,2)]

# 
abb_raw_count <- raw_count_df

# design factor
sampleinfo <- data.frame(row = names(abbr_count),
                         sample = names(abbr_count), 
                         treatment = str_sub(names(abbr_count), start=5),
                         culture = str_sub(names(abbr_count), start=3,end=3),
                         cell    = str_sub(names(abbr_count), start=1,end=1)) %>%
  column_to_rownames('row') 

# Obtain CPMs
countdata <- abb_raw_count %>% as.matrix()
myCPM <- edgeR::cpm(countdata)
# Have a look at the output
head(myCPM)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 77 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)

plot(myCPM[,1],countdata[,1])
plot(myCPM[,3],countdata[,3],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)



# normalization -----------------------------------------------------------
### Unnormalized
unnor_Obj <- edgeR::DGEList(counts.keep)
# have a look at dgeObj
unnor_Obj
# library size
barplot(unnor_Obj$samples$lib.size,names=colnames(y),las=2)
title("Barplot of library sizes")
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(unnor_Obj$samples$lib.size, names=colnames(dgeObj), las=2, 
        ylim = c(0, 80000000), cex.names = 0.7)
title("Barplot of library sizes")

# Get log2 counts per million
unnor_logcounts <- edgeR::cpm(unnor_Obj,log=TRUE)
unnor_counts <- edgeR::cpm(unnor_Obj,log=FALSE)
# Check distributions of samples using boxplots
boxplot(unnor_logcounts, xlab="", ylab="Log2 counts per million",las=2,
        cex.axis=0.7,main="Boxplots of logCPMs (unnormalized)")
abline(h=median(logcounts),col="blue")
# MDS plot
plotMDS(dgeObj)


### normalized
nor_Obj<- edgeR::calcNormFactors(unnor_Obj, method = "TMM")
nor_logcounts <- edgeR::cpm(nor_Obj,log=TRUE,normalized.lib.sizes = T)
nor_counts <- edgeR::cpm(nor_Obj,log=FALSE,normalized.lib.sizes = T)
boxplot(nor_logcounts, xlab="", ylab="Log2 counts per million",las=2,
        main="Boxplots of logCPMs (normalized)")

# before and after
par(mfrow=c(1,2))
boxplot(unnor_logcounts, xlab="", ylab="Log2 counts per million",las=2,
        cex.axis=0.7,main="Boxplots of logCPMs (unnormalized)")
abline(h=median(logcounts),col="blue")
boxplot(nor_logcounts, xlab="", ylab="Log2 counts per million",las=2,
        main="Boxplots of logCPMs (normalized)")
abline(h=median(logcounts),col="blue")
# MD plot
plotMD(unnor_logcounts,column = 11)
abline(h=0,col="grey")
plotMD(nor_logcounts,column =11)
abline(h=0,col="grey")

# get length
counts <- nor_counts %>% as.data.frame() %>% rownames_to_column("GeneID") %>% 
  left_join(length,"GeneID") %>% column_to_rownames("GeneID")
mylength <- counts$length
names(mylength) <- counts$GeneID
my_counts <- counts[,-53]
myTMM <- my_counts %>% as.data.frame()
# saveRDS(myTMM,"myTMM.RDS")
# saveRDS(nor_logcounts,"mylogcount.RDS")

#### 原本的分布
logCPM_table <- log2(edgeR::cpm(unnor_counts) + 1)
col_sample <- c(RColorBrewer::brewer.pal(n = 8, "Set1"), RColorBrewer::brewer.pal(n = 8, "Set2"))
plot(density(logCPM_table[,1]), col = col_sample[1], lwd = 2, las =1,
     xlab = "logCPM", main = "Before filtering")
for (i in seq(2, ncol(logCPM_table))){
  lines(density(logCPM_table[,i]), col = col_sample[i], lwd = 2, las =2)
}


# Remove the low expression genes
### 方法1. 考慮平均表現量
keep.genes <- apply(nor_Obj$counts, 1, mean) > 1
table(keep.genes)
### 方法2. count == 0 的比例
keep.genes2 <- apply(nor_Obj$counts == 0, 1, mean) < 0.3
table(keep.genes2)
### Create a DGEList with keepping genes only
deg_obj_filtered <- nor_Obj[keep.genes,, keep.lib.size = FALSE]
deg_obj_filtered$counts

#### 篩選後的分布
logCPM_table_filtered <- edgeR::cpm(deg_obj_filtered, log = TRUE)
col_sample <- c(RColorBrewer::brewer.pal(n = 8, "Set1"), RColorBrewer::brewer.pal(n = 8, "Set2"))
plot(density(logCPM_table_filtered[,1]), col = col_sample[1], lwd = 2, las =1,
     xlab = "logCPM", main = "Before filtering")
for (i in seq(2, ncol(logCPM_table_filtered))){
  lines(density(logCPM_table_filtered[,i]), col = col_sample[i], lwd = 2, las =2)
}


# heatmap -----------------------------------------------------------------
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(nor_logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcount[select_var,] %>% as.matrix()
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
## Get some nicer colours
mypalette <- RColorBrewer::brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$sample]

# Plot the heatmap
png("heatmap1.png", width = 1500, height = 1000)
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none",
          main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()


# heatmap with annotation
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
png("heatmap2.png", width = 1500, height = 1000)
aheatmap(highly_variable_lcpm, col=rev(morecols(50)), main="Top 500 most variable genes across samples",
         annCol=sampleinfo[, c(2,4)], labCol=sampleinfo$sample, scale="row")
dev.off()

