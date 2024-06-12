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
name_df <- openxlsx::read.xlsx("/Users/benson/Documents/project/RNA-seq1-3/data/name.xlsx")
raw_count <- read.csv("/Users/benson/Documents/raw_data/RNA-seq1-3/df_total.txt") %>% 
  column_to_rownames("X")
length <- raw_count[,1:2]
raw_count_df <- raw_count %>% .[,-c(1,2)]
abb_raw_count <- raw_count_df %>% 
  dplyr::select(name_df$group_name) %>% setNames(name_df$abbreviate) %>%
  as.matrix()

### filter CPM --------------------------------------------------------------
### Normalization
raw_obj <- edgeR::DGEList(abb_raw_count)
raw_obj <- edgeR::calcNormFactors(raw_obj, method = "TMM")
raw_CPM <- edgeR::cpm(raw_obj, log = F, normalized.lib.sizes = T)
# Have a look at the output
head(raw_obj$samples)
head(raw_CPM)
# Which values in raw_CPM are greater than 0.5?
thresh <- raw_CPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of abb_raw_count to keep the more highly expressed genes
counts.keep <- abb_raw_count[keep,]
summary(keep)
dim(counts.keep)

# take a look CPM vs raw_counts
plot(raw_CPM[,3],abb_raw_count[,3],ylim=c(0,50),xlim=c(0,3),xlab = "CPM",ylab="raw_counts")
# Add a vertical line at 0.5 CPM
abline(v=0.5)

### density plot before CPM-filtered
raw_logCPM <- edgeR::cpm(raw_obj, log = T)
col_sample <- c(RColorBrewer::brewer.pal(n = 8, "Set1"), RColorBrewer::brewer.pal(n = 8, "Set2"))
plot(density(raw_logCPM[,1]), col = col_sample[1], lwd = 2, las =1,
     xlab = "logCPM", main = "Before filtering")
for (i in seq(2, ncol(raw_logCPM))){
  lines(density(raw_logCPM[,i]), col = col_sample[i], lwd = 2, las =2)
}

### density plot after CPM-filtered
filtered_Obj <- edgeR::DGEList(counts.keep)
filtered_Obj<- edgeR::calcNormFactors(filtered_Obj, method = "TMM")
head(filtered_Obj$samples)
filtered_logCPM <- edgeR::cpm(filtered_Obj, log = T)
col_sample <- c(RColorBrewer::brewer.pal(n = 8, "Set1"), RColorBrewer::brewer.pal(n = 8, "Set2"))
plot(density(filtered_logCPM[,1]), col = col_sample[1], lwd = 2, las =1,
     xlab = "logCPM", main = "Before filtering")
for (i in seq(2, ncol(filtered_logCPM))){
  lines(density(filtered_logCPM[,i]), col = col_sample[i], lwd = 2, las =2)
}


### library size -----------------------------------------------------------
barplot(raw_obj$samples$lib.size, names=name_df$abbreviate, las=2, 
        ylim = c(0, 80000000), cex.names = 0.8)
title("Barplot of library sizes of 52 samples")

### get log2CPM
head(filtered_Obj$samples)
mylogCPM <- edgeR::cpm(filtered_Obj,log=TRUE,normalized.lib.sizes = T)
dim(mylogCPM)
head(mylogCPM)
# saveRDS(mylogCPM,"mylogCPM.RDS")

### get CPM
myCPM <- edgeR::cpm(filtered_Obj,log=FALSE,normalized.lib.sizes = T)
dim(myCPM)
head(myCPM)
# saveRDS(myCPM,"myCPM.RDS")

### boxplot of logCPM
boxplot(mylogCPM, xlab="", ylab="Log2 counts per million",las=2,
        main="Boxplots of logCPM (normalized)")
abline(h=median(mylogCPM),col="blue")


# heatmap -----------------------------------------------------------------
mylogCPM <- mylogCPM %>% as.data.frame()
abbr_sampleinfo <- data.frame(row = names(mylogCPM),
                              sample = names(mylogCPM), 
                              treatment = str_sub(names(mylogCPM), start=3),
                              cell = str_sub(names(mylogCPM), start=1,end=1)) %>%
  column_to_rownames('row') 
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(mylogCPM, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- mylogCPM[select_var,] %>% as.matrix()
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
## Get some nicer colours
mypalette <- RColorBrewer::brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[abbr_sampleinfo$sample]

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
         annCol=abbr_sampleinfo[, c(2,3)], labCol=abbr_sampleinfo$sample, scale="row")
dev.off()
