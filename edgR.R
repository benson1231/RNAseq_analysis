# library -----------------------------------------------------------------
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

# "ip_L_V_L_CON", "ip_L_V_L_DMS", "ip_L_V_L_AZA", "ip_L_V_L_DAC", "ip_L_V_L_AS",
# "ip_L_V_L_CO", "ip_L_V_L_LCD", "ip_L_V_L_HCD", "ip_L_V_L_BAP", "ip_L_V_L_AS_BAP",
# "ip_L_V_L_CO_BAP", "ip_L_V_L_LCD_BAP", "ip_L_V_L_HCD_BAP",
# pre load ----------------------------------------------------------------
group_names <- c("ip_Y_V_S_CON",
                 "ip_Y_V_S_DMS", "ip_Y_V_S_AZA", "ip_Y_V_S_DAC", "ip_Y_V_S_AS", "ip_Y_V_S_CO",
                 "ip_Y_V_S_LCD", "ip_Y_V_S_HCD", "ip_Y_V_S_BAP", "ip_Y_V_S_AS_BAP", "ip_Y_V_S_CO_BAP",
                 "ip_Y_V_S_LCD_BAP", "ip_Y_V_S_HCD_BAP")

# anno --------------------------------------------------------------------
annotation_df <- "/Users/benson/Documents/project/RNA-seq1-3/anno_gene.RDS" %>% 
  readRDS() %>% as.data.frame() %>% select("ENTREZID","SYMBOL") %>% na.omit()

# load data ---------------------------------------------------------------
### create edgeR object
y <- DGEList(mycount_df)  
group <- group_names
group <- factor(group)
y$samples$group <- group
y$genes <- annotation_df
head(y$samples)
### library size
y$samples$lib.size
# we can also adjust the labelling if we want
barplot(y$samples$lib.size/1e06, names=colnames(y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

### Get log2 counts per million
logcounts <- y %>% as.matrix() %>% log2()
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

### Multidimensional scaling plots
plotMDS(y,xlab = "PC1",ylab = "PC2",pch = 16)

sampleinfo <- myfactors
# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
sampleinfo$sample <- factor(sampleinfo$sample)
sampleinfo$treatment <- factor(sampleinfo$treatment)
sampleinfo$culture <- factor(sampleinfo$culture)
sampleinfo$cell <- factor(sampleinfo$cell)
levels(sampleinfo$sample)
levels(sampleinfo$treatment)
levels(sampleinfo$culture)
levels(sampleinfo$cell)

# 'CON' = '#E0E0E0', 'DMS' = '#ADADAD', 'AZA' = '#FFD2D2', 'DAC' = '#FF9797',
# 'AS' = '#FFFF37', 'CO' = '#FF5151',
# 'LCD' = '#0080FF', 'HCD' = '#005AB5', 'BAP' = '#00DB00',
# 'AS_BAP' = '#FFDC35', 'CO_BAP' = '#EA0000',
# 'LCD_BAP' = '#0000E3', 'HCD_BAP' = '#000079'
# clone = c('WT' = '#FF2D2D', 'L858R' = '#FF9224',
#           'DEL19' = '#BE77FF', 'YAP' = '#6F00D2'))
col.sample <- c("black")[sampleinfo$sample]
col.treatment <- c('#FFFF37','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000',
                   '#E0E0E0','#FF9797','#ADADAD','#005AB5','#000079','#0080FF',
                   '#0000E3')[sampleinfo$treatment]
col.culture <- c("#0072E3","darkblue")[sampleinfo$culture]
col.cell <- c('#FF9224','#6F00D2')[sampleinfo$cell]
# Redo the MDS with cell type colouring
plotMDS(y,col=col.sample,xlab = "PC1",ylab = "PC2")
title("sample")

plotMDS(y,col=col.treatment,xlab = "PC1",ylab = "PC2")
legend("topleft",fill=c('#FFFF37','#FFDC35','#FFD2D2','#00DB00','#FF5151','#EA0000',
                        '#E0E0E0','#FF9797','#ADADAD','#005AB5','#000079','#0080FF',
                        '#0000E3'),
       legend=levels(sampleinfo$treatment))
title("Treatment")

plotMDS(y,col=col.culture,xlab = "PC1",ylab = "PC2")
legend("topleft",fill=c("darkblue","#0072E3"),legend=levels(sampleinfo$culture))
title("Protocol(Short-term/Long-term)")

plotMDS(y,col=col.cell,xlab = "PC1",ylab = "PC2")
legend("topleft",fill=c("orange","#6F00D2"),legend=levels(sampleinfo$cell))
title("Cell type")








# 創建一個 26x26 的矩陣，初始化為全零
matrix_data <- matrix(0, nrow = 26, ncol = 26)
# 將對角線上的元素設置為 1
diag(matrix_data) <- 1
# 檢視結果
colnames(matrix_data) <- group_names

par(mfrow=c(1,1))
v <- voom(y,matrix_data,plot = TRUE,weights = 1)
fit <- lmFit(v)
names(fit)
cont.matrix <- makeContrasts(minus=ip_L_V_L_LCD_BAP-ip_L_V_L_HCD_BAP,levels=matrix_data)
fit.cont <- contrasts.fit(fit, cont.matrix)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=fit.cont[,"ip_L_V_L_LCD_BAP"], values = c(-1, 1), hl.col=c("blue","red"))
