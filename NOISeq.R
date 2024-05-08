# library -----------------------------------------------------------------
library(NOISeq)
library(tidyverse)

# input -------------------------------------------------------------------
### Chromosome table
gene_info <- read_table("gene_info.txt")
mychroms <- gene_info %>% column_to_rownames("ENSEMBL") %>% .[1:3]
### length list
mylength <- gene_info$length
names(mylength) <- gene_info$ENSEMBL
### raw count data(no TMM)
raw_counts_df <- read.csv("/Users/benson/Documents/raw_data/RNA-seq1-3/mycounts_total_f.csv") %>% 
  .[,-1] %>% column_to_rownames("ENSEMBL") 
# mycount_df <- "/Users/benson/Documents/project/RNA-seq1-3/mycount_tmm.RDS" %>% 
#   readRDS()
# df <- raw_counts_df %>% select(.,ip_L_V_L_CON,ip_L_V_L_DMS,ip_L_V_L_AZA,ip_L_V_L_DAC,
#                                ip_L_V_L_AS,ip_L_V_L_CO,ip_L_V_L_LCD,ip_L_V_L_HCD,
#                                ip_L_V_L_BAP, ip_L_V_L_AS_BAP, ip_L_V_L_CO_BAP,
#                                ip_L_V_L_LCD_BAP,ip_L_V_L_HCD_BAP,
#                                ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
#                                ip_Y_V_S_AS,ip_Y_V_S_CO,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
#                                ip_Y_V_S_BAP, ip_Y_V_S_AS_BAP, ip_Y_V_S_CO_BAP,
#                                ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)

# design factor
myfactors <- data.frame(row = names(raw_counts_df),
                        sample = names(raw_counts_df), 
                        treatment = str_sub(names(raw_counts_df), start=10),
                        culture = str_sub(names(raw_counts_df), start=8,end=8),
                        cell    = str_sub(names(raw_counts_df), start=4,end=4)) %>%
  column_to_rownames('row') 

# without TMM
mydata_nt <- readData(data=df,chromosome=mychroms, factors=myfactors,length=mylength)
mycountsbio_nt = dat(mydata_nt, factor = NULL, type = "countsbio")
plot_nt <- explo.plot(mycountsbio_nt, samples = NULL, plottype = "boxplot")

# TMM
myTMM = tmm(assayData(mydata_nt)$exprs, long = 1000, lc = 0) 
mydata_tmm <- readData(data=myTMM,chromosome=mychroms, factors=myfactors,length=mylength)
mycountsbio_tmm = dat(mydata_tmm, factor = NULL, type = "countsbio")
plot_tmm <- explo.plot(mycountsbio_tmm, samples = NULL, plottype = "boxplot")

# save data
mycount_df <- myTMM %>% as.data.frame()
# saveRDS(myTMM,"myTMM.RDS")

# noiseq ------------------------------------------------------------------
mynoiseq.tmm <- noiseq(mydata_tmm, factor = "sample", k = 0.5, norm = "n", pnr = 0.2,
                      nss = 5, v = 0.02, lc = 1, replicates = "no",
                      conditions = c("ip_Y_V_S_CO_BAP","ip_Y_V_S_CO")) ## control要放後面
q_value <- 0.9
mynoiseq.deg = degenes(mynoiseq.tmm, q = q_value, M = NULL)
DE.plot(mynoiseq.tmm, q = q_value, graphic = "MD", xlim=c(-10,10))
DE.plot(mynoiseq.tmm, q = q_value, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseq.tmm, chromosomes = c(1, 2), log.scale = TRUE, join = FALSE,
        q = q_value, graphic = "chrom")
DE.plot(mynoiseq.tmm, chromosomes = NULL, q = q_value, graphic = "distr")

mynoiseq.deg = degenes(mynoiseq.tmm, q = q_value, M = NULL)
mynoiseq.deg1 = degenes(mynoiseq.tmm, q = q_value, M = "up")
mynoiseq.deg2 = degenes(mynoiseq.tmm, q = q_value, M = "down")

ddf <- mynoiseq.tmm@results[[1]]  %>% rownames_to_column("ENSEMBL") %>% left_join(gnno_df)
de <- ddf %>% filter(M>4) %>% pull(SYMBOL)


# output function ---------------------------------------------------------
noiseq_group <- function(seq, conditions,factor=factor,filter=c(0,-0)) {
  mresults <- noiseq(seq, factor = factor, k = 0.5, norm = "n", pnr = 0.2,
                     nss = 5, v = 0.02, lc = 1, replicates = "no",
                     conditions = conditions)
  mnoiseq.deg <- degenes(mresults, q = 0, M = NULL) %>% subset(., M >= filter[1] | M <= filter[2] ) 
  all <- mnoiseq.deg %>% rownames_to_column('ENSEMBL') 
  
  #columns(org.Hs.eg.db) # returns list of available keytypes
  merge_list <- inner_join(anno_df, all, by='ENSEMBL')
  da <- paste0(conditions[1],'_',filter[1])
  colorder_1 = colnames(merge_list)
  sty1 = openxlsx::createStyle(numFmt="0.00")
  wb <- openxlsx::createWorkbook()
  hs <- openxlsx::createStyle(wrapText = T, textDecoration = 'BOLD',fgFill = "#FDE64B",halign='right')
  openxlsx::addWorksheet(wb, da)
  openxlsx::writeData(wb, da, merge_list,   startRow = 1, colNames = T, headerStyle = hs, withFilter=T)
  
  openxlsx::freezePane(wb, da ,firstRow = T, firstCol = T)
  openxlsx::addStyle(wb, da , sty1,  row= 1:nrow(merge_list)+1, cols = 1:length(colorder_1)+1, gridExpand = TRUE)
  #openxlsx::conditionalFormatting(wb, i, row= 1:nrow(da)+1, cols = 1:length(colorder_1)+1, 
  #                                style = c("#5b71f1", "#F7F34C"),  type = "colourScale")
  openxlsx::saveWorkbook(wb, paste0('./',da,"_","deg",".xlsx", sep = ""), overwrite  = T)
  message("generating excel Done")
  return(merge_list)
}

# output ------------------------------------------------------------------
# load count df with TMM
myTMM <- "/Users/benson/Documents/project/RNA-seq1-3/myTMM.RDS" %>% 
  readRDS() %>% as.data.frame %>% 
  select(.,ip_L_V_L_CON,ip_L_V_L_DMS,ip_L_V_L_AZA,ip_L_V_L_DAC,
         ip_L_V_L_AS,ip_L_V_L_CO,ip_L_V_L_LCD,ip_L_V_L_HCD,
         ip_L_V_L_BAP, ip_L_V_L_AS_BAP, ip_L_V_L_CO_BAP,
         ip_L_V_L_LCD_BAP,ip_L_V_L_HCD_BAP,
         ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
         ip_Y_V_S_AS,ip_Y_V_S_CO,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
         ip_Y_V_S_BAP, ip_Y_V_S_AS_BAP, ip_Y_V_S_CO_BAP,
         ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
myfactors <- data.frame(row = names(myTMM),
                        sample = names(myTMM), 
                        treatment = str_sub(names(myTMM), start=10),
                        culture = str_sub(names(myTMM), start=8,end=8),
                        cell    = str_sub(names(myTMM), start=4,end=4)) %>%
  column_to_rownames('row') 
anno_df <- "/Users/benson/Documents/project/RNA-seq1-3/anno_gene.RDS" %>% 
  readRDS()
mydata_tmm <- readData(data=myTMM,chromosome=mychroms, factors=myfactors,length=mylength)
mycountsbio_tmm = dat(mydata_tmm, factor = NULL, type = "countsbio")
plot_tmm <- explo.plot(mycountsbio_tmm, samples = NULL, plottype = "boxplot")

# VTN ---------------------------------------------------------------------
setwd("./Vitro/")
Y_AS <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_AS','ip_Y_V_S_CON'),factor='sample', filter=c(0,-0))
Y_AZA <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_AZA','ip_Y_V_S_CON'),factor='sample', filter=c(0,-0))
Y_DAC <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_DAC','ip_Y_V_S_CON'),factor='sample', filter=c(0,-0))
Y_CO <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_CO','ip_Y_V_S_CON'),factor='sample', filter=c(0,-0))
Y_HCD <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_HCD','ip_Y_V_S_CON'),factor='sample', filter=c(0,-0))
Y_LCD <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_LCD','ip_Y_V_S_CON'),factor='sample', filter=c(0,-0))
Y_B<- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_BAP','ip_Y_V_S_DMS'),factor='sample', filter=c(0,-0))

Y_AS_B <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_AS_BAP','ip_Y_V_S_DMS'),factor='sample', filter=c(0,-0))
Y_CO_B <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_CO_BAP','ip_Y_V_S_DMS'),factor='sample', filter=c(0,-0))
Y_HCD_B <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_LCD_BAP','ip_Y_V_S_DMS'),factor='sample', filter=c(0,-0))
Y_LCD_B <- noiseq_group(seq=mydata_tmm, conditions = c('ip_Y_V_S_HCD_BAP','ip_Y_V_S_DMS'),factor='sample', filter=c(0,-0))


L_AS <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_AS','ip_L_V_L_CON'),factor='sample', filter=c(0,-0))
L_AZA <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_AZA','ip_L_V_L_CON'),factor='sample', filter=c(0,-0))
L_DAC <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_DAC','ip_L_V_L_CON'),factor='sample', filter=c(0,-0))
L_CO <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_CO','ip_L_V_L_CON'),factor='sample', filter=c(0,-0))
L_HCD <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_HCD','ip_L_V_L_CON'),factor='sample', filter=c(0,-0))
L_LCD <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_LCD','ip_L_V_L_CON'),factor='sample', filter=c(0,-0))
L_B<- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_BAP','ip_L_V_L_DMS'),factor='sample', filter=c(0,-0))

L_AS_B <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_AS_BAP','ip_L_V_L_DMS'),factor='sample', filter=c(0,-0))
L_CO_B <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_CO_BAP','ip_L_V_L_DMS'),factor='sample', filter=c(0,-0))
L_HCD_B <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_LCD_BAP','ip_L_V_L_DMS'),factor='sample', filter=c(0,-0))
L_LCD_B <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_V_L_HCD_BAP','ip_L_V_L_DMS'),factor='sample', filter=c(0,-0))


# matrigel ----------------------------------------------------------------
setwd("./Matrigel/")
#W_1np <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_1NP','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#W_1np_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_1NP_NPY','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#W_18d <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_18D','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#W_18d_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_18D_NPY','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
W_BAP <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_BAP','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
W_BAP_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_BAP_NPY','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
W_DBA <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_DBA','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
W_DBA_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_DBA_NPY','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))

#W_18D_1NP <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_18D_1NP','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#W_MX <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_MX','ip_W_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#W_DMH <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_DMH','ip_W_M_L_Con_S'),factor='sample', filter=c(0.0,-0.0))
W_NPY <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_NPY','ip_W_M_L_Con_S'),factor='sample', filter=c(0.0,-0.0))
#W_DAC <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_DAC','ip_W_M_L_Con'),factor='sample', filter=c(0.0,-0.0))
#W_AZA <- noiseq_group(seq=mydata_tmm, conditions = c('ip_W_M_L_AZA','ip_W_M_L_Con'),factor='sample', filter=c(0.0,-0.0))

##################
#D_1np <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_1NP','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#D_1np_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_1NP_NPY','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#D_18d <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_18D','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#D_18d_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_18D_NPY','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
D_BAP <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_BAP','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
D_BAP_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_BAP_NPY','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
D_DBA <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_DBA','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
D_DBA_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_DBA_NPY','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))

#D_18D_1NP <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_18D_1NP','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#D_MX <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_MX','ip_D_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#D_DMH <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_DMH','ip_D_M_L_Con_S'),factor='sample', filter=c(0.0,-0.0))
D_NPY <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_NPY','ip_D_M_L_Con_S'),factor='sample', filter=c(0.0,-0.0))
#D_DAC <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_DAC','ip_D_M_L_Con'),factor='sample', filter=c(0.0,-0.0))
#D_AZA <- noiseq_group(seq=mydata_tmm, conditions = c('ip_D_M_L_AZA','ip_D_M_L_Con'),factor='sample', filter=c(0.0,-0.0))

##################
#L_1np <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_1NP','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#L_1np_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_1NP_NPY','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#L_18d <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_18D','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#L_18d_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_18D_NPY','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
L_BAP <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_BAP','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
L_BAP_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_BAP_NPY','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
L_DBA <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_DBA','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
L_DBA_npy <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_DBA_NPY','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))

#L_18D_1NP <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_18D_1NP','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#L_MX <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_MX','ip_L_M_L_DMS_S'),factor='sample', filter=c(0.0,-0.0))
#L_DMH <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_DMH','ip_L_M_L_Con_S'),factor='sample', filter=c(0.0,-0.0))
L_NPY <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_NPY','ip_L_M_L_Con_S'),factor='sample', filter=c(0.0,-0.0))
#L_DAC <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_DAC','ip_L_M_L_Con'),factor='sample', filter=c(0.0,-0.0))
#L_AZA <- noiseq_group(seq=mydata_tmm, conditions = c('ip_L_M_L_AZA','ip_L_M_L_Con'),factor='sample', filter=c(0.0,-0.0))
