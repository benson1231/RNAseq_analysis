library(NOISeq)
# input data
df <- raw_counts_df[,-1] %>% column_to_rownames("ENSEMBL") %>% 
  select(.,ip_L_V_L_CON,ip_L_V_L_DMS,ip_L_V_L_AZA,ip_L_V_L_DAC,
         ip_L_V_L_AS,ip_L_V_L_CO,ip_L_V_L_LCD,ip_L_V_L_HCD,
         ip_L_V_L_BAP, ip_L_V_L_AS_BAP, ip_L_V_L_CO_BAP,
         ip_L_V_L_LCD_BAP,ip_L_V_L_HCD_BAP,
         ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
         ip_Y_V_S_AS,ip_Y_V_S_CO,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
         ip_Y_V_S_BAP, ip_Y_V_S_AS_BAP, ip_Y_V_S_CO_BAP,
         ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP) 
myfactors = data.frame(group=substr(colnames(df), 10,20),
                       clone=substr(colnames(df), 4,4),
                       group_name=colnames(df))
# without TMM
mydata_nt <- readData(data=df,chromosome=mychroms, factors=myfactors)
mycountsbio_nt = dat(mydata_nt, factor = NULL, type = "countsbio")
plot_nt <- explo.plot(mycountsbio_nt, samples = NULL, plottype = "boxplot")

# TMM
myTMM = tmm(assayData(mydata_nt)$exprs, long = 1000, lc = 0) 
mydata_tmm <- readData(data=myTMM,chromosome=mychroms, factors=myfactors)
mycountsbio_tmm = dat(mydata_tmm, factor = NULL, type = "countsbio")
plot_tmm <- explo.plot(mycountsbio_tmm, samples = NULL, plottype = "boxplot")

# save data
mycount_df <- myTMM %>% as.data.frame()
# saveRDS(myTMM,"myTMM.RDS")


# noiseq ------------------------------------------------------------------
mynoiseq.tmm = noiseq(mydata_nt, k = 0.5, norm = "tmm", factor="group_name", 
                      conditions = c("ip_Y_V_S_CO","ip_Y_V_S_CO_BAP"), 
                      lc = 0, replicates = "technical")

mynoiseq.tmm <- noiseq(mydata_tmm, factor = "group_name", k = NULL, norm = "n", pnr = 0.2,
                      nss = 5, v = 0.02, lc = 1, replicates = "no",
                      conditions = c("ip_Y_V_S_CO","ip_Y_V_S_CO_BAP"))
mynoiseq.deg = degenes(mynoiseq.tmm, q = 0.8, M = NULL)
DE.plot(mynoiseq.tmm, q = 0.7, graphic = "MD", xlim=c(-10,10))
DE.plot(mynoiseq.tmm, q = 0.7, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseq.tmm, chromosomes = c(1, 2), log.scale = TRUE, join = FALSE,
        q = 0.8, graphic = "chrom")
DE.plot(mynoiseq.tmm, chromosomes = NULL, q = 0.8, graphic = "distr")

