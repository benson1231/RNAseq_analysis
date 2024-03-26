# library -----------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(ggVennDiagram)

# load raw reads counts ---------------------------------------------------------------
raw_counts_df <- read.csv("/Users/benson/Documents/raw_data/RNA-seq1/mycounts_f.txt")

mycount_df <- raw_counts_df %>% column_to_rownames(var = "ENSEMBL") %>% 
  select(.,ip_L_V_L_CON,ip_L_V_L_DMS,ip_L_V_L_AZA,ip_L_V_L_DAC,
         ip_L_V_L_AS,ip_L_V_L_CO,ip_L_V_L_LCD,ip_L_V_L_HCD,
         ip_L_V_L_BAP, ip_L_V_L_AS_BAP, ip_L_V_L_CO_BAP,
         ip_L_V_L_LCD_BAP,ip_L_V_L_HCD_BAP,
         ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
         ip_Y_V_S_AS,ip_Y_V_S_CO,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
         ip_Y_V_S_BAP, ip_Y_V_S_AS_BAP, ip_Y_V_S_CO_BAP,
         ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP) 

# df <- raw_counts_df[,grepl("ip_L_V_L|ip_Y_V_S", colnames(mycount_df))]

# 抓出含有ensembl/symbol df
# data_path <- "/Users/benson/Documents/raw_data/RNA-seq1/old/raw/Y_L_V"
# gene_df <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_AS_0_deg.xlsx")) %>% 
#   select(ENSEMBL,SYMBOL) 

# annotation data frame
# saveRDS(anno_v45_g_col,"gene_df.RDS")
gene_df <- "/Users/benson/Documents/project/RNA-seq1/gene_df.RDS" %>% readRDS() %>%
  select(ENSEMBL,SYMBOL)

# heatmap with annotation ---------------------------------------------------------
mat_scale <- mycount_df %>% t() %>% scale() %>% t()
col <- colnames(mycount_df)

# agent name
agent <- factor (
  str_replace_all(col, c("ip_L_V_L_|ip_Y_V_S_"='')),
  levels=c('CON','DMS',"AZA","DAC",'AS',"CO","LCD","HCD","BAP",
           "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
# clone name
clone <- factor( 
  str_replace_all(str_sub(col, 4,4), 
                  c("W" = "WT", "L" = "L858R", "D" = "DEL19", "Y" = "YAP")),
  levels=c('WT','L858R',"DEL19","YAP"))
# draw hp
ha <- HeatmapAnnotation(agent = agent, clone = clone,
                        col = list(agent=c('CON'='#E0E0E0', 'DMS'='#ADADAD', 'AZA'='#FFD2D2', 'DAC'='#FF9797',
                                          'AS'='#FFFF37','CO'='#FF5151','LCD'='#0080FF',
                                          'HCD'='#005AB5', 'BAP'='#00DB00', 'AS_BAP'='#FFDC35', 'CO_BAP'= '#EA0000',
                                          'LCD_BAP'= '#0000E3', "HCD_BAP"="#000079"),
                                  clone=c('WT'="#FF2D2D",'L858R'="#FF9224",
                                          "DEL19"= "#66B3FF","YAP"="#2828FF")))

ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                        show_row_names = F, show_column_names = F,
                        name = "Z-score"
                        )

# DEG heatmap -----------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1/metal_anno/"
co <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_BAP_0_deg.xlsx"))
co <- co %>% filter(abs(M)>1)
list <- co$ENSEMBL
gene_df <- co[,c("ENSEMBL","SYMBOL")] 
mycount_df$ENSEMBL <- rownames(mycount_df)
new_df <- mycount_df %>% left_join(.,gene_df, by="ENSEMBL")
mat <- mycount_df %>% filter(.,rownames(mycount_df) %in% list) %>% 
  left_join(.,gene_df, by="ENSEMBL") %>% group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), sum)) %>% 
  column_to_rownames(., var = "SYMBOL") %>% 
  select(ip_L_V_L_CON,ip_L_V_L_DMS,ip_L_V_L_CO,ip_L_V_L_BAP,ip_L_V_L_CO_BAP)

mat_scale <- mat %>% t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()


col <- colnames(mat_scale)
# agent name
agent <- factor (
  str_replace_all(col, c("ip_L_V_L_|ip_Y_V_S_"='')),
  levels=c('CON','DMS',"AZA","DAC",'AS',"CO","LCD","HCD","BAP",
           "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))

# clone name
clone <- factor( 
  
  str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                       "D" = "DEL19", "Y" = "YAP")),
  
  levels=c('WT','L858R',"DEL19","YAP"))

# draw hp
ha <- HeatmapAnnotation(agent = agent, clone = clone,
                        col = list(agent=c('CON'='#E0E0E0', 'DMS'='#ADADAD', 'AZA'='#FFD2D2', 'DAC'='#FF9797',
                                           'AS'='#FFFF37','CO'='#FF5151','LCD'='#0080FF',
                                           'HCD'='#005AB5', 'BAP'='#00DB00', 'AS_BAP'='#FFDC35', 'CO_BAP'= '#EA0000',
                                           'LCD_BAP'= '#0000E3', "HCD_BAP"="#000079"),
                                   clone=c('WT'="#FF2D2D",'L858R'="#FF9224",
                                           "DEL19"= "#66B3FF","YAP"="#2828FF")))

ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                        show_row_names = T, show_column_names = F,
                        name = "Z-score"
)


# draw heatmap function ---------------------------------------------------
### function
draw_heatmap <- function(data_path=data_path,
                         file=file,
                         log_crit=3,
                         groups=groups
){
  data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
  data <- data %>% filter(abs(M)>log_crit)  # filter log2FC criteria
  list <- data$ENSEMBL  # 抓出差異ensembl id
  gene_df <- data[,c("ENSEMBL","SYMBOL")] # 
  mycount_df$ENSEMBL <- rownames(mycount_df)
  new_df <- mycount_df %>% left_join(.,gene_df, by="ENSEMBL")
  mat <- mycount_df %>% filter(.,rownames(mycount_df) %in% list) %>% 
    left_join(.,gene_df, by="ENSEMBL") %>% group_by(SYMBOL) %>%
    summarize(across(where(is.numeric), sum)) %>% 
    column_to_rownames(., var = "SYMBOL")
  
  if(groups == "AS"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AS,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP)
  } else if(groups == "CO"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_CO,ip_Y_V_S_BAP,ip_Y_V_S_CO_BAP)
  } else if(groups == "CD"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
                               ip_Y_V_S_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else if(groups == "BAP"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP
                               ,ip_Y_V_S_CO_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else if(groups == "ALL"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
                               ip_Y_V_S_AS,ip_Y_V_S_CO,
                               ip_Y_V_S_LCD,ip_Y_V_S_HCD,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP,
                               ip_Y_V_S_CO_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else{
    cat(c("<- please type in groups", "\n"))
    break
  }

  cat(c("<- filtering data", "\n"))
  
  mat_scale <- data_mat %>% t() %>% scale(scale = T) %>% t() %>% as.matrix() %>% na.omit()
  cat(c("<- annotation", "\n"))
  
  col <- colnames(mat_scale)
  
  # agent name
  agent <- factor (
    str_replace_all(col, c("ip_L_V_L_|ip_Y_V_S_"='')),
    levels=c('CON','DMS',"AZA","DAC",'AS',"CO","LCD","HCD","BAP",
             "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
  
  # clone name
  clone <- factor( 
    str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                         "D" = "DEL19", "Y" = "YAP")),
    levels=c('WT','L858R',"DEL19","YAP"))
  
  # draw hp
  ha <- HeatmapAnnotation(agent = agent, clone = clone,
                          col = list(agent=c('CON'='#E0E0E0', 'DMS'='#ADADAD', 'AZA'='#FFD2D2', 'DAC'='#FF9797',
                                             'AS'='#FFFF37','CO'='#FF5151','LCD'='#0080FF',
                                             'HCD'='#005AB5', 'BAP'='#00DB00', 'AS_BAP'='#FFDC35', 'CO_BAP'= '#EA0000',
                                             'LCD_BAP'= '#0000E3', "HCD_BAP"="#000079"),
                                     clone=c('WT'="#FF2D2D",'L858R'="#FF9224",
                                             "DEL19"= "#66B3FF","YAP"="#2828FF")))
  cat(c("<- drawing heatmap", "\n"))
  
  ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                          show_row_names = T, show_column_names = F,
                          name = "Z-score"
  )
  
}

### use function
draw_heatmap(data_path = data_path ,file = "ip_Y_V_S_BAP_0_deg.xlsx",groups = "CO",log_crit = 4)

draw_heatmap(data_path = data_path ,file = "ip_Y_V_S_CO_0_deg.xlsx",groups = "ALL",log_crit = 2)
# draw heatmap form a list ------------------------------------------------
draw_from_list <- function(list,
                           groups=groups,
                           id="ENSEMBL",
                           cluster=TRUE,
                           anno=TRUE,
                           title=""# ENSEMBL/SYMBOL
){
  if(id=="ENSEMBL"){
    list <- list  
    mycount_df$ENSEMBL <- rownames(mycount_df)
    new_df <- mycount_df %>% left_join(.,gene_df, by="ENSEMBL")
    mat <- mycount_df %>% filter(.,rownames(mycount_df) %in% list) %>% 
      left_join(.,gene_df, by="ENSEMBL") %>% group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>% na.omit() %>% 
      column_to_rownames(., var = "SYMBOL")
  }else if(id=="SYMBOL"){
    list <- list  
    mycount_df$ENSEMBL <- rownames(mycount_df)
    new_df <- mycount_df %>% 
      left_join(.,gene_df, by="ENSEMBL") %>% 
      group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>% 
      na.omit() %>% 
      column_to_rownames(.,var="SYMBOL") 
    mat <- new_df %>% 
      filter(.,rownames(new_df) %in% list)
  }else{
    cat(c("<- error, check gene id", "\n"))
    break
  }
  
  if(groups == "AS"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AS,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP)
  } else if(groups == "CO"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_CO,ip_Y_V_S_BAP,ip_Y_V_S_CO_BAP)
  } else if(groups == "CD"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_LCD,ip_Y_V_S_HCD,
                               ip_Y_V_S_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else if(groups == "BAP"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP
                               ,ip_Y_V_S_CO_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else if(groups == "ALL"){
    data_mat <- mat %>% select(ip_Y_V_S_CON,ip_Y_V_S_DMS,ip_Y_V_S_AZA,ip_Y_V_S_DAC,
                               ip_Y_V_S_AS,ip_Y_V_S_CO,
                               ip_Y_V_S_LCD,ip_Y_V_S_HCD,ip_Y_V_S_BAP,ip_Y_V_S_AS_BAP,
                               ip_Y_V_S_CO_BAP,ip_Y_V_S_LCD_BAP,ip_Y_V_S_HCD_BAP)
  } else{
    cat(c("<- please type in groups", "\n"))
    break
  }
  cat(c("<- filtering data", "\n"))
  
  mat_scale <- data_mat %>% t() %>% scale() %>% t() %>% as.matrix() %>% na.omit()
  cat(c("<- annotation", "\n"))
  
  col <- colnames(mat_scale)
  
  # agent name
  agent <- factor (
    str_replace_all(col, c("ip_L_V_L_|ip_Y_V_S_"='')),
    levels=c('CON','DMS',"AZA","DAC",'AS',"CO","LCD","HCD","BAP",
             "AS_BAP","CO_BAP","LCD_BAP","HCD_BAP"))
  
  # clone name
  clone <- factor( 
    str_replace_all(str_sub(col, 4,4), c("W" = "WT", "L" = "L858R", 
                                         "D" = "DEL19", "Y" = "YAP")),
    levels=c('WT','L858R',"DEL19","YAP"))
  
  # draw hp
  ha <- HeatmapAnnotation(agent = agent, clone = clone,
                          col = list(agent=c('CON'='#E0E0E0', 'DMS'='#ADADAD', 'AZA'='#FFD2D2', 'DAC'='#FF9797',
                                             'AS'='#FFFF37','CO'='#FF7575','LCD'='#0080FF',
                                             'HCD'='#005AB5', 'BAP'='#00DB00', 'AS_BAP'='#FFD306', 'CO_BAP'= '#FF0000',
                                             'LCD_BAP'= '#0000E3', "HCD_BAP"="#000079"),
                                     clone=c('WT'="#FF2D2D",'L858R'="#FF9224",
                                             "DEL19"= "#66B3FF","YAP"="#2828FF")))
  cat(c("<- drawing heatmap", "\n"))
  
  if(anno==T){
    ComplexHeatmap::Heatmap(mat_scale, top_annotation = ha, cluster_columns = F, 
                            show_row_names = T, show_column_names = F,cluster_rows = cluster,
                            name = "Z-score", row_title = title
    )
  }else{
    ComplexHeatmap::Heatmap(mat_scale, cluster_columns = F, 
                            show_row_names = T, show_column_names = F,cluster_rows = cluster,
                            name = "Z-score", row_title = title
    )
  }
}

list <- c("CYP1A1","CYP1B1","GADD45B")
draw_from_list(list = list,groups = "ALL",id = "SYMBOL")


# germ layer marker gene list --------------------------------------------------------
ecto_list1 <- read.table("./marker_list/List_Gifford_EctoMarkers.txt") %>% as.list() %>% unlist()
endo_list1 <- read.table("./marker_list/List_Gifford_EndoMarkers.txt") %>% as.list() %>% unlist()
meso_list1 <- read.table("./marker_list/List_Gifford_MesoMarkers.txt") %>% as.list() %>% unlist()
pluri_list1 <- read.table("./marker_list/List_Gifford_PluriMarkers.txt") %>% as.list() %>% unlist()

ecto_list2 <- read.table("./marker_list/List_Hutchins_EctoMarkers.txt") %>% as.list() %>% unlist()
endo_list2 <- read.table("./marker_list/List_Hutchins_EndoMarkers.txt") %>% as.list() %>% unlist()
meso_list2 <- read.table("./marker_list/List_Hutchins_MesoMarkers.txt") %>% as.list() %>% unlist()
pluri_list2 <- read.table("./marker_list/List_Hutchins_PluriMarkers.txt") %>% as.list() %>% unlist()

pluri_list3 <- read.table("./marker_list/List_Sperger_PluriMarkers.txt") %>% as.list() %>% unlist()

ecto <- c("CDH9","COL2A1","DMBX1","DRD4","EN1","LMX1A","MAP2","NR2F2","MYO3B",
          "NOS2","NR2F1","NR2F2","PAX3","PAX6","POU4F1","OLFM3","PAPLN","PRKCA",
          "SDC2","SOX1","TRPM8","WNT1","ZBTB16")
endo <- c("AFP","CABP7","CDH20","CLDN1","CPLX2","ELAVL3","EOMES","FOXA1","FOXA2",
          "FOXP2","GATA4","GATA6","HNF1B","HHEX","HNF4A","KLF5","LEFTY1","LEFTY2",
          "NODAL","PHOX2B","POU3F3","PRDM1","RXRG","SOX17","SST","HMP19")
meso <- c("ABCA4","ALOX15","BMP10","CDH5","CDX2","ESM1","FCN3","FOXF1","HAND1",
          "HAND2","HEY1","RGS4","ODAM","NKX2-5","PDGFRA","PLVAP","SNAI2","TBX3",
          "TM4SF1","COLEC10","HOPX","IL6ST")
mesendo <- c("FGF4","GDF3","NPPB","NR5A2","PTHLH","BRACHYURY")
pluri <- c("CXCL5","DNMT3B","HESX1","IDO1","LCK","NANOG","POU5F1","SOX2","TRIM22")

p1 <- draw_from_list(list = ecto,groups = "ALL",
                     id = "SYMBOL",title ="ectoderm")
p2 <- draw_from_list(list = endo,groups = "ALL",
                     id = "SYMBOL",title ="endoderm",anno = F)
p3 <- draw_from_list(list = meso,groups = "ALL",
                     id = "SYMBOL",title ="mesoderm",anno = F)
p4 <- draw_from_list(list = pluri,groups = "ALL",
                     id = "SYMBOL",title ="pluripotency",anno = F)
p1 %v% p2 %v% p3 %v% p4

# from gene card
genelist <- read.csv("/Users/benson/Downloads/GeneCards-SearchResults.csv") %>% 
  select(Gene.Symbol) %>% 
  as.list() %>% unlist()
draw_from_list(list = genelist[1:100],groups = "CO",id = "SYMBOL")

# get DEG list ------------------------------------------------------------
get_deg <- function(data_path=data_path,
                    file=file,
                    log_crit = c(1,-1),
                    dir="all"  # all/up/down
){
  data <- readxl::read_xlsx(file.path(data_path, file)) # 讀檔
  if(dir=="all"){
    data <- data %>% filter(M >log_crit[1]| M < log_crit[2])
  }else if(dir=="up"){
    data <- data %>% filter(M > log_crit[1])
  }else if(dir=="down"){
    data <- data %>% filter(M < log_crit[2])
  }else{
    cat(c("<- error, check direction", "\n"))
  }
  list <- data$ENSEMBL 
  return(list)# 抓出差異ensembl id
}

BAP_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_BAP_0_deg.xlsx",
                   log_crit = c(1,-1),dir = "down")
CO_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_CO_0_deg.xlsx",
                  log_crit = c(1,-1),dir = "down")
CO_BAP_down <- get_deg(data_path = data_path,file = "ip_Y_V_S_CO_BAP_0_deg.xlsx",
                  log_crit = c(1,-1),dir = "down")

draw_from_list(list = CO_BAP_down,groups = "ALL",id = "ENSEMBL")

# venn diagram ------------------------------------------------------------
CO_gene <- list(CO = CO_down,
                BAP = BAP_down,
                CO_BAP =CO_BAP_down)

ggVennDiagram(CO_gene,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#FF2D2D")
ggVennDiagram(CO_gene,label_percent_digit = 1,label_alpha = 0) +
  scale_fill_gradient(low="white",high = "#6A6AFF")

CO_list <- process_region_data(Venn(CO_gene))
draw_from_list(list = CO_list$item[[7]],groups = "CO",id = "ENSEMBL")

# gene list ---------------------------------------------------------------
meta_list <- c("CYP1A1","CYP1B1","CYP1A2","CYP2A6","CYP2B6","CYP2C19","CYP2C8",
               "CYP2C9","CYP2D6","CYP2E1","CYP3A4","CYP3A5","DPYD","GSTM1","NAT1",
               "SLC15A","SLC22A","SLCO1B1","TPMT","UGT1A1","UGT2B15",
               "UGT2B7","ABCB1","ABCC2","ABCG2")
poly_list <- c("POLA1","POLB","POLD1","POLD2","POLD3","POLD4","POLE1","POLE2",
               "POLE3","POLE4","POLZ","REV7","POLH","POLI","POLQ","POLK","POLL",
               "POLM","POLN","PRIMPOL","DNTT")
ber_list <- c("UNG","SMUG1","MBD4","TDG","OGG1","MYH","NTH1","MPG","NEIL1","NEIL2",
              "NEIL3","APEX1","APEX2","LIG3","XRCC1","PNKP","APLF","HMCES")
ddr_list <- c("PARP1","PARP2","PARP3","PARG","PARPBP","MGMT","ALKBH2","ALKBH3",
              "TDP1","TDP2","SPRTN","NUDT1","DUT","RRM2B","PARK7","DNPH1","NUDT15","NUDT18")
mmr_list <- c("MSH2","MSH3","MSH6","MLH1","PMS2","MSH4","MSH5","MLH3","PMS1","PMS2P3","HFM1")
ner_list <- c("XPC","RAD23B","CETN2","RAD23A","XPA","DDB1","DDB2","RPA1","RPA2",
              "RPA3","ERCC3","ERCC2","GTF2H1","GTF2H2","GTF2H3","GTF2H4","GTF2H5","GTF2E2",
              "CDK7","CCNH","MNAT1","ERCC5","ERCC1","ERCC4","LIG1","ERCC8","ERCC6","UVSSA",
              "XAB2","MMS19")
hr_list <- c("RAD51","RAD51B","RAD51D","HELQ","SWI5","SWSAP1","ZSWIM7","SPIDR",
             "PDS5B","DMC1","XRCC2","XRCC3","RAD52","RAD54L","RAD54B","BRCA1",
             "BARD1","ABRAXAS1","PAXIP1","SMC5","SMC6","SHLD1","SHLD2","SHLD3",
             "SEM1","RAD50","MRE11A","NBN","RBBP8","MUS81","EME1","EME2","SLX1A",
             "SLX1B","GEN1")
fa_list <- c("FANCA","FANCB","FANCC","BRCA2","FANCD2","FANCE","FANCF","FANCG","FANCI",
             "BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4","FAAP20","FAAP24","FAAP100",
             "UBE2T")
nhej_list <- c("XRCC6","XRCC5","PRKDC","LIG4","XRCC4","DCLRE1C","NHEJ1")
chrom_list <- c("H2AX","CHAF1A","SETMAR","ATRX")
cc_list <- c("ATM", "ATR", "BUB1", "BUB1B", "BUB3", "CCNA1", "CCNA2", 
             "CCNB1", "CCNB2", "CCNB3", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", 
             "CCNH", "CDC14A", "CDC14B", "CDC16", "CDC20", "CDC23", "CDC25A", "CDC25B",
             "CDC25C", "CDC26", "CDC27", "CDC6", "CDC7", "CDK2", "CDK4", "CDK6", 
             "CDK7", "CDKN1A", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", 
             "CDKN2D", "CHEK1", "CHEK2", "CREBBP", "CUL1","DBF4", "E2F1", "E2F2", 
             "E2F3", "EP300", "ESPL1", "FZR1", "GADD45A", "GADD45B", "GADD45G", 
             "GSK3B", "HDAC1", "HDAC2", "MAD1L1", "MAD2L1", "MAD2L2", "MCM2", "MCM3",
             "MCM4", "MCM5", "MCM6", "MCM7", "MDM2", "PCNA", "PKMYT1", "PLK1", "PRKDC",
             "PTTG1", "PTTG2", "RB1", "RBL1", "RBL2", "RBX1", "SFN", "SKP1", "SKP2", 
             "SMAD2", "SMAD3", "SMAD4", "SMC1A", "SMC1B", "TFDP1", "TGFB1", "TGFB2", 
             "TGFB3", "TP53")

draw_from_list(list = cc_list,groups = "ALL",id = "SYMBOL",anno = F)

p1 <- draw_from_list(list = ner_list,groups = "ALL",
                     id = "SYMBOL",title ="NER")
p2 <- draw_from_list(list = nhej_list,groups = "ALL",
                     id = "SYMBOL",title ="NHEJ",anno = F)
p1 %v% p2
