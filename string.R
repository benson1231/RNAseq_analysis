# library -----------------------------------------------------------------
library(STRINGdb)
# https://rpubs.com/HWH/913747
# load string database
# 9606 for Human, 10090 for mouse
string_db <- STRINGdb$new(version = "11.5", species = 9606, 
                          score_threshold = 200, input_directory="")
class(string_db)

# code block --------------------------------------------------------------
data_path <- "/Users/benson/Documents/raw_data/RNA-seq1-3/V"
as <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_AS_0_deg.xlsx")) %>% 
  select(M,SYMBOL) %>% setNames(c("logFC","gene")) %>% as.data.frame()%>% 
  filter(abs(logFC)>1)
co <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_CO_0_deg.xlsx")) %>% 
  select(M,SYMBOL) %>% setNames(c("logFC","gene")) %>% as.data.frame() %>% 
  filter(abs(logFC)>1)
hcd <- readxl::read_xlsx(file.path(data_path, "ip_Y_V_S_HCD_0_deg.xlsx")) %>% 
  select(M,SYMBOL) %>% setNames(c("logFC","gene")) %>% as.data.frame() %>% 
  filter(abs(logFC)>1)

as_map <- string_db$map(as, "gene", removeUnmappedRows = TRUE)
co_map <- string_db$map(co, "gene", removeUnmappedRows = TRUE)
hcd_map <- string_db$map(hcd, "gene", removeUnmappedRows = TRUE)

data <- as_map %>% mutate(color="red")
data <- co_map %>% mutate(color="blue") %>% rbind(data) %>% arrange(desc(abs(logFC)))
data <- hcd_map %>% mutate(color="yellow") %>% rbind(data) %>% arrange(desc(abs(logFC)))

hits <- data$STRING_id[1:200]

# post payload information to the STRING server
payload_id <- string_db$post_payload(data$STRING_id,
                                     colors=data$color )
# display a STRING network png with the "halo"
string_db$plot_network(hits, payload_id=payload_id)

# venn diagram list with string -------------------------------------------
co <- get_df("ip_Y_V_S_CO_0_deg.xlsx",list=CO_list$item[[1]])
bap <- get_df("ip_Y_V_S_BAP_0_deg.xlsx",list=CO_list$item[[2]])
cobap <- get_df("ip_Y_V_S_CO_BAP_0_deg.xlsx",list=CO_list$item[[3]])

co_map <- string_db$map(co, "gene", removeUnmappedRows = TRUE)
bap_map <- string_db$map(bap, "gene", removeUnmappedRows = TRUE)
cobap_map <- string_db$map(cobap, "gene", removeUnmappedRows = TRUE)

data <- co_map %>% mutate(color="red")
data <- bap_map %>% mutate(color="blue") %>% rbind(data) %>% arrange(desc(abs(logFC)))
data <- cobap_map %>% mutate(color="yellow") %>% rbind(data) %>% arrange(desc(abs(logFC)))

hits <- data$STRING_id[1:100]

# post payload information to the STRING server
payload_id <- string_db$post_payload(data$STRING_id,
                                     colors=data$color )
# display a STRING network png with the "halo"
string_db$plot_network(hits, payload_id=payload_id)


# cluster -----------------------------------------------------------------
clustersList <- string_db$get_clusters(data$STRING_id[1:600])
string_db$plot_network(clustersList[[1]])

