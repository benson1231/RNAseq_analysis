# library -----------------------------------------------------------------
library(readxl)
library(dplyr)
library(progressr)

# DEG number --------------------------------------------------------------
# 初始化新的数据框来存储结果
results_df <- data.frame(file = character(),
                         count_greater_than_1 = numeric(),
                         count_less_than_minus_1 = numeric(),
                         stringsAsFactors = FALSE)

# 遍历 name_df 中的每一个文件
num <- c(3:13,16:26,29:39,42:52)
# 初始化 progressr
handlers(global = TRUE)
handlers("txtprogressbar")

# 遍历文件
with_progress({
  pb <- progressor(along = name_df$file_name[num])
  
  for (file_path in name_df$file_name[num]) {
    # 更新进度条
    pb()
    
    # 读取 Excel 文件
    data <- read_xlsx(file.path(data_path, file_path))
    
    # 过滤 M 列大于 1 的行
    filtered_data_greater_than_1 <- data %>% filter(M > 1) %>% 
      filter(!(SYMBOL %in% c("havana", "ensembl_havana", "havana_tagene"))) %>% 
      group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) 
    count_greater_than_1 <- nrow(filtered_data_greater_than_1)
    
    # 过滤 M 列小于 -1 的行
    filtered_data_less_than_minus_1 <- data %>% filter(M < -1) %>% 
      filter(!(SYMBOL %in% c("havana", "ensembl_havana", "havana_tagene"))) %>% 
      group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) 
    count_less_than_minus_1 <- nrow(filtered_data_less_than_minus_1)
    
    # 将结果添加到 results_df
    results_df <- rbind(results_df, data.frame(file = file_path, 
                                               count_greater_than_1 = count_greater_than_1,
                                               count_less_than_minus_1 = count_less_than_minus_1,
                                               stringsAsFactors = FALSE))
  }
})

# 查看结果
print(results_df)
results_df$group <- name_df$abbreviate[num]
head(results_df)

# 将数据转为长格式
df_long <- results_df %>%
  mutate(count_less_than_minus_1 = -count_less_than_minus_1) %>%
  pivot_longer(cols = c(count_greater_than_1, count_less_than_minus_1),
               names_to = "count_type", values_to = "count")
df_long$group <- factor(df_long$group, levels = results_df$group)

# 绘制柱状图
ggplot(df_long, aes(x = group, y = count, fill = count_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = abs(count)), vjust = ifelse(df_long$count > 0, -0.3, 1.3), position = position_dodge(width = 0.5)) +
  scale_y_continuous(labels = abs) +  # 使 y 轴标签显示为正值
  labs(title = "DEG",
       x = "Group",
       y = "Count",
       fill = "Count Type") +
  scale_fill_manual(values = c("count_greater_than_1" = "red", "count_less_than_minus_1" = "darkblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# DEG to excel ------------------------------------------------------------
# 遍历 name_df 中的每一个文件
num <- c(3:13,16:26,29:39,42:52)
top <- 100
logcrit <- 1
# 初始化 progressr
handlers(global = TRUE)
handlers("txtprogressbar")
# 初始化新的数据框来存储结果
results_DEG <- data.frame(num = 1:top)
# 遍历文件
with_progress({
  pb <- progressor(along = name_df$file_name[num])
  for (file_path in name_df$file_name[num]) {
    # 更新进度条
    pb()
    
    # 读取 Excel 文件
    data <- read_xlsx(file.path(data_path, file_path))
    
    # 过滤 M 列大于 1 的行
    filtered_data_greater_than_1 <- data %>% filter(M > logcrit) %>% 
      filter(!(SYMBOL %in% c("havana", "ensembl_havana", "havana_tagene"))) %>% 
      group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>%
      arrange(desc((M))) %>% head(top)
    
    # 过滤 M 列小于 -1 的行
    filtered_data_less_than_minus_1 <- data %>% filter(M < -logcrit) %>% 
      filter(!(SYMBOL %in% c("havana", "ensembl_havana", "havana_tagene"))) %>% 
      group_by(SYMBOL) %>%
      summarize(across(where(is.numeric), sum)) %>%
      arrange(desc((M))) %>% head(top)
    
    # 创建新的列名
    up_colname <- paste0(abb(file_path,"file"), "_up")
    down_colname <- paste0(abb(file_path,"file"), "_down")
    
    # 将排序后的 M 列添加到 results_DEG 数据框
    results_DEG[[up_colname]] <- filtered_data_greater_than_1$SYMBOL
    results_DEG[[down_colname]] <- filtered_data_less_than_minus_1$SYMBOL
  }
})
print(results_DEG)

writeLines(results_DEG$W_S_CO_up, "test.txt")
