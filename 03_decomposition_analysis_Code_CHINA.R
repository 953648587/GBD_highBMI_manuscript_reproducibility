# ============================================================
# Three-factor decomposition: aging + population growth + epidemiological change
# Using the Das Gupta three-factor decomposition method
# ============================================================

library(tidyverse)

# ==================== 第一步：准备数据 ====================

# 1. 读取人口数据
pop_raw <- read_csv("C:/Users/95364/Desktop/FEN/population_detailed.csv", 
                    show_col_types = FALSE)

# 2. 筛选中国、20岁以上、合并性别
age_levels <- c("20-24 years", "25-29 years", "30-34 years", "35-39 years",
                "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                "80-84 years", "85-89 years", "90-94 years", "95+ years")

pop_china <- pop_raw %>%
  filter(location_name == "China",
         year %in% c(1990, 2023),
         age_name %in% age_levels) %>%
  group_by(year, age_name) %>%
  summarise(population = sum(val, na.rm = TRUE), .groups = "drop")

# 3. 读取疾病数据
raw_gout <- read_csv("C:/Users/95364/Desktop/FEN/g20_bmi_gout_raw_1990-2023.csv", 
                     show_col_types = FALSE) %>% mutate(cause = "Gout")
raw_oa <- read_csv("C:/Users/95364/Desktop/FEN/g20_bmi_oa_raw_1990-2023.csv", 
                   show_col_types = FALSE) %>% mutate(cause = "Osteoarthritis")
raw_lbp <- read_csv("C:/Users/95364/Desktop/FEN/g20_bmi_lbp_raw_1990-2023.csv", 
                    show_col_types = FALSE) %>% mutate(cause = "Low back pain")

all_data <- bind_rows(raw_gout, raw_oa, raw_lbp)

# 4. 构建分解数据（去重）
decomp_raw <- all_data %>%
  filter(location_name == "China",
         sex_name == "Both",
         metric_name == "Number",
         year %in% c(1990, 2023),
         age_name %in% age_levels) %>%
  select(location_name, cause, age_name, year, val) %>%
  distinct()

# 5. 合并人口数据
decomp_data <- decomp_raw %>%
  left_join(pop_china, by = c("year", "age_name"))

# 6. 检查
cat("人口缺失检查：\n")
decomp_data %>% filter(is.na(population)) %>% nrow() %>% print()

# ==================== 第二步：三因素分解函数 ====================

decompose_three_factors <- function(df) {
  
  # 提取两年数据并按年龄排序
  data_1990 <- df %>% filter(year == 1990) %>% arrange(age_name)
  data_2023 <- df %>% filter(year == 2023) %>% arrange(age_name)
  
  # 检查年龄组是否一致
  if (!all(data_1990$age_name == data_2023$age_name)) stop("年龄组不一致")
  
  # 提取人口和DALYs
  p1 <- data_1990$population
  p2 <- data_2023$population
  y1 <- data_1990$val
  y2 <- data_2023$val
  
  # 计算率
  r1 <- y1 / p1
  r2 <- y2 / p2
  
  # 总人口
  P1 <- sum(p1)
  P2 <- sum(p2)
  
  # 年龄结构（各年龄组人口比例）
  c1 <- p1 / P1
  c2 <- p2 / P2
  
  # ========== 三因素分解公式 ==========
  
  # 1. 人口增长效应：总人口变化带来的影响
  pop_growth <- (P2 - P1) * (sum(c1 * r1) + sum(c2 * r2)) / 2
  
  # 2. 老龄化效应：年龄结构变化带来的影响
  aging <- P1 * (sum(c2 * r1) - sum(c1 * r1)) + 
    P2 * (sum(c2 * r2) - sum(c1 * r2)) / 2
  
  # 3. 流行病学效应：各年龄组率变化带来的影响
  epi <- P1 * (sum(c1 * r2) - sum(c1 * r1)) + 
    P2 * (sum(c2 * r2) - sum(c2 * r1)) / 2
  
  # 总变化
  total_change <- sum(y2) - sum(y1)
  
  # 微小差异按比例分配
  sum_factors <- pop_growth + aging + epi
  if (abs(sum_factors - total_change) > 1e-6 && sum_factors != 0) {
    scaling <- total_change / sum_factors
    pop_growth <- pop_growth * scaling
    aging <- aging * scaling
    epi <- epi * scaling
  }
  
  tibble(
    aging = aging,
    population_growth = pop_growth,
    epidemiological = epi,
    total = total_change
  )
}

# ==================== 第三步：运行三因素分解 ====================

diseases <- unique(decomp_data$cause)
results <- list()

cat("\n开始三因素分解...\n")
for (d in diseases) {
  df_sub <- decomp_data %>% filter(cause == d)
  tryCatch({
    res <- decompose_three_factors(df_sub)
    res$cause <- d
    results[[d]] <- res
  }, error = function(e) {
    cat("分解失败 -", d, ":", e$message, "\n")
  })
}

decomp_all <- bind_rows(results)

# ==================== 第四步：计算百分比 ====================

decomp_pct <- decomp_all %>%
  mutate(
    aging_pct = aging / total * 100,
    pop_growth_pct = population_growth / total * 100,
    epi_pct = epidemiological / total * 100
  ) %>%
  select(cause, aging, population_growth, epidemiological, total,
         aging_pct, pop_growth_pct, epi_pct)

cat("\n========================================\n")
cat("中国三种疾病三因素分解结果\n")
cat("========================================\n")
print(decomp_pct)

# ==================== 第五步：保存结果 ====================

write_csv(decomp_pct, "C:/Users/95364/Desktop/FEN/decomposition_3factors1.csv")

cat("\n✅ 三因素分解完成！结果已保存到：\n")
cat("   - decomposition_3factors1.csv\n")