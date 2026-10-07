# ============================================================
# Description: Three-factor decomposition analysis of high BMI-
#   attributable DALYs for gout, osteoarthritis, and low back pain
#   in G20 countries, 1990-2023.
# Data source: GBD 2023 Results Tool (GHDx)
#   http://ghdx.healthdata.org/gbd-results-tool
#   Download date: January 2026
#
# R version: 4.5.2
# Required packages: tidyverse
#
# NOTE: Change the data_path below to your own working directory
# ============================================================

library(tidyverse)

# ==================== SETUP ====================

# Data path (change to your own directory)
data_path <- "C:/Users/95364/Desktop/FEN/"

# G20 sovereign states
g20_countries <- c(
  "China", "United States of America", "India", "Brazil", 
  "Russian Federation", "Mexico", "Japan", "Germany", 
  "United Kingdom", "France", "Italy", "Canada", 
  "Australia", "Republic of Korea", "Indonesia", "Saudi Arabia", 
  "Türkiye", "South Africa", "Argentina"
)

# Age groups (20 years and above)
age_levels <- c("20-24 years", "25-29 years", "30-34 years", "35-39 years",
                "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                "80-84 years", "85-89 years", "90-94 years", "95+ years")


# ==================== STEP 1: READ DATA ====================

# 1.1 Population data
pop_raw <- read_csv(paste0(data_path, "population_detailed.csv"), 
                    show_col_types = FALSE)

# 1.2 Disease burden data
raw_gout <- read_csv(paste0(data_path, "g20_bmi_gout_raw_1990-2023.csv"), 
                     show_col_types = FALSE) %>% mutate(cause = "Gout")

raw_oa <- read_csv(paste0(data_path, "g20_bmi_oa_raw_1990-2023.csv"), 
                   show_col_types = FALSE) %>% mutate(cause = "Osteoarthritis")

raw_lbp <- read_csv(paste0(data_path, "g20_bmi_lbp_raw_1990-2023.csv"), 
                    show_col_types = FALSE) %>% mutate(cause = "Low back pain")

all_data <- bind_rows(raw_gout, raw_oa, raw_lbp)

cat("STEP 1: Data loaded\n")
cat("  all_data rows:", nrow(all_data), "\n")


# ==================== STEP 2: PREPARE POPULATION DATA ====================

pop_g20 <- pop_raw %>%
  filter(location_name %in% g20_countries,
         year %in% c(1990, 2023),
         age_name %in% age_levels) %>%
  group_by(location_name, year, age_name) %>%
  summarise(population = sum(val, na.rm = TRUE), .groups = "drop")

cat("\nSTEP 2: Population data prepared\n")
cat("  Countries in population data:", 
    length(unique(pop_g20$location_name)), "\n")


# ==================== STEP 3: PREPARE DECOMPOSITION DATA ====================

decomp_raw <- all_data %>%
  filter(location_name %in% g20_countries,
         sex_name == "Both",
         metric_name == "Number",
         year %in% c(1990, 2023),
         age_name %in% age_levels) %>%
  select(location_name, cause, age_name, year, val) %>%
  distinct()

decomp_data <- decomp_raw %>%
  left_join(pop_g20, by = c("location_name", "year", "age_name"))

cat("\nSTEP 3: Decomposition data prepared\n")
cat("  decomp_data rows:", nrow(decomp_data), "\n")

# Check for missing population data
missing_check <- decomp_data %>% filter(is.na(population))
cat("  Missing population rows:", nrow(missing_check), "\n")
if (nrow(missing_check) > 0) {
  cat("  Countries with missing data:\n")
  print(unique(missing_check$location_name))
}


# ==================== STEP 4: DECOMPOSITION FUNCTION ====================

decompose_three_factors <- function(df) {
  
  # Extract data for the two years, sorted by age
  data_1990 <- df %>% filter(year == 1990) %>% arrange(age_name)
  data_2023 <- df %>% filter(year == 2023) %>% arrange(age_name)
  
  # Check age group consistency
  if (!all(data_1990$age_name == data_2023$age_name)) stop("Age groups not consistent")
  
  # Extract population and DALYs
  p1 <- data_1990$population
  p2 <- data_2023$population
  y1 <- data_1990$val
  y2 <- data_2023$val
  
  # Calculate age-specific rates
  r1 <- y1 / p1
  r2 <- y2 / p2
  
  # Total population
  P1 <- sum(p1)
  P2 <- sum(p2)
  
  # Age structure (proportion in each age group)
  c1 <- p1 / P1
  c2 <- p2 / P2
  
  # ========== Das Gupta three-factor decomposition ==========
  
  # 1. Population growth effect
  pop_growth <- (P2 - P1) * (sum(c1 * r1) + sum(c2 * r2)) / 2
  
  # 2. Aging effect
  aging <- P1 * (sum(c2 * r1) - sum(c1 * r1)) + 
    P2 * (sum(c2 * r2) - sum(c1 * r2)) / 2
  
  # 3. Epidemiological change effect
  epi <- P1 * (sum(c1 * r2) - sum(c1 * r1)) + 
    P2 * (sum(c2 * r2) - sum(c2 * r1)) / 2
  
  # Total change
  total_change <- sum(y2) - sum(y1)
  
  # Proportional adjustment for minor discrepancies
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

cat("\nSTEP 4: Decomposition function defined\n")


# ==================== STEP 5: RUN DECOMPOSITION ====================

cat("\nSTEP 5: Running three-factor decomposition...\n")

decomp_results <- decomp_data %>%
  group_by(location_name, cause) %>%
  nest() %>%
  mutate(
    result = map(data, ~{
      tryCatch({
        decompose_three_factors(.x)
      }, error = function(e) {
        tibble(aging = NA, population_growth = NA, 
               epidemiological = NA, total = NA)
      })
    })
  ) %>%
  unnest(result) %>%
  select(location = location_name, cause, aging, population_growth, 
         epidemiological, total) %>%
  ungroup()

cat("  Completed:", nrow(decomp_results), "country-disease combinations\n")


# ==================== STEP 6: CALCULATE PERCENTAGES ====================

decomp_pct <- decomp_results %>%
  mutate(
    aging_pct = aging / total * 100,
    pop_growth_pct = population_growth / total * 100,
    epi_pct = epidemiological / total * 100
  )

cat("\nSTEP 6: Percentages calculated\n")


# ==================== STEP 7: VIEW RESULTS ====================

# 7.1 All G20 countries
cat("\n========================================\n")
cat("Three-factor decomposition: G20 countries\n")
cat("========================================\n")
print(decomp_pct, n = 60)

# 7.2 China only
cat("\n========================================\n")
cat("Three-factor decomposition: China\n")
cat("========================================\n")
china_result <- decomp_pct %>%
  filter(location == "China") %>%
  select(cause, aging, population_growth, epidemiological, total,
         aging_pct, pop_growth_pct, epi_pct)
print(china_result)


# ==================== STEP 8: SAVE RESULTS ====================

write_csv(decomp_pct, paste0(data_path, "decomposition_3factors_G20.csv"))
write_csv(china_result, paste0(data_path, "decomposition_3factors_China.csv"))

cat("\n========================================\n")
cat("Analysis complete! Files saved:\n")
cat("  -", paste0(data_path, "decomposition_3factors_G20.csv"), "\n")
cat("  -", paste0(data_path, "decomposition_3factors_China.csv"), "\n")
cat("========================================\n")