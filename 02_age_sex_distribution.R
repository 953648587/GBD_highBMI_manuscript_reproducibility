# ============================================================
# Age- and sex-specific DALY rates attributable to high BMI, 2023
# Analysis of gout, low back pain, and osteoarthritis in five
# representative countries (China, Japan, Mexico, the United
# Kingdom, and the United States of America)
#
# Input files:
#   g20_bmi_gout_raw_1990-2023.csv
#   g20_bmi_lbp_raw_1990-2023.csv
#   g20_bmi_oa_raw_1990-2023.csv
#
# Output files:
#   Figure5_reproduced.png
#   peak_summary.csv
#   sex_ratio_by_age.csv
# ============================================================

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 2. Read raw data
# ------------------------------------------------------------
base_dir <- "C:/Users/95364/Desktop/joint"

gout_raw <- read_csv(file.path(base_dir, "g20_bmi_gout_raw_1990-2023.csv"))
lbp_raw  <- read_csv(file.path(base_dir, "g20_bmi_lbp_raw_1990-2023.csv"))
oa_raw   <- read_csv(file.path(base_dir, "g20_bmi_oa_raw_1990-2023.csv"))

# ------------------------------------------------------------
# 3. Standardize column names
#    Raw GBD column names: location_name, sex_name, age_name,
#    year, val, measure_name, metric_name, cause_name
# ------------------------------------------------------------
standardize_cols <- function(df) {
  df %>%
    rename(
      location = location_name,
      sex      = sex_name,
      age      = age_name
    )
}

gout_raw <- standardize_cols(gout_raw) %>% mutate(cause = "Gout")
lbp_raw  <- standardize_cols(lbp_raw)  %>% mutate(cause = "Low back pain")
oa_raw   <- standardize_cols(oa_raw)   %>% mutate(cause = "Osteoarthritis")

df_all <- bind_rows(gout_raw, lbp_raw, oa_raw)

# ------------------------------------------------------------
# 4. Define filtering criteria
# ------------------------------------------------------------
target_countries <- c(
  "China", "Japan", "Mexico",
  "United Kingdom", "United States of America"
)

target_causes <- c("Gout", "Low back pain", "Osteoarthritis")

target_ages <- c(
  "20-24 years", "25-29 years", "30-34 years", "35-39 years",
  "40-44 years", "45-49 years", "50-54 years", "55-59 years",
  "60-64 years", "65-69 years", "70-74 years", "75-79 years",
  "80-84 years", "85+ years"
)

# ------------------------------------------------------------
# 5. Filter and harmonize
#    Only DALY rate (per 100,000) is retained.
# ------------------------------------------------------------
df_plot <- df_all %>%
  filter(
    location %in% target_countries,
    year == 2023,
    cause %in% target_causes,
    sex %in% c("Male", "Female"),
    age %in% target_ages,
    measure_name == "DALYs (Disability-Adjusted Life Years)",
    metric_name  == "Rate"
  ) %>%
  mutate(
    location = factor(location, levels = target_countries),
    age      = factor(age, levels = target_ages),
    cause    = factor(cause, levels = target_causes)
  )

cat("Rows retained after filtering:", nrow(df_plot), "\n")
cat("Countries:", paste(unique(df_plot$location), collapse = ", "), "\n")
cat("Causes:",    paste(unique(df_plot$cause),    collapse = ", "), "\n")
cat("Measure:",   unique(df_plot$measure_name), "\n")
cat("Metric:",    unique(df_plot$metric_name),  "\n")

# ------------------------------------------------------------
# 6. Peak age group and peak DALY rate
#    Defined as the age group with the highest age-specific
#    DALY rate within each country x cause x sex stratum.
# ------------------------------------------------------------
peak_summary <- df_plot %>%
  group_by(location, cause, sex) %>%
  slice_max(order_by = val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    country        = location,
    cause,
    sex,
    peak_age_group = age,
    peak_daly_rate = val
  )

print(peak_summary, n = 30)

# ------------------------------------------------------------
# 7. Male-to-female ratio by age group
# ------------------------------------------------------------
sex_ratio_by_age <- df_plot %>%
  select(location, cause, age, sex, val) %>%
  pivot_wider(names_from = sex, values_from = val) %>%
  mutate(ratio_male_female = Male / Female) %>%
  arrange(location, cause, age)

print(sex_ratio_by_age, n = 20)

# ------------------------------------------------------------
# 8. Figure 5: age- and sex-specific DALY rates, 2023
#    Population-pyramid layout with male rates plotted to the
#    left (negative) and female rates to the right (positive).
#    x-axis scales are free across countries.
# ------------------------------------------------------------
df_pyramid <- df_plot %>%
  mutate(val_signed = ifelse(sex == "Male", -val, val))

p <- ggplot(df_pyramid, aes(x = val_signed, y = age, fill = sex)) +
  geom_col(width = 0.8) +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.3) +
  facet_grid(cause ~ location, scales = "free_x") +
  scale_fill_manual(values = c("Male" = "#F4A300", "Female" = "#1F77B4")) +
  scale_x_continuous(
    labels = function(x) abs(x),
    name   = "Age-specific DALY rate per 100,000"
  ) +
  labs(
    title = "Age- and sex-specific DALY rates attributable to high BMI, 2023",
    y     = "Age group",
    fill  = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "gray90", color = NA),
    strip.text         = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    axis.text.y        = element_text(size = 8)
  )

print(p)

ggsave(
  filename = file.path(base_dir, "Figure5_reproduced.png"),
  plot     = p,
  width    = 14,
  height   = 10,
  dpi      = 300
)

# ------------------------------------------------------------
# 9. Export summary tables
# ------------------------------------------------------------
write_csv(
  peak_summary,
  file.path(base_dir, "peak_summary.csv")
)

write_csv(
  sex_ratio_by_age,
  file.path(base_dir, "sex_ratio_by_age.csv")
)