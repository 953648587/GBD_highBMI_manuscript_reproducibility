# ============================================================
# 01a_figure1.R
# Figure 1: Age-standardized DALY rates attributable to high BMI
#           in five G20 countries, 1990-2023
# ============================================================

library(tidyverse)

data_dir <- "C:/Users/95364/Desktop/trend_analysis"

# ---------- 1. Load and filter ASR ----------
files <- c(
  Gout           = file.path(data_dir, "g20_bmi_gout_raw_1990-2023.csv"),
  Low_back_pain  = file.path(data_dir, "g20_bmi_lbp_raw_1990-2023.csv"),
  Osteoarthritis = file.path(data_dir, "g20_bmi_oa_raw_1990-2023.csv")
)

load_asr <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  df_f <- df[
    df[[10]] == "Age-standardized" &
      df[[8]]  == "Both" &
      df[[4]]  == "DALYs (Disability-Adjusted Life Years)" &
      df[[16]] == "Rate",
  ]
  df_out <- df_f[, c(5, 6, 17, 18, 19, 20)]
  names(df_out) <- c("location_id", "location_name", "year",
                     "val", "upper", "lower")
  as_tibble(df_out)
}

asr_data <- list(
  Gout           = load_asr(files[["Gout"]]),
  Low_back_pain  = load_asr(files[["Low_back_pain"]]),
  Osteoarthritis = load_asr(files[["Osteoarthritis"]])
)

# ---------- 2. Select 5 representative countries ----------
representative <- c(
  "China", "Japan", "Mexico",
  "United Kingdom", "United States of America"
)

plot_data <- bind_rows(
  asr_data$Gout           %>% mutate(disease = "Gout"),
  asr_data$Low_back_pain  %>% mutate(disease = "Low back pain"),
  asr_data$Osteoarthritis %>% mutate(disease = "Osteoarthritis")
) %>%
  filter(location_name %in% representative) %>%
  mutate(
    disease       = factor(disease,
                           levels = c("Gout", "Low back pain", "Osteoarthritis")),
    location_name = factor(location_name, levels = representative)
  )

# ---------- 3. Colour mapping ----------
color_map <- c(
  "China"                    = "#D62728",
  "Japan"                    = "#2CA02C",
  "Mexico"                   = "#FF7F0E",
  "United Kingdom"           = "#9467BD",
  "United States of America" = "#1F77B4"
)

# ---------- 4. Figure 1 ----------
fig1 <- ggplot(plot_data, aes(x = year, y = val,
                              color = location_name,
                              fill  = location_name)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ disease, scales = "free_y", ncol = 3) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values  = color_map) +
  labs(
    title = "Age-standardized DALY rates attributable to high BMI in five G20 countries, 1990–2023",
    x     = "Year",
    y     = "ASR per 100,000",
    color = "Country",
    fill  = "Country"
  ) +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(fig1)

# ---------- 5. Save ----------
ggsave(file.path(data_dir, "Figure1.png"), fig1,
       width = 12, height = 5, dpi = 300, bg = "white")
ggsave(file.path(data_dir, "Figure1.pdf"), fig1,
       width = 12, height = 5, device = "pdf", bg = "white")

cat("\nFigure 1 saved: Figure1.png / Figure1.pdf\n")