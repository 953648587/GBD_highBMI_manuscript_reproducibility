# ============================================================
# Full reproducible script (DALYs-only version)
#  1) Joinpoint analysis (segmented log-linear regression, R `segmented`)
#  2) Fixed absolute threshold (75th percentile of G20 2023 ASR, DALYs only)
#  3) Monte Carlo projection for China (2024–2050)
#  4) Year reaching threshold under three scenarios
#
# Output files (saved to C:/Users/95364/Desktop/joint/):
#   - additional_file_10_new.csv
#   - abs_threshold_75.csv
#   - china_forecast_final.csv
#   - china_convergence_final.csv
# ============================================================

# ---------- 0. Packages ----------
library(tidyverse)
library(segmented)

# ---------- 1. Paths ----------
data_dir <- "C:/Users/95364/Desktop/joint"

files <- c(
  Gout           = file.path(data_dir, "g20_bmi_gout_raw_1990-2023.csv"),
  Low_back_pain  = file.path(data_dir, "g20_bmi_lbp_raw_1990-2023.csv"),
  Osteoarthritis = file.path(data_dir, "g20_bmi_oa_raw_1990-2023.csv")
)

# ---------- 2. Load and filter ASR (DALYs only) ----------
# Filters:
#   age_name     == "Age-standardized"
#   sex_name     == "Both"
#   metric_name  == "Rate"
#   measure_name == "DALYs (Disability-Adjusted Life Years)"   <- NEW
# Columns used by position:
#   6 = location_name, 17 = year, 18 = val, 19 = lower, 20 = upper

load_asr <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  df_filtered <- df[
    df[[10]] == "Age-standardized" &
      df[[8]]  == "Both" &
      df[[16]] == "Rate" &
      df[[4]]  == "DALYs (Disability-Adjusted Life Years)",
  ]
  df_out <- df_filtered[, c(6, 17, 18, 19, 20)]
  names(df_out) <- c("location_name", "year", "val", "lower", "upper")
  df_out <- df_out[order(df_out$location_name, df_out$year), ]
  as_tibble(df_out)
}

asr_data <- list(
  Gout           = load_asr(files[["Gout"]]),
  Low_back_pain  = load_asr(files[["Low_back_pain"]]),
  Osteoarthritis = load_asr(files[["Osteoarthritis"]])
)

# Sanity check: each country should have 34 rows per disease
cat("\n===== 加载后行数 =====\n")
cat("Gout:           ", nrow(asr_data$Gout), " (expect 19 × 34 = 646)\n")
cat("Low back pain:  ", nrow(asr_data$Low_back_pain), "\n")
cat("Osteoarthritis: ", nrow(asr_data$Osteoarthritis), "\n")

# ---------- 3. G20 median (19 sovereign countries) ----------
g20_countries <- c(
  "Argentina", "Australia", "Brazil", "Canada", "China", "France",
  "Germany", "India", "Indonesia", "Italy", "Japan", "Mexico",
  "Russian Federation", "Saudi Arabia", "South Africa",
  "Republic of Korea", "Türkiye", "United Kingdom",
  "United States of America"
)

make_g20 <- function(df) {
  sub <- df[df$location_name %in% g20_countries, ]
  sub %>%
    group_by(year) %>%
    summarise(
      val   = median(val,   na.rm = TRUE),
      lower = median(lower, na.rm = TRUE),
      upper = median(upper, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(location_name = "G20")
}

asr_data <- lapply(asr_data, function(df) bind_rows(df, make_g20(df)))

# ---------- 4. Joinpoint regression ----------
run_segmented <- function(df, series_label, K_max = 5) {
  d <- df[, c("year", "val")]
  d <- d[order(d$year), ]
  d <- d[!is.na(d$val), ]
  d <- as.data.frame(d)
  
  if (nrow(d) < 10) {
    return(list(apc = NULL, seg_fit = NULL, d = d))
  }
  
  fit <- lm(log(val) ~ year, data = d)
  
  seg_fit <- tryCatch(
    suppressWarnings(segmented(fit, seg.Z = ~year, npsi = K_max)),
    error = function(e) NULL
  )
  
  if (is.null(seg_fit)) {
    apc <- 100 * (exp(coef(fit)[2]) - 1)
    return(list(
      apc = tibble(series = series_label, segment = "1990-2023",
                   APC = apc, CI_lower = NA, CI_upper = NA),
      seg_fit = NULL, d = d
    ))
  }
  
  psi_est     <- sort(seg_fit$psi[, "Est."])
  year_breaks <- c(min(d$year), round(psi_est), max(d$year) + 1)
  seg_labels  <- paste0(year_breaks[-length(year_breaks)], "-",
                        year_breaks[-1] - 1)
  
  slope_df <- as.data.frame(slope(seg_fit, APC = TRUE)$year)
  
  apc_tab <- tibble(
    series   = series_label,
    segment  = seg_labels[seq_len(nrow(slope_df))],
    APC      = slope_df[[1]],
    CI_lower = slope_df[[2]],
    CI_upper = slope_df[[3]]
  )
  
  return(list(apc = apc_tab, seg_fit = seg_fit, d = d))
}

series_list <- list(
  "China-Gout"           = list(loc = "China", dis = "Gout"),
  "China-Low_back_pain"  = list(loc = "China", dis = "Low_back_pain"),
  "China-Osteoarthritis" = list(loc = "China", dis = "Osteoarthritis"),
  "G20-Gout"             = list(loc = "G20",   dis = "Gout"),
  "G20-Low_back_pain"    = list(loc = "G20",   dis = "Low_back_pain"),
  "G20-Osteoarthritis"   = list(loc = "G20",   dis = "Osteoarthritis")
)

seg_results <- lapply(names(series_list), function(label) {
  loc <- series_list[[label]]$loc
  dis <- series_list[[label]]$dis
  df  <- asr_data[[dis]]
  df  <- df[df$location_name == loc, ]
  run_segmented(df, label)
})
names(seg_results) <- names(series_list)

joinpoint_results <- map_dfr(seg_results, "apc")
write_csv(joinpoint_results,
          file.path(data_dir, "additional_file_10_new.csv"))

# ---------- 5. Fixed absolute threshold (DALYs only) ----------
get_threshold_75 <- function(df) {
  d2023 <- df[df$location_name %in% g20_countries & df$year == 2023, ]
  quantile(d2023$val, 0.75, na.rm = TRUE)
}

abs_threshold <- tibble(
  series       = c("China-Gout", "China-Low_back_pain", "China-Osteoarthritis"),
  disease      = c("Gout", "Low back pain", "Osteoarthritis"),
  threshold_75 = c(
    get_threshold_75(asr_data$Gout),
    get_threshold_75(asr_data$Low_back_pain),
    get_threshold_75(asr_data$Osteoarthritis)
  )
)
write_csv(abs_threshold, file.path(data_dir, "abs_threshold_75.csv"))

# ---------- 6. Monte Carlo projection for China ----------
last_segment <- joinpoint_results %>%
  group_by(series) %>%
  slice_tail(n = 1) %>%
  ungroup()

baseline_2023 <- tibble(
  series   = c("China-Gout", "China-Low_back_pain", "China-Osteoarthritis"),
  asr_2023 = c(7.34, 42.16, 48.51),
  ui_low   = c(4.22, 6.26, 0.01),
  ui_high  = c(11.39, 86.30, 129.57),
  dist     = c("lnorm", "lnorm", "truncnorm")
)

sim_input <- baseline_2023 %>%
  left_join(last_segment, by = "series") %>%
  left_join(abs_threshold[, c("series", "threshold_75")], by = "series")

# truncated normal sampler (base R only)
rtruncnorm_manual <- function(n, mean, sd) {
  a <- -mean / sd
  u <- runif(n)
  p_lower <- pnorm(a)
  x_std <- qnorm(p_lower + u * (1 - p_lower))
  mean + sd * x_std
}

simulate_one <- function(asr_mean, asr_low, asr_high,
                         apc_mean, apc_low, apc_high,
                         dist_type, n = 1000) {
  set.seed(2024)
  years_out <- 2024:2050
  
  if (dist_type == "lnorm") {
    log_se  <- (log(asr_high) - log(asr_low)) / (2 * 1.96)
    asr_sim <- rlnorm(n, meanlog = log(asr_mean), sdlog = log_se)
  } else {
    asr_sd  <- (asr_high - asr_low) / (2 * 1.96)
    asr_sim <- rtruncnorm_manual(n, mean = asr_mean, sd = asr_sd)
  }
  
  apc_se  <- (apc_high - apc_low) / (2 * 1.96)
  apc_sim <- rnorm(n, mean = apc_mean, sd = apc_se)
  
  future_matrix <- sapply(years_out, function(y) {
    asr_sim * ((1 + apc_sim / 100) ^ (y - 2023))
  })
  future_matrix <- pmin(pmax(future_matrix, 0), 500)
  
  tibble(
    year     = years_out,
    lower_95 = apply(future_matrix, 2, quantile, 0.025, na.rm = TRUE),
    median   = apply(future_matrix, 2, quantile, 0.5,   na.rm = TRUE),
    upper_95 = apply(future_matrix, 2, quantile, 0.975, na.rm = TRUE)
  )
}

china_forecast <- sim_input %>%
  rowwise() %>%
  mutate(sim = list(simulate_one(
    asr_mean  = asr_2023, asr_low = ui_low, asr_high = ui_high,
    apc_mean  = APC, apc_low = CI_lower, apc_high = CI_upper,
    dist_type = dist
  ))) %>%
  unnest(sim) %>%
  as.data.frame()

china_forecast <- china_forecast[, c("series", "year",
                                     "lower_95", "median", "upper_95")]
write_csv(china_forecast, file.path(data_dir, "china_forecast_final.csv"))

# ---------- 7. Year reaching fixed threshold (3 scenarios) ----------
convergence <- china_forecast %>%
  left_join(abs_threshold[, c("series", "threshold_75")], by = "series") %>%
  group_by(series) %>%
  summarise(
    threshold_75      = first(threshold_75),
    year_reach_median = ifelse(any(median   >= threshold_75),
                               min(year[median   >= threshold_75]), NA),
    year_reach_lower  = ifelse(any(lower_95 >= threshold_75),
                               min(year[lower_95 >= threshold_75]), NA),
    year_reach_upper  = ifelse(any(upper_95 >= threshold_75),
                               min(year[upper_95 >= threshold_75]), NA),
    .groups = "drop"
  )
write_csv(convergence, file.path(data_dir, "china_convergence_final.csv"))

# ---------- 8. Console output ----------
cat("\n===== Joinpoint (Additional file 10) =====\n")
print(joinpoint_results, n = 50)

cat("\n===== Fixed threshold (75th percentile, DALYs only) =====\n")
print(abs_threshold)

cat("\n===== China 2050 forecast =====\n")
print(china_forecast[china_forecast$year == 2050, ])

cat("\n===== Convergence (3 scenarios) =====\n")
print(convergence)

cat("\nDone. Outputs saved to:", data_dir, "\n")