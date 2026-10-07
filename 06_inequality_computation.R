# ============================================================
# Recompute SII, CI, and their significance tests
# from raw GBD 2023 data (G20, 2023)
#
# Fixes:
#   1. Correct population file: population_detailed.csv (19 countries)
#   2. Sum population across age groups AND sexes
#   3. Türkiye unified
# ============================================================

setwd("C:/Users/95364/Desktop/inequality_analysis")
cat("Working directory:", getwd(), "\n\n")

library(tidyverse)

# ------------------------------------------------------------
# 0. G20 sovereign member states (19 countries)
# ------------------------------------------------------------

g20_countries <- c(
  "Argentina", "Australia", "Brazil", "Canada", "China",
  "France", "Germany", "India", "Indonesia", "Italy",
  "Japan", "Mexico", "Republic of Korea", "Russian Federation",
  "Saudi Arabia", "South Africa", "Türkiye", "United Kingdom",
  "United States of America"
)

fix_turkiye <- function(x) {
  x <- gsub("T.*rkiye", "Türkiye", x)
  x <- gsub("T.*ky.*",   "Türkiye", x)
  x
}

# ------------------------------------------------------------
# 1. Extract 2023 ASR
# ------------------------------------------------------------

extract_asr_2023 <- function(path, disease_label) {
  df <- read.csv(path, fileEncoding = "UTF-8")
  df[[6]] <- fix_turkiye(df[[6]])
  d <- df[df[[6]] %in% g20_countries &
            df[[17]] == 2023 &
            df[[16]] == "Rate" &
            df[[10]] == "Age-standardized" &
            df[[8]]  == "Both" &
            df[[4]]  == "DALYs (Disability-Adjusted Life Years)", ]
  out <- d[, c(6, 18)]
  names(out) <- c("location_name", "asr_2023")
  out$disease <- disease_label
  out <- unique(out)
  out[order(out$location_name), ]
}

asr_gout <- extract_asr_2023("g20_bmi_gout_raw_1990-2023.csv", "Gout")
asr_lbp  <- extract_asr_2023("g20_bmi_lbp_raw_1990-2023.csv",  "Low back pain")
asr_oa   <- extract_asr_2023("g20_bmi_oa_raw_1990-2023.csv",  "Osteoarthritis")

cat("=== ASR country count ===\n")
cat("Gout:", nrow(asr_gout), " LBP:", nrow(asr_lbp),
    " OA:", nrow(asr_oa), "\n")

# ------------------------------------------------------------
# 2. Population 2023 (use population_detailed.csv, sum age + sex)
# ------------------------------------------------------------

popu <- read.csv("population_detailed.csv", fileEncoding = "UTF-8")
popu$location_name <- fix_turkiye(popu$location_name)

popu_2023_raw <- popu[popu$year == 2023 &
                        popu$location_name %in% g20_countries, ]

popu_2023 <- aggregate(val ~ location_name,
                       data = popu_2023_raw,
                       FUN = sum, na.rm = TRUE)
names(popu_2023) <- c("location_name", "population")

cat("\n=== Population 2023 ===\n")
print(popu_2023, row.names = FALSE)
cat("Rows:", nrow(popu_2023), "\n")

# ------------------------------------------------------------
# 3. SDI 2023 (with duplicate Turkiye removed)
# ------------------------------------------------------------

sdi <- read.csv("GBD_2023_SDI_1950_2023.csv", fileEncoding = "UTF-8")
sdi$location_name <- fix_turkiye(sdi$location_name)

sdi_2023 <- sdi[sdi$year_id == 2023 &
                  sdi$covariate_name_short == "sdi" &
                  sdi$sex == "Both" &
                  sdi$age_group_name == "All Ages", ]

# Remove duplicate Turkiye (keep location_id = 155)
sdi_2023 <- sdi_2023[!(sdi_2023$location_name == "Türkiye" &
                         sdi_2023$location_id != 155), ]

sdi_2023 <- sdi_2023[sdi_2023$location_name %in% g20_countries,
                     c("location_name", "mean_value")]
names(sdi_2023) <- c("location_name", "sdi_2023")
sdi_2023 <- unique(sdi_2023)

cat("\n=== SDI 2023 ===\n")
print(sdi_2023, row.names = FALSE)
cat("Rows:", nrow(sdi_2023), "\n")

# ------------------------------------------------------------
# 4. Merge
# ------------------------------------------------------------

merge_data <- function(asr_dat, popu_dat, sdi_dat) {
  d <- merge(asr_dat, popu_dat, by = "location_name")
  d <- merge(d, sdi_dat, by = "location_name")
  d <- d[order(d$sdi_2023), ]
  d$mid_rank <- (rank(d$sdi_2023) - 0.5) / nrow(d)
  d
}

dat_gout <- merge_data(asr_gout, popu_2023, sdi_2023)
dat_lbp  <- merge_data(asr_lbp,  popu_2023, sdi_2023)
dat_oa   <- merge_data(asr_oa,   popu_2023, sdi_2023)

cat("\n=== Merged rows ===\n")
cat("Gout:", nrow(dat_gout), " LBP:", nrow(dat_lbp),
    " OA:", nrow(dat_oa), "\n")

# ------------------------------------------------------------
# 5. SII
# ------------------------------------------------------------

compute_sii <- function(dat) {
  fit <- lm(asr_2023 ~ mid_rank, data = dat, weights = population)
  s   <- summary(fit)$coefficients
  list(SII = s["mid_rank", "Estimate"],
       SII_lower = confint(fit)["mid_rank", 1],
       SII_upper = confint(fit)["mid_rank", 2],
       SII_t = s["mid_rank", "t value"],
       SII_p = s["mid_rank", "Pr(>|t|)"])
}

sii_gout <- compute_sii(dat_gout)
sii_lbp  <- compute_sii(dat_lbp)
sii_oa   <- compute_sii(dat_oa)

# ------------------------------------------------------------
# 6. CI
# ------------------------------------------------------------

compute_ci <- function(dat) {
  d <- dat[order(dat$sdi_2023), ]
  d$cum_pop <- cumsum(d$population) / sum(d$population)
  d$burden  <- d$asr_2023 * d$population
  d$cum_burden <- cumsum(d$burden) / sum(d$burden)
  cp <- c(0, d$cum_pop); cb <- c(0, d$cum_burden)
  auc <- sum(diff(cp) * (head(cb, -1) + tail(cb, -1)) / 2)
  ci  <- 1 - 2 * auc
  list(ci = ci, cum_pop = cp, cum_burden = cb)
}

ci_gout <- compute_ci(dat_gout)
ci_lbp  <- compute_ci(dat_lbp)
ci_oa   <- compute_ci(dat_oa)

# ------------------------------------------------------------
# 7. Bootstrap CI
# ------------------------------------------------------------

bootstrap_ci <- function(dat, n_boot = 1000) {
  set.seed(20260126)
  n <- nrow(dat)
  ci_boot <- numeric(n_boot)
  for (i in 1:n_boot) {
    idx <- sample(1:n, n, replace = TRUE)
    d_boot <- dat[idx, ]
    d_boot <- d_boot[order(d_boot$sdi_2023), ]
    cp <- cumsum(d_boot$population) / sum(d_boot$population)
    cb <- cumsum(d_boot$asr_2023 * d_boot$population) /
      sum(d_boot$asr_2023 * d_boot$population)
    cp <- c(0, cp); cb <- c(0, cb)
    auc <- sum(diff(cp) * (head(cb, -1) + tail(cb, -1)) / 2)
    ci_boot[i] <- 1 - 2 * auc
  }
  list(ci = quantile(ci_boot, c(0.025, 0.975)),
       p  = 2 * min(mean(ci_boot <= 0), mean(ci_boot >= 0)))
}

boot_gout <- bootstrap_ci(dat_gout)
boot_lbp  <- bootstrap_ci(dat_lbp)
boot_oa   <- bootstrap_ci(dat_oa)

# ------------------------------------------------------------
# 8. Assemble results
# ------------------------------------------------------------

results <- data.frame(
  disease   = c("Gout", "Low back pain", "Osteoarthritis"),
  SII       = c(sii_gout$SII,       sii_lbp$SII,       sii_oa$SII),
  SII_lower = c(sii_gout$SII_lower, sii_lbp$SII_lower, sii_oa$SII_lower),
  SII_upper = c(sii_gout$SII_upper, sii_lbp$SII_upper, sii_oa$SII_upper),
  SII_t     = c(sii_gout$SII_t,     sii_lbp$SII_t,     sii_oa$SII_t),
  SII_p     = c(sii_gout$SII_p,     sii_lbp$SII_p,     sii_oa$SII_p),
  CI        = c(ci_gout$ci,         ci_lbp$ci,         ci_oa$ci),
  CI_lower  = c(boot_gout$ci[1],    boot_lbp$ci[1],    boot_oa$ci[1]),
  CI_upper  = c(boot_gout$ci[2],    boot_lbp$ci[2],    boot_oa$ci[2]),
  CI_p      = c(boot_gout$p,        boot_lbp$p,        boot_oa$p)
)

cat("\n=== Final SII/CI results (19 countries) ===\n")
print(results, row.names = FALSE, digits = 4)

write.csv(results, "SII_CI_results_full.csv", row.names = FALSE)
cat("\nSaved: SII_CI_results_full.csv\n")