# ============================================================
# Figure 5. SII and CI for high-BMI-attributable
# musculoskeletal diseases in G20, 2023
#
# Redrawn with recalculated SII and CI values
# Includes blue regression line label
# ============================================================

setwd("C:/Users/95364/Desktop/inequality_analysis")
cat("Working directory:", getwd(), "\n\n")

library(tidyverse)
library(ggrepel)
library(patchwork)

# ------------------------------------------------------------
# 1. G20 countries
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
# 2. Extract 2023 ASR
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

# ------------------------------------------------------------
# 3. Population 2023
# ------------------------------------------------------------

popu <- read.csv("population_detailed.csv", fileEncoding = "UTF-8")
popu$location_name <- fix_turkiye(popu$location_name)

popu_2023_raw <- popu[popu$year == 2023 &
                        popu$location_name %in% g20_countries, ]
popu_2023 <- aggregate(val ~ location_name,
                       data = popu_2023_raw,
                       FUN = sum, na.rm = TRUE)
names(popu_2023) <- c("location_name", "population")

# ------------------------------------------------------------
# 4. SDI 2023
# ------------------------------------------------------------

sdi <- read.csv("GBD_2023_SDI_1950_2023.csv", fileEncoding = "UTF-8")
sdi$location_name <- fix_turkiye(sdi$location_name)

sdi_2023 <- sdi[sdi$year_id == 2023 &
                  sdi$covariate_name_short == "sdi" &
                  sdi$sex == "Both" &
                  sdi$age_group_name == "All Ages", ]
sdi_2023 <- sdi_2023[!(sdi_2023$location_name == "Türkiye" &
                         sdi_2023$location_id != 155), ]
sdi_2023 <- sdi_2023[sdi_2023$location_name %in% g20_countries,
                     c("location_name", "mean_value")]
names(sdi_2023) <- c("location_name", "sdi_2023")
sdi_2023 <- unique(sdi_2023)

# ------------------------------------------------------------
# 5. Merge
# ------------------------------------------------------------

merge_data <- function(asr_dat) {
  d <- merge(asr_dat, popu_2023, by = "location_name")
  d <- merge(d, sdi_2023, by = "location_name")
  d <- d[order(d$sdi_2023), ]
  d$mid_rank <- (rank(d$sdi_2023) - 0.5) / nrow(d)
  d
}

dat_gout <- merge_data(asr_gout)
dat_lbp  <- merge_data(asr_lbp)
dat_oa   <- merge_data(asr_oa)

# ------------------------------------------------------------
# 6. Compute SII (with fitted line)
# ------------------------------------------------------------

compute_sii <- function(dat) {
  fit <- lm(asr_2023 ~ mid_rank, data = dat, weights = population)
  dat$fitted <- fitted(fit)
  s <- summary(fit)$coefficients
  sii    <- s["mid_rank", "Estimate"]
  ci <- confint(fit)["mid_rank", ]
  list(data = dat, SII = sii, SII_lower = ci[1], SII_upper = ci[2])
}

sii_gout <- compute_sii(dat_gout)
sii_lbp  <- compute_sii(dat_lbp)
sii_oa   <- compute_sii(dat_oa)

# ------------------------------------------------------------
# 7. Compute CI
# ------------------------------------------------------------

compute_ci <- function(dat) {
  d <- dat[order(dat$sdi_2023), ]
  d$cum_pop <- cumsum(d$population) / sum(d$population)
  d$burden  <- d$asr_2023 * d$population
  d$cum_burden <- cumsum(d$burden) / sum(d$burden)
  cp <- c(0, d$cum_pop); cb <- c(0, d$cum_burden)
  auc <- sum(diff(cp) * (head(cb, -1) + tail(cb, -1)) / 2)
  ci  <- 1 - 2 * auc
  list(ci = ci, cum_pop = cp, cum_burden = cb, data = d)
}

ci_gout <- compute_ci(dat_gout)
ci_lbp  <- compute_ci(dat_lbp)
ci_oa   <- compute_ci(dat_oa)

# ------------------------------------------------------------
# 8. Bootstrap CI
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
# 9. Summary table
# ------------------------------------------------------------

summary_tab <- data.frame(
  disease   = c("Gout", "Low back pain", "Osteoarthritis"),
  SII       = c(sii_gout$SII, sii_lbp$SII, sii_oa$SII),
  SII_lower = c(sii_gout$SII_lower, sii_lbp$SII_lower, sii_oa$SII_lower),
  SII_upper = c(sii_gout$SII_upper, sii_lbp$SII_upper, sii_oa$SII_upper),
  CI        = c(ci_gout$ci, ci_lbp$ci, ci_oa$ci),
  CI_lower  = c(boot_gout$ci[1], boot_lbp$ci[1], boot_oa$ci[1]),
  CI_upper  = c(boot_gout$ci[2], boot_lbp$ci[2], boot_oa$ci[2]),
  CI_p      = c(boot_gout$p, boot_lbp$p, boot_oa$p)
)

cat("\n=== Figure 5 SII/CI results ===\n")
print(summary_tab, row.names = FALSE, digits = 4)

write.csv(summary_tab, "Figure5_SII_CI_summary.csv", row.names = FALSE)

# ------------------------------------------------------------
# 10. Prepare plotting data
# ------------------------------------------------------------

diseases <- c("Gout", "Low back pain", "Osteoarthritis")

scatter_data <- rbind(sii_gout$data, sii_lbp$data, sii_oa$data)
scatter_data$disease <- factor(scatter_data$disease, levels = diseases)

curve_data <- data.frame(
  disease    = rep(diseases, times = c(length(ci_gout$cum_pop),
                                       length(ci_lbp$cum_pop),
                                       length(ci_oa$cum_pop))),
  cum_pop    = c(ci_gout$cum_pop, ci_lbp$cum_pop, ci_oa$cum_pop),
  cum_burden = c(ci_gout$cum_burden, ci_lbp$cum_burden, ci_oa$cum_burden)
)
curve_data$disease <- factor(curve_data$disease, levels = diseases)

# SII / CI labels
sii_labels <- summary_tab %>%
  mutate(
    disease = factor(disease, levels = diseases),
    label   = paste0("SII = ", sprintf("%.2f", SII),
                     " (", sprintf("%.2f", SII_lower),
                     ", ", sprintf("%.2f", SII_upper), ")")
  )

ci_labels <- summary_tab %>%
  mutate(
    disease = factor(disease, levels = diseases),
    label   = paste0("CI = ", sprintf("%.3f", CI))
  )

# China highlight
scatter_data$is_china <- factor(
  ifelse(scatter_data$location_name == "China", "China", "Other"),
  levels = c("Other", "China"))

# ------------------------------------------------------------
# 11. Left panels: SII scatter (with regression line label)
# ------------------------------------------------------------

p_left <- ggplot(scatter_data, aes(x = mid_rank, y = asr_2023)) +
  geom_point(aes(size = population, color = is_china), alpha = 0.7) +
  geom_line(aes(y = fitted), color = "#2C7FB8", linewidth = 0.8) +
  geom_text_repel(aes(label = location_name), size = 2.3,
                  box.padding = 0.3, point.padding = 0.3,
                  max.overlaps = 20) +
  geom_text(data = sii_labels,
            aes(x = 0.02, y = Inf, label = label),
            hjust = 0, vjust = 1.4, size = 3, inherit.aes = FALSE) +
  # --- 蓝色回归线标注 ---
  annotate("text", x = 0.55, y = -Inf,
           label = "Weighted regression line (slope = SII)",
           color = "#2C7FB8", size = 2.4, hjust = 0.5, vjust = -1) +
  facet_wrap(~ disease, scales = "free_y", ncol = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_size_continuous(range = c(1, 6), guide = "none") +
  scale_color_manual(values = c("Other" = "grey40", "China" = "red"),
                     guide = "none") +
  labs(x = "Midpoint rank of SDI",
       y = "DALY rate (per 100,000)") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 11, face = "bold"),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------
# 12. Right panels: concentration curves
# ------------------------------------------------------------

p_right <- ggplot(curve_data, aes(x = cum_pop, y = cum_burden)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  geom_line(color = "darkgreen", linewidth = 0.9) +
  geom_text(data = ci_labels,
            aes(x = 0.02, y = 0.98, label = label),
            hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ disease, scales = "fixed", ncol = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Cumulative population proportion",
       y = "Cumulative burden proportion") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 11, face = "bold"),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------
# 13. Combine
# ------------------------------------------------------------

final_fig <- p_left + p_right +
  plot_layout(widths = c(1.4, 1)) +
  plot_annotation(
    title = "Health inequality of high-BMI-attributable burden within G20",
    theme = theme(plot.title = element_text(hjust = 0.5,
                                            size = 13, face = "bold"))
  )

print(final_fig)

# ------------------------------------------------------------
# 14. Save
# ------------------------------------------------------------

ggsave("Figure5_inequality_new.png", final_fig,
       width = 11, height = 12, dpi = 300, bg = "white")
ggsave("Figure5_inequality_new.pdf", final_fig,
       width = 11, height = 12, device = "pdf", bg = "white")

cat("\nSaved: Figure5_inequality_new.png / .pdf\n")