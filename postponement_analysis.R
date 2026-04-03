# =============================================================================
# POSTPONEMENT CONTRIBUTION TO TFR DECLINE: CROSS-COUNTRY ANALYSIS
#
# Methodology (following San 2025, https://x.com/san_muly):
#
#   For each birth cohort C, define:
#     F(C)    = sum_a f(a, C+a)          cohort total fertility (quantum)
#     S0(a)   = f(a, 1951+a) / F(1951)  age distribution of the 1951 cohort
#
#   Counterfactual TFR in year t:
#     CTFR(t) = sum_a  F(t-a) * S0(a)
#
#   This asks: "What would period TFR be if every cohort spread its births
#   across ages exactly as the 1951 birth cohort did, while keeping its own
#   total quantum of children?"
#
#   Timing (postponement) effect:
#     TE(t)   = TFR(t) - CTFR(t)
#   Negative TE means postponement is depressing period TFR below the quantum.
#
#   Postponement share of TFR decline (comparing early 1970s to ~2018):
#     PS = [TE(early) - TE(recent)] / [TFR(early) - TFR(recent)]
#
# Data:  Our World in Data / UN WPP
#   fertility-rate-by-age-group.csv   age-specific rates per 1,000 women
#   period-average-age-of-mothers.csv mean age at childbearing
#
# Author: [Author name]
# Date:   April 2026
# =============================================================================

library(tidyverse)
library(ggrepel)

# ---- 0. Paths ----------------------------------------------------------------
data_path <- "."   # folder with CSVs; adjust if needed

# ---- 1. Load and reshape age-specific fertility ------------------------------
asfr_raw <- read_csv(file.path(data_path, "fertility-rate-by-age-group.csv"),
                     show_col_types = FALSE)

# Rename columns to something clean
# Expected columns: Entity, Code, Year, f12, f17, f22, f27, f32, f37, f42, f47, f52
asfr_raw <- asfr_raw %>%
  rename_with(tolower) %>%
  rename(country = entity)

# Drop regional/world aggregates — keep only 3-letter ISO codes
asfr_raw <- asfr_raw %>%
  filter(!is.na(code), nchar(code) == 3)

# Pivot to long format; age = midpoint of the 5-year age group
asfr_long <- asfr_raw %>%
  pivot_longer(cols = starts_with("f"),
               names_to  = "age",
               values_to = "rate_per_1000") %>%
  mutate(age = as.integer(str_remove(age, "f")))

# Convert to TFR contribution: rate per 1000 women × 5-yr width / 1000
asfr_long <- asfr_long %>%
  mutate(tfr_contrib = rate_per_1000 * 5 / 1000)

# Focus on prime reproductive ages 20–44 (midpoints 22, 27, 32, 37, 42)
asfr_long <- asfr_long %>%
  filter(age >= 22, age <= 42)

# ---- 2. Cohort alignment ----------------------------------------------------
# birth_year = calendar year − midpoint age
# Since age ∈ {22,27,32,37,42} (all ≡ 2 mod 5), birth_year ≡ year−2 mod 5.
# Keeping birth_year mod 5 == 1 selects cohorts 1931, 1936, 1941, 1946, 1951,
# 1956, … each appearing exactly once per period year that is ≡ 3 mod 5.
# This avoids overlap between adjacent 5-year age groups.

asfr_cohort <- asfr_long %>%
  mutate(birth_year = year - age) %>%
  filter(birth_year %% 5 == 1)

# ---- 3. Cohort total fertility F(C) -----------------------------------------
cohort_totals <- asfr_cohort %>%
  group_by(country, code, birth_year) %>%
  summarise(F_cohort  = sum(tfr_contrib, na.rm = TRUE),
            n_ages    = n(),
            .groups   = "drop")

asfr_cohort <- asfr_cohort %>%
  left_join(cohort_totals, by = c("country", "code", "birth_year"))

# ---- 4. Reference age distribution S0(a) from the 1951 birth cohort ---------
ref_dist <- asfr_cohort %>%
  filter(birth_year == 1951) %>%
  group_by(country, code) %>%
  mutate(S0 = tfr_contrib / F_cohort) %>%
  select(country, code, age, S0)

# Keep only countries with a complete 1951 reference (all 5 age groups present)
ref_complete <- ref_dist %>%
  group_by(country, code) %>%
  summarise(n_ref = n(), .groups = "drop") %>%
  filter(n_ref >= 5)

asfr_cohort <- asfr_cohort %>%
  semi_join(ref_complete, by = c("country", "code")) %>%
  left_join(ref_dist, by = c("country", "code", "age"))

# ---- 5. Counterfactual births -----------------------------------------------
# cf_births(a, t) = F(t−a) × S0(a)
asfr_cohort <- asfr_cohort %>%
  mutate(cf_births = F_cohort * S0)

# ---- 6. Collapse to period year ---------------------------------------------
tfr_ts <- asfr_cohort %>%
  group_by(country, code, year) %>%
  summarise(
    tfr_actual        = sum(tfr_contrib, na.rm = TRUE),
    tfr_counterfactual = sum(cf_births,  na.rm = TRUE),
    n_ages_observed   = n(),
    .groups = "drop"
  ) %>%
  # Only keep period-years with complete age coverage (all 5 groups)
  filter(n_ages_observed >= 5) %>%
  mutate(timing_effect = tfr_actual - tfr_counterfactual)

# ---- 7. Cross-country postponement summary (early 1970s vs. 2015–21) --------
early  <- tfr_ts %>%
  filter(year >= 1970, year <= 1977) %>%
  group_by(country, code) %>%
  summarise(tfr_early  = mean(tfr_actual,         na.rm = TRUE),
            ctfr_early = mean(tfr_counterfactual,  na.rm = TRUE),
            te_early   = mean(timing_effect,        na.rm = TRUE),
            .groups = "drop")

recent <- tfr_ts %>%
  filter(year >= 2015, year <= 2021) %>%
  group_by(country, code) %>%
  summarise(tfr_recent  = mean(tfr_actual,         na.rm = TRUE),
            ctfr_recent = mean(tfr_counterfactual,  na.rm = TRUE),
            te_recent   = mean(timing_effect,        na.rm = TRUE),
            .groups = "drop")

summary_df <- inner_join(early, recent, by = c("country", "code")) %>%
  mutate(
    delta_tfr          = tfr_early  - tfr_recent,   # >0 if TFR fell
    delta_ctfr         = ctfr_early - ctfr_recent,  # quantum decline
    timing_contrib     = te_early   - te_recent,    # >0 if more postponed
    postponement_share = if_else(delta_tfr > 0.05,
                                 timing_contrib / delta_tfr,
                                 NA_real_)
  )

# World region labels (simplified)
region_lookup <- tribble(
  ~code, ~region,
  # East Asia & Pacific
  "CHN","East Asia & Pacific","JPN","East Asia & Pacific","KOR","East Asia & Pacific",
  "SGP","East Asia & Pacific","AUS","East Asia & Pacific","NZL","East Asia & Pacific",
  "VNM","East Asia & Pacific","THA","East Asia & Pacific","IDN","East Asia & Pacific",
  "PHL","East Asia & Pacific","MYS","East Asia & Pacific","MNG","East Asia & Pacific",
  "KHM","East Asia & Pacific","MMR","East Asia & Pacific","LAO","East Asia & Pacific",
  "PNG","East Asia & Pacific","FJI","East Asia & Pacific","TWN","East Asia & Pacific",
  # Europe
  "DEU","Europe","FRA","Europe","GBR","Europe","ITA","Europe","ESP","Europe",
  "PRT","Europe","BEL","Europe","NLD","Europe","AUT","Europe","CHE","Europe",
  "SWE","Europe","NOR","Europe","DNK","Europe","FIN","Europe","IRL","Europe",
  "GRC","Europe","POL","Europe","CZE","Europe","HUN","Europe","SVK","Europe",
  "SVN","Europe","HRV","Europe","BGR","Europe","ROU","Europe","LTU","Europe",
  "LVA","Europe","EST","Europe","BLR","Europe","UKR","Europe","MDA","Europe",
  "SRB","Europe","BIH","Europe","MKD","Europe","ALB","Europe","MNE","Europe",
  "LUX","Europe","ISL","Europe","MLT","Europe","CYP","Europe","RUS","Europe",
  # Latin America
  "BRA","Latin America","MEX","Latin America","COL","Latin America","ARG","Latin America",
  "PER","Latin America","VEN","Latin America","CHL","Latin America","ECU","Latin America",
  "BOL","Latin America","PRY","Latin America","URY","Latin America","PAN","Latin America",
  "CRI","Latin America","GTM","Latin America","HND","Latin America","SLV","Latin America",
  "NIC","Latin America","DOM","Latin America","CUB","Latin America","HTI","Latin America",
  "JAM","Latin America","TTO","Latin America","GUY","Latin America","SUR","Latin America",
  # North America
  "USA","North America","CAN","North America",
  # South Asia
  "IND","South Asia","PAK","South Asia","BGD","South Asia","NPL","South Asia",
  "LKA","South Asia","AFG","South Asia","MDV","South Asia","BTN","South Asia",
  # Middle East & North Africa
  "TUR","Middle East & N. Africa","IRN","Middle East & N. Africa",
  "EGY","Middle East & N. Africa","DZA","Middle East & N. Africa",
  "MAR","Middle East & N. Africa","TUN","Middle East & N. Africa",
  "SAU","Middle East & N. Africa","IRQ","Middle East & N. Africa",
  "SYR","Middle East & N. Africa","JOR","Middle East & N. Africa",
  "LBN","Middle East & N. Africa","ISR","Middle East & N. Africa",
  "YEM","Middle East & N. Africa","OMN","Middle East & N. Africa",
  "ARE","Middle East & N. Africa","KWT","Middle East & N. Africa",
  "QAT","Middle East & N. Africa","BHR","Middle East & N. Africa",
  "LBY","Middle East & N. Africa"
)

summary_df <- summary_df %>%
  left_join(region_lookup, by = "code") %>%
  mutate(region = replace_na(region, "Sub-Saharan Africa"))

write_csv(summary_df, file.path(data_path, "postponement_summary.csv"))

# =============================================================================
# FIGURES
# =============================================================================

region_colors <- c(
  "East Asia & Pacific"    = "#1f78b4",
  "Europe"                 = "#e31a1c",
  "Latin America"          = "#33a02c",
  "North America"          = "#ff7f00",
  "South Asia"             = "#6a3d9a",
  "Middle East & N. Africa"= "#b15928",
  "Sub-Saharan Africa"     = "#a6a6a6"
)

# ---- Figure 1: Actual vs. Counterfactual TFR for selected countries ----------
illustrative <- c("Chile", "Spain", "South Korea", "Sweden", "Israel",
                  "United States", "Japan", "Brazil")

fig1_data <- tfr_ts %>%
  filter(country %in% illustrative) %>%
  pivot_longer(c(tfr_actual, tfr_counterfactual),
               names_to = "series", values_to = "tfr") %>%
  mutate(series = recode(series,
    tfr_actual         = "Actual TFR",
    tfr_counterfactual = "Counterfactual TFR\n(1951 cohort timing)"))

fig1 <- ggplot(fig1_data, aes(x = year, y = tfr, color = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~country, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Actual TFR" = "#1f4e79",
                                "Counterfactual TFR\n(1951 cohort timing)" = "#c00000")) +
  scale_linetype_manual(values = c("Actual TFR" = "solid",
                                   "Counterfactual TFR\n(1951 cohort timing)" = "dashed")) +
  labs(
    title    = "Actual and Counterfactual TFR by Year",
    subtitle = "Counterfactual holds the 1951 birth cohort's age distribution of fertility constant",
    x = "Year", y = "Total Fertility Rate",
    color = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        panel.grid.minor = element_blank())

ggsave(file.path(data_path, "fig1_actual_vs_cf.pdf"),
       fig1, width = 14, height = 8)
ggsave(file.path(data_path, "fig1_actual_vs_cf.png"),
       fig1, width = 14, height = 8, dpi = 300)

message("Figure 1 saved.")

# ---- Figure 2: Cohort age profiles (births by birth year and age) -----------
# Replicates the core chart from the Twitter thread for selected countries
cohort_fig_countries <- c("Israel", "Chile", "South Korea", "Spain")

cohort_fig_data <- asfr_cohort %>%
  filter(country %in% cohort_fig_countries,
         age >= 22, age <= 42,
         birth_year >= 1910) %>%
  mutate(age_label = paste0("Age ", age - 2, " to ", age + 2))

fig2 <- ggplot(cohort_fig_data,
               aes(x = birth_year, y = tfr_contrib,
                   color = age_label, group = age_label)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~country, ncol = 2, scales = "free_y") +
  labs(
    title    = "Births by Birth Year and Age Group",
    subtitle = "Each line shows the contribution to TFR from women of that age group, by birth cohort",
    x = "Birth year of mother", y = "TFR contribution (births per woman)",
    color = "Age group"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey90"),
        panel.grid.minor = element_blank())

ggsave(file.path(data_path, "fig2_cohort_age_profiles.pdf"),
       fig2, width = 12, height = 8)
ggsave(file.path(data_path, "fig2_cohort_age_profiles.png"),
       fig2, width = 12, height = 8, dpi = 300)

message("Figure 2 saved.")

# ---- Figure 3: Cross-country scatter ----------------------------------------
scatter_data <- summary_df %>%
  filter(!is.na(postponement_share),
         postponement_share > -0.5,
         postponement_share < 2,
         delta_tfr > 0.1)

fig3 <- ggplot(scatter_data,
               aes(x = delta_tfr, y = postponement_share,
                   color = region, label = country)) +
  geom_hline(yintercept = c(0, 0.5, 1), linetype = "dotted", color = "grey60") +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              linetype = "dashed", linewidth = 0.7, alpha = 0.15,
              inherit.aes = FALSE,
              aes(x = delta_tfr, y = postponement_share)) +
  geom_text_repel(size = 2.8, max.overlaps = 20, segment.color = "grey60") +
  scale_color_manual(values = region_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Postponement Contribution to TFR Decline",
    subtitle = "Share of total TFR decline (early 1970s to 2015–21) attributable to changes in birth timing",
    x = "Total TFR decline (early 1970s to 2015–21)",
    y = "Share of decline due to postponement",
    color = "Region"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave(file.path(data_path, "fig3_postponement_scatter.pdf"),
       fig3, width = 12, height = 8)
ggsave(file.path(data_path, "fig3_postponement_scatter.png"),
       fig3, width = 12, height = 8, dpi = 300)

message("Figure 3 saved.")

# ---- Figure 4: Regional average timing effect over time ---------------------
tfr_ts_region <- tfr_ts %>%
  left_join(region_lookup, by = "code") %>%
  mutate(region = replace_na(region, "Sub-Saharan Africa")) %>%
  filter(region %in% c("East Asia & Pacific", "Europe", "Latin America",
                        "North America", "Sub-Saharan Africa")) %>%
  group_by(region, year) %>%
  summarise(timing_effect = mean(timing_effect, na.rm = TRUE), .groups = "drop")

fig4 <- ggplot(tfr_ts_region, aes(x = year, y = timing_effect,
                                   color = region, group = region)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = region_colors) +
  labs(
    title    = "Regional Average Timing Effect on TFR Over Time",
    subtitle = "Timing effect = Actual TFR − Counterfactual TFR; negative = postponement depresses TFR",
    x = "Year", y = "Timing effect (TFR points)",
    color = "Region"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave(file.path(data_path, "fig4_timing_effect_regions.pdf"),
       fig4, width = 10, height = 6)
ggsave(file.path(data_path, "fig4_timing_effect_regions.png"),
       fig4, width = 10, height = 6, dpi = 300)

message("Figure 4 saved.")

# ---- Figure 5: Bar chart — top and bottom countries by postponement share ---
bar_data <- summary_df %>%
  filter(!is.na(postponement_share),
         delta_tfr > 0.1,
         postponement_share > -0.5,
         postponement_share < 2) %>%
  arrange(postponement_share) %>%
  mutate(rank_asc = row_number(),
         rank_desc = n() - row_number() + 1) %>%
  filter(rank_asc <= 15 | rank_desc <= 15) %>%
  mutate(group = if_else(rank_desc <= 15, "High postponement share", "Low postponement share"))

fig5 <- ggplot(bar_data,
               aes(x = reorder(country, postponement_share),
                   y = postponement_share, fill = group)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("High postponement share" = "#1f4e79",
                                "Low postponement share"  = "#c00000")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Countries by Share of TFR Decline Attributed to Postponement",
    subtitle = "Top and bottom 15 countries (minimum TFR decline of 0.1)",
    x = NULL, y = "Postponement share of TFR decline", fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(data_path, "fig5_postponement_bar.pdf"),
       fig5, width = 10, height = 8)
ggsave(file.path(data_path, "fig5_postponement_bar.png"),
       fig5, width = 10, height = 8, dpi = 300)

message("Figure 5 saved.")

# ---- Figure 6: Israel decomposition (to match Twitter chart) ----------------
israel_ts <- tfr_ts %>% filter(country == "Israel")

fig6 <- ggplot(israel_ts, aes(x = year)) +
  geom_line(aes(y = tfr_actual,         color = "Actual TFR"),         linewidth = 1) +
  geom_line(aes(y = tfr_counterfactual, color = "Counterfactual TFR"), linewidth = 1,
            linetype = "dashed") +
  scale_color_manual(values = c("Actual TFR" = "#1f4e79",
                                 "Counterfactual TFR" = "#c00000")) +
  labs(
    title    = "Actual and Counterfactual TFR by Year: Israel",
    subtitle = "Counterfactual uses the 1951 birth cohort's age distribution of fertility",
    x = "Year", y = "Total Fertility Rate", color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(data_path, "fig6_israel_decomposition.pdf"),
       fig6, width = 9, height = 6)
ggsave(file.path(data_path, "fig6_israel_decomposition.png"),
       fig6, width = 9, height = 6, dpi = 300)

message("Figure 6 (Israel) saved.")

# =============================================================================
# REGRESSION ANALYSIS
# =============================================================================

# Load mean age at childbearing data
mac_raw <- read_csv(file.path(data_path, "period-average-age-of-mothers.csv"),
                    show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  rename(country = entity,
         mean_age = 4)   # 4th column is the MAC variable

mac_early <- mac_raw %>%
  filter(year >= 1970, year <= 1977) %>%
  group_by(country, code) %>%
  summarise(mean_age_early = mean(mean_age, na.rm = TRUE), .groups = "drop")

mac_recent <- mac_raw %>%
  filter(year >= 2015, year <= 2021) %>%
  group_by(country, code) %>%
  summarise(mean_age_recent = mean(mean_age, na.rm = TRUE), .groups = "drop")

reg_data <- summary_df %>%
  left_join(mac_early,  by = c("country", "code")) %>%
  left_join(mac_recent, by = c("country", "code")) %>%
  mutate(delta_mean_age = mean_age_recent - mean_age_early) %>%
  filter(!is.na(postponement_share), !is.na(delta_mean_age),
         delta_tfr > 0.1,
         postponement_share > -0.5, postponement_share < 2)

# OLS models
m1 <- lm(postponement_share ~ delta_mean_age,              data = reg_data)
m2 <- lm(postponement_share ~ delta_mean_age + mean_age_early + delta_tfr,
                                                             data = reg_data)
m3 <- lm(postponement_share ~ delta_mean_age + mean_age_early + delta_tfr + region,
                                                             data = reg_data)

# Print summaries
cat("\n=== Model 1 ===\n"); print(summary(m1))
cat("\n=== Model 2 ===\n"); print(summary(m2))
cat("\n=== Model 3 (with region FE) ===\n"); print(summary(m3))

# Save regression output as text
sink(file.path(data_path, "regression_output.txt"))
cat("=== Model 1: Postponement share ~ Change in mean age ===\n")
print(summary(m1))
cat("\n=== Model 2: + Early mean age + Total TFR decline ===\n")
print(summary(m2))
cat("\n=== Model 3: + Region fixed effects ===\n")
print(summary(m3))
sink()

message("Regression output saved.")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n--- Summary statistics (countries with TFR decline > 0.1) ---\n")
summary_df %>%
  filter(delta_tfr > 0.1, !is.na(postponement_share)) %>%
  select(tfr_early, tfr_recent, delta_tfr, ctfr_early, ctfr_recent,
         delta_ctfr, timing_contrib, postponement_share) %>%
  summary() %>%
  print()

cat("\n--- Top 20 countries by postponement share ---\n")
summary_df %>%
  filter(delta_tfr > 0.05, !is.na(postponement_share)) %>%
  arrange(desc(postponement_share)) %>%
  select(country, region, tfr_early, tfr_recent, delta_tfr,
         timing_contrib, postponement_share) %>%
  slice_head(n = 20) %>%
  print(n = 20)

cat("\n--- Bottom 20 countries by postponement share ---\n")
summary_df %>%
  filter(delta_tfr > 0.05, !is.na(postponement_share)) %>%
  arrange(postponement_share) %>%
  select(country, region, tfr_early, tfr_recent, delta_tfr,
         timing_contrib, postponement_share) %>%
  slice_head(n = 20) %>%
  print(n = 20)

message("\nAll done. Output files written to: ", data_path)
