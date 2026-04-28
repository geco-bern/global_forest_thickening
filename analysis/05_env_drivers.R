# Fixed effects plot -----------------------------------------------------------
# This script analyses the main environmental drivers affecting the STL changes

## Load packages ---------------------------------------------------------------
library(readr)
library(dplyr)
library(tibble)
library(lubridate)
library(lme4)
library(lmerTest)
library(here)
library(broom)
library(kableExtra)
library(modelsummary)
library(broom.mixed)
library(forcats)
library(ggplot2)

## Load data -----------
# note: always use (use_slopefilter = TRUE). Results are consistent for MI and
# Ndep against alternative filter options, but not for other factors.
df <- read_rds(here("data/inputs/df_unm_withfilters.rds")) |> 

  # # Re-define all bnp plots as management_cat == 2 (primary), see email Rupert Seidl, 16.04.2026
  # mutate(management_cat = ifelse(dataset == "bnp", 2, management_cat)) |>
  
  # # filter level: no disturbed, no QMD shifts, no outlier slope
  # filter(ndisturbed == 0) #, !badqmdbin, !badslope)

  # use only "primary" for global model
  filter(
    ndisturbed == 0,
    !badqmdbin,
    !badslope,

    # Pristine/primary/old-growth/protected
    management_cat == 2
  )

# Downsample for even distribution.
# Perform a stratified random sampling, where strata are defined by all
# combinations of three levels in tavg, ai, and ndep. The number of
# sampled points per stratum should be equal to the total number of
# points in the sparsest stratum.

# # Create 3-level strata for each variable
# df_strat <- df %>%
#   ungroup() |>
#   select(plotID, tavg, ai, ndep) |>
#   distinct(plotID, .keep_all = TRUE) |>
#   mutate(
#     tavg_strata = ntile(tavg, 3),
#     ai_strata = ntile(ai, 3),
#     ndep_strata = ntile(ndep, 3)
#   )

# # Find the size of the smallest stratum
# min_n <- df_strat %>%
#   count(tavg_strata, ai_strata, ndep_strata) %>%
#   summarise(min_n = quantile(n, 0.25)) %>%
#   pull(min_n)

# # Perform stratified sampling
# set.seed(123)
# df_sampled <- df_strat %>%
#   group_by(tavg_strata, ai_strata, ndep_strata) %>%
#   slice_sample(n = min_n) %>%
#   ungroup()

# # Check balance
# df_sampled %>%
#   count(tavg_strata, ai_strata, ndep_strata) |> 
#   View()

# # subset original data frame to sampled plots
# df <- df |> 
#   filter(plotID %in% df_sampled$plotID)

# # xxx try
# df <- read_rds(here("data/data_unm_slopefilter.rds"))

# LMM with lmer() --------------------------------------------------------------
## Fit model -------------------------------------------------------------------
# with all environmental factors and their interaction with time as predictors
### complete model --------
mod_lmer_env_complete <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    scale(year) * scale(ORGC) +
    scale(year) * scale(PBR) +
    scale(year) * scale(CNrt) +
    (1 | plotID) + (1 | species),
  data = df,
  na.action = "na.exclude"
)

summary(mod_lmer_env_complete)

write_rds(mod_lmer_env_complete, file = here("data/mod_lmer_env_complete.rds"))

### no soil variables ----------
mod_lmer_env_nosoil <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    (1 | plotID) + (1 | species),
  data = df,
  na.action = "na.exclude"
)

summary(mod_lmer_env_nosoil)

write_rds(mod_lmer_env_nosoil, file = here("data/mod_lmer_env_nosoil.rds"))


### no PBR ---------
# fit model without PBR - this turns out as the best one
mod_lmer_env_nopbr <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    scale(year) * scale(ORGC) +
    # scale(year) * scale(PBR) + # excluded because neither fixed nor interactive effect is significant
    scale(year) * scale(CNrt) +
    (1 | plotID) + (1 | species),
  data = df,
  na.action = "na.exclude"
)

write_rds(mod_lmer_env_nopbr, file = here("data/mod_lmer_env_nopbr.rds"))

### no PBR and ORGC -----------
# fit model without PBR and ORGC
mod_lmer_env_nopbr_noorgc <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    # scale(year) * scale(ORGC) +
    # scale(year) * scale(PBR) + # excluded because neither fixed nor interactive effect is significant
    scale(year) * scale(CNrt) +
    (1 | plotID) + (1 | species),
  data = df,
  na.action = "na.exclude"
)

write_rds(mod_lmer_env_nopbr_noorgc, file = here("data/mod_lmer_env_nopbr_noorgc.rds"))

### no PBR and CNrt ------------
# fit model without PBR and CNrt
mod_lmer_env_nopbr_nocn <- lmer(
  logDensity ~ scale(logQMD) +
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    scale(year) * scale(ORGC) +
    # scale(year) * scale(PBR) + # excluded because neither fixed nor interactive effect is significant
    # scale(year) * scale(CNrt) +
    (1 | plotID) + (1 | species),
  data = df,
  na.action = "na.exclude"
)

write_rds(mod_lmer_env_nopbr_nocn, file = here("data/mod_lmer_env_nopbr_nocn.rds"))

### complete model with interactions on QMD ----------------
mod_lmer_env_complete_interactions <- lmer(
  logDensity ~
    scale(logQMD) +

    # year × environment interactions
    scale(year) * scale(tavg) +
    scale(year) * scale(ai) +
    scale(year) * scale(ndep) +
    # scale(year) * scale(ORGC) +
    # scale(year) * scale(PBR) +
    # scale(year) * scale(CNrt) +

    # logQMD × environment interactions
    scale(logQMD) * scale(tavg) +
    scale(logQMD) * scale(ai) +
    scale(logQMD) * scale(ndep) +
    # scale(logQMD) * scale(ORGC) +
    # scale(logQMD) * scale(PBR) +
    # scale(logQMD) * scale(CNrt) +

    (1 | plotID) + (1 | species),
  data = df,
  na.action = "na.exclude"
)

summary(mod_lmer_env_complete_interactions)

## Write summary of model fits to latex table
options("modelsummary_format_numeric_latex" = "plain")

# no interactions between slope (QMD) and environmental factors
modelsummary(
  list(
    "Complete" = mod_lmer_env_complete,
    "No soil" = mod_lmer_env_nosoil,
    "Complete with interactions" = mod_lmer_env_complete_interactions
  ),
  output = here(paste0("manuscript/tables/mods_env.tex")),
  # estimate  = "p.value",
  estimate = "{estimate}{stars}",
  # estimate  = "{estimate} [{conf.low}, {conf.high}], {p.value}",
  title = "Regression Results",
  statistic = "conf.int",
  coef_omit = "Intercept"
)

# with interactions between slope (QMD) and environmental factors
modelsummary(
  list(
    "Complete interactions" = mod_lmer_env_complete_interactions
  ),
  output = here(paste0("manuscript/tables/mods_env_int_", lab_filter, ".tex")),
  # estimate  = "p.value",
  estimate = "{estimate}{stars}",
  # estimate  = "{estimate} [{conf.low}, {conf.high}], {p.value}",
  title = "Regression Results",
  statistic = "conf.int",
  coef_omit = "Intercept"
)

# Warning: this is not dynamic
mod_lmer_env <- mod_lmer_env_nosoil  # based on slope filter

## Visualise fixed effects -----------------------------------------------------
tl <- attr(terms(mod_lmer_env), "term.labels")
main_effects <- tl[ !grepl(":", tl) ]   # remove interaction terms
len_main_effects <- length(main_effects)

df_coef_plot <- broom.mixed::tidy(
  mod_lmer_env,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.90
) |>
  mutate(
    eff = ifelse(
      row_number() %in% 1:(len_main_effects + 1),
      "Main effect",
      "Interaction terms"
    ),
    eff = as_factor(eff)
  ) |>
  rename(var = term) |>
  mutate(
    varnew = ifelse(var == "scale(year)", "year", var),
    varnew = ifelse(var == "scale(tavg)", "MAT", varnew),
    varnew = ifelse(var == "scale(ai)", "MI", varnew),
    varnew = ifelse(var == "scale(ndep)", "Ndep", varnew),
    varnew = ifelse(var == "scale(ORGC)", "ORGC", varnew),
    varnew = ifelse(var == "scale(PBR)", "PBR", varnew),
    varnew = ifelse(var == "scale(CNrt)", "C:N", varnew),
    varnew = ifelse(var == "scale(year):scale(tavg)", "MAT", varnew),
    varnew = ifelse(var == "scale(year):scale(ai)", "MI", varnew),
    varnew = ifelse(var == "scale(year):scale(ndep)", "Ndep", varnew),
    varnew = ifelse(var == "scale(year):scale(ORGC)", "ORGC", varnew),
    varnew = ifelse(var == "scale(year):scale(PBR)", "PBR", varnew),
    varnew = ifelse(var == "scale(year):scale(CNrt)", "C:N", varnew),
    varnew = as_factor(varnew),
    varnew = fct_relevel(varnew, c("C:N", "Ndep", "ORGC", "MI", "MAT"))
  ) |>
  filter(
    varnew == "MAT" |
      varnew == "MI" |
      varnew == "Ndep" |
      varnew == "PBR" |
      varnew == "ORGC" |
      varnew == "C:N"
  ) |>
  mutate(
    varnew = forcats::fct_relevel(
      varnew,
      "PBR",
      "ORGC",
      "C:N",
      "Ndep",
      "MI",
      "MAT"
    )
  )

## Save model object and coefficients table ------------------------------------
saveRDS(mod_lmer_env, file = here("data/mod_lmer_env.rds"))
saveRDS(df_coef_plot, file = here("data/df_coef_plot.rds"))

## Plot effects ----------------------------------------------------------------
df_coef_plot<- read_rds(here("data/df_coef_plot.rds"))

fig2a <- df_coef_plot |>
  filter(eff == "Main effect") |>
  ggplot() +
  geom_point(
    aes(
      x = varnew,
      y = estimate
    ),
    size = 2
  ) +
  geom_errorbar(
    aes(
      x = varnew,
      ymin = conf.low,
      ymax = conf.high
    ),
    width = 0
  ) +
  labs(x = "", y = "Coefficient (scaled)", title = "Main effects") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_classic() +
  scale_x_discrete(
    labels = c(
      "MAT" = "Mean\n temperature",
      "MI" = "Moisture\n index",
      "ORGC" = "Organic\n carbon",
      "Ndep" = "Nitrogen\n deposition",
      "PBR" = "Phosporus\n availability"
    )
  ) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
  ) +
  coord_flip()

fig2b <- df_coef_plot |>
  filter(eff == "Interaction terms") |>
  ggplot() +
  geom_point(
    aes(
      x = varnew,
      y = estimate
    ),
    size = 2
  ) +
  geom_errorbar(
    aes(
      x = varnew,
      ymin = conf.low,
      ymax = conf.high
    ),
    width = 0
  ) +
  labs(x = "", y = "Coefficient (scaled)", title = expression(paste("Interactions with ", italic("year")))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_classic() +
  scale_x_discrete(
    labels = NULL
  ) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
  ) +
  coord_flip()

cowplot::plot_grid(
  fig2a,
  fig2b,
  labels = letters[1:2],
  rel_widths = c(1, 0.83)
)

ggsave(
  here(paste0("manuscript/figures/coefs_env.png")),
  width = 6,
  height = 3,
  dpi = 500
  )

ggsave(
  here(paste0("manuscript/figures/coefs_env.pdf")),
  width = 6,
  height = 3
)

# # LMM with lqmm() --------------------------------------------------------------
#
# # data_unm <- readRDS(here("data/inputs/data_unm.rds")) |>
# #   tidyr::drop_na(
# #     year, logQMD, tavg, ai, ndep, ORGC, PBR, CNrt, lai
# #   ) |>
# #   mutate(
# #     year_sc = scale(year),
# #     logQMD_sc = scale(logQMD),
# #     tavg_sc = scale(tavg),
# #     ai_sc = scale(ai),
# #     ndep_sc = scale(ndep),
# #     ORGC_sc = scale(ORGC),
# #     PBR_sc = scale(PBR),
# #     CNrt_sc = scale(CNrt),
# #     lai_sc = scale(lai)
# #   )
# #
# # visdat::vis_miss(data_unm, warn_large_data = FALSE)
# #
# # mod_lqmm_env <- lqmm(
# #   logDensity ~ logQMD_sc +
# #     year_sc * tavg_sc +
# #     year_sc * ai_sc +
# #     year_sc * ndep_sc +
# #     year_sc * ORGC_sc +
# #     year_sc * PBR_sc +
# #     year_sc * CNrt_sc,
# #   random = ~1,
# #   group = plotID,
# #   tau = 0.9,
# #   data = data_unm,
# #   type = "normal",
# #   control = lqmmControl(
# #     LP_max_iter = 5000, # inner loop iterations
# #     LP_tol_ll   = 1e-05, # inner loop tolerance
# #     startQR     = TRUE
# #   )
# # )
#
#
# data_unm <- readRDS(here("data/inputs/data_unm.rds"))
#
# data_unm <- data_unm |>
#   mutate(
#     year_sc = scale(year),
#     logQMD_sc = scale(logQMD),
#     tavg_sc = scale(tavg),
#     ai_sc = scale(ai),
#     ndep_sc = scale(ndep),
#     orgc_sc = scale(ORGC),
#     pbr_sc = scale(PBR)
#   )
#
# ### Identify disturbed plots ---------------------------------------------------
# source(here("R/identify_disturbed_plots.R"))
# source(here("R/get_breaks.R"))
#
# data_unm <- data_unm |>
#   identify_disturbed_plots()
#
# breaks <- get_breaks(data_unm$year)
#
# df_disturbed <- data_unm |>
#   mutate(year_bin = cut(
#     year,
#     breaks = breaks,
#     labels = breaks[1:length(breaks) - 1] + 2.5,
#     include.lowest = TRUE
#   )) |>
#   group_by(year_bin) |>
#   summarise(
#     nplots = length(unique(plotID)),
#     ndisturbed = sum(disturbed, na.rm = TRUE)
#   ) |>
#   filter(nplots > 30) |>
#   mutate(fdisturbed = ndisturbed / nplots) |>
#   mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))
#
# ### Plot disturbed plots -------------------------------------------------------
# gg_fdisturbed <- df_disturbed |>
#   ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
#   geom_point() +
#   geom_smooth(method = "lm", color = "red") +
#   theme_classic() +
#   labs(
#     x = "Year"
#   ) +
#   scale_y_continuous(
#     name = expression(logit(Fraction ~ disturbed)),
#     sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
#   )
#
# gg_fdisturbed
#
# # Additional filter: remove plots with no change in ln(N)
# data_unm <- data_unm |>
#   group_by(plotID) |>
#   mutate(var_logdensity = diff(range(logDensity))) |>
#   filter(var_logdensity > 0.001)
#
# # remove disturbed plots
# data_unm_including_disturbed <- data_unm
# data_unm <- data_unm |>
#   filter(ndisturbed == 0)
#
# ### LQMM fit -------------------------------------------------------------------
# set.seed(123)
# fit_lqmm <- lqmm(
#   logDensity ~ logQMD_sc + year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm,
#   type = "normal",
#   control = lqmmControl(startQR = TRUE)
# )
#
# # with all environmental factors and their interaction with time as predictors
# set.seed(123)
# mod_lqmm_env <- lqmm(
#   logDensity ~ logQMD_sc +
#     year_sc * tavg_sc +
#     year_sc * ai_sc +
#     year_sc * ndep_sc +
#     year_sc * orgc_sc +
#     year_sc * pbr_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.70,
#   data = data_unm |>
#     drop_na(),
#   type = "normal",
#   control = list(
#     LP_max_iter = 5000,
#     LP_tol_ll   = 1e-05,
#     startQR     = TRUE
#   )
# )
#
# summary(mod_lqmm_env)
#
# out <- summary(mod_lqmm_env)
#
# plot_lqmm_bybiome(
#   data_unm,
#   mod_lqmm_env
# )
