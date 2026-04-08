# Analyse effects of filters defined in 02_filter_selfthinning.R on the self-thinning slope,
# fitted with quantile regression.
#
# MANAGEMENT INFO:
# df_all contains information of management:
# - management_since_census1_yrs (before called years_last_management)
# - management_cat = 1,2 or 3.
# - management = 0 (indicating no management between census)
# 
# Among management_cat:
# 1. Recorded history
#    The time since the last management before the first census is known and recorded.
#    the number of years since the last management before the first census, or
#    the calendar year of the last management intervention.
# 2. Pristine/primary/old-growth/protected
#    The time since the last management before the first census is unknown. The forest was considered pristine, primary, old-growth, protected, or minimally used at the time of the first census.
# 3. Unrecorded history
#    The time since the last management before the first census is unknown. The forest may have a history of substantial past wood harvesting and/or visible signs of forest use before the first census.

library(tidyverse)
library(lqmm)
library(purrr)
library(forcats)

df_unm <- read_rds(here("data/inputs/df_unm_withfilters.rds")) |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD)
  )
  
  # # xxx test
  # filter(dataset_major %in% c("fia_us", "wytham", "rainfor_plots", "euforia", "uholka"))

# Define filter stages
list_df_filtered <- list(

  # all data
  "raw" = df_unm,

  # no disturbance-affecte plots
  "no_disturbance" = df_unm |>
    filter(!disturbed),

  # no ingrowth-affected plots
  "no_ingrowth" = df_unm |>
    filter(!disturbed, !ingrowth),

  # no data from years with shifted QMD distribution
  "no_badqmd" = df_unm |>
    filter(!disturbed, !ingrowth, !badqmdbin),

  # no data from plots with outlying self-thinning slope
  "no_badslope" = df_unm |>
    filter(!disturbed, !ingrowth, !badqmdbin, !badslope),

  # no data from plots without management history and management <30 years prior to first census
  "no_mgmt_30" = df_unm |>
    filter(
      !disturbed, !ingrowth, !badqmdbin, !badslope, 

      # recorded history and management > 30 years prior to first census OR pristine
      (management_cat == 1 & management_since_census1_yrs >= 30) | management_cat == 2
    ),

  # no data from plots without management history and management <100 years prior to first census
  "no_mgmt_100" = df_unm |>
    filter(
      !disturbed, !ingrowth, !badqmdbin, !badslope, 

      # recorded history and management > 30 years prior to first census OR pristine
      (management_cat == 1 & management_since_census1_yrs >= 100) | management_cat == 2
    ),

  # only old-growth
  "pristine" = df_unm |> 
    filter(
      !disturbed, !ingrowth, !badqmdbin, !badslope, 

      # Pristine/primary/old-growth/protected
      management_cat == 2
    )
)

# Fit model for each filter stage
fit_model <- function(df) {
  set.seed(123)

  fit <- lqmm(
    logDensity ~ logQMD_sc + year_sc,
    random = ~1,
    group = plotID,
    tau = 0.9,
    data = df,
    type = "normal",
    control = lqmmControl(
      LP_max_iter = 200, # xxx TEST: revert back to 500,
      LP_tol_ll   = 1e-05,
      startQR     = TRUE
    )
  )

  return(fit)
}

# Extract coefficient + CI
extract_coef <- function(fit) {
  cf <- summary(fit)$tTable

  tibble(
    term = rownames(cf),
    estimate = cf[, "Value"],
    se = cf[, "Std. Error"],
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  ) |>
    filter(term == "year_sc")
}

# Analyse sensitivity of fit (year coefficient) subject to filters,
# separately for each (major) dataset.
df_filtereffects <- imap_dfr(list_df_filtered, \(df, step_name) {
  df |>
    group_by(dataset_major) |>
    nest() |>
    mutate(
      n = map_int(data, nrow),
      res = map(data, \(d) {
        d |>
          fit_model() |>
          extract_coef()
      })
    ) |>
    dplyr::select(-data) |>
    unnest(res) |>
    mutate(step = step_name)
})

# Order steps properly
df_filtereffects <- df_filtereffects |>
  mutate(
    step = fct_relevel(
      step,
      "pristine",
      "no_mgmt_100",
      "no_mgmt_30",
      "no_badslope",
      "no_badqmd",
      "no_ingrowth",
      "no_disturbance",
      "raw",
    )
  )

write_rds(df_filtereffects, here("data/df_filtereffects.rds"))

# Plot
ggplot(df_filtereffects, aes(x = step, y = estimate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  # geom_text(aes(label = paste0("n=", n)), hjust = -0.9, size = 3) +
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = "Filtering step",
    y = "Coefficient (year)"
  ) +
  # Add observation counts to the right of points
  geom_text(
    aes(label = n), 
    y = Inf,
    hjust = 1.1,    # nudge to the right
    size = 3
  ) +      # small font
  facet_wrap(~dataset_major, scales = "free_y", ncol = 3) +
  theme_classic() +
  coord_flip()

# ggsave(
#   here("fig/filtereffects_TEST.pdf"),
#   width = 12,
#   height = 10
# )

ggsave(
  here("fig/filtereffects.pdf"),
  width = 12,
  height = 15
)
