# Analyse effects of filters defined in 02_filter_selfthinning.R on the self-thinning slope,
# fitted with quantile regression.

library(tidyverse)
library(lqmm)
library(purrr)
library(forcats)

data_unm <- read_rds(here("data/inputs/data_unm_withfilters.rds")) |>
  mutate(
    year_sc = scale(year),
    logQMD_sc = scale(logQMD)
  ) |> 
  
  # xxx test
  filter(dataset_major %in% c("nfi_spain", "aus_plots"))

# Define filter stages
data_steps <- list(
  "raw" = data_unm,

  "no_disturbance" = data_unm |>
    filter(!disturbed),

  "no_ingrowth" = data_unm |>
    filter(!disturbed, !ingrowth),

  "no_badqmd" = data_unm |>
    filter(!disturbed, !ingrowth, !badqmdbin),

  "no_badslope" = data_unm |>
    filter(!disturbed, !ingrowth, !badqmdbin, !badslope)
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
      LP_max_iter = 500,
      LP_tol_ll   = 1e-05,
      startQR     = TRUE
    )
  )

  return(fit)
}

# Run for all filter stages, separately for each (major) dataset

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

df_filtereffects <- imap_dfr(data_steps, ~{
  .x |>
    group_by(dataset_major) |>
    group_modify(~{
      fit <- fit_model(.x)
      extract_coef(fit)
    }) |>
    mutate(step = .y)
})

# Order steps properly
df_filtereffects <- df_filtereffects |>
  mutate(
    step = fct_relevel(
      step,
      "raw",
      "no_disturbance",
      "no_ingrowth",
      "no_badqmd",
      "no_badslope"
    )
  )

# Add sample size
n_obs <- map_dfr(data_steps, ~tibble(n = nrow(.x)), .id = "step")

df_filtereffects <- df_filtereffects |>
  left_join(n_obs, by = "step")

write_rds(df_filtereffects, here("data/df_filtereffects.rds"))

# Plot
ggplot(df_filtereffects, aes(x = step, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  # geom_text(aes(label = paste0("n=", n)), hjust = -0.9, size = 3) +
  labs(
    x = "Filtering step",
    y = "Coefficient (year)"
  ) +
  facet_wrap(~dataset_major) +
  theme_classic() +
  coord_flip()

