# This script analyses the changes in the STLs over time (calendar year) by biome.

# For each biome
#   Analyse disturbance trends
#   Histogram of data availability (to be moved to other script)
#   Fit quantile regression on data after applying filters up to bad-slope (with and without interaction)
#   Fit quantile regression for each QMD bin and for each filter level.
#   Bootstrap quantile regression

# Load packages ----------------------------------------------------------------
# library(renv)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
# library(rFIA)
# library(lme4)
# library(lmerTest)
# library(ggeffects)
# library(effects)
# library(sjPlot)
# library(measurements)
library(lqmm)
# library(ggforce)
# library(MuMIn)
# library(DescTools)
library(here)
library(viridis)
library(purrr)
library(rsample)
library(cowplot)

# Load functions ---------------------------------------------------------------
source(here("R/functions.R"))
source(here("R/plot_stl_bybiome.R"))
source(here("R/identify_disturbed_plots.R"))
source(here("R/get_breaks.R"))
source(here("R/plot_lqmm_bybiome.R"))
source(here("R/calc_lqmm_byqmdbin.R"))
source(here("R/create_table_latex.R"))
source(here("R/process_cite_lines.R"))
source(here("R/wrap_fit_lqmm.R"))
source(here("R/calc_percent_change.R"))
source(here("R/plot_disturbed.R"))
source(here("R/fit_model.R"))
source(here("R/extract_coef.R"))

## Define workflow -------------
analyse_biome <- function(df, biome_name){

  # Analyse disturbance trend (before applying disturbed-filter)
  breaks <- get_breaks(df$year)
  gg_fdisturbed <- plot_disturbed(df, biome_name, breaks)

  # Define filter stages --------------
  list_df_filtered <- list(

    # all data
    "raw" = df,

    # no disturbance-affecte plots
    "no_disturbance" = df |>
      filter(ndisturbed == 0),

    # no data from years with shifted QMD distribution
    "no_badqmd" = df |>
      filter(ndisturbed == 0, !badqmdbin), 

    # no data from plots with outlying self-thinning slope
    "no_badslope" = df |>
      filter(ndisturbed == 0, !badqmdbin, !badslope),

    # no data from plots without management history and management <30 years prior to first census
    "no_mgmt_30" = df |>
      filter(
        ndisturbed == 0, !badqmdbin, !badslope,

        # recorded history and management > 30 years prior to first census OR pristine
        (management_cat == 1 & management_since_census1_yrs >= 30) | management_cat == 2
      ),

    # no data from plots without management history and management <100 years prior to first census
    "no_mgmt_100" = df |>
      filter(
        ndisturbed == 0, !badqmdbin, !badslope,

        # recorded history and management > 30 years prior to first census OR pristine
        (management_cat == 1 & management_since_census1_yrs >= 100) | management_cat == 2
      ),

    # only old-growth
    "primary" = df |> 
      filter(
        ndisturbed == 0, !badqmdbin, !badslope,

        # Pristine/primary/old-growth/protected
        management_cat == 2
      )
  )

  # Analyse sensitivity of fit ---------------
  # (year coefficient) subject filter levels
  df_filtereffects <- imap_dfr(list_df_filtered, \(df, step_name) {
    df |>
      ungroup() |> 
      nest() |>
      mutate(
        n = map_int(data, nrow),

        # store both model + coefficients
        res = map(data, \(d) {
          model <- fit_model(d, lqmm = TRUE)
          coefs <- extract_coef(model, lqmm = TRUE)

          list(
            model = model,
            coefs = coefs
          )
        })
      ) |>
      mutate(
        fit_lqmm = purrr::map(res, "model"),
        coef_lqmm = purrr::map(res, "coefs")
      ) |>     # creates columns: model, coefs
      tidyr::unnest(coef_lqmm) |>         # expand coefficient table
      dplyr::select(-data) |>
      mutate(step = step_name)
  })

  # Create plot --------------
  # of coefficient 'year' for all filter levels
  gg_coef_filters <- df_filtereffects |>
    mutate(
      step = fct_relevel(
        step,
        "primary",
        "no_mgmt_100",
        "no_mgmt_30",
        "no_badslope",
        "no_badqmd",
        "no_disturbance",
        "raw",
      )
    ) |> 
    ggplot(aes(x = step, y = estimate)) +
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
      hjust = 1.1,
      size = 3
    ) + 
    theme_classic() +
    ylim(-0.5, 0.5) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # Create classic LQMM plot -------------
  # for filter level no_badslope
  fit_lqmm <- df_filtereffects |> 
      filter(step == "no_badslope") |> 
      pull(fit_lqmm)
  
  gg_lqmm <- plot_lqmm_bybiome(
    list_df_filtered$no_badslope,
    fit_lqmm[[1]]
  )

  # LQMM by QMD bin ----------------
  # for filter no_badslope
  df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(list_df_filtered$no_badslope)
  gg_lqmm_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df)

  return(
    list(
      gg_fdisturbed = gg_fdisturbed,
      df_filtereffects = df_filtereffects,
      gg_lqmm = gg_lqmm,
      gg_coef_filters = gg_coef_filters,
      gg_lqmm_byqmdbin = gg_lqmm_byqmdbin
    )
  )

}

## Run workflow ---------------------
### Read and process data ------------
df <- read_rds(here("data/inputs/df_unm_withfilters.rds"))

### Apply workflow by biome -------------
df_biome_analysis <- df |> 
  group_by(biome_major) |> 
  nest() |> 

  # avoid tropical/subtropical coniferous
  filter(!is.na(biome_major)) |> 
    
  # # xxx test
  # ungroup() |> 
  # tail(1) |> 
  
  # apply workflow by biome
  mutate(out = purrr::map2(data, biome_major, ~analyse_biome(df = .x, biome_name = .y)))

### Save results -------------
# (includes fitted models and ggplot objects for each major biome)
# warning: large file, approx 6 GB
write_rds(
  df_biome_analysis,
  file = here("data/df_biome_analysis.rds")
)
 
# ## Plots -------------------------
# # construct panel
# left_panel <- cowplot::plot_grid(
#   gg_lqmm,
#   gg_lqmm_byqmdbin,
#   ncol = 1,
#   rel_heights = c(1, 0.6),
#   align = "v",
#   # labels = c("",  "g"),
#   label_y = 1.1
# )

# panel <- cowplot::plot_grid(
#   left_panel,
#   gg_coef_filters,
#   rel_widths = c(1, 0.5)
# )


# XXXXXXXXXXXXXXXXXX

# OLD CODE BELOW ----------------------------------------------------------

# # load data
# data_unm <- readRDS(here("data/inputs/df_unm.rds"))

# # optionally subset
# do_subset_primary <- FALSE  # <- manually adjust here

# if (do_subset_primary){
#   suffix_subset <- "_SUBSET"

#   table_s1 <- read_csv(here("data/table_s1_exported.csv"))
#   vec_datasets_subset <- table_s1 |>
#     filter(`Can be considered primary/old growth forest?` == "YES") |>
#     pull(Dataset)

#   data_unm <- data_unm |>
#     filter(dataset %in% vec_datasets_subset)

# } else {
#   suffix_subset <- ""
# }


# ### Plot disturbed plots -------------------------------------------------------
# breaks <- get_breaks(data_unm_biome$year)

# gg_fdisturbed_biome1 <- plot_disturbed(
#   data_unm_biome,
#   "Tropical & Subtropical Moist Broadleaf Forests"
#   )

# ### Remove disturbed plots -----------------------------------------------------
# data_unm_biome_including_disturbed <- data_unm_biome
# data_unm_biome <- data_unm_biome |>
#   filter(ndisturbed == 0)

# write_rds(
#   data_unm_biome,
#   file = here(paste0("data/data_unm_undist_biome1", suffix_subset, ".rds"))
# )

# ### LQMM fit -------------------------------------------------------------------
# set.seed(123)
# fit_lqmm <- lqmm(
#   logDensity ~ logQMD_sc + year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.9, # c(0.75, 0.90),
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(
#     LP_max_iter = 500, # inner loop iterations
#     LP_tol_ll   = 1e-05, # inner loop tolerance
#     startQR     = TRUE
#   )
# )
# # summary(fit_lqmm)

# write_rds(fit_lqmm, file = here(paste0("data/outputs/fit_lqmm_biome1", suffix_subset, ".rds")))

# ### LQMM fit with interaction --------------------------------------------------
# set.seed(123)
# fit_lqmm_int <- lqmm(
#   logDensity ~ logQMD_sc * year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.9, # c(0.75, 0.90),
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(
#     LP_max_iter = 500, # inner loop iterations
#     LP_tol_ll   = 1e-05, # inner loop tolerance
#     startQR     = TRUE
#   )
# )
# # summary(fit_lqmm_int)

# write_rds(fit_lqmm_int, file = here(paste0("data/outputs/fit_lqmm_int_biome1", suffix_subset, ".rds")))

# #### STL shift -----------------------------------------------------------------
# # Opt 2: predicted logDensity at two years differing by one year (in scaled units)
# percent_change <- calc_percent_change(data_unm_biome, fit_lqmm)

# #### Bootstrapping LQMM fit ----------------------------------------------------
# # boot_data <- rsample::bootstraps(
# #   data_unm_biome %>%
# #     group_by(plotID),
# #   times = 500,
# #   apparent = FALSE
# # )
# #
# # # Apply model to each bootstrap sample
# # boot_results <- boot_data %>%
# #   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
# #   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
# #   unnest(coefs) |>
# #   dplyr::select(-splits)
# #
# # write_rds(boot_results, file = here(paste0("data/boot_results_biome1", suffix_subset, ".rds")))
# #
# # df_boot <- boot_results |>
# #   filter(term == "year") |>
# #   mutate(biome = "Tropical & Subtropical Moist Broadleaf Forests")
# #
# # boot_results <- read_rds(file = here(paste0("data/boot_results_biome1", suffix_subset, ".rds")))

# ### Plot STL from LQMM ---------------------------------------------------------
# gg_lqmm_biome1 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm,
#   name = bquote(bold("a") ~ ~"Tropical Moist Broadleaf Forests")
# )
# gg_lqmm_biome1

# ### Plot STL-interaction from LQMM ---------------------------------------------
# gg_lqmm_int_biome1 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm_int,
#   name = bquote(bold("a") ~ ~"Tropical Moist Broadleaf Forests")
# )
# gg_lqmm_int_biome1

# ### Within QMD bins ------------------------------------------------------------
# # Test whether upward shift of 90% quantile is significant within logQMD-bins
# # returns data frame with pval indicating significance level of a positive
# # effect of year.
# df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
# df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
#   data_unm_biome_including_disturbed,
#   breaks = df_lqmm_byqmdbin$breaks
# )

# # Build the plot to access internal structure
# gg_lqmm_biome1_byqmdbin <- plot_lqmm_byqmdbin(
#   df_lqmm_byqmdbin$df,
#   df_lqmm_byqmdbin_including_disturbed$df
#   )

# gg_lqmm_biome1_both <- cowplot::plot_grid(
#   gg_lqmm_biome1,
#   gg_lqmm_biome1_byqmdbin,
#   ncol = 1,
#   rel_heights = c(1, 0.6),
#   align = "v",
#   # labels = c("",  "g"),
#   label_y = 1.1
# )
# gg_lqmm_biome1_both


# ### Plot disturbed plots -------------------------------------------------------
# breaks <- get_breaks(data_unm_biome$year)

# gg_fdisturbed_biome2 <- plot_disturbed(
#   data_unm_biome,
#   "Tropical Dry Broadleaf Forests"
# )

# gg_fdisturbed_biome2

# ### Remove disturbed plots -----------------------------------------------------
# data_unm_biome_including_disturbed <- data_unm_biome
# data_unm_biome <- data_unm_biome |>
#   filter(ndisturbed == 0)

# write_rds(
#   data_unm_biome,
#   file = here(paste0("data/data_unm_undist_biome2", suffix_subset, ".rds"))
# )

# ### LQMM fit -------------------------------------------------------------------
# set.seed(123)
# fit_lqmm <- lqmm(
#   logDensity ~ logQMD_sc + year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(
#     LP_max_iter = 500, # inner loop iterations
#     LP_tol_ll   = 1e-05, # inner loop tolerance
#     startQR     = TRUE
#   )
# )
# # summary(fit_lqmm)

# write_rds(fit_lqmm, file = here(paste0("data/outputs/fit_lqmm_biome2", suffix_subset, ".rds")))

# ### LQMM interaction fit -------------------------------------------------------
# set.seed(123)
# fit_lqmm_int <- lqmm(
#   logDensity ~ logQMD_sc * year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(
#     LP_max_iter = 500, # inner loop iterations
#     LP_tol_ll   = 1e-05, # inner loop tolerance
#     startQR     = TRUE
#   )
# )
# # summary(fit_lqmm_int)

# write_rds(fit_lqmm_int, file = here(paste0("data/outputs/fit_lqmm_int_biome2", suffix_subset, ".rds")))

# #### Bootstrapping LQMM fit -----------------------------------------------------
# # boot_data <- rsample::bootstraps(
# #   data_unm_biome %>%
# #     group_by(plotID),
# #   times = 500,
# #   apparent = FALSE
# # )
# #
# # # Apply model to each bootstrap sample
# # boot_results <- boot_data %>%
# #   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
# #   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
# #   unnest(coefs) |>
# #   dplyr::select(-splits)
# #
# # write_rds(boot_results, file = here(paste0("data/boot_results_biome2", suffix_subset, ".rds")))
# #
# # df_boot <- bind_rows(
# #   df_boot,
# #   boot_results |>
# #     filter(term == "year") |>
# #     mutate(biome = "Tropical Dry Broadleaf Forests")
# # )

# ### Plot STL from LQMM ---------------------------------------------------------
# gg_lqmm_biome2 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm,
#   name = bquote(bold("b") ~ ~"Tropical Dry Broadleaf Forests")
# )
# gg_lqmm_biome2

# ### Plot STL-interaction from LQMM ---------------------------------------------
# gg_lqmm_int_biome2 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm_int,
#   name = bquote(bold("b") ~ ~"Tropical Dry Broadleaf Forests")
# )
# gg_lqmm_int_biome2

# ### Within QMD bins ----------------------------------------
# # Test whether upward shift of 90% quantile is significant within logQMD-bins
# # returns data frame with pval indicating significance level of a positive
# # effect of year.
# df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
# df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
#   data_unm_biome_including_disturbed,
#   breaks = df_lqmm_byqmdbin$breaks
# )

# # Build the plot to access internal structure
# gg_lqmm_biome2_byqmdbin <- plot_lqmm_byqmdbin(
#   df_lqmm_byqmdbin$df,
#   df_lqmm_byqmdbin_including_disturbed$df
#   )

# gg_lqmm_biome2_both <- cowplot::plot_grid(
#   gg_lqmm_biome2,
#   gg_lqmm_biome2_byqmdbin,
#   ncol = 1,
#   rel_heights = c(1, 0.6),
#   align = "v",
#   # labels = c("",  "h"),
#   label_y = 1.1
# )
# gg_lqmm_biome2_both


# ### Plot disturbed plots -------------------------------------------------------
# breaks <- get_breaks(data_unm_biome$year)

# gg_fdisturbed_biome4 <- plot_disturbed(
#   data_unm_biome,
#   "Temperate Broadleaf & Mixed Forests"
# )

# gg_fdisturbed_biome4

# ### Remove disturbed plots -----------------------------------------------------
# data_unm_biome_including_disturbed <- data_unm_biome
# data_unm_biome <- data_unm_biome |>
#   filter(ndisturbed == 0)

# write_rds(
#   data_unm_biome,
#   file = here(paste0("data/data_unm_undist_biome4", suffix_subset, ".rds"))
#   )

# ### LQMM fit -------------------------------------------------------------------
# set.seed(123)
# fit_lqmm <- lqmm(
#   logDensity ~ logQMD_sc + year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(startQR = TRUE)
# )
# # summary(fit_lqmm)

# write_rds(fit_lqmm, file = here(paste0("data/outputs/fit_lqmm_biome4", suffix_subset, ".rds")))

# ### LQMM with interaction fit -------------------------------------------------------------------
# set.seed(123)
# fit_lqmm_int <- lqmm(
#   logDensity ~ logQMD_sc * year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(startQR = TRUE)
# )
# # summary(fit_lqmm_int)

# write_rds(fit_lqmm_int, file = here(paste0("data/outputs/fit_lqmm_int_biome4", suffix_subset, ".rds")))

# #### Bootstrapping LQMM fit -----------------------------------------------------
# # boot_data <- rsample::bootstraps(
# #   data_unm_biome %>%
# #     group_by(plotID),
# #   times = 500,
# #   apparent = FALSE
# # )
# #
# # # Apply model to each bootstrap sample
# # boot_results <- boot_data %>%
# #   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
# #   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
# #   unnest(coefs) |>
# #   dplyr::select(-splits)
# #
# # write_rds(boot_results, file = here(paste0("data/boot_results_biome4", suffix_subset, ".rds")))
# #
# # df_boot <- bind_rows(
# #   df_boot,
# #   boot_results |>
# #     filter(term == "year") |>
# #     mutate(biome = "Temperate Broadleaf & Mixed Forests")
# # )

# ### Plot STL from LQMM ---------------------------------------------------------
# gg_lqmm_biome4 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm,
#   name = bquote(bold("c") ~ ~"Temperate Broadleaf & Mixed Forests")
# )
# gg_lqmm_biome4

# gg_lqmm_biome4 + 
#   labs(title = NULL, subtitle = NULL)

# ggsave(here("fig/gg_lqmm_biome4.pdf"), width = 4, height = 3)

# ### Plot STL from LQMM with interaction ----------------------------------------
# gg_lqmm_int_biome4 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm_int,
#   name = bquote(bold("c") ~ ~"Temperate Broadleaf & Mixed Forests")
# )
# gg_lqmm_int_biome4

# ### Within QMD bins ----------------------------------------
# # Test whether upward shift of 90% quantile is significant within logQMD-bins
# # returns data frame with pval indicating significance level of a positive
# # effect of year.
# df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
# df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
#   data_unm_biome_including_disturbed,
#   breaks = df_lqmm_byqmdbin$breaks
# )

# # Build the plot to access internal structure
# gg_lqmm_biome4_byqmdbin <- plot_lqmm_byqmdbin(
#   df_lqmm_byqmdbin$df,
#   df_lqmm_byqmdbin_including_disturbed$df
#   )

# gg_lqmm_biome4_both <- cowplot::plot_grid(
#   gg_lqmm_biome4,
#   gg_lqmm_biome4_byqmdbin,
#   ncol = 1,
#   rel_heights = c(1, 0.6),
#   # labels = c("",  "i"),
#   align = "v",
#   label_y = 1.1
# )
# gg_lqmm_biome4_both


# ### Plot disturbed plots -------------------------------------------------------
# breaks <- get_breaks(data_unm_biome$year)

# gg_fdisturbed_biome5 <- plot_disturbed(
#   data_unm_biome,
#   "Temperate Conifer Forests"
# )

# gg_fdisturbed_biome5

# ### Remove disturbed plots -----------------------------------------------------
# data_unm_biome_including_disturbed <- data_unm_biome
# data_unm_biome <- data_unm_biome |>
#   filter(ndisturbed == 0)

# write_rds(
#   data_unm_biome,
#   file = here(paste0("data/data_unm_undist_biome5", suffix_subset, ".rds"))
# )

# ### LQMM fit -------------------------------------------------------------------
# set.seed(123)
# fit_lqmm <- lqmm(
#   logDensity ~ logQMD_sc + year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(startQR = TRUE)
# )
# # summary(fit_lqmm)

# write_rds(fit_lqmm, file = here(paste0("data/outputs/fit_lqmm_biome5", suffix_subset, ".rds")))

# ### LQMM-interaction fit -------------------------------------------------------
# set.seed(123)
# fit_lqmm_int <- lqmm(
#   logDensity ~ logQMD_sc * year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(startQR = TRUE)
# )
# # summary(fit_lqmm_int)

# write_rds(fit_lqmm_int, file = here(paste0("data/outputs/fit_lqmm_int_biome5", suffix_subset, ".rds")))

# #### Bootstrapping LQMM fit -----------------------------------------------------
# # boot_data <- rsample::bootstraps(
# #   data_unm_biome %>%
# #     group_by(plotID),
# #   times = 500,
# #   apparent = FALSE
# # )
# #
# # # Apply model to each bootstrap sample
# # boot_results <- boot_data %>%
# #   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
# #   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
# #   unnest(coefs) |>
# #   dplyr::select(-splits)
# #
# # write_rds(boot_results, file = here(paste0("data/boot_results_biome5", suffix_subset, ".rds")))
# #
# # df_boot <- bind_rows(
# #   df_boot,
# #   boot_results |>
# #     filter(term == "year") |>
# #     mutate(biome = "Temperate Conifer Forests")
# # )

# ### Plot STL from LQMM ---------------------------------------------------------
# gg_lqmm_biome5 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm,
#   name = bquote(bold("d") ~ ~"Temperate Conifer Forest")
# )
# gg_lqmm_biome5

# ### Plot STL-interaction from LQMM ---------------------------------------------
# gg_lqmm_int_biome5 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm_int,
#   name = bquote(bold("d") ~ ~"Temperate Conifer Forest")
# )
# gg_lqmm_int_biome5

# ### Within QMD bins ----------------------------------------
# # Test whether upward shift of 90% quantile is significant within logQMD-bins
# # returns data frame with pval indicating significance level of a positive
# # effect of year.
# df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
# df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
#   data_unm_biome_including_disturbed,
#   breaks = df_lqmm_byqmdbin$breaks
# )

# # Build the plot to access internal structure
# gg_lqmm_biome5_byqmdbin <- plot_lqmm_byqmdbin(
#   df_lqmm_byqmdbin$df,
#   df_lqmm_byqmdbin_including_disturbed$df
#   )

# gg_lqmm_biome5_both <- cowplot::plot_grid(
#   gg_lqmm_biome5,
#   gg_lqmm_biome5_byqmdbin,
#   ncol = 1,
#   rel_heights = c(1, 0.6),
#   # labels = c("",  "j"),
#   align = "v",
#   label_y = 1.1
# )
# gg_lqmm_biome5_both


# ### Plot disturbed plots -------------------------------------------------------
# breaks <- get_breaks(data_unm_biome$year)

# gg_fdisturbed_biome6 <- plot_disturbed(
#   data_unm_biome,
#   "Boreal Forests/Taiga"
# )

# gg_fdisturbed_biome6

# ### Remove disturbed plots -----------------------------------------------------
# data_unm_biome_including_disturbed <- data_unm_biome
# data_unm_biome <- data_unm_biome |>
#   filter(ndisturbed == 0)

# write_rds(
#   data_unm_biome,
#   file = here(paste0("data/data_unm_undist_biome6", suffix_subset, ".rds"))
# )

# ### LQMM fit -------------------------------------------------------------------
# set.seed(123)
# fit_lqmm <- lqmm(
#   logDensity ~ logQMD_sc + year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(
#     LP_max_iter = 500, # increase max iterations
#     LP_tol_ll = 1e-4, # relax tolerance slightly (default is 1e-5)
#     startQR = TRUE # good to keep this TRUE
#   )
# )
# # summary(fit_lqmm)

# write_rds(fit_lqmm, file = here(paste0("data/outputs/fit_lqmm_biome6", suffix_subset, ".rds")))

# ### LQMM-interaction fit -------------------------------------------------------
# set.seed(123)
# fit_lqmm_int <- lqmm(
#   logDensity ~ logQMD_sc * year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(
#     LP_max_iter = 500, # increase max iterations
#     LP_tol_ll = 1e-4, # relax tolerance slightly (default is 1e-5)
#     startQR = TRUE # good to keep this TRUE
#   )
# )
# # summary(fit_lqmm)

# write_rds(fit_lqmm_int, file = here(paste0("data/outputs/fit_lqmm_int_biome6", suffix_subset, ".rds")))

# #### Bootstrapping LQMM fit -----------------------------------------------------
# # boot_data <- rsample::bootstraps(
# #   data_unm_biome %>%
# #     group_by(plotID),
# #   times = 500,
# #   apparent = FALSE
# # )
# #
# # # Apply model to each bootstrap sample
# # boot_results <- boot_data %>%
# #   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
# #   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
# #   unnest(coefs) |>
# #   dplyr::select(-splits)
# #
# # write_rds(boot_results, file = here(paste0("data/boot_results_biome6", suffix_subset, ".rds")))
# #
# # df_boot <- bind_rows(
# #   df_boot,
# #   boot_results |>
# #     filter(term == "year") |>
# #     mutate(biome = "Boreal Forests/Taiga")
# # )

# ### Plot STL from LQMM ---------------------------------------------------------
# gg_lqmm_biome6 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm,
#   name = bquote(bold("e") ~ ~"Boreal Forests/Taiga")
# )
# gg_lqmm_biome6

# ### Plot STL-interaction from LQMM ---------------------------------------------
# gg_lqmm_int_biome6 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm_int,
#   name = bquote(bold("e") ~ ~"Boreal Forests/Taiga")
# )
# gg_lqmm_int_biome6

# ### Within QMD bins ----------------------------------------
# # Test whether upward shift of 90% quantile is significant within logQMD-bins
# # returns data frame with pval indicating significance level of a positive
# # effect of year.
# df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
# df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
#   data_unm_biome_including_disturbed,
#   breaks = df_lqmm_byqmdbin$breaks
# )

# # Build the plot to access internal structure
# gg_lqmm_biome6_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df, df_lqmm_byqmdbin_including_disturbed$df)

# gg_lqmm_biome6_both <- cowplot::plot_grid(
#   gg_lqmm_biome6,
#   gg_lqmm_biome6_byqmdbin,
#   ncol = 1,
#   rel_heights = c(1, 0.6),
#   # labels = c("",  "k"),
#   align = "v",
#   label_y = 1.1
# )
# gg_lqmm_biome6_both


# ### Plot disturbed plots -------------------------------------------------------
# breaks <- get_breaks(data_unm_biome$year)

# gg_fdisturbed_biome12 <- plot_disturbed(
#   data_unm_biome,
#   "Mediterranean Forests"
# )

# gg_fdisturbed_biome12

# ### Remove disturbed plots -----------------------------------------------------
# data_unm_biome_including_disturbed <- data_unm_biome
# data_unm_biome <- data_unm_biome |>
#   filter(ndisturbed == 0)

# write_rds(
#   data_unm_biome,
#   file = here(paste0("data/data_unm_undist_biome12", suffix_subset, ".rds"))
# )

# ### LQMM fit -------------------------------------------------------------------
# set.seed(123)
# fit_lqmm <- lqmm(
#   logDensity ~ logQMD_sc + year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(startQR = TRUE)
# )
# # summary(fit_lqmm)

# write_rds(fit_lqmm, file = here(paste0("data/outputs/fit_lqmm_biome12", suffix_subset, ".rds")))

# ### LQMM-interaction fit -------------------------------------------------------
# set.seed(123)
# fit_lqmm_int <- lqmm(
#   logDensity ~ logQMD_sc * year_sc,
#   random = ~1,
#   group = plotID,
#   tau = 0.90,
#   data = data_unm_biome,
#   type = "normal",
#   control = lqmmControl(startQR = TRUE)
# )
# # summary(fit_lqmm_int)

# write_rds(fit_lqmm_int, file = here(paste0("data/outputs/fit_lqmm_int_biome12", suffix_subset, ".rds")))

# #### Bootstrapping LQMM fit -----------------------------------------------------
# # boot_data <- rsample::bootstraps(
# #   data_unm_biome %>%
# #     group_by(plotID),
# #   times = 500,
# #   apparent = FALSE
# # )
# #
# # # Apply model to each bootstrap sample
# # boot_results <- boot_data %>%
# #   mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
# #   filter(!map_lgl(coefs, is.null)) %>% # drop failed fits
# #   unnest(coefs) |>
# #   dplyr::select(-splits)
# #
# # write_rds(boot_results, file = here(paste0("data/boot_results_biome12", suffix_subset, ".rds")))
# #
# # df_boot <- bind_rows(
# #   df_boot,
# #   boot_results |>
# #     filter(term == "year") |>
# #     mutate(biome = "Mediterranean Forests")
# # )

# ### Plot STL from LQMM ---------------------------------------------------------
# gg_lqmm_biome12 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm,
#   name = bquote(bold("f") ~ ~"Mediterranean Forests")
# )
# gg_lqmm_biome12

# ### Plot STL-interaction from LQMM ---------------------------------------------
# gg_lqmm_int_biome12 <- plot_lqmm_bybiome(
#   data_unm_biome,
#   fit_lqmm_int,
#   name = bquote(bold("f") ~ ~"Mediterranean Forests")
# )
# gg_lqmm_int_biome12

# ### Within QMD bins ----------------------------------------
# # Test whether upward shift of 90% quantile is significant within logQMD-bins
# # returns data frame with pval indicating significance level of a positive
# # effect of year.
# df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
# df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
#   data_unm_biome_including_disturbed,
#   breaks = df_lqmm_byqmdbin$breaks
# )

# # Build the plot to access internal structure
# gg_lqmm_biome12_byqmdbin <- plot_lqmm_byqmdbin(df_lqmm_byqmdbin$df, df_lqmm_byqmdbin_including_disturbed$df)

# gg_lqmm_biome12_both <- cowplot::plot_grid(
#   gg_lqmm_biome12,
#   gg_lqmm_biome12_byqmdbin,
#   ncol = 1,
#   rel_heights = c(1, 0.6),
#   # labels = c("",  "l"),
#   align = "v",
#   label_y = 1.1
# )
# gg_lqmm_biome12_both

# Publication figures and tables  ----------------------------------------------

## Figure 1 --------------------------------------------------------------------
### STL and dots combined ------------------------------------------------------
legend <- get_legend(
  gg_lqmm_biome1 +
    theme(legend.position = "right")
)

# Arrange the 9 plots in a 3x3 grid
fig1_lqmm <- cowplot::plot_grid(
  gg_lqmm_biome1_both,
  gg_lqmm_biome2_both,
  gg_lqmm_biome4_both,
  gg_lqmm_biome5_both,
  gg_lqmm_biome6_both,
  gg_lqmm_biome12_both,
  ncol = 3
)

fig1_lqmm

# Combine grid and legend
fig1_lqmm <- cowplot::plot_grid(
  fig1_lqmm,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm.pdf"),
  plot = fig1_lqmm,
  width = 11,
  height = 12
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm.png"),
  plot = fig1_lqmm,
  width = 11,
  height = 12
)

### Fig 1 only STL ------------------
fig1_lqmm_onlystl <- cowplot::plot_grid(
  gg_lqmm_biome1,
  gg_lqmm_biome2,
  gg_lqmm_biome4,
  gg_lqmm_biome5,
  gg_lqmm_biome6,
  gg_lqmm_biome12,
  ncol = 3
)

fig1_lqmm_onlystl <- cowplot::plot_grid(
  fig1_lqmm_onlystl,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm_onlystl.pdf"),
  plot = fig1_lqmm_onlystl,
  width = 11,
  height = 8
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm_onlystl.png"),
  plot = fig1_lqmm_onlystl,
  width = 11,
  height = 8
)

### Fig 1 ALTERNATIVE only STL with interactions  ------------------
fig1_lqmm_int_onlystl <- cowplot::plot_grid(
  gg_lqmm_int_biome1,
  gg_lqmm_int_biome2,
  gg_lqmm_int_biome4,
  gg_lqmm_int_biome5,
  gg_lqmm_int_biome6,
  gg_lqmm_int_biome12,
  ncol = 3
)

fig1_lqmm_int_onlystl <- cowplot::plot_grid(
  fig1_lqmm_int_onlystl,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm_int_onlystl.pdf"),
  plot = fig1_lqmm_int_onlystl,
  width = 11,
  height = 8
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm_int_onlystl.png"),
  plot = fig1_lqmm_int_onlystl,
  width = 11,
  height = 8
)


### Fig 1 only dots ------------------
fig1_lqmm_onlydots <- cowplot::plot_grid(
  gg_lqmm_biome1_byqmdbin,
  gg_lqmm_biome2_byqmdbin,
  gg_lqmm_biome4_byqmdbin,
  gg_lqmm_biome5_byqmdbin,
  gg_lqmm_biome6_byqmdbin,
  gg_lqmm_biome12_byqmdbin,
  ncol = 3,
  labels = letters[1:6]
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm_onlydots.pdf"),
  plot = fig1_lqmm_onlydots,
  width = 9,
  height = 5
)

ggsave(
  filename = here("manuscript/figures/fig1_lqmm_onlystl.png"),
  plot = fig1_lqmm_onlydots,
  width = 9,
  height = 5
)



### Number: years covered -------------
df_tmp <- read_rds(here("data/data_unm_undist_biome1", suffix_subset, ".rds"))) |>
  bind_rows(
    read_rds(here("data/data_unm_undist_biome2", suffix_subset, ".rds")))
  ) |>
  bind_rows(
    read_rds(here("data/data_unm_undist_biome4", suffix_subset, ".rds")))
  ) |>
  bind_rows(
    read_rds(here("data/data_unm_undist_biome5", suffix_subset, ".rds")))
  ) |>
  bind_rows(
    read_rds(here("data/data_unm_undist_biome6", suffix_subset, ".rds")))
  ) |>
  bind_rows(
    read_rds(here("data/data_unm_undist_biome12", suffix_subset, ".rds")))
  )

df_tmp |>
  ggplot(aes(x = year, y = after_stat(density))) +
  geom_histogram() +
  facet_wrap(~biome)

# total years covered
min(df_tmp$year)
max(df_tmp$year)

## SI Figure: Length of repeated observations ----------------------------------
df_len <- data_unm |>
  group_by(plotID, biome, biomeID) |>
  summarise(start = min(year), end = max(year)) |>
  mutate(len = end - start)

# numbers for paper
df_len |>
  ungroup() |>
  group_by(biomeID, biome) |>
  summarise(
    len_median = median(len),
    len_mean = mean(len)
  )

df_len |>
  ungroup() |>
  summarise(
    len_median = median(len),
    len_mean = mean(len)
  )

df_len |>
  ggplot(aes(x = len, color = biome, fill = biome)) +
  geom_density(adjust = 3, alpha = 0.5) +
  scale_fill_manual(
    values = c(
      "Boreal Forests/Taiga"                                       = "dodgerblue4",
      "Mediterranean Forests, Woodlands & Scrub"                   = "orangered3",
      "Temperate Broadleaf & Mixed Forests"                        = "darkgreen",
      "Temperate Conifer Forests"                                  = "lightseagreen",
      "Tropical & Subtropical Dry Broadleaf Forests"               = "goldenrod4",
      "Tropical & Subtropical Moist Broadleaf Forests"             = "springgreen3"
    ),
    na.value = NA,
    breaks = ~ .x[!is.na(.x)],
    name = ""
  ) +
  scale_color_manual(
    values = c(
      "Boreal Forests/Taiga"                                       = "dodgerblue4",
      "Mediterranean Forests, Woodlands & Scrub"                   = "orangered3",
      "Temperate Broadleaf & Mixed Forests"                        = "darkgreen",
      "Temperate Conifer Forests"                                  = "lightseagreen",
      "Tropical & Subtropical Dry Broadleaf Forests"               = "goldenrod4",
      "Tropical & Subtropical Moist Broadleaf Forests"             = "springgreen3"
    ),
    na.value = NA,
    breaks = ~ .x[!is.na(.x)],
    name = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = "Length (years)",
    y = "Density"
  )

ggsave(
  filename = here("manuscript/figures/distribution_length", suffix_subset, ".pdf")),
  width = 10,
  height = 5
)

## SI Figure: Bootstrapped percent change of N per year ------------------------
# write_rds(df_boot, file = here(paste0("data/df_boot", suffix_subset, ".rds")))
df_boot <- read_rds(here("data/df_boot", suffix_subset, ".rds")))

df_boot |>
  mutate(percent_change = 100*(exp(estimate) - 1)) |>
  ggplot(aes(x = percent_change, group = biome, color = biome, fill = biome)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(
    values = c(
      "Boreal Forests/Taiga"                                       = "dodgerblue4",
      "Mediterranean Forests, Woodlands & Scrub"                   = "orangered3",
      "Temperate Broadleaf & Mixed Forests"                        = "darkgreen",
      "Temperate Conifer Forests"                                  = "lightseagreen",
      "Tropical & Subtropical Dry Broadleaf Forests"               = "goldenrod4",
      "Tropical & Subtropical Moist Broadleaf Forests"             = "springgreen3"
    ),
    na.value = NA,
    breaks = ~ .x[!is.na(.x)],
    name = ""
  ) +
  scale_color_manual(
    values = c(
      "Boreal Forests/Taiga"                                       = "dodgerblue4",
      "Mediterranean Forests, Woodlands & Scrub"                   = "orangered3",
      "Temperate Broadleaf & Mixed Forests"                        = "darkgreen",
      "Temperate Conifer Forests"                                  = "lightseagreen",
      "Tropical & Subtropical Dry Broadleaf Forests"               = "goldenrod4",
      "Tropical & Subtropical Moist Broadleaf Forests"             = "springgreen3"
    ),
    na.value = NA,
    breaks = ~ .x[!is.na(.x)],
    name = ""
  ) +
  theme_classic() +
  labs(
    x = expression(paste("Change in density (%"^{-yr}, ")")),
    y = "Density",
    color = "",
    fill = ""
    ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = here("manuscript/figures/distribution_percent_change", suffix_subset, ".pdf")),
  width = 8,
  height = 4
)

# table of bootstrapped estimates for coefficient and percent change
summary_stats <- df_boot |>
  mutate(percent_change = 100*(exp(estimate) - 1)) |>
  ungroup() |>
  group_by(biome) %>%
  summarise(
    estimate_mean = mean(estimate),
    estimate_sd = sd(estimate),
    estimate_ci_low = quantile(estimate, 0.025),
    estimate_ci_high = quantile(estimate, 0.975),

    percent_change_mean = mean(percent_change),
    percent_change_sd = sd(percent_change),
    percent_change_ci_low = quantile(percent_change, 0.025),
    percent_change_ci_high = quantile(percent_change, 0.975),

    .groups = "drop"
  )

write_rds(summary_stats, file = here(paste0("data/summary_stats.csv"))

create_table_latex(
  summary_stats |>
    select(
      Biome = biome,
      Mean = percent_change_mean,
      SD = percent_change_sd
      ),
    caption = "Percentage change of forest stand density.",
    filn = here("manuscript/tables/table_percentage_change", suffix_subset, ".tex"))
    # align = c("p{0.1cm}", "p{5cm}", "p{7cm}")
    )

## SI Figure: Disturbed plots --------------------------------------------------
cowplot::plot_grid(
  gg_fdisturbed_biome1,
  gg_fdisturbed_biome2,
  gg_fdisturbed_biome4,
  gg_fdisturbed_biome5,
  gg_fdisturbed_biome6,
  gg_fdisturbed_biome12,
  ncol = 2,
  labels = letters
)

ggsave(
  filename = here("manuscript/figures/fdisturbed", suffix_subset, ".pdf")),
  width = 9,
  height = 9
)

## Table: datasets
df_datasets <- data_unm |>
  group_by(dataset) |>
  summarise(count = n()) |>
  arrange(-count) |>
  rename(Dataset = dataset, N = count) |>
  mutate(Description = "", Reference = "")

create_table_latex(
  df_datasets,
  caption = "Constituent forest dataset sizes and descriptions.",
  filn = here("manuscript/tables/datasets", suffix_subset, ".tex"))
  # align = c("p{0.1cm}", "p{5cm}", "p{7cm}")
)
