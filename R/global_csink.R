# Functions for the environmental upscaling in 08_globalcsink.R.

# Fit the STL model within each bootstrap sample -------------------------------
# The model terms follow the selection results from 07_env_drivers.R.
fit_stl_byboot <- function(df) {
  vars_to_scale <- c("logQMD", "year", "tavg", "ai", "ndep", "ORGC", "CNrt", "PBR")

  rec <- recipe(
    logDensity ~ logQMD + year + tavg + ai + ndep + ORGC + CNrt + PBR + dataset + plotID + species,
    data = df |>
      drop_na(
        all_of(
          c(
            "logDensity",
            "logQMD",
            "year",
            "tavg",
            "ai",
            "ndep",
            "ORGC",
            "CNrt",
            "PBR",
            "dataset",
            "plotID",
            "species"
          )
        )
      )
  ) |>
    step_center(all_of(vars_to_scale)) |>
    step_scale(all_of(vars_to_scale))

  # fit the recipe (stores the training means and SDs)
  rec_prep <- prep(rec)

  # get the preprocessed training data
  df_scaled <- bake(rec_prep, new_data = NULL)

  # Fit the model using the preprocessing estimated from this training sample.
  model <- lmer(
    logDensity ~ logQMD +
      year * tavg +
      year * ai +
      year * ndep +
      year * ORGC +
      year * CNrt +
      year * PBR +
      (1 | plotID) + (1 | species),
    data = df_scaled
  )

  # Keep the recipe with its fitted model for consistent prediction scaling.
  model_bundle <- list(recipe = rec_prep, model = model)

  return(model_bundle)
}

# Fit the biomass relationship within each bootstrap sample --------------------
fit_biomass_byboot <- function(df) {
  # biomass data is in kg-DM/ha
  df <- df |>
    drop_na(
      all_of(
        c(
          "biomass",
          "NQMD2",
          # "dataset",
          "plotID"
        )
      )
    )
  if (nrow(df) == 0) {
    return(NA)
  } else {
    model <- lmer(
      biomass ~ NQMD2 + 0 + (1 | plotID),
      data = df
    )
    return(model)
  }
}

# Predict density change at fixed QMD ------------------------------------------
predict_dn <- function(bundle_stl, vec_qmd, df) {
  # sample QMD from empirical distribution and add dummy random factor levels
  # (ignored for prediction but required to avoid error)
  df <- df |>
    mutate(
      qmd = sample(vec_qmd, nrow(df), replace = TRUE),
      year = 2000,
      dataset = "dummy_dataset",
      plotID = "dummy_plot",
      species = "dummy_species"
    ) |>
    mutate(logQMD = log(qmd))

  # apply same preprocessing as for the STL data
  new_df_scaled <- bake(bundle_stl$recipe, new_data = df)

  # predict N, given QMD, year, and environmental factors for all forest gridcells globally
  df$logDensity_0 <- predict(
    bundle_stl$model,
    newdata = new_df_scaled,
    re.form = NA
  ) # Predict ignoring random effects

  df <- df |>
    mutate(year = 2001)

  # apply same preprocessing as for the STL data
  new_df_scaled <- bake(bundle_stl$recipe, new_data = df)

  # predict N, given QMD, year, and environmental factors for all forest gridcells globally
  df$logDensity_1 <- predict(
    bundle_stl$model,
    newdata = new_df_scaled,
    re.form = NA
  ) # Predict ignoring random effects

  df <- df |>
    mutate(dn = exp(logDensity_0) * expm1(logDensity_1 - logDensity_0)) |>
    mutate(dnqmd2 = dn * qmd^2) |>
    select(lon_i, lat_i, dn, dnqmd2)

  return(df)
}

# Convert the density change to biomass change ---------------------------------
predict_db <- function(df, model) {
  if (identical(model, NA)) {
    df <- df |>
      mutate(db = NA)
  } else {
    # add dummy random factor levels
    # (ignored for prediction but required to avoid error)
    df <- df |>
      mutate(
        dataset = "dummy_dataset",
        plotID = "dummy_plot"
      )

    # Apply the through-origin biomass model directly to the change in N * QMD^2.
    df$db <- predict(
      model,
      newdata = df |>
        rename(NQMD2 = dnqmd2), # predict one with difference rather than twice with absolute
      re.form = NA # Predict ignoring random effects
    )

    df <- df |>
      select(-dataset, -plotID)
  }

  return(df)
}

# Integrate the C sink over forest-covered area --------------------------------
calc_csink_global <- function(df, df_info) {
  # biomass data is in kg-DM/ha
  df |>
    left_join(
      df_info,
      by = join_by(lon_i, lat_i)
    ) |>

    # convert from kg-DM ha-1 on forest land to gC (integral across gridcell)
    mutate(db_gc_per_gridcell = db * 1e3 * 0.5 * fcf * area_ha) |>

    # sum across all gridcells globally and convert to PgC
    summarise(db_pgc = 1e-15 * sum_csink_available(db_gc_per_gridcell)) |>
    pull(db_pgc)
}

calc_csink_bylat <- function(df, df_info) {
  # biomass data is in kg-DM/ha
  df |>
    left_join(
      df_info,
      by = join_by(lon_i, lat_i)
    ) |>

    # convert from kg-DM ha-1 on forest land to gC (integral across gridcell)
    mutate(db_gc_per_gridcell = db * 1e3 * 0.5 * fcf * area_ha) |>

    # sum across all gridcells by latitudinal band and convert to TgC
    group_by(lat_i) |>
    summarise(db_tgc = 1e-12 * sum_csink_available(db_gc_per_gridcell))
}
# An entirely missing prediction is not a zero sink. Partial sums retain the
# existing available-cell convention; coverage must be assessed separately.
sum_csink_available <- function(x) {
  if (any(is.infinite(x))) {
    stop("Non-finite carbon sink: inspect extrapolated predictions.")
  }

  if (all(is.na(x))) {
    return(NA_real_)
  }

  sum(x, na.rm = TRUE)
}
