# Run from the project root: Rscript tests/test_fia_management.R
suppressPackageStartupMessages(library(dplyr))
source("R/prepare_fia_management.R")

ids <- c("oldgrowth", "known_age", "prior_cut", "unknown", "later_age", "later_oldgrowth",
         "planted_before", "undated_prior_cut", "baseline_year_cut", "older_history",
         "mixed_known", "mixed_unknown", "mixed_oldgrowth", "treated", "undated_treatment",
         "new_stumps", "new_planting", "nonreserve", "unknown_reserve", "unknown_followup",
         "missing_census", "missing_baseline", "mixed_managed", paste0("sentinel", 1:5))
stands <- tibble(plotID = rep(ids, each = 3), year = rep(c(2000, 2005, 2010), length(ids)),
                 density = 100, QMD = 50, management_cat = 2,
                 management_since_census1_yrs = 999)
conditions <- stands |>
  transmute(plotID, INVYR = year, CONDID = 1L, COND_STATUS_CD = 1L,
            RESERVCD = 1L, STDORGCD = 0L, STDAGE = NA_real_,
            STND_COND_CD_PNWRS = NA_integer_, STUMP_CD_PNWRS = "N",
            TRTCD1 = 0L, TRTCD2 = NA_integer_, TRTCD3 = 0L,
            TRTYR1 = NA_real_, TRTYR2 = NA_real_, TRTYR3 = NA_real_)
set_values <- function(id, field, values) {
  conditions[[field]][conditions$plotID == id] <<- values
}
set_values("oldgrowth", "STND_COND_CD_PNWRS", c(7L, NA, NA))
set_values("oldgrowth", "STDAGE", c(300, 305, 310))
set_values("known_age", "STDAGE", c(60, 61, 62))
set_values("prior_cut", "STDAGE", 80)
set_values("prior_cut", "STUMP_CD_PNWRS", "Y")
set_values("prior_cut", "TRTCD1", c(10L, 0L, 0L))
set_values("prior_cut", "TRTYR1", c(1990, NA, NA))
set_values("later_age", "STDAGE", c(NA, 70, 75))
set_values("later_oldgrowth", "STND_COND_CD_PNWRS", c(6L, 7L, 7L))
set_values("planted_before", "STDORGCD", 1L)
set_values("planted_before", "STDAGE", 35)
set_values("undated_prior_cut", "TRTCD1", c(10L, 0L, 0L))
set_values("baseline_year_cut", "TRTCD1", c(10L, 0L, 0L))
set_values("baseline_year_cut", "TRTYR1", c(2000, NA, NA))
set_values("mixed_known", "STDAGE", 60)
set_values("mixed_unknown", "STDAGE", 60)
set_values("mixed_oldgrowth", "STND_COND_CD_PNWRS", c(7L, NA, NA))
set_values("mixed_oldgrowth", "STDAGE", 300)
set_values("treated", "TRTCD2", c(NA, 10L, NA))
set_values("treated", "TRTYR2", c(NA, 2003, NA))
set_values("undated_treatment", "TRTCD1", c(0L, 0L, 10L))
set_values("new_stumps", "STUMP_CD_PNWRS", c("N", NA, "Y"))
set_values("new_planting", "STDORGCD", c(0L, 0L, 1L))
set_values("nonreserve", "RESERVCD", c(1L, 1L, 0L))
set_values("unknown_reserve", "RESERVCD", NA_integer_)
set_values("unknown_followup", "TRTCD1", c(0L, NA, 0L))
for (i in 1:5) set_values(paste0("sentinel", i), "STDAGE", c(0, -999, 998, 999, 9999)[i])
conditions <- conditions |>
  filter(!(plotID == "missing_census" & INVYR == 2005),
         !(plotID == "missing_baseline" & INVYR == 2000))
extra_known <- conditions |> filter(plotID == "mixed_known") |> mutate(CONDID = 2L, STDAGE = 25)
extra_unknown <- conditions |> filter(plotID == "mixed_unknown") |> mutate(CONDID = 2L, STDAGE = NA_real_)
extra_young <- conditions |> filter(plotID == "mixed_oldgrowth") |>
  mutate(CONDID = 2L, STND_COND_CD_PNWRS = 6L, STDAGE = 90)
extra_managed <- conditions |> filter(plotID == "mixed_managed") |>
  mutate(CONDID = 2L, TRTCD1 = if_else(INVYR == 2010, 10L, 0L))
earlier <- conditions |> filter(plotID == "older_history", INVYR == 2000) |>
  mutate(INVYR = 1990, RESERVCD = 0L, TRTCD1 = 10L, TRTYR1 = 1985)
future <- conditions |> filter(plotID == "oldgrowth", INVYR == 2010) |>
  mutate(INVYR = 2020, TRTCD1 = 10L, TRTYR1 = 2019)
conditions <- bind_rows(conditions, extra_known, extra_unknown, extra_young,
                        extra_managed, earlier, future)

result <- prepare_fia_management(stands, conditions)
expected <- tibble(
  plotID = c("oldgrowth", "known_age", "prior_cut", "unknown", "later_age", "later_oldgrowth",
             "planted_before", "undated_prior_cut", "baseline_year_cut", "older_history",
             "mixed_known", "mixed_unknown", "mixed_oldgrowth", paste0("sentinel", 1:5)),
  management_cat = c(2L, 1L, 1L, 3L, 3L, 3L, 1L, 3L, 1L, 1L, 1L, 3L, 1L, rep(3L, 5)),
  management_since_census1_yrs = c(NA, 60, 10, NA, NA, NA, 35, NA, 0, 15, 25, NA, 90, rep(NA, 5))
)
stopifnot(
  identical(result |> distinct(plotID, management_cat, management_since_census1_yrs), expected),
  nrow(result) == nrow(expected) * 3L, all(result$management == 0),
  all(result$density == 100), all(result$QMD == 50),
  identical(result$year, rep(c(2000, 2005, 2010), nrow(expected))),
  identical(result, prepare_fia_management(stands, conditions[nrow(conditions):1, ])),
  identical(result, prepare_fia_management(result, conditions)),
  nrow(prepare_fia_management(stands, conditions[FALSE, ])) == 0L,
  nrow(prepare_fia_management(stands[FALSE, ], conditions)) == 0L
)

# A treatment reported later but explicitly dated before the baseline is prior
# history. Placeholder dates must not be interpreted as real calendar years.
reported_prior <- conditions |> filter(plotID == "unknown") |>
  mutate(TRTCD1 = if_else(INVYR == 2010, 10L, 0L),
         TRTYR1 = if_else(INVYR == 2010, 1995, NA_real_))
prior_result <- prepare_fia_management(filter(stands, plotID == "unknown"), reported_prior)
stopifnot(all(prior_result$management_cat == 1L),
          all(prior_result$management_since_census1_yrs == 5))
invalid_date <- conditions |> filter(plotID == "undated_prior_cut") |> mutate(TRTYR1 = 9999)
stopifnot(all(prepare_fia_management(filter(stands, plotID == "undated_prior_cut"),
                                     invalid_date)$management_cat == 3L))
# Genuine ages greater than 999 are possible; exclude sentinel codes, not a
# blanket age range. Missing pre-census treatment information is allowed.
ancient <- conditions |> filter(plotID == "unknown") |>
  mutate(STDAGE = 1000, TRTCD1 = if_else(INVYR == 2000, NA_integer_, 0L))
stopifnot(all(prepare_fia_management(filter(stands, plotID == "unknown"),
                                     ancient)$management_since_census1_yrs == 1000))
cat("FIA between-census management and three-category history checks passed.\n")
