# Compile stand records for all datasets, before environmental extraction.

# Load packages ----
library(tidyverse)

# Load functions ----
source(here::here("R/functions.R"))
source(here::here("R/extract_plot_environment.R"))
source(here::here("R/plot_logn_logqmd_dataset.R"))
source(here::here("R/deduplicate_parent_plots.R"))

# Collect data ----------------
# Collect and combine data from standard-formatted files.
## NFI Spain ----
df_nfi_spain <- readRDS(here::here("data/inputs/df_nfi_spa_v2026.rds"))

## NFI Sweeden ----
df_nfi_sweeden <- readRDS(here::here("data/inputs/df_nfi_swe_v2026.rds"))

## FIA US ----
df_fia_us <- readRDS(here::here("data/inputs/df_fia_us_v2026.rds"))

## NFI Switzerland ----
df_nfi_switzerland <- readRDS(here::here("data/inputs/df_nfi_swi_v2026.rds"))

## NFI Norway ----
df_nfi_norway <- readRDS(here::here("data/inputs/df_nfi_nor_v2026.rds"))

## EFM Switzerland ----
df_efm_swi <- readRDS(here::here("data/inputs/df_efm_swi_v2026.rds"))

## Uholka plot ----
df_uholka <- readRDS(here::here("data/inputs/df_uholka_v2026.rds"))

## Greece sites ----
df_fep_gre <- readRDS(here::here("data/inputs/df_fep_gre_v2026.rds"))

## France sites ----
df_inrae_lessem <- readRDS(here::here("data/inputs/df_inrae_lessem_v2026.rds"))

## EuFoRia plots ----

### bnp ----
df_bnp <- readRDS(here::here("data/inputs/df_bnp_v2026.rds"))

### czu ----
df_czu <- readRDS(here::here("data/inputs/df_czu_v2026.rds"))

### fvabw ----
df_fvabw <- readRDS(here::here("data/inputs/df_fvabw_v2026.rds"))

### iberbas ----
df_iberbas <- readRDS(here::here("data/inputs/df_iberbas_v2026.rds"))

### unitbv ----
df_unitbv <- readRDS(here::here("data/inputs/df_unitbv_v2026.rds"))

### lwf ----
df_lwf <- readRDS(here::here("data/inputs/df_lwf_v2026.rds"))

### npvbw ----
df_npvbw <- readRDS(here::here("data/inputs/df_npvbw_v2026.rds"))

### nfr ----
df_nfr_swi <- readRDS(here::here("data/inputs/df_nfr_swi_v2026.rds"))

### nwfva ----
df_nwfva <- readRDS(here::here("data/inputs/df_nwfva_v2026.rds"))

### tuzvo ----
df_tuzvo <- readRDS(here::here("data/inputs/df_tuzvo_v2026.rds"))

### ul ----
df_ul <- readRDS(here::here("data/inputs/df_ul_v2026.rds"))

### unito ----
df_unito <- readRDS(here::here("data/inputs/df_unito_v2026.rds"))

### urk ----
df_urk <- readRDS(here::here("data/inputs/df_urk_v2026.rds"))

### wuls ----
df_wuls <- readRDS(here::here("data/inputs/df_wuls_v2026.rds"))

## ForestGEO ----

### Luquillo ----
df_luquillo <- readRDS(here::here("data/inputs/df_luquillo_v2026.rds"))

### BCI ----
df_bci <- readRDS(here::here("data/inputs/df_bci_v2026.rds"))

### SCBI ----
df_scbi <- readRDS(here::here("data/inputs/df_scbi_v2026.rds"))

### Palanam ----
df_palanam <- readRDS(here::here("data/inputs/df_palanam_v2026.rds"))

### SERC ----
df_serc <- readRDS(here::here("data/inputs/df_serc_v2026.rds"))

### Wytham woods ----
df_wytham <- readRDS(here::here("data/inputs/df_wytham_v2026.rds"))

### Pasoh woods ----
df_pasoh <- readRDS(here::here("data/inputs/df_pasoh_v2026.rds"))

### Mudumalai ----
df_mudumalai <- readRDS(here::here("data/inputs/df_mudumalai_v2026.rds"))

## Forest plots ----
df_forestplots <- readRDS(here::here("data/inputs/df_forestplots_v2026.rds"))

## Australian plots ----
df_aus <- readRDS(here::here("data/inputs/df_aus_v2026.rds"))

## RAINFOR plots ----
df_rainfor <- readRDS(here::here("data/inputs/df_rainfor_v2026.rds"))

## Join all -----
df_all <- bind_rows(
  df_nfi_spain,
  df_nfi_sweeden,
  df_fia_us,
  df_nfi_switzerland,
  df_nfi_norway,
  df_efm_swi,
  df_uholka,
  df_fep_gre,
  df_inrae_lessem,
  df_bnp,
  df_czu,
  df_fvabw,
  df_iberbas,
  df_unitbv,
  df_lwf,
  df_npvbw,
  df_nfr_swi,
  df_nwfva,
  df_tuzvo,
  df_ul,
  df_unito,
  df_urk,
  df_wuls,
  # df_luquillo,  # not included because QMD of all plots lie far outside the range of qmd of remaining plots (except Pasoh, which is also excluded)
  # df_bci,  # not included. Comment in `00_raw_data.R`: finally not added to the global dataset due to difficulty to exclude management (communication with  Salomón Aguilar, Research Technician)
  df_scbi,
  # df_palanam,  # not included. None of the 1,401 raw quadrat–census groups spans more than one year. The forest data doesn't exhibit a self-thinning pattern, potentially due to recent disturbances (Typhoons in 2003 and 2010, see https://forestgeo.si.edu/sites/philippines/palanan).
  df_serc,
  df_wytham,
  # df_pasoh,  # not included because QMD of all plots lie far outside the range of qmd of remaining plots (except Luquillo, which is also excluded)
  df_mudumalai,
  df_forestplots,
  df_aus,
  df_rainfor
)

# Harmonise dataset names
df_all <- df_all |>
  mutate(dataset = harmonise_dataset_names(dataset))

# Remove duplicate parent plots before adding dataset prefixes to plot IDs.
# Matching leading codes (e.g. ACJ-01) identify source representations of the
# same parent, including all its subplots. Keep the longest census span;
# ties favour more census years, then the latest census, then dataset name.
# Different codes within 10 km count as duplicates only when parent-level BA
# and density agree within 1% in every shared census year (at least one needed).
# Missing comparisons remain unconfirmed. Subplot values are area-weighted.
parent_plot_dedup <- deduplicate_parent_plots(
  df_all,
  review_distance_km = 10,
  measurement_tolerance_pct = 1
)
df_all <- parent_plot_dedup$data
dedup_audit_dir <- here::here("analysis/diagnostics/parent_plot_duplicates")
dir.create(
  dedup_audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
write_csv(
  parent_plot_dedup$code_matches,
  file.path(dedup_audit_dir, "code_matches.csv")
)
write_csv(
  parent_plot_dedup$geographic_candidates,
  file.path(dedup_audit_dir, "geographic_candidates.csv")
)
write_csv(
  parent_plot_dedup$geographic_year_comparisons,
  file.path(dedup_audit_dir, "geographic_year_comparisons.csv")
)

df_all <- df_all |>
  # unite plotID and dataset to ensure there are not matches form different sources
  rename(plotIDD = plotID) |>
  unite(plotID, dataset, plotIDD, sep = "_", remove = FALSE) |>
  select(-plotIDD) |>

  # Remove repeated census rows within the retained source and subplot.
  distinct(country, plotID, year, .keep_all = TRUE)

# Apply basic filters ------
# Apply filters relevant for all analyses
df_unm <- df_all |>

  # Beni: this is an important step: not to remove this a priori.
  # # filter for min qmd
  # filter(QMD >= 10) |>

  # filter for unmanaged plots (between censi)
  filter(management == 0) |>

  # Filter by minimum 3 censi
  group_by(dataset, plotID) |>
  mutate(n_census_unm = n()) |>
  filter(n_census_unm >= 3) |>

  # Beni: since we don't remove disturbed and bad-slope anymore, also don't remove non-declining densities
  # # remove plots with no change in ln(N)
  # mutate(var_logdensity = diff(range(log(density)))) |>
  # filter(var_logdensity > 0.01) |>

  # clean up
  ungroup() |>
  relocate(n_census_unm, .after = n_census)

# Add biome and environmental covariates ------
# Extract global environmental fields for all plot locations ----
# This section can also be rerun independently using the saved compilation.

# Edit individual paths here if the global fields are stored elsewhere.
environment_paths <- plot_environment_paths()

# GSDE/Shangguan is required on the server. Set FALSE only when deliberately
# omitting phosphorus; PBR will then be NA. LAI is no longer extracted.
include_phosphorus <- TRUE

df_unm <- extract_plot_environment(
  df_unm,
  paths = environment_paths,
  include_phosphorus = include_phosphorus
)

# some biome information is implausibly ascribed from global maps. Overwrite with NA, to later fill with dominant biome from the same dataset.
df_unm <- df_unm |>
  mutate(
    biome = ifelse(
      biomeID %in%
        c(
          3, # Tropical & Subtropical Coniferous Forests
          13, # Deserts & Xeric Shrublands
          9, # Flooded Grasslands & Savannas
          8, # Temperate Grasslands, Savannas & Shrublands
          14, # Mangroves
          7, # Tropical & Subtropical Grasslands, Savannas & Shrublands
          10 # Montane Grasslands & Shrublands
        ),
        NA, 
        biome
      ),
    biomeID = ifelse(
      biomeID %in%
        c(
          3, # Tropical & Subtropical Coniferous Forests
          13, # Deserts & Xeric Shrublands
          9, # Flooded Grasslands & Savannas
          8, # Temperate Grasslands, Savannas & Shrublands
          14, # Mangroves
          7, # Tropical & Subtropical Grasslands, Savannas & Shrublands
          10 # Montane Grasslands & Shrublands
        ),
        NA, 
        biomeID
      )
    )


# Fill missing biome assignments from the most frequent biome in the dataset
# (counting distinct plots), before deriving biome_major from biomeID.
df_unm <- fill_missing_biomes(df_unm)

df_unm |> 
  select(biome, biomeID) |> 
  distinct()

saveRDS(df_unm, file = here::here("data/inputs/df_unm_step_environment.rds"))

# Construct 'major biome' --------
# df_unm <- readRDS(here::here("data/inputs/df_unm_step_environment.rds"))
df_unm <- df_unm |> 

  # Re-group into "major biomes"
  mutate(
    biome_major = ifelse(
      biomeID %in% c(1, 2),
      "Tropical & Subtropical Broadleaf Forests",
      # biomeID == 1,
      # "Tropical & Subtropical Moist Broadleaf Forests",
      # ifelse(
      #   biomeID == 2,
      #   "Tropical & Subtropical Dry Broadleaf Forests",
      ifelse(
        biomeID %in% c(4, 5),
        "Temperate Forests",
        ifelse(
          biomeID == 6,
          "Boreal Forests/Taiga",
          ifelse(
            biomeID == 12,
            "Mediterranean Forests",
            NA
          )
        )
      )
    )
  )

df_unm |>
  select(biome, biomeID, biome_major) |>
  distinct()

# Construct 'major dataset' ----
df_unm <- df_unm |>
  mutate(
    dataset_major = ifelse(
      dataset %in%
        c(
          # "bnp", # too big, keeping separate
          "czu",
          "fvabw",
          "iberbas",
          "unitbv",
          "lwf",
          "nbw",
          "nwfva",
          "tuzvo",
          "ul",
          "unito",
          "urk"
          # "wuls" # too big, keeping separate
        ),
      "euforia",
      ifelse(
        dataset %in%
          c(
            "wytham",
            "uholka",
            "nfr_swi",
            "efm_swi",
            "fep_gre_stand",
            "inrae_lessem"
          ),
        "other_europe",
        dataset
      )
    )
  )

  # # old
# mutate(
#   dataset_major = ifelse(
#     dataset %in% c("bnp", "czu", "fvabw", "iberbas", "unitbv", "lwf", "nbw", "nfr", "nwfva", "tuzvo", "ul", "unito", "urk", "wuls", "uholka"),
#     "euforia",
#     ifelse(
#       dataset %in% c("luquillo", "bci", "scbi", "palanan", "serc", "pasoh", "mudumalai"),
#       "forestgeo",
#       ifelse(
#         dataset == "forestplots",
#         "forestplots",
#         ifelse(
#           dataset == "rainfor",
#           "rainfor",
#           dataset
#         )
#       )
#     )
#   )
# )

# additionally exclude tropical data from fia_us
df_unm <- df_unm |>
  filter(!(biome_major == "Tropical & Subtropical Broadleaf Forests" & dataset == "fia_us"))

# Write file ----
saveRDS(df_unm, file = here::here("data/inputs/df_unm_v2026.rds"))

# Plots by dataset ----
# df_unm <- readRDS(here::here("data/inputs/df_unm_v2026.rds"))
plots_logn_logqmd <- split(df_unm, df_unm$dataset, drop = TRUE) |>
  purrr::map(plot_logn_logqmd_dataset)

# Divide the dataset plots evenly between two multi-panel plots.
n_first <- ceiling(length(plots_logn_logqmd) / 2)
plots_logn_logqmd_1 <- head(plots_logn_logqmd, n_first)
plots_logn_logqmd_2 <- tail(plots_logn_logqmd, length(plots_logn_logqmd) - n_first)

gg_logn_logqmd_1 <- cowplot::plot_grid(
  plotlist = plots_logn_logqmd_1,
  ncol = 3,
  labels = letters[1:length(plots_logn_logqmd_1)]
)

gg_logn_logqmd_2 <- cowplot::plot_grid(
  plotlist = plots_logn_logqmd_2,
  ncol = 3,
  labels = letters[1:length(plots_logn_logqmd_2)]
)

ggsave(
  here::here("manuscript/figures/logn_logqmd_1.pdf"), 
  plot = gg_logn_logqmd_1,
  width = 14,
  height = 20
)

ggsave(
  here::here("manuscript/figures/logn_logqmd_2.pdf"),
  plot = gg_logn_logqmd_2,
  width = 14,
  height = 20
)

df_unm |>
  ggplot(aes(
    x = forcats::fct_reorder(biome_major, QMD, .fun = median, .desc = TRUE),
    y = QMD
  )) +
  geom_boxplot(outliers = FALSE) +
  labs(x = NULL, y = "QMD (cm)") +
  theme_classic() +
  coord_flip()

df_unm |>
  group_by(plotID) |>
  filter(year == first(year)) |>
  ggplot(aes(
    x = tidytext::reorder_within(
      dataset,
      QMD,
      biome_major,
      fun = function(x) -median(x, na.rm = TRUE)
    ),
    y = QMD
  )) +
  geom_boxplot(outliers = FALSE) +
  # geom_density() +
  tidytext::scale_x_reordered() +
  scale_y_log10() +
  facet_wrap(~biome_major, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "QMD (cm)") +
  theme_classic() +
  coord_flip()

# /////////////////////////////////////////////////////////////////////
# Legacy code --------------
# # Next run the plot-extraction section in 00_collect_global_drivers.R,
# # then apply the common analysis filters in 02_create_filters.R.

# # Historical analyses retained for reference; superseded by 02_create_filters.R.
# # They require environmental and log columns that this compilation no longer adds.
# run_legacy_analyses <- FALSE
# if (run_legacy_analyses) {
# # Write unmanaged data ----
# # NEW APPROACH:
# # Management categories 1,2,3 depending on the knowledge of the forest use history.

# # THIS WAS THE INITIAL APPROACH:
# # Select only unmanaged forests. These filters are now applied in 02_create_filters.R
# df_unm <- data_unm_fc(df_all)
# df_unm
# df_unm |>
#   nrow()

# plot_stl(df_unm)
# plot_map(df_unm)

# saveRDS(df_unm, file = here::here("data/inputs/df_unm.rds"))

# # Figures for SI -------------------------


# # Generate filtered by biomes ----
# # The unmanaged data is divided by biomes, and then for each biome we apply the filter of the upper quantile

# # read data written just above
# df_unm <- read_rds(here::here("data/inputs/df_unm.rds"))

# ## Biome 1 ----
# # Tropical & Subtropical Moist Broadleaf Forests
# df_unm_biome1 <- df_unm |>
#   filter(biomeID == 1)
# df_unm_biome1 |>
#   nrow()

# plot_stl(df_unm_biome1)

# # filter data with upper quantile
# df_fil75_biome1 <- df_filter75_fc(df_unm_biome1)
# plot_stl(df_fil75_biome1)
# saveRDS(df_fil75_biome1, file = here::here("data/inputs/df_fil75_biome1.rds"))

# df_fil90_biome1 <- df_filter90_fc(df_unm_biome1)
# plot_stl(df_fil90_biome1)
# saveRDS(df_fil90_biome1, file = here::here("data/inputs/df_fil90_biome1.rds"))

# df_fil55_biome1 <- df_filter55_fc(df_unm_biome1)
# plot_stl(df_fil55_biome1)
# saveRDS(df_fil55_biome1, file = here::here("data/inputs/df_fil55_biome1.rds"))

# ## Biome 2 ----
# # Tropical & Subtropical Dry Broadleaf Forests Forest
# df_unm_biome2 <- df_unm |> 
#   filter(biomeID==2)
# df_unm_biome2 |> 
#   nrow()
# plot_stl(df_unm_biome2)

# # filter data with upper quantile
# df_fil75_biome2 <- df_filter75_fc(df_unm_biome2)
# plot_stl(df_fil75_biome2)
# saveRDS(df_fil75_biome2, file = here::here("data/inputs/df_fil75_biome2.rds"))

# df_fil90_biome2 <- df_filter90_fc(df_unm_biome2)
# plot_stl(df_fil90_biome2)
# saveRDS(df_fil90_biome2, file = here::here("data/inputs/df_fil90_biome2.rds"))

# df_fil55_biome2 <- df_filter55_fc(df_unm_biome2)
# plot_stl(df_fil55_biome2)
# saveRDS(df_fil55_biome2, file = here::here("data/inputs/df_fil55_biome2.rds"))

# ## Biome 3 ----
# # Tropical & Subtropical Coniferous Forests Forest
# # NOTE: There are no enough data

# ## Biome 4 ----
# # Temperate Broadleaf & Mixed Forests Forest
# df_unm_biome4 <- df_unm |>
#   filter(biomeID == 4)
# plot_stl(df_unm_biome4)

# # filter data with upper quantile
# df_fil75_biome4 <- df_filter75_fc(df_unm_biome4)
# plot_stl(df_fil75_biome4)
# saveRDS(df_fil75_biome4, file = here::here("data/inputs/df_fil75_biome4.rds"))

# df_fil90_biome4 <- df_filter90_fc(df_unm_biome4)
# plot_stl(df_fil90_biome4)
# saveRDS(df_fil90_biome4, file = here::here("data/inputs/df_fil90_biome4.rds"))

# df_fil55_biome4 <- df_filter55_fc(df_unm_biome4)
# plot_stl(df_fil55_biome4)
# saveRDS(df_fil55_biome4, file = here::here("data/inputs/df_fil55_biome4.rds"))

# ## Biome 5 ----
# # Temperate Conifer Forests Forest
# df_unm_biome5 <- df_unm |>
#   filter(biomeID == 5)
# plot_stl(df_unm_biome5)

# # filter data with upper quantile
# df_fil75_biome5 <- df_filter75_fc(df_unm_biome5)
# plot_stl(df_fil75_biome5)
# saveRDS(df_fil75_biome5, file = here::here("data/inputs/df_fil75_biome5.rds"))

# df_fil90_biome5 <- df_filter90_fc(df_unm_biome5)
# plot_stl(df_fil90_biome5)
# saveRDS(df_fil90_biome5, file = here::here("data/inputs/df_fil90_biome5.rds"))

# df_fil55_biome5 <- df_filter55_fc(df_unm_biome5)
# plot_stl(df_fil55_biome5)
# saveRDS(df_fil55_biome5, file = here::here("data/inputs/df_fil55_biome5.rds"))

# ## Biome 6 ----
# # Boreal Forests/Taiga Forest
# df_unm_biome6 <- df_unm |>
#   filter(biomeID == 6)
# plot_stl(df_unm_biome6)

# # filter data with upper quantile
# df_fil75_biome6 <- df_filter75_fc(df_unm_biome6)
# plot_stl(df_fil75_biome6)
# saveRDS(df_fil75_biome6, file = here::here("data/inputs/df_fil75_biome6.rds"))

# df_fil90_biome6 <- df_filter90_fc(df_unm_biome6)
# plot_stl(df_fil90_biome6)
# saveRDS(df_fil90_biome6, file = here::here("data/inputs/df_fil90_biome6.rds"))

# df_fil55_biome6 <- df_filter55_fc(df_unm_biome6)
# plot_stl(df_fil55_biome6)
# saveRDS(df_fil55_biome6, file = here::here("data/inputs/df_fil55_biome6.rds"))

# ## Biome 12 ----
# # Mediterranean Forests
# df_unm_biome12 <- df_unm |>
#   filter(biomeID == 12)
# plot_stl(df_unm_biome12)

# # filter data with upper quantile
# df_fil75_biome12 <- df_filter75_fc(df_unm_biome12)
# plot_stl(df_fil75_biome12)
# saveRDS(df_fil75_biome12, file = here::here("data/inputs/df_fil75_biome12.rds"))

# df_fil90_biome12 <- df_filter90_fc(df_unm_biome12)
# plot_stl(df_fil90_biome12)
# saveRDS(df_fil90_biome12, file = here::here("data/inputs/df_fil90_biome126.rds"))

# df_fil55_biome12 <- df_filter55_fc(df_unm_biome12)
# plot_stl(df_fil55_biome12)
# saveRDS(df_fil55_biome12, file = here::here("data/inputs/df_fil55_biome12.rds"))

# # Write filtered data ----
# # join filtered datasets
# df_fil75_biomes <- bind_rows(
#   df_fil75_biome1,
#   df_fil75_biome2,
#   df_fil75_biome4,
#   df_fil75_biome5,
#   df_fil75_biome6,
#   df_fil75_biome12
# )
# saveRDS(df_fil75_biomes, file = here::here("data/inputs/df_fil75_biomes.rds"))

# df_fil90_biomes <- bind_rows(
#   df_fil90_biome1,
#   df_fil90_biome2,
#   df_fil90_biome4,
#   df_fil90_biome5,
#   df_fil90_biome6,
#   df_fil90_biome12
# )
# saveRDS(df_fil90_biomes, file = here::here("data/inputs/df_fil90_biomes.rds"))

# df_fil55_biomes <- bind_rows(
#   df_fil55_biome1,
#   df_fil55_biome2,
#   df_fil55_biome4,
#   df_fil55_biome5,
#   df_fil55_biome6,
#   df_fil55_biome12
# )
# saveRDS(df_fil55_biomes, file = here::here("data/inputs/df_fil55_biomes.rds"))

# # Plot filtered data --------
# df_fil75_biomes <- read_rds(here::here("data/inputs/df_fil75_biomes.rds"))

# df_fil75_biomes |>
#   drop_na(years_since_management) |>
#   distinct(dataset)

# # Check DBH distributions
# ggplot() +
#   geom_density(data = df_fil75_biomes, aes(QMD, ..density.., col = factor(biome), fill = factor(biome)), alpha = 0.7) +
#   theme_classic() +
#   # color palette from okabeito_colors()
#   scale_fill_manual(name = "Biome", values = c(
#     "Tropical & Subtropical Moist Broadleaf Forests" = "#E69F00", 
#     "Tropical & Subtropical Dry Broadleaf Forests" = "#56B4E9",
#     "Temperate Broadleaf & Mixed Forests" = "#009E73",
#     "Temperate Conifer Forests" = "#F5C710",
#     "Boreal Forests/Taiga" = "#0072B2",
#     "Mediterranean Forests, Woodlands & Scrub" = "#D55E00"
#   )) +
#   scale_color_manual(name = "Biome", values = c(
#     "Tropical & Subtropical Moist Broadleaf Forests" = "#E69F00",
#     "Tropical & Subtropical Dry Broadleaf Forests" = "#56B4E9",
#     "Temperate Broadleaf & Mixed Forests" = "#009E73",
#     "Temperate Conifer Forests" = "#F5C710",
#     "Boreal Forests/Taiga" = "#0072B2",
#     "Mediterranean Forests, Woodlands & Scrub" = "#D55E00"
#   )) +
#   theme(legend.position = "bottom") +
#   scale_x_continuous(limits = c(0, 80)) 

# # Table S1 ----
# data <- readRDS(here::here("/data/inputs/df_unm.rds"))
# data <- readRDS(here::here("/data/inputs/df_fil75_biomes.rds"))

# # Unmanaged years
# data |>
#   filter(years_since_management != 999) |>
#   summarise(mean=mean(years_since_management, na.rm=T),
#             min=min(years_since_management, na.rm=T),
#             max=max(years_since_management, na.rm=T))
# # plot size
# data |>
#   summarise(mean=mean(plotsize, na.rm=T),
#             min=min(plotsize, na.rm=T),
#             max=max(plotsize, na.rm=T))

# # Mean year for first and last census per plot
# data %>%
#   group_by(plotID) %>%
#   filter(census == 1) %>%
#   ungroup() |>
#   summarise(mean=mean(year, na.rm=T))

# data %>%
#   group_by(plotID) %>%
#   filter(census == max(census, na.rm = TRUE)) %>%
#   ungroup() |>
#   summarise(mean=mean(year, na.rm=T))

# species <- as.data.frame(sort(unique(df_unm$species)))
# length(unique(df_unm$plotID))
# df_unm
# summary(data)

# df_unm |>
#   filter(plotID == "nwfva_tree_03-025_2")

# # Table at the plot level
# table_s1 <- data %>% group_by(country,dataset) %>% 
#   summarise(n_plots=n_distinct(plotID),
#             min_lon = min(lon),
#             max_lon = max(lon),
#             min_lat = min(lat),
#             max_lat = max(lat),
#             min_year = min(year),
#             max_year = max(year),
#             timespan=max_year-min_year,
#             n_census_mean = as.integer(mean(n_census,na.rm=T)),
#             interval=as.integer(mean(period,na.rm=T)),
#             interval_sd =sd(period,na.rm=T),
#             qmd_mean=mean(QMD,na.rm=T),
#             QMD_sd=sd(QMD,na.rm=T),
#             density_mean=mean(density,na.rm=T),
#             density_sd=sd(density,na.rm=T),
#             biomass_mean=mean(biomass,na.rm=T),
#             biomass_sd=sd(biomass,na.rm=T))
# table_s1
# write.csv(table_s1, file = here::here("/data/outputs/tableS1.csv"))

# # Table S2 ----
# # Table at the biome level
# table_s2 <- data %>% 
#   group_by(biome) %>% 
#   summarise(n_plots=n_distinct(plotID),
#             min_year = min(year),
#             max_year = max(year),
#             timespan=mean(period,na.rm=T),
#             timespan_sd =sd(period,na.rm=T),
#             qmd_mean=mean(QMD,na.rm=T),
#             QMD_sd=sd(QMD,na.rm=T),
#             density_mean=mean(density,na.rm=T),
#             density_sd=sd(density,na.rm=T))
# table_s2
# write.csv(table_s2, file = here::here("/data/outputs/tableS2.csv"))

# }
