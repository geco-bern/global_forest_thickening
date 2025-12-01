library(lqmm)
library(tidyverse)
library(MASS)

# load data
data_fil_biomes <- readRDS(here("data/data_fil_biomes.rds"))

get_samples_biomass_change_bybiome <- function(biome_number, data_fil_biomes, n_sim_biomass){
  
  df_biome <- data_fil_biomes |>
    filter(biomeID == biome_number) |>
    mutate(NQMD2 = density * QMD^2)
  
  # Get samples of coefficients estimated from bootstrapping
  boot_results <- read_rds(here::here(paste0("data/boot_results_biome", biome_number, ".rds")))
  
  # Fit biomass model (again)
  fit_biomass = lmer(
    biomass ~ NQMD2 + 0 + (1|plotID), 
    data = df_biome , 
    na.action = "na.exclude"
  )
  
  fixef_means_biomass <- fixef(fit_biomass)
  vcov_matrix_biomass <- vcov(fit_biomass)  # variance-covariance matrix of fixed effects
  
  # Simulate 1000 draws of the coefficients from a multivariate normal distribution
  get_df_by_biomass_sample <- function(
    idx,
    df_biome,
    boot_results, 
    fixef_means_biomass,
    vcov_matrix_biomass,
    n_sim_plots
  ){
    
    # extract sampled coefficients of the STL-model from bootstrapping results
    coef_samples_selfthinning <- boot_results |> 
      pivot_wider(
        names_from = "term",
        values_from = "estimate"
      ) |> 
      dplyr::select(-id) |> 
      as.matrix()
    
    # generate a single sample of coefficients of the biomass model from variance-covariance matrix
    coef_samples_biomass <- MASS::mvrnorm(
      n = 1, 
      mu = fixef_means_biomass, 
      Sigma = vcov_matrix_biomass
    )
    
    # New observation (must match predictor names)
    new_data_t0 <- df_biome |> 
      dplyr::select(logDensity, logQMD, NQMD2, density, QMD, plotID, dataset, species) |> 
      slice_sample(n = n_sim_plots) |> 
      mutate(year = 2000) |> 
      mutate(idx_sample_plot = 1:n())
    
    new_data_t1 <- new_data_t0 |> 
      mutate(year = 2001)
    
    # Construct model matrix for new data
    X_new_t0 <- model.matrix( ~ logQMD + year, data = new_data_t0)
    X_new_t1 <- model.matrix( ~ logQMD + year, data = new_data_t1)
    
    # Simulate predicted values
    pred_sim_t0 <- coef_samples_selfthinning %*% t(X_new_t0)  # matrix multiply
    pred_sim_t1 <- coef_samples_selfthinning %*% t(X_new_t1)  # matrix multiply
    
    # construct data frame
    df <- pred_sim_t0 |> 
      as_tibble() %>%
      rowid_to_column(var = "idx_sample_coef") %>%
      pivot_longer(
        cols = -idx_sample_coef,
        names_to = "idx_sample_plot",
        names_transform = list(column = as.integer),
        values_to = "pred_sim_t0"
      ) |> 
      left_join(
        pred_sim_t1 |> 
          as_tibble() %>%
          rowid_to_column(var = "idx_sample_coef") %>%
          pivot_longer(
            cols = -idx_sample_coef,
            names_to = "idx_sample_plot",
            names_transform = list(column = as.integer),
            values_to = "pred_sim_t1"
          ),
        by = c("idx_sample_coef", "idx_sample_plot")
      ) |> 
      mutate(
        # Back-transform
        pred_density_t0 = exp(pred_sim_t0),
        pred_density_t1 = exp(pred_sim_t1)
      ) |> 
      mutate(
        # change in number of trees per ha in one year
        diff_density = pred_density_t1 - pred_density_t0
      ) |> 
      mutate(idx_sample_plot = as.integer(idx_sample_plot)) |> 
      
      # add data for biomass change estimate dB = b * diff_density * QMD2
      left_join(
        new_data_t0,
        by = "idx_sample_plot"
      ) |> 
      
      # biomass change estimate dB = b * diff_density * QMD2
      mutate(
        dB = c(coef_samples_biomass) * diff_density * QMD^2
      ) |> 
      mutate(
        idx_sample_biomass = idx
      )
    return(df)
  }
  
  df <- purrr::map_dfr(
    1:n_sim_biomass,
    ~get_df_by_biomass_sample(
      .,
      df_biome = data_fil_biome1,
      boot_results, 
      fixef_means_biomass,
      vcov_matrix_biomass,
      n_sim_plots = 30
    )
  ) |> 
    mutate(biome_number = biome_number)
  
  return(df)
}


# get_samples_biomass_change_bybiome(
#   biome_number = 1, 
#   data_fil_biomes, 
#   n_sim_biomass = 30
# )

df <- purrr::map_dfr(
  as.list(c(1, 4, 5, 6, 12)),
  ~get_samples_biomass_change_bybiome(
    biome_number = ., 
    data_fil_biomes, 
    n_sim_biomass = 30
  )
)

ggplot() +
  geom_density(
    aes(
      dB, 
      group = as.factor(biome_number), 
      fill =  as.factor(biome_number),
      color = as.factor(biome_number)), 
    data = df,
    alpha = 0.5
  ) +
  khroma::scale_color_okabeito() +
  khroma::scale_fill_okabeito() +
  theme_classic()

ggsave(here::here("manuscript/figures/csink_bybiome.pdf"), width = 6, height = 4)

#---------------------------OLD BELOW-------------------------------------------

# # This script evaluates the mature forests C sink for each forest  biome represented in the data
# 
# library(tidyverse)
# 
# # load data ----
# data_fil_biome <- readRDS(here::here("data/data_fil_biomes.rds"))
# 
# data_fil_biome <- data_fil_biome |>
#   mutate(NQMD2 = density * QMD^2)  # new variable for modelling biomass 
# 
# ## Biome 1 ---------------------------------------------------------------------
# # Tropical & Subtropical Moist Broadleaf Forests
# df_sub <- data_fil_biome |> 
#   filter(biomeID == 1)
# 
# fit_selfthinning = lmer(
#   log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), 
#   data = df_sub, 
#   na.action = "na.exclude"
# )
# 
# fit_biomass = lmer(
#   biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
#   data = df_sub , 
#   na.action = "na.exclude"
# )
# a_mean <- summary(fit_biomass)$coefficient[1,1]
# a_sd <- summary(fit_biomass)$coefficient[1,2]
# 
# # Original method: considering uncertainty in coefficient relating N*QMD^2 to biomass
# filn <- here::here("data/out_csink_biome1_ORI.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_ORI(df_sub,a_mean,a_sd))
#   saveRDS(out, here::here("data/out_csink_biome1_ORI.rds"))
# } 
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# # New method: considering residual error (prediction error) in biomass estimate, given N*QMD^2
# filn <- here::here("data/out_csink_biome1_NEW.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_NEW(df_sub))
#   out
#   saveRDS(out, here::here("data/out_csink_biome1_NEW.rds"))
# }
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# ## Biome 4 ----
# # Temperate Broadleaf & Mixed Forests 
# df_sub <- data_fil_biome |> 
#   filter(biomeID == 4)
# 
# df_sub |>
#   ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + geom_point() + 
#   coord_obs_pred()
# 
# fit_selfthinning = lmer(
#   log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), 
#   data = df_sub, 
#   na.action = "na.exclude"
# )
# 
# fit_biomass = lmer(
#   biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
#   data = df_sub , 
#   na.action = "na.exclude"
# )
# a_mean <- summary(fit_biomass)$coefficient[1,1]
# a_sd <- summary(fit_biomass)$coefficient[1,2]
# 
# # Original method: considering uncertainty in coefficient relating N*QMD^2 to biomass
# filn <- here::here("data/out_csink_biome4_ORI.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_ORI(df_sub,a_mean,a_sd))
#   saveRDS(out, here::here("data/out_csink_biome4_ORI.rds"))
# } 
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# # New method: considering residual error (prediction error) in biomass estimate, given N*QMD^2
# filn <- here::here("data/out_csink_biome4_NEW.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_NEW(df_sub))
#   out
#   saveRDS(out, here::here("data/out_csink_biome4_NEW.rds"))
# }
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# ## Biome 5 ----
# # Temperate Conifer Forests  
# df_sub <- data_fil_biome |> 
#   filter(biomeID == 5)
# 
# df_sub |>
#   ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + geom_point() + 
#   coord_obs_pred()
# 
# fit_selfthinning = lmer(
#   log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), 
#   data = df_sub, 
#   na.action = "na.exclude"
# )
# 
# fit_biomass = lmer(
#   biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
#   data = df_sub , 
#   na.action = "na.exclude"
# )
# a_mean <- summary(fit_biomass)$coefficient[1,1]
# a_sd <- summary(fit_biomass)$coefficient[1,2]
# 
# # Original method: considering uncertainty in coefficient relating N*QMD^2 to biomass
# filn <- here::here("data/out_csink_biome5_ORI.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_ORI(df_sub,a_mean,a_sd))
#   saveRDS(out, here::here("data/out_csink_biome5_ORI.rds"))
# } 
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# # New method: considering residual error (prediction error) in biomass estimate, given N*QMD^2
# filn <- here::here("data/out_csink_biome5_NEW.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_NEW(df_sub))
#   out
#   saveRDS(out, here::here("data/out_csink_biome5_NEW.rds"))
# }
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# ## Biome 6 ----
# # Boreal Forests/Taiga 
# df_sub <- data_fil_biome |> 
#   filter(biomeID == 6)
# 
# df_sub |>
#   ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + geom_point() + 
#   coord_obs_pred()
# 
# fit_selfthinning = lmer(
#   log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), 
#   data = df_sub, 
#   na.action = "na.exclude"
# )
# 
# fit_biomass = lmer(
#   biomass ~ NQMD2 + 0 + (1|plotID), 
#   data = df_sub , 
#   na.action = "na.exclude"
# )
# a_mean <- summary(fit_biomass)$coefficient[1,1]
# a_sd <- summary(fit_biomass)$coefficient[1,2]
# 
# # Original method: considering uncertainty in coefficient relating N*QMD^2 to biomass
# filn <- here::here("data/out_csink_biome6_ORI.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_ORI(df_sub,a_mean,a_sd))
#   saveRDS(out, here::here("data/out_csink_biome6_ORI.rds"))
# } 
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# # New method: considering residual error (prediction error) in biomass estimate, given N*QMD^2
# filn <- here::here("data/out_csink_biome6_NEW.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_NEW(df_sub))
#   out
#   saveRDS(out, here::here("data/out_csink_biome6_NEW.rds"))
# }
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# ## Biome 12 ----
# # Mediterranean Forests
# df_sub <- data_fil_biome |> 
#   filter(biomeID == 12)
# 
# df_sub |>
#   ggplot(aes(y = biomass, x = density*QMD^2, col=dataset)) + geom_point() + 
#   coord_obs_pred()
# 
# fit_selfthinning = lmer(
#   log(density) ~ logQMD + year + (1|dataset/plotID) + (1|species), 
#   data = df_sub, 
#   na.action = "na.exclude"
# )
# 
# fit_biomass = lmer(
#   biomass ~ NQMD2 + 0 + (1|dataset/plotID), 
#   data = df_sub, 
#   na.action = "na.exclude",
#   control = lmerControl(optimizer = "bobyqa")
# )
# a_mean <- summary(fit_biomass)$coefficient[1,1]
# a_sd <- summary(fit_biomass)$coefficient[1,2]
# 
# # Original method: considering uncertainty in coefficient relating N*QMD^2 to biomass
# filn <- here::here("data/out_csink_biome12_ORI.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_ORI(df_sub,a_mean,a_sd))
#   saveRDS(out, here::here("data/out_csink_biome12_ORI.rds"))
# } 
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# # New method: considering residual error (prediction error) in biomass estimate, given N*QMD^2
# filn <- here::here("data/out_csink_biome12_NEW.rds")
# if (!file.exists(filn)){
#   out <- purrr::map_dfr(
#     as.list(seq(1e5)), #1e5
#     ~csink_NEW(df_sub))
#   out
#   saveRDS(out, here::here("data/out_csink_biome12_NEW.rds"))
# }
# 
# out |> 
#   ggplot(aes(dB_Mg_ha, ..density..)) +
#   geom_density() +
#   theme_classic() +
#   labs(title = "Biomass change, dB per ha")
# 
# # Total mature C sink ----
# # dB per yr
# 
# # Get forest fraction cover for each biome
# 
# # load biomes vector from the WWF Ecoregions data
# v_biomes <- terra::vect(here::here("data/open_data/wwf/wwf_terr_ecos.shp"))
# values(v_biomes)
# # Check the CRS
# crs_info <- crs(v_biomes)
# 
# # load modis fraction forest cover raster
# r_fcf <- terra::rast(here::here("data/open_data/modis/MODIS_ForestCoverFraction.nc"))
# names(r_fcf)
# # select only the forestcoverfraction
# r_fcf <- r_fcf[[1]]
# plot(r_fcf)
# 
# # Ensure CRS Alignment
# # Reproject polygon to raster CRS, if needed
# v_biomes <- project(v_biomes, crs(r_fcf))  
# 
# # Calculate area of the grid in projected units
# v_biomes$area_m2 <- expanse(v_biomes, unit = "m")
# values(v_biomes)
# 
# # Extract raster values for each polygon
# extracted <- extract(r_fcf, v_biomes, fun = NULL, na.rm = TRUE, touches = TRUE)
# 
# # Combine extracted values with polygon attributes
# v_biomes$fcf <- extracted$forestcoverfraction  # Assuming 'layer' contains the raster values
# values(v_biomes)
# 
# # View combined polygon attributes
# biomes_fcf <- as.data.frame(v_biomes) |>
#   mutate(fcf_perc = fcf*1e-2,
#          forestcover_area_ha = area_m2*1e-4*fcf_perc) |>
#   group_by(BIOME) |>
#   summarise(total_forestcover_area_ha = sum(forestcover_area_ha, na.rm = T)) |>
#   mutate(forestcover_scinot = (format(total_forestcover_area_ha, scientific = TRUE)))
# 
# # 5. Taken the dB estimates per biome, multiple this by the total forest area within this biome 
# # to get the distribution of the biome-level total mature forest C sink.
# # Mg = 1e-9 Petagrams. 1PG of C = 1e15 grams
# 
# out_csink_biome1 <- #readRDS(here::here("data/out_csink_biome1_ORI.rds")) |>
#   readRDS(here::here("data/out_csink_biome1_NEW.rds")) |>
#   mutate(biome = "Tropical Moist Broadleaf Forests") |>
#   mutate(dB_Mg_yr = dB_Mg_ha * (biomes_fcf |> filter(BIOME==1))$total_forestcover_area_ha) 
# 
# out_csink_biome4 <- #readRDS(here::here("data/out_csink_biome4_ORI.rds")) |>
#   readRDS(here::here("data/out_csink_biome4_NEW.rds")) |>
#   mutate(biome = "Temperate Broadleaf & Mixed Forests") |>
#   mutate(dB_Mg_yr = dB_Mg_ha * (biomes_fcf |> filter(BIOME==4))$total_forestcover_area_ha) 
# 
# out_csink_biome5 <- #readRDS(here::here("data/out_csink_biome5_ORI.rds")) |>
#   readRDS(here::here("data/out_csink_biome5_NEW.rds")) |>
#   mutate(biome = "Temperate Conifer Forests") |>
#   mutate(dB_Mg_yr = dB_Mg_ha * (biomes_fcf |> filter(BIOME==5))$total_forestcover_area_ha) 
# 
# out_csink_biome6 <- #readRDS(here::here("data/out_csink_biome6_ORI.rds")) |>
#   readRDS(here::here("data/out_csink_biome6_NEW.rds")) |>
#   mutate(biome = "Boreal Forests/Taiga") |>
#   mutate(dB_Mg_yr = dB_Mg_ha * (biomes_fcf |> filter(BIOME==6))$total_forestcover_area_ha) 
# 
# out_csink_biome12 <- #readRDS(here::here("data/out_csink_biome12_ORI.rds")) |>
#   readRDS(here::here("data/out_csink_biome12_NEW.rds")) |>
#   mutate(biome = "Mediterranean forests") |>
#   mutate(dB_Mg_yr = dB_Mg_ha * (biomes_fcf |> filter(BIOME==12))$total_forestcover_area_ha) 
# 
# # All biomes plots ----
# 
# out_csink <- out_csink_biome1 |>
#   bind_rows(out_csink_biome4) |>
#   bind_rows(out_csink_biome5) |>
#   bind_rows(out_csink_biome6) |>
#   bind_rows(out_csink_biome12) |>
#   mutate(dC_Mg_ha = dB_Mg_ha*0.5,
#          dB_Pg_yr = dB_Mg_yr*1e-9,
#          dC_Pg_yr = dB_Pg_yr*0.5) |>
#   mutate(biome = factor(biome, 
#                         levels =c("Tropical Moist Broadleaf Forests",
#                                   "Temperate Broadleaf & Mixed Forests",
#                                   "Temperate Conifer Forests",
#                                   "Boreal Forests/Taiga",
#                                   "Mediterranean forests")))
# 
# # C sink changes ----
# dC_Mg_ha <- out_csink |> 
#   ggplot(aes(dC_Mg_ha, ..density.., col=factor(biome), fill=factor(biome))) +
#   geom_density() +
#   theme_classic() +
#   labs(#title = "Forest carbon per ha and year",
#     x = expression(paste("Mg C ", ha^-1, " ",yr^-1)), y = "Density") +
#   scale_fill_okabe_ito(name = "Biome", alpha = .7) +
#   scale_color_okabe_ito(name = "Biome", alpha = .7) +
#   scale_x_continuous(limits = c(0,2)) +  
#   theme(legend.position = "bottom")
# dC_Mg_ha
# 
# # Figure 3b ----
# dC_Pg_yr <- out_csink |> 
#   ggplot(aes(dC_Pg_yr, ..density.., col=biome, fill=biome)) +
#   geom_density() +
#   theme_classic() +
#   labs(#title = "Forest carbon per biome and year",
#     x = expression(paste("Pg C ", yr^-1)), y = "Density") +
#   scale_fill_okabe_ito(name = "Biome", alpha = .7) +
#   scale_color_okabe_ito(name = "Biome", alpha = .7) +
#   scale_x_continuous(limits = c(0,1)) +  
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),axis.title = element_text(size = 12),
#         axis.text.y = element_text(hjust = 0.5),
#         legend.text = element_text(size = 9),legend.title = element_text(size = 10),
#         plot.title = element_text(size = 12),
#         legend.key = element_rect(fill = NA, color = NA),
#         legend.position = "none", #c(0.65, 0.75),
#         legend.direction="vertical",
#         legend.box = "horizontal",
#         legend.margin = margin(2, 2, 2, 2),
#         legend.key.size = unit(.5, 'cm'),
#         legend.box.margin = margin(1, 1, 1, 1)) 
# dC_Pg_yr
# 
# # Weighted mean of Mg per ha ----
# 
# # Sum across the same row_id from each group
# biomes_fcf_observed <- biomes_fcf |> 
#   filter(BIOME==1|BIOME==4|BIOME==5|BIOME==6|BIOME==12)
# 
# out_csink_wm <- out_csink |>
#   group_by(biome) |>
#   mutate(row_id = row_number()) |>
#   ungroup() |>
#   group_by(row_id) |>
#   summarise(dB_Mg_ha_allbiomes = sum(dB_Mg_yr)/sum(biomes_fcf_observed$total_forestcover_area_ha), #weighted mean 
#             dB_Pg_yr_allbiomes = sum(dB_Pg_yr)) |>
#   mutate(dC_Mg_ha_allbiomes = dB_Mg_ha_allbiomes*0.5,
#          dC_Pg_yr_allbiomes = dB_Pg_yr_allbiomes*0.5)
# 
# # add the agreggated all biomes as new raws for plotting
# out_csink <- out_csink |>
#   bind_rows(data.frame(dC_Mg_ha = out_csink_wm$dC_Mg_ha_allbiomes, biome = "Observed biomes")) |>
#   mutate(biome = factor(biome, 
#                         levels =c("Observed biomes",
#                                   "Tropical Moist Broadleaf Forests",
#                                   "Temperate Broadleaf & Mixed Forests",
#                                   "Temperate Conifer Forests",
#                                   "Boreal Forests/Taiga",
#                                   "Mediterranean forests")))
# 
# # Figure 3a ----
# dC_Mg_ha_allbiomes <- ggplot() +
#   geom_density(data = out_csink, aes(dC_Mg_ha, ..density.., col=factor(biome), fill=factor(biome)), alpha= 0.7) +
#   theme_classic() +
#   labs(#title = "Forest carbon per ha and year",
#     x = expression(paste("Mg C ", ha^-1, " ",yr^-1)), y = "Density") +
#   #scale_fill_okabe_ito(name = "Biome", alpha = .7) +
#   #scale_color_okabe_ito(name = "Biome", alpha = .7) +
#   scale_fill_manual(name = "Biome", values = c("Observed biomes" = "#999999",
#                                                "Tropical Moist Broadleaf Forests" = "#E69F00", 
#                                                "Temperate Broadleaf & Mixed Forests" = "#56B4E9",
#                                                "Temperate Conifer Forests" = "#009E73", 
#                                                "Boreal Forests/Taiga" = "#F5C710", 
#                                                "Mediterranean forests" = "#0072B2")) +
#   scale_color_manual(name = "Biome", values = c("Observed biomes" = "#999999",
#                                                 "Tropical Moist Broadleaf Forests" = "#E69F00", 
#                                                 "Temperate Broadleaf & Mixed Forests" = "#56B4E9",
#                                                 "Temperate Conifer Forests" = "#009E73", 
#                                                 "Boreal Forests/Taiga" = "#F5C710", 
#                                                 "Mediterranean forests" = "#0072B2")) +
#   scale_x_continuous(limits = c(0,2)) +  
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),axis.title = element_text(size = 12),
#         axis.text.y = element_text(hjust = 0.5),
#         legend.text = element_text(size = 10),legend.title = element_text(size = 10),
#         plot.title = element_text(size = 12),
#         legend.key = element_rect(fill = NA, color = NA),
#         legend.position = c(0.65, 0.75),
#         legend.direction="vertical",
#         legend.box = "horizontal",
#         legend.margin = margin(2, 2, 2, 2),
#         legend.key.size = unit(.5, 'cm'),
#         legend.box.margin = margin(1, 1, 1, 1))
# dC_Mg_ha_allbiomes
# 
# # Figure 3c ----
# dC_Pg_yr_allbiomes <- out_csink_wm |> 
#   ggplot(aes(dC_Pg_yr_allbiomes, ..density..)) +
#   geom_density(color="#999999",fill="#999999", alpha =0.7) +
#   theme_classic() +
#   labs(#title = "Total forest carbon per year",
#     x = expression(paste("Pg C ", yr^-1)), y = "Density") +
#   scale_x_continuous(limits = c(0.6,1.6)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),axis.title = element_text(size = 12),
#         axis.text.y = element_text(hjust = 0.5)) 
# dC_Pg_yr_allbiomes
# 
# # Figure 3abc ----
# fig3abc <- dC_Mg_ha_allbiomes + dC_Pg_yr + dC_Pg_yr_allbiomes +
#   plot_layout(ncol = 2)
# fig3abc
# 
# # Mean estimates (for paper)
# out_csink |>
#   group_by(biome) |>
#   summarise(dC_Mg_ha_mean = mean(dC_Mg_ha),
#             dC_Mg_ha_sd = sd(dC_Mg_ha),
#             dC_Pg_yr_mean = mean(dC_Pg_yr),
#             dC_Pg_yr_sd = sd(dC_Pg_yr))
