# This script analyses the changes in the STLs over time (calendar year) by biome.

# Load packages ----------------------------------------------------------------
#library(renv)
library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(rFIA)
library(patchwork)
library(terra)
library(sf)
library(lme4)
library(lmerTest)
library(ggeffects)
library(effects)
library(sjPlot)
library(measurements)
library(sp)
library(lqmm)
library(ggforce)
library(MuMIn)
# library(ingestr)
library(DescTools)
library(here)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(purrr)
library(rsample)

# Load functions ---------------------------------------------------------------
# files <- list.files(path = "./R", pattern = "\\.R$", full.names = TRUE)
# walk(files, source)
source(here("R/functions.R"))
source(here("R/plot_stl_bybiome.R"))
source(here("R/identify_disturbed_plots.R"))
source(here("R/get_breaks.R"))
source(here("R/plot_lqmm_bybiome.R"))
source(here("R/calc_lqmm_byqmdbin.R"))

# load data
data_fil_biomes <- readRDS(here("data/data_fil_biomes.rds"))
data_fil_biomes <- readRDS(here("data/75perc/data_fil_biomes.rds"))

plot_map_fil <-  plot_map(data_fil_biomes)
plot_map_fil

# get temporally fixed self-thinning line
plot_stl_fil <- plot_stl(data_fil_biomes)
plot_stl_fil

# STL LMM without interactions--------------------------------------------------
## Biome 1: Tropical & Subtropical Moist Broadleaf Forests  ---------------------

# waldo::compare(
#   readRDS(here::here("data/data_fil_biome1.rds")), 
#   data_fil_biomes |>
#     filter(biomeID == 1)
#   )

data_fil_biome1 <- data_fil_biomes |>
  filter(biomeID == 1)

### Linear mixed effects model --------------------------------------------------
#### Fit model ------------------------------------------------------------------
mod_lmm_biome1 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome1, 
  na.action = "na.exclude"
  )
summary(mod_lmm_biome1)

write_rds(mod_lmm_biome1, file = here::here("data/mod_lmm_biome1.rds"))

#### Plot self-thinning line ----------------------------------------------------
gg_stl_biome1 <- plot_stl_bybiome(
  data_fil_biome1, 
  #data_unm |> filter(biomeID==1),
  mod_lmm_biome1, 
  name = bquote(bold("a") ~~ "Tropical Moist Broadleaf Forests"), 
  years = c(1985, 2000, 2015) 
) 
gg_stl_biome1

#### Data over years ------------------------------------------------------------
gg_hist_year_biome1 <- ggplot(data_fil_biome1, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic() +
  labs(
    title = bquote(bold("a") ~~ "Tropical Moist Broadleaf Forests"),
    x = "Year", 
    y = "Number of invenories"
    )
gg_hist_year_biome1

#### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome1, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
  )

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
    ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

#### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome1)
# print(out)
# r.squaredGLMM(mod_lmm_biome1)
# 
# plot(allEffects(mod_lmm_biome1))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome1,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
#   ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome1$year))
# years

## Biome 2: Tropical & Subtropical Dry Broadleaf Forests  ----------------------

data_fil_biome2 <- data_fil_biomes |>
  filter(biomeID == 2)

### Linear mixed effects model --------------------------------------------------
#### Fit model ------------------------------------------------------------------
mod_lmm_biome2 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome2, 
  na.action = "na.exclude"
)
summary(mod_lmm_biome2)

write_rds(mod_lmm_biome2, file = here::here("data/mod_lmm_biome2.rds"))

#### Plot self-thinning line ----------------------------------------------------
gg_stl_biome2 <- plot_stl_bybiome(
  data_fil_biome2, 
  mod_lmm_biome2, 
  name = bquote(bold("a") ~~ "Tropical Moist Broadleaf Forests"), 
  years = c(1985, 2000, 2015)
)
gg_stl_biome2

ggplot(data_fil_biome2) + geom_point(aes(x=logQMD, y=logDensity, color = plotID))

data_fil_biome2 <- data_fil_biome2 |>
  filter(!(plotID %in% c("mudumalai_7_1","mudumalai_7_2", "mudumalai_7_3", "mudumalai_7_5", "mudumalai_7_6" )))

#### Data over years ------------------------------------------------------------
gg_hist_year_biome2 <- ggplot(data_fil_biome2, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic() +
  labs(
    title = bquote(bold("a") ~~ "Tropical Dry Broadleaf Forests"),
    x = "Year", 
    y = "Number of invenories"
  )
gg_hist_year_biome2

#### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome2, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

## Biome 4: Temperate Broadleaf & Mixed Forests  --------------------------------

data_fil_biome4 <- data_fil_biomes |>
  filter(biomeID == 4)

### Linear mixed effects model --------------------------------------------------
#### Fit model ------------------------------------------------------------------
mod_lmm_biome4 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome4, 
  na.action = "na.exclude"
)

write_rds(mod_lmm_biome4, file = here::here("data/mod_lmm_biome4.rds"))

#### Plot self-thinning line ----------------------------------------------------
gg_stl_biome4 <- plot_stl_bybiome(
  data_fil_biome4, 
  mod_lmm_biome4, 
  name = bquote(bold("b") ~~ "Temperate Broadleaf & Mixed Forests"), 
  years = c(1985, 2000, 2015)
)

#### Data over years ------------------------------------------------------------
gg_hist_year_biome4 <- ggplot(data_fil_biome4, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic() +
  labs(
    title = bquote(bold("b") ~~ "Temperate Broadleaf & Mixed Forests"), 
    x = "Year", 
    y = "Number of invenories"
    )

gg_hist_year_biome4

#### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome4, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

#### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome4)
# print(out)
# r.squaredGLMM(mod_lmm_biome4)
# 
# plot(allEffects(mod_lmm_biome4))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome4,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
# ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome4$year))
# years

## Biome 5: Temperate Conifer Forests Forest  --------------------------------
data_fil_biome5 <- data_fil_biomes |>
  filter(biomeID == 5)

### Linear mixed effects model --------------------------------------------------
#### Fit model ------------------------------------------------------------------
mod_lmm_biome5 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome5, 
  na.action = "na.exclude"
)

write_rds(mod_lmm_biome5, file = here::here("data/mod_lmm_biome5.rds"))

#### Plot self-thinning line ----------------------------------------------------
gg_stl_biome5 <- plot_stl_bybiome(
  data_fil_biome5, 
  mod_lmm_biome5, 
  name = bquote(bold("c") ~~ "Temperate Conifer Forests Forest"), 
  years = c(1985, 2000, 2015)
)

#### Data over years ------------------------------------------------------------
gg_hist_year_biome5 <- ggplot(data_fil_biome5, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic() +
  labs(
    title = bquote(bold("c") ~~ "Temperate Conifer Forests Forest"), 
    x = "Year", 
    y = "Number of invenories"
    )

gg_hist_year_biome5

#### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome5, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

#### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome5)
# print(out)
# r.squaredGLMM(mod_lmm_biome5)
# 
# plot(allEffects(mod_lmm_biome5))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome5,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
# ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome2$year))
# years

## Biome 6: Boreal Forests/Taiga Forest  --------------------------------
data_fil_biome6 <- data_fil_biomes |>
  filter(biomeID == 6)

### Linear mixed effects model --------------------------------------------------
#### Fit model ------------------------------------------------------------------
mod_lmm_biome6 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome6, 
  na.action = "na.exclude"
)

write_rds(mod_lmm_biome6, file = here::here("data/mod_lmm_biome6.rds"))


#### Plot self-thinning line ----------------------------------------------------
gg_stl_biome6 <- plot_stl_bybiome(
  data_fil_biome6, 
  mod_lmm_biome6, 
  name = bquote(bold("d") ~~ "Boreal Forests/Taiga Forest"), 
  years = c(1985, 2000, 2015)
)

#### Data over years ------------------------------------------------------------
gg_hist_year_biome6 <- ggplot(data_fil_biome6, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic() +
  labs(
    title = bquote(bold("d") ~~ "Boreal Forests/Taiga Forest"), 
    x = "Year", 
    y = "Number of invenories"
    )

gg_hist_year_biome6

#### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome6, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

#### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome6)
# print(out)
# r.squaredGLMM(mod_lmm_biome6)
# 
# plot(allEffects(mod_lmm_biome6))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome6,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
# ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome2$year))
# years


## Biome 12: Mediterranean Forests ----------------------------------------------
data_fil_biome12 <- data_fil_biomes |>
  filter(biomeID == 12)

### Linear mixed effects model --------------------------------------------------
#### Fit model ------------------------------------------------------------------
mod_lmm_biome12 = lmer(
  logDensity ~ scale(logQMD) + 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), # + (1|years_since_management_bins), 
  data = data_fil_biome12, 
  na.action = "na.exclude"
)

write_rds(mod_lmm_biome12, file = here::here("data/mod_lmm_biome12.rds"))


#### Plot self-thinning line ----------------------------------------------------
gg_stl_biome12 <- plot_stl_bybiome(
  data_fil_biome12, 
  mod_lmm_biome12, 
  name = bquote(bold("e") ~~ "Mediterranean Forests"), 
  years = c(1985, 2000, 2015)
)

#### Data over years ------------------------------------------------------------
gg_hist_year_biome12 <- ggplot(data_fil_biome12, aes(x = year)) + 
  geom_histogram(color = "black", fill = "grey70", bins = 12) + 
  theme_classic() +
  labs(
    title = bquote(bold("e") ~~ "Mediterranean Forests"), 
    x = "Year", 
    y = "Number of invenories"
    )

gg_hist_year_biome12

#### STL shift ------------------------------------------------------------------
# Mean percent increase in N per year
pred <- ggpredict(
  mod_lmm_biome12, 
  terms = c("logQMD","year [1994, 1995]"), 
  full.data = TRUE
)

change_n <- pred|> 
  as_tibble() |>
  mutate(predicted_trans = exp(predicted)) |> 
  select(x, predicted_trans, group) |> 
  pivot_wider(
    values_from = predicted_trans, 
    names_from = group, 
    names_prefix = "y"
  ) |> 
  ungroup() |> 
  mutate(upSTL = y1995 - y1994) |> 
  mutate(percent = upSTL * 100 / y1994) |> 
  summarise(percent = mean(percent))

change_n

#### Various fit info -----------------------------------------------------------
# out <- summary(mod_lmm_biome12)
# print(out)
# r.squaredGLMM(mod_lmm_biome12)
# 
# plot(allEffects(mod_lmm_biome12))
# 
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate = V1, pvalue = V2) |>
#   mutate(estimate = round(estimate, 3),
#          pvalue = signif(pvalue, digits = 3),
#          pvalue = ifelse(pvalue < 0.001, "< 0.001 ***", pvalue))
# 
# # alternative plot
# plot_model(
#   mod_lmm_biome12,
#   type = "pred",
#   show.data = TRUE, 
#   dot.size = 1.0, 
#   colors = c("#21918c","#fde725", "#440154"),
#   terms = c("logQMD", "year [1985, 2000, 2015]")
# ) + 
#   theme_classic() 
# 
# # years covered
# years <- as.integer(summary(data_fil_biome2$year))
# years
# 


## Special case of BCI ----------------------------------------------------------
# 
# data_bci <- readRDS(here::here("data/inputs/data_bci.rds"))
# data_bci
# plot_stl(data_bci)
# plot_map(data_bci)
# 
# data_bci_1_3 <- data_bci |> 
#   filter(census <= 3)
# 
# data_bci_1_3 <- data_bci |> 
#   filter(census ==1 |
#          census ==2|
#          census ==6)
# 
# 
# data_bci_4_8 <- data_bci |> 
#   filter(census > 3)
# 
# data_bci_1_3_unm <- data_unm_fc(data_bci_1_3)
# data_bci_1_3_fil <- data_filter_fc(data_bci_1_3_unm) 
# data_bci_4_8_unm <- data_unm_fc(data_bci_4_8)
# data_bci_4_8_fil <- data_filter_fc(data_bci_4_8_unm) 
# 
# ggplot() + 
#   geom_point(data = data_bci_1_3_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1.5, col="black",inherit.aes = FALSE) 
# 
# data_bci_1_3_fil <- data_bci_1_3_fil |>
#   filter(logDensity >=7)
# 
# mod_lmm_BCI = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|plotID) + (1|species), # + (1|years_since_management_bins),
#                 data = data_bci_1_3_fil, na.action = "na.exclude")
# out <- summary(mod_lmm_BCI)
# out
# plot(allEffects(mod_lmm_BCI))
# plot_model(mod_lmm_BCI,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate=V1, pvalue=V2) |>
#   mutate(estimate = round(estimate,3),
#          pvalue = signif(pvalue,digits=3),
#          pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
# figbci_1_3 <- plot_model(mod_lmm_BCI,type = "pred",show.data=F ,dot.size=1.0,
#                          terms = c("logQMD", "year")) + theme_classic() + 
#   theme(text=element_text(size=10), plot.title=element_text(size=10)) + 
#   geom_point(data = data_bci_1_3_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
#   labs(title = "BCI - census 1-3", 
#        caption = paste0("Year estimate = ", caption$estimate, '\n', "Year p-value = ", caption$pvalue))
# figbci_1_3
# 
# mod_lmm_BCI = lmer(logDensity ~ scale(logQMD) + scale(year) + (1|plotID) + (1|species), # + (1|years_since_management_bins),
#                 data = data_bci_4_8_fil, na.action = "na.exclude")
# out <- summary(mod_lmm_BCI)
# out
# plot(allEffects(mod_lmm_BCI))
# plot_model(mod_lmm_BCI,type = "pred",show.data=TRUE, dot.size=1.0, terms = c("logQMD","year"))
# caption <- out$coefficients[3,1] |>
#   cbind(out$coefficients[3,5]) |>
#   as_tibble() |>
#   rename(estimate=V1, pvalue=V2) |>
#   mutate(estimate = round(estimate,3),
#          pvalue = signif(pvalue,digits=3),
#          pvalue=ifelse(pvalue<0.001,"<0.001 ***",pvalue))
# figbci_4_8 <- plot_model(mod_lmm_BCI,type = "pred",show.data=F ,dot.size=1.0,
#                          terms = c("logQMD", "year")) + theme_classic() + 
#   theme(text=element_text(size=10), plot.title=element_text(size=10)) + 
#   geom_point(data = data_bci_4_8_fil, aes(x = logQMD, y = logDensity), alpha=0.5, size = 1,col="darkgrey", shape = 16, inherit.aes = FALSE) + 
#   labs(title = "BCI - census 4-8", 
#        caption = paste0("Year estimate = ", caption$estimate, '\n', "Year p-value = ", caption$pvalue))
# figbci_4_8

# STL LMM with interactions  ---------------------------------------------------

## Biome 1: Tropical & Subtropical Moist Broadleaf Forests  --------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome1 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome1, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome1)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome1 <- plot_stl_bybiome(
  data_fil_biome1, 
  mod_lmm_int_biome1, 
  name = bquote(bold("a") ~~ "Tropical Moist Broadleaf Forests"), 
  years = c(1985, 2000, 2015), 
  interactions = TRUE
)

gg_qmd_int_biome1 <- plot_model(
  mod_lmm_int_biome1,
  type = "pred",
  show.data=TRUE, 
  dot.size = 1.0, 
  colors = c("#21918c","#fde725", "#440154"),
  terms = c("year", "logQMD[2.5,3,3.5]")) + 
  theme_classic() + 
  labs(title = bquote(bold("a") ~~ "Tropical Moist Broadleaf Forests"))


## Biome 4: Temperate Broadleaf & Mixed Forests Forest  ------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome4 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome4, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome4)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome4 <- plot_stl_bybiome(
  data_fil_biome4, 
  mod_lmm_int_biome4, 
  name = bquote(bold("b") ~~ "Temperate Broadleaf & Mixed Forests"), 
  years = c(1985, 2000, 2015), 
  interactions = TRUE
)

gg_qmd_int_biome4 <- plot_model(
  mod_lmm_int_biome4,
  type = "pred",
  show.data = TRUE, 
  dot.size = 1.0, 
  colors = c("#21918c", "#fde725", "#440154"),
  terms = c("year","logQMD[2.5,3,3.5]")) + 
  theme_classic() + 
  labs(title = bquote(bold("b") ~~ "Temperate Broadleaf & Mixed Forests"))


## Biome 5: Temperate Conifer Forests    ---------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome5 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome5, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome5)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome5 <- plot_stl_bybiome(
  data_fil_biome5, 
  mod_lmm_int_biome5, 
  name = bquote(bold("c") ~~ "Temperate Conifer Forests Forest"), 
  years = c(1985, 2000, 2015), 
  interactions = TRUE
)

gg_qmd_int_biome5 <- plot_model(
  mod_lmm_int_biome5,
  type = "pred",
  show.data = TRUE, 
  dot.size = 1.0, 
  colors = c("#21918c", "#fde725", "#440154"),
  terms = c("year","logQMD[2.5,3,3.5]")) + 
  theme_classic() + 
  labs(title = bquote(bold("c") ~~ "Temperate Conifer Forests Forest"))


## Biome 6: Boreal Forests/Taiga  ----------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome6 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome6, 
  control = lmerControl(optimizer = "bobyqa"),
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome6)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome6 <- plot_stl_bybiome(
  data_fil_biome6, 
  mod_lmm_int_biome6, 
  name = bquote(bold("d") ~~ "Boreal Forests/Taiga Forest"), 
  years = c(1985, 2000, 2015), 
  interactions = TRUE
)

gg_qmd_int_biome6 <- plot_model(
  mod_lmm_int_biome6,
  type = "pred",
  show.data = TRUE, 
  dot.size = 1.0, 
  colors = c("#21918c", "#fde725", "#440154"),
  terms = c("year","logQMD[2.5,3,3.5]")) + 
  theme_classic() + 
  labs(title = bquote(bold("d") ~~ "Boreal Forests/Taiga Forest"))


## Biome 12: Mediterranean Forests  --------------------------------------------
### Fit model ------------------------------------------------------------------
mod_lmm_int_biome12 = lmer(
  logDensity ~ scale(logQMD) * 
    scale(year) + 
    (1|dataset/plotID) + 
    (1|species), 
  data = data_fil_biome12, 
  na.action = "na.exclude"
)
summary(mod_lmm_int_biome12)

### Plot self-thinning line ----------------------------------------------------
gg_stl_int_biome12 <- plot_stl_bybiome(
  data_fil_biome12, 
  mod_lmm_int_biome12, 
  name = bquote(bold("e") ~~ "Mediterranean Forests"), 
  years = c(1985, 2000, 2015), 
  interactions = TRUE
)

gg_qmd_int_biome12 <- plot_model(
  mod_lmm_int_biome12,
  type = "pred",
  show.data = TRUE, 
  dot.size = 1.0, 
  colors = c("#21918c", "#fde725", "#440154"),
  terms = c("year","logQMD[2.5,3,3.5]")) + 
  theme_classic() + 
  labs(title = bquote(bold("e") ~~ "Mediterranean Forests"))

# Quantile regression ----------------------------------------------------------
data_unm <- readRDS(here::here("data/inputs/agg/data_unm.rds"))

## Biome 1 Tropical & Subtropical Moist Broadleaf Forests ----------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 1) |>
  mutate(year_sc = scale(year),
         logQMD_sc = scale(logQMD))

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots) |> 
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

### Plot disturbed plots -------------------------------------------------------
gg_fdisturbed_biome1 <- df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",  
    title = bquote(bold("a") ~~ "Tropical & Subtropical Moist Broadleaf Forests")
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  # ggplot(aes(var_logdensity, after_stat(count))) + geom_histogram()
  filter(var_logdensity > 0.001)

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.9, #c(0.70, 0.90),
  data = data_unm_biome,
  type = "normal",
  #control = lqmmControl(startQR = TRUE)
  control = lqmmControl(
    LP_max_iter = 2000,   # inner loop iterations
    LP_tol_ll   = 1e-06,  # inner loop tolerance
    startQR     = FALSE
    )
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/fit_lqmm_biome1.rds"))

#### STL shift ------------------------------------------------------------------
# Estimated change in N per unit increase in year

out <- summary(fit_lqmm)

# Extract model coefficient
beta_year_sc <- out$tTable[c("year_sc"), "Value"]

# SD of the original (unscaled) year variable
sd_year <- sd(data_unm_biome$year, na.rm = TRUE)

# Change in logDensity per calendar year
real_coef_year <- beta_year_sc * sd_year

# Convert to % change in tree density per year
percent_change_per_year <- (exp(real_coef_year) - 1) * 100

# Opt 2: predicted logDensity at two years differing by one year (in scaled units)

# Pick a fixed value for logQMD_sc (e.g., mean)
mean_logQMD_sc <- mean(data_unm_biome$logQMD_sc, na.rm = TRUE)
sd_year <- sd(data_unm_biome$year, na.rm = TRUE)

# Predict at year_sc = 0 and year_sc = 1 / sd_year (since one calendar year corresponds to 1/sd_year in scaled units)
newdata1 <- data.frame(logQMD_sc = mean_logQMD_sc, year_sc = 0, plotID = data_unm_biome$plotID[1])
newdata2 <- data.frame(logQMD_sc = mean_logQMD_sc, year_sc = 1 / sd_year, plotID = data_unm_biome$plotID[1])

pred1 <- predict(fit_lqmm, newdata = newdata1)
pred2 <- predict(fit_lqmm, newdata = newdata2)

delta_logDensity <- pred2 - pred1
percent_change <- (exp(delta_logDensity) - 1) * 100

### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>% 
    group_by(plotID), 
  times = 5000, 
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>%  # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome1.rds"))

boot_results |>
  ggplot(aes(estimate)) +
  geom_density() +
  facet_wrap(~term, scales = "free", ncol = 1)

# summarise across bootstraps
summary_stats <- boot_results %>%
  group_by(term) %>%
  summarise(
    estimate = mean(estimate),
    std.error = sd(estimate),
    ci_low = quantile(estimate, 0.025),
    ci_high = quantile(estimate, 0.975),
    .groups = "drop"
  )

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome1 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm, 
  name = bquote(bold("a") ~~ "Tropical Moist Broadleaf Forests")
  )
gg_lqmm_biome1

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
  )

# Build the plot to access internal structure
gg_lqmm_biome1_byqmdbin <- ggplot() +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    size = 1.5,
    color = "grey"
  ) +  
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    width = 0,
    color = "grey"
    ) +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin$df,
    size = 1.5
    ) +
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin$df,
    width = 0
    ) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = expression(ln(QMD)),
    y = expression(italic(beta)(year))
  ) +
  scale_x_continuous(limits = c(2.4, 4.5), breaks = seq(3,4,1)) +
  scale_y_continuous(breaks = seq(-0.02,0.02,0.02))
  
gg_lqmm_biome1_both <- cowplot::plot_grid(
  gg_lqmm_biome1,
  gg_lqmm_biome1_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.4),
  align = "v",
  #labels = c("",  "g"),
  label_y = 1.1
)

gg_lqmm_biome1_both

## Biome 2 Tropical & Subtropical Dry Broadleaf Forests ------------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 2) |>
  mutate(year_sc = scale(year),
         logQMD_sc = scale(logQMD))

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots) |> 
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

### Plot disturbed plots -------------------------------------------------------
gg_fdisturbed_biome2 <- df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",  
    title = bquote(bold("a") ~~ "Tropical & Subtropical Dry Broadleaf Forests")
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  # ggplot(aes(var_logdensity, after_stat(count))) + geom_histogram()
  filter(var_logdensity > 0.001)

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(startQR = TRUE)
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/fit_lqmm_biome2.rds"))

#### STL shift ------------------------------------------------------------------
# Estimated change in N per unit increase in year



### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>% 
    group_by(plotID), 
  times = 5000, 
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>%  # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome1.rds"))

boot_results |>
  ggplot(aes(estimate)) +
  geom_density() +
  facet_wrap(~term, scales = "free", ncol = 1)

# summarise across bootstraps
summary_stats <- boot_results %>%
  group_by(term) %>%
  summarise(
    estimate = mean(estimate),
    std.error = sd(estimate),
    ci_low = quantile(estimate, 0.025),
    ci_high = quantile(estimate, 0.975),
    .groups = "drop"
  )

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome2 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm, 
  name = bquote(bold("b") ~~ "Tropical Dry Broadleaf Forests")
)
gg_lqmm_biome2

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
)

# Build the plot to access internal structure
gg_lqmm_biome2_byqmdbin <- ggplot() +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    size = 1.5,
    color = "grey"
  ) +  
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
    ), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    width = 0,
    color = "grey"
  ) +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin$df,
    size = 1.5
  ) +
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
    ), 
    data = df_lqmm_byqmdbin$df,
    width = 0
  ) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = expression(ln(QMD)),
    y = expression(italic(beta)(year))
  ) +
  scale_x_continuous(limits = c(2.4, 4.5), breaks = seq(3,4,1)) +
  scale_y_continuous(breaks = seq(-0.01,0,0.01))

gg_lqmm_biome2_both <- cowplot::plot_grid(
  gg_lqmm_biome2,
  gg_lqmm_biome2_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.4),
  align = "v",
  #labels = c("",  "h"),
  label_y = 1.1
)

gg_lqmm_biome2_both

## Biome 4 Temperate Broadleaf & Mixed Forests ---------------------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 4) |>
  mutate(year_sc = scale(year),
         logQMD_sc = scale(logQMD))

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots) |> 
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

### Plot disturbed plots -------------------------------------------------------
gg_fdisturbed_biome4 <- df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",  
    title = bquote(bold("b") ~~ "Temperate Broadleaf & Mixed Forests")
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  #control = list(LP_max_iter = 1000, LP_tol_ll = 5e-5),
  control = lqmmControl(startQR = TRUE)
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/fit_lqmm_biome4.rds"))

### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>% 
    group_by(plotID), 
  times = 5000, 
  apparent = FALSE
  )

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>%  # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome4.rds"))

boot_results |>
  ggplot(aes(estimate)) +
  geom_density() +
  facet_wrap(~term, scales = "free", ncol = 1)

# summarise across bootstraps
summary_stats <- boot_results %>%
  group_by(term) %>%
  summarise(
    estimate = mean(estimate),
    std.error = sd(estimate),
    ci_low = quantile(estimate, 0.025),
    ci_high = quantile(estimate, 0.975),
    .groups = "drop"
  )

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome4 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm, 
  name = bquote(bold("c") ~~ "Temperate Broadleaf & Mixed Forests")
)
gg_lqmm_biome4

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
  )

# Build the plot to access internal structure
gg_lqmm_biome4_byqmdbin <- ggplot() +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    size = 1.5,
    color = "grey"
  ) +  
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    width = 0,
    color = "grey"
    ) +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin$df,
    size = 1.5
    ) +
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin$df,
    width = 0
    ) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = expression(ln(QMD)),
    y = expression(italic(beta)(year))
  ) +
  scale_x_continuous(limits = c(2.4, 4.5), breaks = seq(3,4,1)) +
  scale_y_continuous(limits = c(-0.02, 0.04), breaks = seq(-0.02,0.02,0.02))

gg_lqmm_biome4_both <- cowplot::plot_grid(
  gg_lqmm_biome4,
  gg_lqmm_biome4_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.4),
  #labels = c("",  "i"),
  align = "v",
  label_y = 1.1
)

gg_lqmm_biome4_both

## Biome 5  Temperate Conifer Forests Forest -----------------------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 5) |>
  mutate(year_sc = scale(year),
         logQMD_sc = scale(logQMD))

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots) |> 
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

### Plot disturbed plots -------------------------------------------------------
gg_fdisturbed_biome5 <- df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",  
    title = bquote(bold("c") ~~ "Temperate Conifer Forests Forest")
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(startQR = TRUE)
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/fit_lqmm_biome5.rds"))

### Bootstrapping LQMM fit -----------------------------------------------------
# create bootstraps
boot_data <- rsample::bootstraps(
  data_unm_biome %>% 
    group_by(plotID), 
  times = 5000, 
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>%  # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome5.rds"))

boot_results |>
  ggplot(aes(estimate)) +
  geom_density() +
  facet_wrap(~term, scales = "free", ncol = 1)

# summarise across bootstraps
summary_stats <- boot_results %>%
  group_by(term) %>%
  summarise(
    estimate = mean(estimate),
    std.error = sd(estimate),
    ci_low = quantile(estimate, 0.025),
    ci_high = quantile(estimate, 0.975),
    .groups = "drop"
  )

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome5 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm, 
  name = bquote(bold("d") ~~ "Temperate Conifer Forest"))

gg_lqmm_biome5  

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
  )

# Build the plot to access internal structure
gg_lqmm_biome5_byqmdbin <- ggplot() +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    size = 1.5,
    color = "grey"
  ) +  
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    width = 0,
    color = "grey"
    ) +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin$df,
    size = 1.5
    ) +
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin$df,
    width = 0
    ) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = expression(ln(QMD)),
    y = expression(italic(beta)(year))
  ) +
  scale_x_continuous(limits = c(2.4, 4.5), breaks = seq(3,4,1)) +
  scale_y_continuous(limits = c(-0.02, 0.05), breaks = seq(-0.04,0.04,0.02))
  
gg_lqmm_biome5_both <- cowplot::plot_grid(
  gg_lqmm_biome5,
  gg_lqmm_biome5_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.4),
  #labels = c("",  "j"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome5_both

## Biome 6 Boreal Forests/Taiga ------------------------------------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 6) |>
  mutate(year_sc = scale(year),
         logQMD_sc = scale(logQMD))

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots) |> 
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

### Plot disturbed plots -------------------------------------------------------
gg_fdisturbed_biome6 <- df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",  
    title = bquote(bold("d") ~~ "Boreal Forests/Taiga")
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(
    LP_max_iter = 1000,     # increase max iterations
    LP_tol_ll = 1e-4,       # relax tolerance slightly (default is 1e-5)
    startQR = TRUE          # good to keep this TRUE
)
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/fit_lqmm_biome6.rds"))

### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>% 
    group_by(plotID), 
  times = 5000, 
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>%  # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome6.rds"))

boot_results |>
  ggplot(aes(estimate)) +
  geom_density() +
  facet_wrap(~term, scales = "free", ncol = 1)

# summarise across bootstraps
summary_stats <- boot_results %>%
  group_by(term) %>%
  summarise(
    estimate = mean(estimate),
    std.error = sd(estimate),
    ci_low = quantile(estimate, 0.025),
    ci_high = quantile(estimate, 0.975),
    .groups = "drop"
  )

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome6 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm, 
  name = bquote(bold("e") ~~ "Boreal Forests/Taiga")
)

gg_lqmm_biome6

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
  )

# Build the plot to access internal structure
gg_lqmm_biome6_byqmdbin <- ggplot() +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    size = 1.5,
    color = "grey"
  ) +  
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    width = 0,
    color = "grey"
    ) +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin$df,
    size = 1.5
    ) +
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin$df,
    width = 0
    ) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = expression(ln(QMD)),
    y = expression(italic(beta)(year))
  ) +
  scale_x_continuous(limits = c(2.4, 4.5), breaks = seq(3,4,1)) +
  scale_y_continuous(limits = c(-0.03, 0.03), breaks = seq(-0.04,0.04,0.02))
  
gg_lqmm_biome6_both <- cowplot::plot_grid(
  gg_lqmm_biome6,
  gg_lqmm_biome6_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.4),
  #labels = c("",  "k"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome6_both

## Biome 12 Mediterranean Forests ----------------------
data_unm_biome <- data_unm |> 
  filter(biomeID == 12) |>
  mutate(year_sc = scale(year),
         logQMD_sc = scale(logQMD))

### Identify disturbed plots ---------------------------------------------------
data_unm_biome <- data_unm_biome |> 
  identify_disturbed_plots()

breaks <- get_breaks(data_unm_biome$year)

df_disturbed <- data_unm_biome |> 
  mutate(year_bin = cut(
    year, 
    breaks = breaks, 
    labels = breaks[1:length(breaks)-1] + 2.5,
    include.lowest = TRUE
  )) |> 
  group_by(year_bin) |> 
  summarise(
    nplots = length(unique(plotID)),
    ndisturbed = sum(disturbed, na.rm = TRUE)
  ) |> 
  mutate(fdisturbed = ndisturbed / nplots) |> 
  mutate(fdisturbed_logit = log(fdisturbed / (1 - fdisturbed)))

### Plot disturbed plots -------------------------------------------------------
gg_fdisturbed_biome12 <- df_disturbed |> 
  ggplot(aes(as.numeric(as.character(year_bin)), fdisturbed_logit)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_classic() +
  labs(
    x = "Year",  
    title = bquote(bold("d") ~~ "Mediterranean Forests")
  ) +
  scale_y_continuous(
    name = expression(logit(Fraction ~ disturbed)),
    sec.axis = sec_axis(~ plogis(.), name = "Fraction disturbed")
  )

# Additional filter: remove plots with no change in ln(N)
data_unm_biome <- data_unm_biome |>
  group_by(plotID) |>
  mutate(var_logdensity = diff(range(logDensity))) |>
  filter(var_logdensity > 0.001)

# remove disturbed plots
data_unm_biome_including_disturbed <- data_unm_biome
data_unm_biome <- data_unm_biome |> 
  filter(ndisturbed == 0)

### LQMM fit -------------------------------------------------------------------
set.seed(123)
fit_lqmm <- lqmm(
  logDensity ~ logQMD_sc + year_sc,
  random = ~1,
  group = plotID,
  tau = 0.90,
  data = data_unm_biome,
  type = "normal",
  control = lqmmControl(startQR = TRUE)
)
summary(fit_lqmm)

write_rds(fit_lqmm, file = here::here("data/fit_lqmm_biome12.rds"))

#### STL shift ------------------------------------------------------------------
# Opt 1: estimated change in N per unit increase in year

out <- summary(fit_lqmm)

sd(data_unm_biome$year_sc)

# Extract model coefficient
beta_year_sc <- out$tTable[c("year_sc"), "Value"]

# SD of the original (unscaled) year variable
sd_year <- sd(data_unm_biome$year, na.rm = TRUE)

# Change in logDensity per calendar year
real_coef_year <- beta_year_sc * sd_year

# Convert to % change in tree density per year
percent_change_per_year <- (exp(real_coef_year) - 1) * 100

# Opt 2: predicted logDensity at two years differing by one year (in scaled units)

# Pick a fixed value for logQMD_sc (e.g., mean)
mean_logQMD_sc <- mean(data_unm_biome$logQMD_sc, na.rm = TRUE)
sd_year <- sd(data_unm_biome$year, na.rm = TRUE)

# Predict at year_sc = 0 and year_sc = 1 / sd_year (since one calendar year corresponds to 1/sd_year in scaled units)
newdata1 <- data.frame(logQMD_sc = mean_logQMD_sc, year_sc = 0, plotID = data_unm_biome$plotID[1])
newdata2 <- data.frame(logQMD_sc = mean_logQMD_sc, year_sc = 1 / sd_year, plotID = data_unm_biome$plotID[1])

pred1 <- predict(fit_lqmm, newdata = newdata1)
pred2 <- predict(fit_lqmm, newdata = newdata2)

delta_logDensity <- pred2 - pred1
percent_change <- (exp(delta_logDensity) - 1) * 100

### Bootstrapping LQMM fit -----------------------------------------------------
boot_data <- rsample::bootstraps(
  data_unm_biome %>% 
    group_by(plotID), 
  times = 5000, 
  apparent = FALSE
)

# Apply model to each bootstrap sample
boot_results <- boot_data %>%
  mutate(coefs = map(splits, wrap_fit_lqmm)) %>%
  filter(!map_lgl(coefs, is.null)) %>%  # drop failed fits
  unnest(coefs) |>
  dplyr::select(-splits)

write_rds(boot_results, file = here::here("data/boot_results_biome12.rds"))

boot_results |>
  ggplot(aes(estimate)) +
  geom_density() +
  facet_wrap(~term, scales = "free", ncol = 1)

# summarise across bootstraps
summary_stats <- boot_results %>%
  group_by(term) %>%
  summarise(
    estimate = mean(estimate),
    std.error = sd(estimate),
    ci_low = quantile(estimate, 0.025),
    ci_high = quantile(estimate, 0.975),
    .groups = "drop"
  )

### Plot STL from LQMM ---------------------------------------------------------
gg_lqmm_biome12 <- plot_lqmm_bybiome(
  data_unm_biome,
  fit_lqmm, 
  name = bquote(bold("f") ~~ "Mediterranean Forests")
)
gg_lqmm_biome12

### Within QMD bins ----------------------------------------
# Test whether upward shift of 90% quantile is significant within logQMD-bins
# returns data frame with pval indicating significance level of a positive
# effect of year.
df_lqmm_byqmdbin <- calc_lqmm_byqmdbin(data_unm_biome)
df_lqmm_byqmdbin_including_disturbed <- calc_lqmm_byqmdbin(
  data_unm_biome_including_disturbed,
  breaks = df_lqmm_byqmdbin$breaks
  )

# Build the plot to access internal structure
gg_lqmm_biome12_byqmdbin <- ggplot() +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    size = 1.5,
    color = "grey"
  ) +  
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin_including_disturbed$df,
    width = 0,
    color = "grey"
    ) +
  geom_point(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year), 
    data = df_lqmm_byqmdbin$df,
    size = 1.5
    ) +
  geom_errorbar(
    aes(
      as.numeric(as.character(bin_lqmm)), 
      coef_year,
      ymin = coef_year_lower, 
      ymax = coef_year_upper
      ), 
    data = df_lqmm_byqmdbin$df,
    width = 0
    ) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = expression(ln(QMD)),
    y = expression(italic(beta)(year))
  ) +
  scale_x_continuous(limits = c(2.4, 4.5), breaks = seq(3,4,1)) +
  scale_y_continuous(limits = c(-0.03, 0.03), breaks = seq(-0.04,0.04,0.02))
  
gg_lqmm_biome12_both <- cowplot::plot_grid(
  gg_lqmm_biome12,
  gg_lqmm_biome12_byqmdbin,
  ncol = 1,
  rel_heights = c(1, 0.4),
  #labels = c("",  "l"),
  align = "v",
  label_y = 1.1
)
gg_lqmm_biome12_both

# Publication figures ----------------------------------------------------------
## Figure 1 --------------------------------------------------------------------
### No interactions ------------------------------------------------------------
legend <- get_legend(
  gg_stl_biome1 + 
    theme(legend.position = "bottom")
)

fig1 <- cowplot::plot_grid(
  gg_stl_biome1, 
  gg_stl_biome4, 
  gg_stl_biome5, 
  gg_stl_biome6, 
  gg_stl_biome12,
  legend,
  ncol = 3
)

ggsave(
  filename = here::here("manuscript/figures/fig1.pdf"),
  plot = fig1,
  width = 11, 
  height = 7.5
)

ggsave(
  filename = here::here("manuscript/figures/fig1.png"),
  plot = fig1,
  width = 11, 
  height = 7.5
)

### With interactions ----------------------------------------------------------
legend <- get_legend(
  gg_stl_int_biome1 + 
    theme(legend.position = "right")
)

fig1_int <- cowplot::plot_grid(
  gg_stl_int_biome1, 
  gg_stl_int_biome4, 
  gg_stl_int_biome5, 
  gg_stl_int_biome6, 
  gg_stl_int_biome12,
  legend,
  ncol = 3
)

ggsave(
  filename = here::here("manuscript/figures/fig1_int.pdf"),
  plot = fig1_int,
  width = 11, 
  height = 7.5
)

ggsave(
  filename = here::here("manuscript/figures/fig1_int.png"),
  plot = fig1_int,
  width = 11, 
  height = 7.5
)


### Quantile regression ----------------------------------------------------------
library(cowplot)
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
  filename = here::here("manuscript/figures/fig1_lqmm_.pdf"),
  plot = fig1_lqmm,
  width = 11, 
  height = 10
)

ggsave(
  filename = here::here("manuscript/figures/fig1_lqmm_.png"),
  plot = fig1_lqmm,
  width = 11, 
  height = 10
)

## SI Figure histogram over years ----------------------------------------------
fig_hist_year <- cowplot::plot_grid(
  gg_hist_year_biome1, 
  gg_hist_year_biome4, 
  gg_hist_year_biome5, 
  gg_hist_year_biome6, 
  gg_hist_year_biome12,
  ncol = 3
)

ggsave(
  filename = here::here("manuscript/figures/fig_hist_year.pdf"),
  plot = fig_hist_year,
  width = 11, 
  height = 6
)
