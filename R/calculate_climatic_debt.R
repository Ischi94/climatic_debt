library(here)
library(tidyverse)
library(rioja)


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                        "spp_by_depth.rds")) %>% 
  # right format for WAPLS
  mutate(rel_abund = rel_abund * 100)

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 



# estimate species temperature optima -------------------------------------


# select training data for WA-PLS
# select bins with average temperature based on quantiles
dat_temp_quant <- dat_mean_temp %>% 
  summarize(temp_quant = quantile(temp_ym_0m)) %>% 
  pull(temp_quant) 

dat_bins <- dat_mean_temp %>% 
  filter(between(temp_ym_0m, dat_temp_quant[2], dat_temp_quant[4]))

dat_train <- dat_spp %>% 
  filter(bin %in% unique(dat_bins$bin))


# get the species data in the right format
dat_train_spec <- dat_train %>% 
  mutate(rel_abund = rel_abund * 100) %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292)) %>% 
  arrange(bin, core_uniq) %>% 
  select(-c(core_uniq, bin)) %>% 
  as.data.frame() 
  

# same with temperature data
dat_train_temp <- dat_train %>% 
  distinct(core_uniq, bin, temp_surface, temp_depth) %>% 
  filter(!(core_uniq == "124_767B" & bin == 292)) %>% 
  arrange(bin, core_uniq) %>% 
  pull(temp_surface) 


# fit the weighted average partial least squares model (WAPLS)
mod_wapls <- WAPLS(dat_train_spec, dat_train_temp)

# cross-validate model
mod_cv <- crossval(mod_wapls, cv.method = "loo")

# How many components to use based on RMSE via loo
nr_comp <- performance(mod_cv)$crossval[, 1] %>% 
  which.min(.)







# using surface water temperature
dat_spp_surf <- dat_spp %>% 
  pluck("temp_ym_0m") %>% 
  pluck("sp") %>% 
  as_tibble() 

# convert cell numbers to coordinates
rEmpt <- raster(ncols=288, nrows=144, xmn=-180, xmx=180, ymn=-90, ymx=90)
dat_spp_surf[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, dat_spp_surf$cell)

# calculate climate debt per cell and per bin
dat_debt <- dat_spp_surf %>% 
  group_by(species) %>%
  summarise(sti = mean(temp_ym)) %>% 
  full_join(dat_spp_surf) %>% 
  arrange(cell) %>% 
  group_by(bin, cell, centroid_long, centroid_lat) %>% 
  summarise(cti = mean(sti), 
            act_temp = mean(temp_ym)) %>% 
  ungroup() %>% 
  select(cell, centroid_long, centroid_lat, bin, act_temp, cti) %>% 
  mutate(temp_lag = act_temp - cti, 
         cell = as.factor(cell))


# save data
write_rds(dat_debt, 
          here("data", 
               "cleaned_debt_raster.rds"))


# visualise over time -----------------------------------------------------


# quick plot
dat_debt %>%
  mutate(dist_equ = abs(centroid_lat)) %>%
  ggplot(aes(bin, temp_lag)) +
  geom_line(aes(group = cell, colour = dist_equ),
            alpha = 0.4) +
  geom_smooth(colour = "coral3") +
  scale_x_reverse() +
  theme_minimal()




