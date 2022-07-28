library(here)
library(tidyverse)
library(rioja)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                        "spp_by_depth.rds")) %>% 
  # right format for WAPLS
  mutate(rel_abund = rel_abund * 100) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292))

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 



# fit WAPLS -------------------------------------


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
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) %>% 
  arrange(bin, core_uniq) %>% 
  select(-c(core_uniq, bin)) %>% 
  as.data.frame() 
  

# same with temperature data
dat_train_temp <- dat_train %>% 
  distinct(core_uniq, bin, temp_surface, temp_depth) %>% 
  arrange(bin, core_uniq) %>% 
  pull(temp_surface) 

# fit the weighted average partial least squares model (WAPLS)
mod_wapls <- WAPLS(dat_train_spec, dat_train_temp)

# cross-validate model
mod_cv <- crossval(mod_wapls, cv.method = "loo")

# How many components to use based on RMSE via loo
nr_comp <- performance(mod_cv)$crossval[, 1] %>% 
  which.min(.)

# visualize decision
# root mean squared error
plot_rmse <- tibble(model_rmse = performance(mod_cv)$object[, 1],
                    crossval_rmse = performance(mod_cv)$crossval[, 1],
                    nr_components = 1:5) %>% 
  pivot_longer(cols = contains("rmse")) %>% 
 ggplot(aes(nr_components, value, colour = name)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Components", 
       y = "RMSE") +
  scale_color_manual(name = NULL, 
                     values = c(colour_coral, "grey20"), 
                     labels = c("Cross-validated", "Model")) +
  theme(legend.position = c(0.8, 0.8))

# r-squared
plot_rsq <- tibble(model_rsq = performance(mod_cv)$object[, 2],
                   crossval_rsq = performance(mod_cv)$crossval[, 2],
                   nr_components = 1:5) %>% 
  pivot_longer(cols = contains("rsq")) %>% 
  ggplot(aes(nr_components, value, fill = name)) +
  geom_hline(yintercept = 0.9112843, linetype = "dotted", 
             colour = "grey70") +
  geom_col(position = "dodge") +
  coord_cartesian(ylim = c(0.88, 0.92)) +
  labs(x = "Number of Components", 
       y = "R-squared") +
  scale_fill_manual(name = NULL, 
                     values = alpha(c(colour_coral, "grey20"), 0.7), 
                     labels = c("Cross-validated", "Model")) +
  theme(legend.position = "none")

# combine 
plot_final <- plot_rmse + plot_rsq +
  plot_annotation(tag_levels = "a") 

# save
ggsave(plot_final, filename = here("figures",
                                   "supplemental",
                                   "wapls_summary.png"), 
         width = image_width, height = image_height, units = image_units, 
         bg = "white", device = ragg::agg_png)


# predict community temperature index -------------------------------------


# predict temperature on whole dataset using WAPLS
# bring data in right format
dat_spp_full <- dat_spp %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) 

# predict community temperature index (cti)
dat_pred <- dat_spp_full %>% select(-c(core_uniq, bin)) %>% 
  as.data.frame() %>% 
  predict(mod_wapls, newdata = ., npls = nr_comp,
          nboot = 1000) %>% 
  pluck("fit") %>% 
  .[, nr_comp] %>% 
  as_tibble() %>% 
  bind_cols(dat_spp_full) %>% 
  select(core_uniq, bin, cti = value)

# add estimated temperature from aogcm's
dat_debt <- dat_spp %>% 
  distinct(core_uniq, bin, 
           temp_surface, 
           pal.lat, pal.long) %>% 
  left_join(dat_pred) %>% 
  # calculate offset between cti
  mutate(climatic_debt = temp_surface - cti) 


  
# save data
write_rds(dat_debt, 
          here("data", 
               "cleaned_debt_wapls.rds"))


# visualise over time -----------------------------------------------------


# quick plot
dat_debt %>%
  mutate(dist_equ = abs(pal.lat)) %>%
  ggplot(aes(bin, climatic_debt)) +
  geom_line(aes(group = core_uniq, colour = dist_equ),
            alpha = 0.4) +
  geom_smooth(colour = "coral3") +
  scale_x_reverse() 




# based on cti ------------------------------------------------------------



# # convert cell numbers to coordinates
# rEmpt <- raster(ncols=288, nrows=144, xmn=-180, xmx=180, ymn=-90, ymx=90)
# dat_spp_surf[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, dat_spp_surf$cell)
# 
# # calculate climate debt per cell and per bin
# dat_debt <- dat_spp_surf %>% 
#   group_by(species) %>%
#   summarise(sti = mean(temp_ym)) %>% 
#   full_join(dat_spp_surf) %>% 
#   arrange(cell) %>% 
#   group_by(bin, cell, centroid_long, centroid_lat) %>% 
#   summarise(cti = mean(sti), 
#             act_temp = mean(temp_ym)) %>% 
#   ungroup() %>% 
#   select(cell, centroid_long, centroid_lat, bin, act_temp, cti) %>% 
#   mutate(temp_lag = act_temp - cti, 
#          cell = as.factor(cell))
# 
# 
# # save data
# write_rds(dat_debt, 
#           here("data", 
#                "cleaned_debt_raster.rds"))







