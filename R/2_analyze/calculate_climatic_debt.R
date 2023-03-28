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


# get the species data in the right format
dat_train_spec <- dat_spp %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) %>% 
  arrange(bin, core_uniq) %>% 
  select(-c(core_uniq, bin)) %>% 
  as.data.frame() 
  

# same with temperature data
dat_train_temp <- dat_spp %>% 
  distinct(core_uniq, bin, temp_surface, temp_depth) %>% 
  arrange(bin, core_uniq) %>% 
  pull(temp_surface) 

# fit the weighted average partial least squares model (WAPLS)
mod_wapls <- WAPLS(dat_train_spec, dat_train_temp, npls = 10)

# cross-validate model
mod_cv <- crossval(mod_wapls, cv.method = "loo")

# How many components to use based on RMSE via loo
nr_comp <- performance(mod_cv)$crossval[, 1] %>% 
  which.min(.)

# visualize decision
# root mean squared error
plot_rmse <- tibble(model_rmse = performance(mod_cv)$object[, 1],
                    crossval_rmse = performance(mod_cv)$crossval[, 1],
                    nr_components = 1:10) %>% 
  pivot_longer(cols = contains("rmse")) %>% 
 ggplot(aes(nr_components, value, colour = name)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Components", 
       y = "RMSE") +
  scale_color_manual(name = NULL, 
                     values = c(colour_coral, "grey20"), 
                     labels = c("Cross-validated", "Model")) +
  scale_x_continuous(breaks = 1:10) +
  theme(legend.position = c(0.8, 0.8))

# r-squared
plot_rsq <- tibble(model_rsq = performance(mod_cv)$object[, 2],
                   crossval_rsq = performance(mod_cv)$crossval[, 2],
                   nr_components = 1:10) %>% 
  pivot_longer(cols = contains("rsq")) %>% 
  ggplot(aes(nr_components, value, fill = name)) +
  geom_hline(yintercept = performance(mod_cv)$crossval[7, 2], 
             linetype = "dotted", 
             colour = "grey70") +
  geom_col(position = "dodge") +
  coord_cartesian(ylim = c(0.86, 0.91)) +
  labs(x = "Number of Components", 
       y = "R-squared") +
  scale_fill_manual(name = NULL, 
                     values = alpha(c(colour_coral, "grey20"), 0.7), 
                     labels = c("Cross-validated", "Model")) +
  scale_x_continuous(breaks = 1:10) +
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
dat_pred <- dat_spp_full %>% 
  select(-c(core_uniq, bin)) %>% 
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




# extract species niches  -----------------------------------


plot_niche <- dat_spp %>% 
  group_by(species) %>% 
  mutate(mean_temp = mean(temp_surface)) %>% 
  ungroup() %>% 
  ggplot(aes(temp_surface, group = species)) +
  geom_density(aes(fill = mean_temp), 
               alpha = 0.4, 
               colour = "grey60", 
               linewidth = 0.1) +
  geom_hline(yintercept = 0, colour = "white", linewidth = 1) +
  annotate(geom = "point",
           x = c(9, 13, 17), y = -0.012,
           shape = 21,
           size = 5, 
           colour = "grey30") +
  annotate(geom = "text", 
           x = c(9, 13, 17), y = -0.011, 
           label = c("1", "2", "3"),
           size = 10/.pt, 
           colour = "grey30") +
  theme(legend.position = "none") +
  labs(x = "Temperature [°C]", 
       y = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_gradient2(low = "#105157ff",
                       mid = "#8f7270ff",
                       high = colour_lavender, 
                       midpoint = 24,
                       n.breaks = 10) 

# save plot
ggsave(plot_niche, filename = here("figures", 
                                   "fig4_niches.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)
