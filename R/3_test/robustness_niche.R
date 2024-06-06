library(here)
library(tidyverse)
library(patchwork)
library(rioja)

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

# debt data
dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) 


# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 

# niches ------------------------------------------------------------------

dat_niche <- dat_spp %>% 
  group_by(species) %>% 
  summarise(mean_st = mean(temp_surface), 
            sd_st = sd(temp_surface), 
            mean_dt = mean(temp_depth), 
            sd_dt = sd(temp_depth)) %>% 
  arrange(mean_st) 


# compression effect ------------------------------------------------------

plot_cti <- dat_debt %>% 
  ggplot(aes(temp_surface, cti)) +
  geom_point(shape = 21, 
             fill = "grey90", 
             colour = "grey20", 
             alpha = 0.8) + 
  geom_line(data = tibble(temp_surface = seq(0, 30,
                                              length.out = 20),
                          cti = seq(0, 30,
                                            length.out = 20)),
            colour = "orange") +
  labs(y = "Bio-indicated temperature [°C]", 
       x = "Mean annual surface temperature [°C]")

plot_dev <- dat_debt %>%
  ggplot(aes(temp_surface, climatic_debt)) +
  geom_point(shape = 21, 
             fill = "grey90", 
             colour = "grey20", 
             alpha = 0.8) + 
  geom_hline(yintercept = 0, 
             colour = "orange") +
  labs(y = "Thermal deviance [°C]", 
       x = "Mean annual surface temperature [°C]")





# modern analogue ---------------------------------------------------------

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
mod_mat <- MAT(dat_train_spec, dat_train_temp)

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
  predict(mod_mat, newdata = .) %>% 
  pluck("fit") %>% 
  .[, 2] %>%
  as_tibble() %>% 
  bind_cols(dat_spp_full) %>% 
  select(core_uniq, bin, cti_mat = value)

# merge
plot_debt <- dat_debt %>% 
  full_join(dat_pred) %>% 
  ggplot(aes(cti_mat, cti)) +
  geom_point(shape = 21, 
             fill = "grey90", 
             colour = "grey20", 
             alpha = 0.8) +
  geom_line(data = tibble(cti_mat = seq(0, 33,
                                        length.out = 20),
                          cti = seq(0, 33,
                                    length.out = 20)),
            colour = "orange") +
  labs(y = "WAPLS", 
       x = "Modern analogue technique")

plot_cti/plot_dev/plot_debt + plot_annotation(tag_levels = "a")
