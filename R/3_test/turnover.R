library(here)
library(tidyverse)
library(vegan)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                         "spp_by_depth.rds")) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292))


dat_proxy <- read_rds(here("data", 
                           "proxy_temperature.rds")) %>% 
  # nas are encode as -999
  filter(temp != -999)



# preprocess data ---------------------------------------------------------


# get the species data in the right format
dat_spec <- dat_spp %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  nest_by(bin, zone) %>% 
  ungroup() %>% 
  mutate(data = map(data, 
                    ~ pivot_wider(.x, 
                                  id_cols = core_uniq, 
                                  names_from = species, 
                                  values_from = rel_abund, 
                                  values_fn = mean, 
                                  values_fill = 0) %>% 
                      select(-core_uniq) %>% 
                      colSums() %>% 
                      enframe() %>% 
                      mutate(rel_abund = value/sum(value)) %>% 
                      select(-value) %>% 
                      pivot_wider(names_from = name, values_from = rel_abund, 
                                  values_fill = 0))) %>% 
  unnest(cols = c(data)) %>% 
  replace(is.na(.), 0) 
  

# calculate turnover ------------------------------------------------------


# convert to matrix and calculate turnover
dat_mat <- dat_spec %>% 
  select(- c(bin, zone)) %>% 
  vegdist(method = "chisq") %>% 
  as.matrix()  
  
# add the turnover to the existing data
dat_turnover <- c(dat_mat[row(dat_mat) == col(dat_mat) + 1], 0) %>% 
  as_tibble_col(column_name = "turnover") %>% 
  # join with zones and bin
  bind_cols(dat_spec %>% 
              select(bin, zone)) 

# add temperature
dat_final <-  dat_spp %>%
  # add latitudinal zones
  mutate(
    abs_lat = abs(pal.lat),
    zone = case_when(
      abs_lat >= 60 ~ "High",
      between(abs_lat, 30, 60) ~ "Mid",
      between(abs_lat, 0, 30) ~ "Low"
    )
  ) %>%
  # preprocess temperature data
  group_by(zone, bin) %>%
  summarise(mean_temp = mean(temp_surface)) %>%
  mutate(temp_change = mean_temp - lead(mean_temp)) %>%
  ungroup() %>% 
  # merge with turnover dataframe
  right_join(dat_turnover) %>% 
  drop_na(temp_change)
  

# model turnover over zone and temperature change
dat_model <- dat_final %>% 
  mutate(temp = if_else(temp_change >= 0, "warm", "cool")) %>% 
  group_by(zone, temp) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(mod_gam = map(data, 
                       ~ lm(turnover ~ temp_change, 
                             data = .x)), 
         dat_pred = if_else(temp == "warm",
                            list(tibble(temp_change = seq(0, 2, by = 0.1))),
                            list(tibble(temp_change = seq(0, -2, by = -0.1)))), 
         pred_turnover = map2(mod_gam, 
                              dat_pred, 
                              ~ predict(.x, 
                                        newdata = .y, 
                                        se.fit = TRUE)), 
         pred_points = map2(mod_gam,
                            data,
                            ~ predict(.x,
                                      newdata = .y)),
         pred_mean = map(pred_turnover, 
                         ~ pluck(.x, "fit")), 
         pred_se = map(pred_turnover, 
                       ~ pluck(.x, "se.fit")))  
  

# visualise
plot_turnover <- dat_model %>% 
  unnest(c(dat_pred, pred_mean, pred_se)) %>% 
  mutate(pred_low = pred_mean - 1.96 * pred_se, 
         pred_high = pred_mean + 1.96 * pred_se) %>% 
  mutate(zone = factor(zone, levels = c("High", 
                                        "Mid", 
                                        "Low"))) %>% 
  ggplot(aes(temp_change, pred_mean)) +
  geom_vline(xintercept = 0, colour = "grey80", 
             linetype = "dotted") +
  geom_point(aes(y = pred_points, 
                 colour = zone), 
             data = dat_model %>% 
               unnest(c(data, pred_points)) %>% 
               filter(between(temp_change, -2, 2)) %>% 
               mutate(zone = factor(zone, levels = c("High", 
                                                     "Mid", 
                                                     "Low"))), 
             alpha = 0.1) +
  geom_line(aes(colour = zone, 
                group = interaction(zone, temp))) +
  geom_ribbon(aes(ymin = pred_low, 
                  ymax = pred_high, 
                  group = interaction(zone, temp), 
                  fill = zone), 
              alpha = 0.1) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = 1:3, 
                     limits = c(0.9, 3)) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                     colour_brown,
                                     colour_green), 0.4)) +
  scale_colour_manual(values = c(colour_lavender,
                                 colour_brown,
                                 colour_green)) +
  labs(y = "Compositional Turnover", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  facet_wrap(~ zone) +
  theme(legend.position = "none", 
        axis.ticks = element_blank())
  





# same for diversity ------------------------------------------------------

# calculate species diversity per bin and zone
dat_diversity <- dat_spp %>% 
  # add latitudinal zones
  mutate(
    abs_lat = abs(pal.lat),
    zone = case_when(
      abs_lat >= 60 ~ "High",
      between(abs_lat, 30, 60) ~ "Mid",
      between(abs_lat, 0, 30) ~ "Low"
    )
  ) %>% 
  # ca
  group_by(bin, zone) %>% 
  distinct(species) %>% 
  count(name = "sp_div") %>% 
  ungroup() %>% 
  # add temperature
  left_join(dat_spp %>%
              # add latitudinal zones
              mutate(
                abs_lat = abs(pal.lat),
                zone = case_when(
                  abs_lat >= 60 ~ "High",
                  between(abs_lat, 30, 60) ~ "Mid",
                  between(abs_lat, 0, 30) ~ "Low"
                )
              ) %>%
              # preprocess temperature data
              group_by(zone, bin) %>%
              summarise(mean_temp = mean(temp_surface)) %>%
              mutate(temp_change = mean_temp - lead(mean_temp)) %>%
              ungroup())


dat_diversity %>% 
  ggplot(aes(temp_change, sp_div)) +
  geom_smooth(colour = zone, 
              group = )





# save
ggsave(plot_turnover, filename = here("figures",
                                      "supplemental",
                                      "turnover.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)


# same for proxy data -----------------------------------------------------

# join with dataset
dat_final_proxy <- dat_proxy %>% 
  rename(bin = age) %>% 
  filter(bin %in% dat_spec$bin) %>% 
  mutate(temp_change = temp - lead(temp)) %>% 
  add_column(turnover = dat_diag[1:15]) 

# visualise
plot_turnover_proxy <- dat_final_proxy %>%
  mutate(temp = if_else(temp_change >= 0, "warm", "cool")) %>% 
  drop_na(temp) %>% 
  ggplot(aes(abs(temp_change), turnover)) +
  geom_hline(yintercept = 0, colour = "grey80", 
             linetype = "dotted") +
  geom_point() +
  geom_smooth(method = "lm", 
              formula = y ~ x + poly(x, 2),
              colour = colour_coral, 
              fill = colour_coral, 
              alpha = 0.2) +
  labs(y = "Compositional Turnover", 
       x = expression(paste("Absolute " ,Delta, " Temperature [°C]"))) +
  theme(legend.position = "none")

# save
ggsave(plot_turnover_proxy, filename = here("figures",
                                      "supplemental",
                                      "turnover_proxy.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)

