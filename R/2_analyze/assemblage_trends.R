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
  slice(rep(1:n(), each = 3)) %>% 
  mutate(temp = if_else(temp_change >= 0, "warm", "cool")) %>% 
  group_by(zone, temp) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(mod_gam = map(data, 
                       ~ lm(turnover ~ temp_change, 
                             data = .x)), 
         dat_pred = if_else(temp == "warm",
                            list(tibble(temp_change = seq(0, 3, by = 0.1))),
                            list(tibble(temp_change = seq(0, -3, by = -0.1)))), 
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
  group_by(zone) %>% 
  mutate(data = map(data, ~ filter(.x, between(.x$temp_change, -3, 3))), 
         max_cooling = map_dbl(data, ~ min(.x$temp_change)), 
         max_warming =  map_dbl(data, ~ max(.x$temp_change)), 
         line_range = if_else(between(temp_change, 
                                      max_cooling, 
                                      max_warming), 
                              1, 
                              0)) %>% 
  ungroup() %>% 
  ggplot(aes(temp_change, pred_mean)) +
  geom_vline(xintercept = 0, colour = "grey80", 
             linetype = "dotted") +
  geom_point(aes(y = pred_points, 
                 colour = zone), 
             data = dat_model %>% 
               unnest(c(data, pred_points)) %>% 
               filter(between(temp_change, -3, 3)) %>% 
               mutate(pred_points = pred_points + rnorm(n = n(), 
                                                        0, 0.2)) %>% 
               mutate(zone = factor(zone, levels = c("High", 
                                                     "Mid", 
                                                     "Low"))), 
             alpha = 0.1) +
  geom_line(aes(colour = zone, 
                alpha = line_range, 
                group = interaction(zone, temp))) +
  scale_alpha_identity() +
  geom_ribbon(aes(ymin = pred_low, 
                  ymax = pred_high, 
                  group = interaction(zone, temp), 
                  fill = zone), 
              alpha = 0.1) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
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
  




# same for extirpation ----------------------------------------------------

# calculate local extinction (extirpation) per bin and zone
dat_ext <- dat_spec %>%
  mutate(across(-c(bin, zone), 
                ~ if_else(.x > 0, 1, 0) %>% 
                  as.integer)) %>% 
  pivot_wider(names_from = zone, 
              values_from = -c(bin, zone)) %>% 
  mutate(across(-bin, ~ .x - lead(.x))) %>% 
  pivot_longer(cols = -bin, 
               names_to = c("species", "zone"),
               names_sep = "_") %>% 
  mutate(ext_signal = if_else(value == -1, 1, 0)) %>% 
  # merge with temperature
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
              group_by(zone, bin) %>%
              summarise(mean_temp = mean(temp_surface)) %>%
              mutate(temp_change = mean_temp - lead(mean_temp)) %>%
              ungroup() %>% 
              mutate(temp = if_else(temp_change >= 0, "warm", "cool")))
  

# model extirpation risk
dat_ext_model <- dat_ext %>% 
  group_by(zone, temp) %>% 
  nest() %>% 
  ungroup() %>% 
  drop_na(temp) %>% 
  mutate(mod_risk = map(data, 
                        ~ glm(ext_signal ~ temp_change, 
                              family = "binomial", 
                              data = .x)), 
         risk_prob = map2(mod_risk, 
                          data, 
                          ~ predict(.x, .y, 
                                    type = "response"))) %>% 
  unnest(risk_prob) %>% 
  arrange(zone) %>% 
  mutate(risk_prob = risk_prob * 100) %>% 
  mutate(zone = factor(zone, levels = c("High", 
                                        "Mid", 
                                        "Low")))
  
# visualise
plot_ext <- dat_ext_model %>%
  ggplot(aes(temp, risk_prob)) +
  geom_boxplot(aes(colour = zone), 
               outlier.alpha = 0.04, 
               outlier.stroke = 0) +
  coord_cartesian(ylim = c(0, 6)) +
  scale_colour_manual(values = alpha(c(colour_lavender,
                                 colour_brown,
                                 colour_green), 0.6)) +
  scale_x_discrete(labels = c("Cooling", "Warming")) +
  labs(y = "Extirpation Probability [%]", 
       x = NULL) +
  theme(legend.position = "none", 
        axis.ticks = element_blank()) +
  facet_wrap(~ zone)
  







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
  # calculate species richness
  group_by(bin, zone) %>% 
  distinct(species) %>% 
  count(name = "sp_rich") %>% 
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


# model diversity over zone and temperature change
dat_model_div <- dat_diversity %>% 
  mutate(temp = if_else(temp_change >= 0, "warm", "cool")) %>% 
  group_by(zone, temp) %>% 
  nest() %>% 
  drop_na(temp) %>% 
  ungroup() %>% 
  mutate(mod_gam = map(data, 
                       ~ lm(sp_rich ~ temp_change, 
                            data = .x)), 
         dat_pred = if_else(temp == "warm",
                            list(tibble(temp_change = seq(0, 3, by = 0.1))),
                            list(tibble(temp_change = seq(0, -3, by = -0.1)))), 
         pred_div = map2(mod_gam, 
                              dat_pred, 
                              ~ predict(.x, 
                                        newdata = .y, 
                                        se.fit = TRUE)), 
         pred_mean = map(pred_div, 
                         ~ pluck(.x, "fit")), 
         pred_se = map(pred_div, 
                       ~ pluck(.x, "se.fit")))  


# visualise
plot_richness <- dat_model_div %>%
  unnest(c(dat_pred, pred_mean, pred_se)) %>% 
  mutate(pred_low = pred_mean - 1.96 * pred_se, 
         pred_high = pred_mean + 1.96 * pred_se) %>% 
  mutate(zone = factor(zone, levels = c("High", 
                                        "Mid", 
                                        "Low"))) %>% 
  group_by(zone) %>% 
  mutate(data = map(data, ~ filter(.x, between(.x$temp_change, -3, 3))), 
         max_cooling = map_dbl(data, ~ min(.x$temp_change)), 
         max_warming =  map_dbl(data, ~ max(.x$temp_change)), 
         line_range = if_else(between(temp_change, 
                                      max_cooling, 
                                      max_warming), 
                              1, 
                              0)) %>% 
  ungroup() %>% 
  ggplot(aes(temp_change, pred_mean)) +
  geom_vline(xintercept = 0, colour = "grey80", 
             linetype = "dotted") +
  geom_point(aes(y = sp_rich, 
                 colour = zone), 
             data = dat_diversity %>% 
               mutate(zone = factor(zone, levels = c("High", 
                                                     "Mid", 
                                                     "Low"))), 
             alpha = 0.1) +
  geom_line(aes(colour = zone, 
                alpha = line_range,
                group = interaction(zone, temp))) +
  geom_ribbon(aes(ymin = pred_low, 
                  ymax = pred_high, 
                  group = interaction(zone, temp), 
                  fill = zone), 
              alpha = 0.1) +
  scale_alpha_identity() +
  scale_x_continuous(breaks = c(-2, 0, 2), 
                     limits = c(-3, 3)) +
  scale_y_continuous(limits = c(0, 35)) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                     colour_brown,
                                     colour_green), 1)) +
  scale_colour_manual(values = c(colour_lavender,
                                 colour_brown,
                                 colour_green)) +
  labs(y = "Species Richness", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  facet_wrap(~ zone) +
  theme(legend.position = "none", 
        axis.ticks = element_blank())






# patch together and save -------------------------------------------------

# combine plots
plot_final <- plot_turnover / plot_richness / plot_ext +
  plot_annotation(tag_levels = "a")


# save
ggsave(plot_final, filename = here("figures",
                                   "fig3_assemblage_trends.png"), 
       width = image_width, height = image_height*1.5, units = image_units, 
       bg = "white", device = ragg::agg_png)




