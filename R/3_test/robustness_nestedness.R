library(here)
library(tidyverse)
library(betapart)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                         "spp_by_depth.rds")) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292))



# preprocess data ---------------------------------------------------------


# get the species data in the right format 
dat_nest <- dat_spp %>%
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
  replace(is.na(.), 0) %>% 
  mutate(across(-c(bin, zone), ~ if_else(.x > 0, 1, .x))) %>% 
  group_by(zone) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(nested = map(data, 
                      ~ .x %>% 
                        select(-bin) %>% 
                        as.data.frame() %>% 
                        beta.pair(.) %>% 
                        pluck("beta.sne") %>%
                        as.matrix() %>% 
                        .[-1, ] %>% 
                        diag()), 
         bin = map(data, 
                   ~ pull(.x, bin) %>% 
                     .[-1])) %>% 
  unnest(c(nested, bin)) %>% 
  select(-data) %>% 
  mutate(bin = bin - 8)
  
# reformat and plot
dat_spp %>%
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
  full_join(dat_nest) %>% 
  mutate(short_term = if_else(temp_change >= 0, "warming", "cooling")) %>% 
  group_by(zone, short_term) %>% 
  nest() %>% 
  ungroup() %>% 
  drop_na(short_term) %>% 
  mutate(mod_risk = map(data, 
                        ~ glm(nested ~ temp_change, 
                              family = binomial(link = "logit"), 
                              data = .x)), 
         dat_pred = if_else(short_term == "warming",
                            list(tibble(temp_change = seq(0, 3, by = 0.1))),
                            list(tibble(temp_change = seq(0, -3, by = -0.1)))), 
         pred_div = map2(mod_risk, 
                         dat_pred, 
                         ~ predict(.x, 
                                   newdata = .y, 
                                   se.fit = TRUE)), 
         pred_mean = map(pred_div, 
                         ~ pluck(.x, "fit")), 
         pred_se = map(pred_div, 
                       ~ pluck(.x, "se.fit"))) %>%  
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
  geom_line(aes(colour = zone, 
                alpha = line_range,
                group = interaction(zone, short_term))) +
  geom_ribbon(aes(ymin = pred_low, 
                  ymax = pred_high, 
                  group = interaction(zone, short_term), 
                  fill = zone), 
              alpha = 0.1) +
  scale_alpha_identity() +
  scale_x_continuous(breaks = c(-2, 0, 2), 
                     limits = c(-3, 3)) +
  scale_y_continuous(breaks = c(-15, -2.5, 10), 
                     labels = c("0", "0.5", "1")) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                     colour_brown,
                                     colour_green), 1)) +
  scale_colour_manual(values = c(colour_lavender,
                                 colour_brown,
                                 colour_green)) +
  labs(y = "Nestedness", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  facet_wrap(~ zone) +
  theme(legend.position = "none", 
        axis.ticks = element_blank())

