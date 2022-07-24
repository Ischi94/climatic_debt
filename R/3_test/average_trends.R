library(here)
library(tidyverse)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))




# load data ---------------------------------------------------------------

# debt data
dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  # convert temperature in temperature anomaly
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) / 
           sd(temp_surface, na.rm = TRUE))



# average rate of temperature change --------------------------------------

# calculate average rate of temp change per zone
dat_temp_time <- dat_debt %>% 
  group_by(zone, bin) %>% 
  summarise(temp = mean(temp_surface)) %>% 
  pivot_wider(names_from = zone, values_from = temp) %>% 
  ungroup() %>% 
  # calculate short-term changes
  mutate(across(High:Mid, ~ .x - lead(.x, default = mean(.x)))) %>% 
  pivot_longer(cols = -bin, names_to = "zone", 
               values_to = "temp_rate") %>% 
  drop_na() %>% 
  mutate(temp_type = if_else(temp_rate > 0, "Warming", "Cooling"), 
         zone = factor(zone, levels = c("High", "Mid", "Low"))) 

# plot
plot_temp <- dat_temp_time %>% 
  ggplot(aes(zone, temp_rate, fill = temp_type)) +
  geom_hline(yintercept = 0, colour = colour_grey, 
             linetype = "dashed") +
  geom_boxplot(colour = "grey20") +
  geom_label(aes(label = median_temp), 
             data = dat_temp_time %>%
               group_by(zone, temp_type) %>%
               summarise(median_temp = median(temp_rate), 
                         sd_clim = sd(temp_rate)) %>% 
               mutate(temp_rate = median_temp, 
                      median_temp = round(median_temp, 2) %>% as.character()), 
             label.size = 0, 
             fill = "white", 
             colour = "grey20",
             alpha = 0.8,
             size = 10/.pt,
             label.padding = unit(0.1, "lines"),
             position = position_nudge(x = c(-0.2, 0.2))) +
  scale_fill_manual(values = alpha(c("steelblue", 
                                     "darkred"), 0.2), 
                    name = NULL) +
  scale_y_continuous(breaks = seq(-4, 6, 2)) +
  labs(x = "Latitude", 
       y = "Rate of temperature\nchange [°C/8ka]") +
  theme(legend.position = c(0.8, 0.8))


# average climatic lag --------------------------------------

# calculate average lag
dat_lag_time <- dat_debt %>% 
  group_by(zone, bin) %>% 
  summarise(lag = mean(climatic_debt)) %>% 
  pivot_wider(names_from = zone, values_from = lag) %>% 
  ungroup() %>% 
  # calculate short-term changes
  mutate(across(High:Mid, ~ .x - lead(.x, default = mean(.x)))) %>% 
  pivot_longer(cols = -bin, names_to = "zone", 
               values_to = "lag") %>% 
  drop_na() %>% 
  mutate(zone = factor(zone, levels = c("High", "Mid", "Low"))) %>% 
  left_join(dat_temp_time) 
  
plot_lag <- dat_lag_time %>% 
  ggplot(aes(zone, lag, fill = temp_type)) +
  geom_hline(yintercept = 0, colour = colour_grey, 
             linetype = "dashed") +
  geom_boxplot(colour = "grey20") +
  geom_label(aes(label = median_lag), 
             data = dat_lag_time %>%
               group_by(zone, temp_type) %>%
               summarise(median_lag = median(lag), 
                         sd_clim = sd(lag)) %>% 
               mutate(lag = median_lag, 
                      median_lag = round(median_lag, 2) %>% as.character()), 
             label.size = 0, 
             fill = "white", 
             colour = "grey20",
             alpha = 0.8,
             size = 10/.pt,
             label.padding = unit(0.1, "lines"),
             position = position_nudge(x = c(-0.2, 0.2))) +
  scale_fill_manual(values = alpha(c("steelblue", 
                                     "darkred"), 0.2), 
                    name = NULL) +
  labs(x = "Latitude", 
       y = "Average Global\nClimatic Lag [°C/8ka]") +
  scale_y_continuous(breaks = seq(-4, 4, 2)) +
  theme(legend.position = c(0.8, 0.8))
  

# save visualisations -----------------------------------------------------

ggsave(plot_temp, filename = here("figures",
                                 "supplemental",
                                 "average_temperature.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)

ggsave(plot_lag, filename = here("figures",
                                 "supplemental",
                                 "average_lag.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)
