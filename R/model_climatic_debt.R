library(here)
library(tidyverse)
library(patchwork)
library(mgcv)
library(raster)

# load data ---------------------------------------------------------------

dat_debt <- read_rds(here("data", 
                 "cleaned_debt_raster.rds"))

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 


# world map outlines
world <- map_data("world")


# fit model ---------------------------------------------------------------

# this might take a few days depending on the computational power 
model_surf <- gam(
  temp_lag ~ s(centroid_long, centroid_lat) + s(bin) +
    ti(centroid_long, centroid_lat, bin, d = c(2, 1)), 
  data = dat_debt,
  family = gaussian(),
  method = "REML", 
  control = list(nthreads = 8) # use 8 physical cores
)



# model predictions -------------------------------------------------------

# create equally sampled grid over space and time
dat_pred <-  expand.grid(
  centroid_lat = seq(-80,  80, length = 10),
  centroid_long = seq(-180, 180, length = 10),
  bin = seq(700, 4, by = -8)) 

# add predictions to actual data
dat_pred <- dat_pred %>% 
  mutate(pred_debt = predict(model_surf,
                             dat_pred, type = "response")) %>% 
  as_tibble()

# predict at higher resolution at glacial maxima and mimima

# enter glacial maxima and mimima manually
glacial <- list("min" = c(28, 156, 268, 356, 436, 532), 
                "max" = c(124, 212, 332, 412, 492, 588))

dat_pred_glacial <- expand.grid(
  centroid_lat = seq(-80,  80, length = 50),
  centroid_long = seq(-180, 180, length = 50), 
  bin = c(glacial[["min"]], glacial[["max"]])) %>% 
  mutate(pred_debt = predict(model_surf,
                             ., type = "response")) %>% 
  full_join(glacial %>% 
              as_tibble() %>% 
              pivot_longer(cols = everything(), 
                           values_to = "bin", 
                           names_to = "glacial")) %>% 
  as_tibble()


# visualisation -----------------------------------------------------------


# separate into low, mid, and high latitude
dat_pred_lat <- dat_pred %>% 
  as_tibble() %>% 
  mutate(abs_lat = abs(centroid_lat), 
         zone = case_when(
           abs_lat >= 60 ~ "high",
           between(abs_lat, 30, 60) ~ "mid", 
           between(abs_lat, 0, 30) ~ "low"))  
  
# overall model trend
plot_ov <- dat_pred_lat %>% 
  group_by(bin) %>% 
  summarise(y = quantile(pred_debt, c(0.25, 0.5, 0.75)), q = c("lower", "median", "upper")) %>%
  pivot_wider(values_from = y, names_from = q) %>% 
  ggplot(aes(bin, median, 
             ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.3) +
  geom_line(colour = "coral", lwd = 2) +
  labs(x = "Age [ka]", 
       y = "Climatic Debt [°C]") +
  scale_x_reverse() +
  theme_minimal() +
  theme(panel.grid = element_blank())



# global temperature through time
plot_temp <- dat_mean_temp %>% 
  ggplot(aes(bin, temp_ym_0m)) +
  geom_line() +
  geom_point(data = dat_mean_temp %>% 
               filter(bin %in% glacial[["min"]]), 
             fill = "steelblue", shape = 21, 
             size = 3) +
  geom_point(data = dat_mean_temp %>% 
               filter(bin %in% glacial[["max"]]), 
             fill = "firebrick", shape = 21, 
             size = 3) +
  labs(x = "Age [ka]", 
       y = "Average Global\nTemperature [°C]") +
  scale_x_reverse() +
  coord_cartesian(xlim = c(700, 0)) +
  theme_minimal() +
  theme(panel.grid = element_blank())

 
# plot_glacial
dat_pred_glacial %>% 
  select(-bin) %>% 
  pivot_wider(values_from = pred_debt, 
              names_from = glacial, 
              values_fn = list) %>% 
  mutate(debt_diff = map2(min, max, ~ .x - .y)) %>% 
  select(-c(min, max)) %>% 
  unnest(debt_diff) %>% 
  ggplot() +
  scale_x_continuous(name = '', limits = c(-180, 180), expand = c(0,0)) +
  scale_y_continuous(name = '', limits = c(-90, 90),   expand = c(0,0)) +
  geom_raster(aes(x = centroid_long, y = centroid_lat, fill = debt_diff), 
              na.rm = TRUE) +
  geom_map(aes(map_id = region), 
           data = world, map = world, fill = "grey40") + 
  scale_fill_gradient2(low = "steelblue", mid = "grey85", 
                       high = "darkred", midpoint = 0) +
  coord_quickmap() +
  theme_minimal() +
  theme(panel.grid = element_blank())







ggplot(aes(bin, pred_debt, group = interaction(centroid_lat, 
                                                 centroid_long),
             colour = zone)) +
  geom_hline(yintercept = 0) +
  geom_line(alpha = 0.3) +
  scale_x_reverse() +
  theme_minimal() +
  facet_wrap(~zone) +
  theme(legend.position = "none", 
        panel.grid = element_blank())


dat_pred %>% 
  as_tibble() %>% 
  mutate(abs_lat = abs(centroid_lat), 
         zone = case_when(
           abs_lat >= 60 ~ "pole",
           between(abs_lat, 30, 60) ~ "mid", 
           between(abs_lat, 0, 30) ~ "equator"
         )) %>% 
  group_by(zone) %>%
  summarise(mean(pred_debt), 
            sd(pred_debt)) 
