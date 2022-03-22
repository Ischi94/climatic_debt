library(here)
library(tidyverse)
library(patchwork)
library(mgcv)


# load data ---------------------------------------------------------------

dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds"))

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 


# world map outlines
world <- map_data("world")


# fit model ---------------------------------------------------------------

# spatio-temporal model
model_surf <- gam(
  climatic_debt ~ s(pal.long, pal.lat) + s(bin) +
    ti(pal.long, pal.lat, bin, d = c(2, 1)), 
  data = dat_debt,
  family = gaussian(),
  method = "REML", 
  control = list(nthreads = 8) # use 8 physical cores
)



# model predictions -------------------------------------------------------

# create equally sampled grid over space and time
dat_pred <-  expand.grid(
  pal.lat = seq(-80,  80, length = 10),
  pal.long = seq(-180, 180, length = 10),
  bin = seq(700, 4, by = -8)) 

# add predictions to actual data
dat_pred <- dat_pred %>% 
  mutate(pred_debt = predict(model_surf,
                             dat_pred, type = "response")) %>% 
  as_tibble()

# split into background and unusual temperatures based on quantiles
dat_temp_quant <- dat_mean_temp %>% 
  summarize(temp_quant = quantile(temp_ym_0m)) %>% 
  pull(temp_quant) 

dat_pred <- dat_mean_temp %>% 
  mutate(temp_type = case_when(
    temp_ym_0m <= dat_temp_quant[2] ~ "cool", 
    between(temp_ym_0m, 
            dat_temp_quant[2], 
            dat_temp_quant[4]) ~ "background", 
    temp_ym_0m >= dat_temp_quant[4] ~ "warm")) %>% 
  select(-temp_ym_0m) %>% 
  right_join(dat_pred)
  




# visualisation -----------------------------------------------------------


# separate into low, mid, and high latitude
dat_pred_lat <- dat_pred %>% 
  as_tibble() %>% 
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "high",
           between(abs_lat, 30, 60) ~ "mid", 
           between(abs_lat, 0, 30) ~ "low"))  
  
# overall model trend
plot_ov <- dat_pred_lat %>% 
  group_by(bin) %>% 
  summarise(mean_cl_boot(pred_debt)) %>%
  select(bin, mean = y, lower = ymin, upper = ymax)  %>% 
  ggplot(aes(bin, mean, 
             ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.3) +
  geom_line(colour = "coral", lwd = 1.3) +
  labs(x = "Age [ka]", 
       y = "Climatic Debt [°C]") +
  scale_x_reverse() +
  theme_minimal() +
  theme(panel.grid = element_blank())



# global temperature through time
plot_temp <- dat_mean_temp %>% 
  ggplot(aes(bin, temp_ym_0m)) +
  geom_line() +
  labs(x = "Age [ka]", 
       y = "Average Global\nTemperature [°C]") +
  scale_x_reverse() +
  coord_cartesian(xlim = c(700, 0)) +
  theme_minimal() +
  theme(panel.grid = element_blank())


plot_ov/plot_temp +
  plot_layout(heights = c(3, 1))


dat_pred_zone <- dat_pred %>% 
  as_tibble() %>% 
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "pole",
           between(abs_lat, 30, 60) ~ "mid", 
           between(abs_lat, 0, 30) ~ "equator"
         )) 

dat_pred_zone %>% 
  ggplot(aes(bin, pred_debt, 
             group = interaction(pal.lat,
                                 pal.long),
             colour = zone)) +
  geom_hline(yintercept = 0) +
  geom_line(alpha = 0.3) +
  scale_x_reverse() +
  facet_wrap(~ zone) +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank())

dat_pred_zone %>% 
  group_by(bin, zone) %>% 
  summarise(mean_cl_boot(pred_debt)) %>% 
  ggplot(aes(bin, y, 
             colour = zone)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = ymin,
                  ymax = ymax, 
                  colour = NULL, 
                  fill = zone), 
              alpha = 0.4) +
  geom_line(alpha = 0.8) +
  scale_x_reverse() +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank())

dat_pred_zone %>% 
  group_by(zone) %>% 
  summarise(sd(pred_debt))


# background versus maxima
dat_pred  %>% 
  # group_by(temp_type) %>% 
  summarise(mean_cl_boot(pred_debt))



dat_pred_glacial %>% 
  select(-bin) %>% 
  group_by(glacial) %>% 
  summarise(mean_cl_boot(pred_debt))
















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
  geom_raster(aes(x = pal.long, y = pal.lat, fill = debt_diff), 
              na.rm = TRUE) +
  geom_map(aes(map_id = region), 
           data = world, map = world, fill = "grey40") + 
  scale_fill_gradient2(low = "darkblue", mid = "white", 
                       high = "darkred") +
  coord_quickmap() +
  theme_minimal() +
  theme(panel.grid = element_blank())







plot4 <- dat_debt %>%
  group_by(bin) %>% 
  summarise(mean_cl_boot(climatic_debt)) %>% 
  select(bin, mean = y, lower = ymin, upper = ymax)  %>% 
  ggplot(aes(bin, mean)) +
  geom_hline(yintercept = 0) +
  geom_line(colour = "coral", lwd = 1.3) +
  labs(x = "Age [ka]", 
       y = "Climatic Debt [°C]") +
  scale_x_reverse() +
  theme_minimal() +
  theme(panel.grid = element_blank())
  

plot4/plot_temp +
  plot_layout(heights = c(2, 1))
