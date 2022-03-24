library(here)
library(tidyverse)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))




# load data ---------------------------------------------------------------

dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  # cnovert temperature in temperature anomaly
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) / 
           sd(temp_surface, na.rm = TRUE))

# world map outlines
world <- map_data("world")




# transform to geographic debt -------------------------------------

# calculate distance to equator of each observation
dat_dist_eq <- dat_debt  %>% 
  # convert from decimal degrees to radians
  mutate(across(c(pal.lat, pal.long), ~ .x*pi/180)) %>% 
  # calculate distance to equator in km / the earth mean radius is 6371 km 
  # this the great circle distance based on spherical law of cosines
  mutate(dist_eq = acos(sin(0)*sin(pal.lat) + cos(0)*cos(pal.lat) * cos(pal.long-pal.long)) * 6371)



# calculate the rate of change in community composition in response to climate
# change (either warming or cooling) in (°C/8ka)
dat_cti_per_time <- dat_debt %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         cti_change = cti - lead(cti, default = mean(cti))) %>% 
  drop_na(short_term) %>% 
  group_by(short_term) %>% 
  summarise(mean_cl_boot(cti_change))


# calculate rate of change in CTI in kilometers from the equator to the poles in (°C/km)
dat_dist_eq %>% 
  ggplot(aes(y = cti, x = dist_eq)) +
  geom_point() + 
  geom_smooth(method = "lm")

cti_change <- lm(cti ~ dist_eq, data = dat_dist_eq)$coefficients[[2]]


# calculate the ratio between the temporal change and the spatial gradient in CTI 
# (°C/8ka) / (°C/km) = km/8ka
# # this is the velocity of cti
dat_cti_velocity <- dat_cti_per_time %>% 
  mutate(across(y:ymax, ~ .x/cti_change))


# repeat for temperature
# calculate the rate of change in temperature in response to climate
# change (either warming or cooling) in (°C/8ka)
dat_temp_per_time <- dat_debt %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                       default = mean(temp_anom))) %>% 
  drop_na(short_term) %>% 
  group_by(short_term) %>% 
  summarise(mean_cl_boot(temp_change))


# calculate rate of change in CTI in kilometers from the equator to the poles in (°C/km)
dat_dist_eq %>% 
  ggplot(aes(y = temp_anom, x = dist_eq)) +
  geom_point() + 
  geom_smooth(method = "lm")

temp_change <- lm(temp_anom ~ dist_eq, data = dat_dist_eq)$coefficients[[2]]


# calculate the ratio between the temporal change and the spatial gradient in CTI 
# (°C/8ka) / (°C/km) = km/8ka
# this is the velocity of cti
dat_temp_velocity <- dat_temp_per_time %>% 
  mutate(across(y:ymax, ~ .x/temp_change))


# transform to geographic range shift -------------------------------------


dat_dist_eq <- dat_debt  %>% 
  group_by(bin, zone) %>% 
  summarise(av_anom = mean(temp_anom)) %>% 
  pivot_wider(names_from = zone, values_from = av_anom) %>% 
  ungroup() %>% 
  mutate(across(High:Mid, ~ if_else(.x > lead(.x),
                                    "warming",
                                    "cooling"))) %>%
  full_join(dat_debt %>% select(-zone)) %>% 
  # convert from decimal degrees to radians
  mutate(across(c(pal.lat, pal.long), ~ .x*pi/180)) %>% 
  # calculate distance to equator in km / the earth mean radius is 6371 km 
  # this the great circle distance based on spherical law of cosines
  mutate(dist_eq = acos(sin(0)*sin(pal.lat) + cos(0)*cos(pal.lat) * cos(pal.long-pal.long)) * 6371) %>% 
  # bring into right format
  pivot_longer(cols = High:Mid, 
               names_to = "zone", 
               values_to = "short_term") %>% 
  drop_na(short_term) 



# rate of change in community composition in response to climate change (warming or cooling) (°C per 8ka)
dat_dist_eq %>% 
  ggplot(aes(bin, cti)) +
  geom_point()

cti_change <- lm(cti ~ dist_eq, data = dat_dist_eq)$coefficients[[2]]

# estimate average community temperature index 
cti_mean <- dat_dist_eq %>% 
  group_by(zone) %>% 
  summarise(mean_cl_boot(cti)) %>% 
  rename(mean_cti = y, 
         lwr = ymin, 
         upr = ymax)


# compare latitudinal
dat_dist_eq %>% 
  group_by(temp_change, zone) %>% 
  summarise(mean_cl_boot(cti)) %>% 
  ungroup() %>% 
  # calculate change to average
  full_join(cti_mean) %>% 
  mutate(mean_range = y - mean_cti,
         lwr_range = ymin - lwr,
         upr_range = ymax - upr) %>%
  mutate(across(mean_range:upr_range, ~ .x / cti_change))


# what range is actually needed

dat_dist_eq %>% 
  ggplot(aes(dist_eq, short_term)) +
  geom_point()

# estimate average polward shift based on cti
# the change in cti per 1 km further away from the equator
temp_change <- lm(dist_eq ~ temp_surface, data = dat_dist_eq)$coefficients[[2]]

# estimate average community temperature index 
temp_mean <- dat_dist_eq %>% 
  group_by(zone) %>% 
  summarise(mean_cl_boot(short_term)) %>% 
  rename(mean_temp = y, 
         lwr = ymin, 
         upr = ymax)


# compare latitudinal
dat_dist_eq %>% 
  group_by(temp_change, zone) %>% 
  summarise(mean_cl_boot(short_term)) %>% 
  ungroup() %>% 
  full_join(temp_mean) %>% 
  mutate(mean_range = y - mean_temp,
         lwr_range = ymin - lwr,
         upr_range = ymax - upr) %>%
  mutate(across(mean_range:upr_range, ~ .x * temp_change))


  # # calculate change to average
  # mutate(across(y:ymax, ~ .x - temp_mean), 
  #        across(y:ymax, ~ .x * temp_change)) 



# visualize ---------------------------------------------------------------

# plot world map of samples
plot_map <- dat_debt %>% 
  distinct(pal.lat, pal.long) %>% 
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  ggplot() +
  geom_point(aes(x = pal.long, y = pal.lat,
                 fill = zone),
             na.rm = TRUE,
             shape = 21,
             colour = colour_grey,
             stroke = 0.3,
             size = 2.5) +
  geom_map(aes(map_id = region), 
           data = world, map = world, fill = "grey40") + 
  scale_x_continuous(name = '', limits = c(-180, 180), expand = c(0,0), 
                     labels = NULL, breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(name = '', limits = c(-90, 90),   expand = c(0,0), 
                     labels = NULL) +
  scale_fill_manual(name = "Latitude", 
                    values = alpha(c(colour_lavender,
                                     colour_green,
                                     colour_brown), 0.5)) +
  theme_minimal() +
  coord_map(projection = "mollweide") +
  theme(legend.position = 'bottom')


# plot temperature anomaly versus climatic debt
plot_regr <- dat_debt %>% 
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) / 
           sd(temp_surface, na.rm = TRUE)) %>% 
  ggplot(aes(climatic_debt, temp_anom)) +
  geom_jitter(aes(fill = zone), 
              width = 0, height = 0.1, 
              shape = 21,
              colour = colour_grey,
              stroke = 0.3, 
              size = 2.5) +
  geom_smooth(aes(colour = zone),
              method = "lm", 
              lwd = 1.5, 
              fill = colour_grey) +
  scale_color_manual(values = c(colour_lavender,
                                colour_green,
                                colour_brown)) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                colour_green,
                                colour_brown), 0.5)) +
  theme_minimal() +
  theme(legend.position = "none")









plot3 / (plot1 + plot2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(2, 1, 1))
  