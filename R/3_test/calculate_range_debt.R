library(here)
library(tidyverse)
library(lme4)


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



# prepare trend data
dat_trends <- dat_debt %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)


# calculate range debt ----------------------------------------------------


# calculate distance to equator of each observation
dat_dist_eq <- dat_debt  %>% 
  # convert from decimal degrees to radians
  mutate(across(c(pal.lat, pal.long), ~ .x*pi/180)) %>% 
  # calculate distance to equator in km / the earth mean radius is 6371 km 
  # this the great circle distance based on spherical law of cosines
  mutate(dist_eq = acos(sin(0)*sin(pal.lat) + cos(0)*cos(pal.lat) * cos(pal.long-pal.long)) * 6371)



# calculate the rate of change in community composition in response to climate
# change (either warming or cooling) in (°C/8ka)
# we model it as a mixed effect to account for differential sampling in bins
mod_cti <- dat_trends %>% 
  mutate(cti_change = cti - lead(cti, default = mean(cti))) %>% 
  lmer(cti_change ~ 0 + zone:short_term + (1 | bin),
       data = .)

# extract point estimate and confidence intervals
dat_cti_per_time <- tibble(zone = rep(c("High", "Low", "Mid"), 2), 
                           short_term = c(rep("cooling", 3), rep("warming", 3))) %>% 
  mutate(y = fixef(mod_cti), 
         y_sd = bootMer(mod_cti, fixef, nsim = 1000)$t %>% apply(., 2, sd), 
         ymin = y - 1.96 * y_sd, 
         ymax = y + 1.96 * y_sd) %>% 
  arrange(zone) %>% 
  select(-y_sd)



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
  mutate(across(y:ymax, ~ .x/cti_change)) %>% 
  add_column(type = "cti")


# repeat for temperature
# calculate the rate of change in temperature in response to climate
# change (either warming or cooling) in (°C/8ka)
mod_temp <- lmer(temp_change ~ 0 + zone:short_term + (1 | bin),
                 data = dat_trends)

# extract point estimate and CI's
dat_temp_per_time <- tibble(zone = rep(c("High", "Low", "Mid"), 2), 
                            short_term = c(rep("cooling", 3), rep("warming", 3))) %>% 
  mutate(y = fixef(mod_temp), 
         y_sd = bootMer(mod_temp, fixef, nsim = 1000)$t %>% apply(., 2, sd), 
         ymin = y - 1.96 * y_sd, 
         ymax = y + 1.96 * y_sd) %>% 
  arrange(zone) %>% 
  select(-y_sd)


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
  mutate(across(y:ymax, ~ .x/temp_change)) %>% 
  add_column(type = "temperature")


# combine with cti velocity
dat_velocity <- dat_temp_velocity %>% 
  full_join(dat_cti_velocity) %>% 
  ungroup()

# velocity change aka range debt
dat_velocity_debt <- dat_velocity %>% 
  # calculate difference between what is needed (temperature) and the actual
  # shift (cti)
  pivot_wider(names_from = type, values_from = y:ymax) %>% 
  mutate(range_debt = y_temperature - y_cti, 
         lwr = ymin_temperature - ymin_cti, 
         upr = ymax_temperature - ymax_cti, 
         across(range_debt:upr, ~ round(.x) %>% 
                  as.character()), 
         range_debt_unc = paste0(" [", 
                                 lwr, ", ", 
                                 upr, "]"), 
         new_x = (y_temperature + y_cti)/ 2, 
         new_x_unc = if_else(new_x < 0, 
                             new_x + 250, 
                             new_x + 210), 
         new_x_mean = new_x - 100, 
         new_x_mean = if_else(zone == "Mid" & short_term == "warming", 
                              new_x, 
                              new_x_mean)) %>% 
  select(zone, short_term, range_debt, range_debt_unc, 
         new_x_unc, new_x, new_x_mean) %>% 
  full_join(dat_velocity) %>% 
  mutate(short_term = str_to_title(short_term), 
         zone = factor(zone, levels = c("Low", 
                                        "Mid", 
                                        "High"))) 


# visualize ---------------------------------------------------------------


# plot range debt as emerging from different velocities
plot_velocity <- dat_velocity_debt %>%
  ggplot(aes(y = short_term, x = y, 
             xmin = ymin, xmax = ymax, 
             fill = zone, colour = type)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_label(aes(label = range_debt_unc, x = new_x_unc, 
                 group = zone),
             colour = colour_grey,
             fill = "white", 
             label.size = 0,
             position = position_dodge(width = 1), 
             size = 10/.pt) +
  geom_text(aes(label = range_debt, x = new_x_mean),
            colour = colour_coral,
            alpha = 0.7,
            position = position_dodge(width = 1), 
            size = 10/.pt) +
  geom_linerange(position = position_dodge2(width = 1), 
                 lwd = 0.7) +
  geom_point(position = position_dodge2(width = 1), 
             shape = 21, stroke = 0.5, size = 3) +
  geom_hline(yintercept = 1.5, colour = colour_grey, 
             linetype = "dotted") +
  annotate(geom = "rect", 
           xmin = 500, xmax = 1500, 
           ymin = 1.9, ymax = 2.47, 
           fill = "white", 
           colour = colour_grey, lwd = 0.4) +
  annotate(geom = "text", 
           x = 740, y = 2.53, 
           label = "Interpretation", 
           colour = colour_grey, 
           size = 10/.pt) +
  annotate(geom = "text", 
           x = 1000, y = 2.36, 
           label = "Required Velocity", 
           colour = "grey70",
           size = 10/.pt) +
  annotate(geom = "text", 
           x = 810, y = 2.2, 
           label = "Range Lag", 
           colour = colour_coral, 
           alpha = 0.7, 
           size = 10/.pt) +
  annotate(geom = "text", 
           x = 1250, y = 2.2, 
           label = " [95% CI]", 
           colour = colour_grey, 
           size = 10/.pt) +
  annotate(geom = "text", 
           x = 1000, y = 2.04, 
           label = "Realized Velocity", 
           colour = "grey30",
           size = 10/.pt) +
  annotate(geom = "curve", 
           x = -1720, xend = -1300, 
           y = 1.6, yend = 1, 
           curvature = 0.5, 
           arrow = arrow(length = unit(0.05, "inch"), 
                         ends = "first"), 
           colour = colour_grey, lwd = 0.3) +
  annotate(geom = "text", 
           x = -1000, y = 1, 
           label = "Range Lag\nin km per 8ka", 
           colour = colour_grey, 
           size = 10/.pt) +
  labs(y = NULL, 
       x = "Equatorward Range Velocity [km/8ka]") +
  scale_fill_manual(values = c(colour_green,
                               colour_brown,
                               colour_lavender), 
                    name = NULL) +
  scale_colour_manual(values = c("grey30", "grey70")) +
  scale_x_continuous(breaks = c(-1000, 0, 1000)) +
  guides(colour = "none", 
         fill = guide_legend(override.aes = list(alpha = 0.6))) +
  theme(legend.position = c(-0.05, 0.01), 
        legend.text = element_text(colour = "grey70"))

plot_velocity



# save --------------------------------------------------------------------

# save range velocity plot
ggsave(plot_velocity, filename = here("figures", 
                                      "supplemental",
                                      "range_debt.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
