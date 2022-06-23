library(here)
library(tidyverse)
library(patchwork)
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


# world map outlines
world <- map_data("world")


# calculate range debt ----------------------------------------------------


# prepare trend data
dat_trends <- dat_debt %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)


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
  mutate(range_debt = abs(y_temperature) - abs(y_cti), 
         lwr = abs(ymin_temperature) - abs(ymin_cti), 
         upr = abs(ymax_temperature) - abs(ymax_cti), 
         across(range_debt:upr, ~ round(.x) %>% as.character(.x)), 
         range_debt_unc = paste0(range_debt, 
                                 " [", 
                                 lwr, ", ", 
                                 upr, "]"), 
         new_x = (y_temperature + y_cti)/ 2) %>% 
  select(zone, short_term, range_debt, range_debt_unc, new_x) %>% 
  full_join(dat_velocity) %>% 
  mutate(short_term = str_to_title(short_term), 
         zone = factor(zone, levels = c("Low", 
                                        "Mid", 
                                        "High"))) 

# model climatic debt versus temperature ----------------------------------


# fit overall model 
mod1 <- lmer(climatic_debt ~ temp_change + (1 | bin),
             data = dat_trends)

# create equal grid for inference (accounting for differential sampling in bins)
new_data <- tibble(temp_change = seq(-3, 3, by = 0.1))  

# estimate over this equal grid  
dat_trends_pred <- predict(mod1, newdata = new_data,
                           re.form = ~ 0,
                           allow.new.levels = TRUE) %>%
  as_tibble() %>% 
  add_column(new_data) %>% 
  rename(predicted_debt = value)

# summarize beta coefficient
bootMer(mod1, fixef, nsim = 1000)

  

# latitudinal wise
# split data into latitudinal zones and then apply mixed effect models,
# accounting for sampling
dat_mod <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lmer(climatic_debt ~ temp_change + (1 | bin),
                             data = df)
                      })
  )

# summarize the beta coefficient
dat_mod %>% 
  mutate(fix_eff = map(lm_mod, ~ bootMer(.x, fixef, nsim = 1000)), 
         beta_coef = map(lm_mod, fixef),
         beta_coef = map_dbl(beta_coef, pluck(2)),
         beta_coef_sd = map_dbl(fix_eff, ~ sd(.x$t[, 2])), 
         ci_low = beta_coef - 1.96 * beta_coef_sd, 
         ci_high = beta_coef + 1.96 * beta_coef_sd) %>%
  select(zone, short_term, beta_coef, ci_low, ci_high) %>% 
  write_csv(here("data", 
                 "beta_coefficient_per_latitude.csv"))


new_data_lat <- dat_mod %>%
  mutate(new_data = if_else(short_term == "warming", 
                            list(tibble(temp_change = seq(0, 3, by = 0.1))), 
                            list(tibble(temp_change = seq(-3, 0, by = 0.1)))), 
    predicted_debt = map2(.x = lm_mod,
                          .y = new_data, 
                          .f = ~ predict(.x,
                                         newdata = .y,
                                         re.form = ~ 0,
                                         allow.new.levels = TRUE))) %>% 
  select(-c(data, lm_mod)) %>% 
  unnest(cols = c(predicted_debt, new_data))





# visualize ---------------------------------------------------------------



# plot temperature anomaly versus climatic debt

# in total
plot_trends_total <- dat_trends_pred %>% 
  ggplot(aes(temp_change, predicted_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_line(data = tibble(predicted_debt = seq(-3 - abs(mean(dat_trends_pred$predicted_debt)),
                                               3 - abs(mean(dat_trends_pred$predicted_debt)),
                                               length.out = 20),
                          temp_change = seq(-3, 3,
                                            length.out = 20)),
            colour = "grey40") +
  geom_hline(yintercept = 0 - abs(mean(dat_trends_pred$predicted_debt)), 
             colour = "grey40") +
  geom_line(colour = colour_coral, 
              lwd = 1, 
              alpha = 0.8) +
  annotate(geom = "text", 
           x = 1.4, y = 1.7, 
           label = "No response", 
           angle = 30, 
           colour = "grey20",
           size = 10/.pt) +
  annotate(geom = "text", 
           x = -1.4, y = 0.17, 
           label = "Equilibrium", 
           colour = "grey20",
           size = 10/.pt) +
  labs(y = "Climatic Debt [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.5, 2.5)) 


plot_trends_total


# and per latitudinal zone
plot_trends_lat <- new_data_lat %>% 
  ggplot(aes(temp_change, predicted_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_line(aes(colour = zone, 
                  group = interaction(short_term, zone)), 
              lwd = 1) + 
  scale_color_manual(values = alpha(c(colour_lavender,
                                colour_green,
                                colour_brown), 0.8)) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                      colour_green,
                                      colour_brown), 0.3)) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.5, 2.5)) +
  labs(y = "Climatic Debt [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  theme(legend.position = "none")


plot_trends_lat




# plot range debt as emerging from different velocities
plot_velocity <- dat_velocity_debt %>%
  ggplot(aes(y = short_term, x = y, 
           xmin = ymin, xmax = ymax, 
           fill = zone, colour = type)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_linerange(position = position_dodge2(width = 1), 
                 lwd = 0.7) +
  geom_point(position = position_dodge2(width = 1), 
             shape = 21, stroke = 0.5, size = 3) +
  geom_text(aes(label = range_debt, x = new_x),
            colour = colour_coral,
            alpha = 0.7,
            position = position_dodge(width = 1), 
            size = 10/.pt) +
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
           x = 1000, y = 2.2, 
           label = "Range Debt", 
           colour = colour_coral, 
           alpha = 0.7, 
           size = 10/.pt) +
  annotate(geom = "text", 
           x = 1000, y = 2.04, 
           label = "Realized Velocity", 
           colour = "grey30",
           size = 10/.pt) +
  annotate(geom = "curve", 
           x = -1600, xend = -1300, 
           y = 1.64, yend = 1, 
           curvature = 0.5, 
           arrow = arrow(length = unit(0.05, "inch"), 
                         ends = "first"), 
           colour = colour_grey, lwd = 0.3) +
  annotate(geom = "text", 
           x = -1000, y = 1, 
           label = "Range Debt\nin km per 8ka", 
           colour = colour_grey, 
           size = 10/.pt) +
  labs(y = NULL, 
       x = "Poleward Range Velocity [km/8ka]") +
  scale_fill_manual(values = c(colour_green,
                               colour_brown,
                               colour_lavender)) +
  scale_colour_manual(values = c("grey30", "grey70")) +
  scale_x_continuous(breaks = c(-1000, 0, 1000)) +
  theme(legend.position = "none")


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
  coord_map(projection = "mollweide") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = 'top', 
        panel.grid.major = element_line(colour = "grey80"))



# combine and save --------------------------------------------------------

# patch plots together
plot_final <- plot_map / (plot_trends_total + plot_trends_lat) +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1, 1))


# save plot_final
ggsave(plot_final, filename = here("figures", 
                                   "fig2_debt_trends.png"), 
       width = image_width, height = image_height*1.5, units = image_units, 
       bg = "white", device = ragg::agg_png)


# save range velocity plot
ggsave(plot_velocity, filename = here("figures", 
                                   "fig3_range_debt.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)


  