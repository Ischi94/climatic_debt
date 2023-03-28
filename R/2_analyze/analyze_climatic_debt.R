library(here)
library(tidyverse)
library(patchwork)
library(mgcv)
library(gratia)


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

# prepare trend data
dat_trends <- dat_debt %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)


# model climatic debt versus temperature ----------------------------------


# fit overall generalized additive model
mod1 <- gam(climatic_debt ~ s(temp_change, bs = "cs"),
           data = dat_trends)

# estimate the slope from the model
derivatives(mod1, 
            newdata = tibble(temp_change = seq(-0.5, 0.5, by = 0.1))) 

# accordingly, the change points from a positive to a negative slope are at
# -0.3 and 0.3

# get predictions
dat_pred <- predict(mod1, se.fit = TRUE) 
dat_pred <- tibble(pred_debt = dat_pred[[1]]) %>% 
  bind_cols(dat_trends)

# latitudinal wise
# split data into latitudinal zones and then apply fixed effect models
dat_mod <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(climatic_debt ~ temp_change,
                           data = df)
                      })
  )

# summarize the beta coefficient and save in csv
dat_mod %>%
  mutate(beta_coef = map_dbl(lm_mod, 
                             ~ coef(.x) %>% pluck(2)), 
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 2),  
         ci_high = map_dbl(ci, pluck, 4)) %>% 
  select(zone, short_term, beta_coef, ci_low, ci_high) %>%
  write_csv(here("data",
                 "beta_coefficient_per_latitude.csv"))

# same for the intercept
dat_mod %>%
  mutate(beta_coef = map_dbl(lm_mod, 
                             ~ coef(.x) %>% pluck(1)), 
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 1),  
         ci_high = map_dbl(ci, pluck, 3)) %>% 
  select(zone, short_term, beta_coef, ci_low, ci_high) %>%
  write_csv(here("data",
                 "beta_intercept_per_latitude.csv"))


new_data_lat <- dat_mod %>%
  mutate(new_data = if_else(short_term == "warming", 
                            list(tibble(temp_change = seq(0, 3, by = 0.1), 
                                        bin = 0)), 
                            list(tibble(temp_change = seq(-3, 0, by = 0.1), 
                                        bin = 0))), 
         predicted_debt = map2(.x = lm_mod,
                               .y = new_data, 
                               .f = ~ predict(.x,
                                              newdata = .y,
                                              interval = "confidence") %>% 
                                 as_tibble)) %>% 
  select(-c(data, lm_mod)) %>% 
  unnest(cols = c(predicted_debt, new_data)) %>% 
  rename(predicted_debt = fit)





# visualize ---------------------------------------------------------------



# plot temperature anomaly versus climatic debt

# in total
plot_trends_total <- dat_trends %>%
  ggplot(aes(temp_change, climatic_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_line(data = tibble(climatic_debt = seq(-3, 3,
                                              length.out = 20),
                          temp_change = seq(-3, 3,
                                            length.out = 20)),
            colour = "grey40") +
  geom_hline(yintercept = 0, 
             colour = "grey40") +
  geom_smooth(colour = colour_coral, 
              linewidth = 0.6, 
              alpha = 0.2,
              fill = colour_coral) +
  annotate(geom = "text", 
           x = 1.9, y = 2.25, 
           label = "No response", 
           angle = 30, 
           colour = "grey20",
           size = 10/.pt) +
  annotate(geom = "text", 
           x = -1.4, y = 0.25, 
           label = "Equilibrium", 
           colour = "grey20",
           size = 10/.pt) +
  annotate(geom = "curve",
           x = -0.25, xend = 0.8,
           y = -0.9, yend = -2.6,
           curvature = 0.3,
           arrow = arrow(length = unit(0.05, "inch"),
                         ends = "first"),
           colour = colour_grey, lwd = 0.3) +
  annotate(geom = "text",
           x = 1.27, y = -2.6,
           label = "GAM",
           colour = colour_grey,
           size = 10/.pt) +
  labs(y = "Climatic Mismatch [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.2, 2.2)) +
  theme(axis.ticks = element_blank())


plot_trends_total


# and per latitudinal zone
plot_trends_lat <- new_data_lat %>% 
  ggplot(aes(temp_change, predicted_debt)) +
  geom_point(aes(temp_change, climatic_debt/3, 
                 colour = zone), 
             data = dat_trends, 
             alpha = 0.1, 
             size = 0.6) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_ribbon(aes(fill = zone, 
                  group = interaction(short_term, zone),
                  ymin = lwr, ymax = upr)) +
  geom_line(aes(colour = zone, 
                group = interaction(short_term, zone)), 
            lwd = 1) + 
  scale_color_manual(values = alpha(c(colour_lavender,
                                      colour_green,
                                      colour_brown), 0.6)) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                     colour_green,
                                     colour_brown), 0.15)) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.5, 2.5)) +
  labs(y = "Climatic Mismatch [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  theme(legend.position = "none") +
  theme(axis.ticks = element_blank())



plot_trends_lat



# plot world map of samples
plot_map <- dat_debt %>% 
  distinct(pal.lat, pal.long) %>% 
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low"), 
         zone = factor(zone, levels = c("High", 
                                        "Mid", 
                                        "Low"))) %>% 
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
                                     colour_brown, 
                                     colour_green), 0.4)) +
  coord_map(projection = "mollweide") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = 'top', 
        panel.grid.major = element_line(colour = "grey80"), 
        axis.ticks = element_blank())

plot_map


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





