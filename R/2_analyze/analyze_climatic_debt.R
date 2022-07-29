library(here)
library(merTools)
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

# prepare trend data
dat_trends <- dat_debt %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)


# model climatic debt versus temperature ----------------------------------


# fit overall model 
mod1 <- lmer(climatic_debt ~ temp_change + (1 | bin),
             data = dat_trends)

# create equal grid for inference (accounting for differential sampling in bins)
new_data <- tibble(temp_change = seq(-3, 3, by = 0.1), 
                   bin = 0)  

# estimate over this equal grid  
dat_trends_pred <- predictInterval(mod1,
                                   newdata = new_data, 
                                   which = "full", 
                                   include.resid.var = FALSE) %>%
  as_tibble() %>% 
  add_column(new_data) %>% 
  rename(predicted_debt = fit)

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

# # summarize the beta coefficient and save in csv
# dat_mod %>%
#   mutate(fix_eff = map(lm_mod, ~ bootMer(.x, fixef, nsim = 1000)),
#          beta_coef = map(lm_mod, fixef),
#          beta_coef = map_dbl(beta_coef, pluck(2)),
#          beta_coef_sd = map_dbl(fix_eff, ~ sd(.x$t[, 2])),
#          ci_low = beta_coef - 1.96 * beta_coef_sd,
#          ci_high = beta_coef + 1.96 * beta_coef_sd) %>%
#   select(zone, short_term, beta_coef, ci_low, ci_high) %>%
#   write_csv(here("data",
#                  "beta_coefficient_per_latitude.csv"))
# 
# # same for the intercept
# dat_mod %>%
#   mutate(fix_eff = map(lm_mod, ~ bootMer(.x, fixef, nsim = 1000)),
#          intercept = map(lm_mod, fixef),
#          intercept = map_dbl(intercept, pluck(1)),
#          intercept_sd = map_dbl(fix_eff, ~ sd(.x$t[, 1])),
#          ci_low = intercept - 1.96 * intercept_sd,
#          ci_high = intercept + 1.96 * intercept_sd) %>%
#   select(zone, short_term, intercept, ci_low, ci_high) %>%
#   write_csv(here("data",
#                  "intercept_per_latitude.csv"))


new_data_lat <- dat_mod %>%
  mutate(new_data = if_else(short_term == "warming", 
                            list(tibble(temp_change = seq(0, 3, by = 0.1), 
                                        bin = 0)), 
                            list(tibble(temp_change = seq(-3, 0, by = 0.1), 
                                        bin = 0))), 
         predicted_debt = map2(.x = lm_mod,
                               .y = new_data, 
                               .f = ~ predictInterval(.x,
                                                      newdata = .y, 
                                                      which = "full", 
                                                      include.resid.var = FALSE))) %>% 
  select(-c(data, lm_mod)) %>% 
  unnest(cols = c(predicted_debt, new_data)) %>% 
  rename(predicted_debt = fit)





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
  geom_ribbon(aes(temp_change, ymin = lwr, ymax = upr), 
              fill = colour_coral,
              alpha = 0.2) +
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
  labs(y = "Climatic Lag [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.5, 2.5)) +
  theme(axis.ticks = element_blank())


plot_trends_total


# and per latitudinal zone
plot_trends_lat <- new_data_lat %>% 
  ggplot(aes(temp_change, predicted_debt)) +
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
                                colour_brown), 0.8)) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                      colour_green,
                                      colour_brown), 0.2)) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.5, 2.5)) +
  labs(y = "Climatic Lag [°C]", 
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
                                     colour_green), 0.5)) +
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





  