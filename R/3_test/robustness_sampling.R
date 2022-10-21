library(here)
library(merTools)
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


# reported beta coefficients based on fixed effects
dat_reported <- read_csv(here("data",
                              "beta_coefficient_per_latitude.csv")) %>% 
  add_column(type = "Fixed Effects")




# random effects ----------------------------------------------------------


# latitudinal wise
# split data into latitudinal zones and then apply mixed effect models,
# accounting for sampling
dat_re <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lmer(climatic_debt ~ temp_change + (1 | bin),
                             data = df)
                      }), 
         beta_coef_boot = map(.x = lm_mod,
                         .f = ~ bootMer(.x, fixef, 
                                        nsim = 1000)), 
         beta_coef = map(.x = beta_coef_boot,
                         .f = ~ pluck(.x, 2) %>%
                           as_tibble() %>%
                           summarise(median_hilow(temp_change)))
  ) %>% 
  select(zone, short_term, beta_coef) %>% 
  unnest(beta_coef) %>% 
  rename(beta_coef = y, ci_low = ymin, ci_high = ymax) %>% 
  add_column(type = "Random Effects")





# omiting low sampled time bins -------------------------------------------

# remove assemblages older than 400ka
dat_omit_old <- dat_trends %>%
  filter(bin <= 400) %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(climatic_debt ~ temp_change,
                           data = df)
                      })
         ) %>% 
  mutate(beta_coef = map_dbl(lm_mod, 
                             ~ coef(.x) %>% pluck(2)), 
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 2),  
         ci_high = map_dbl(ci, pluck, 4)) %>% 
  select(zone, short_term, beta_coef, ci_low, ci_high) %>% 
  add_column(type = "Younger than 400ka")




# omiting highest sampled time bin -------------------------------------------


# remove assemblages from the most recent bin (with the highest sampling
# intensity)
dat_omit_young <- dat_trends %>%
  filter(bin > 4) %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(climatic_debt ~ temp_change,
                           data = df)
                      })
  ) %>% 
  mutate(beta_coef = map_dbl(lm_mod, 
                             ~ coef(.x) %>% pluck(2)), 
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 2),  
         ci_high = map_dbl(ci, pluck, 4)) %>% 
  select(zone, short_term, beta_coef, ci_low, ci_high) %>% 
  add_column(type = "Older than 4ka")





# visualize ---------------------------------------------------------------

# combine dataframes into one
dat_comp <- dat_reported %>% 
  full_join(dat_re) %>% 
  full_join(dat_omit_old) %>%
  full_join(dat_omit_young) 



# plot
plot_comparison <- dat_comp %>% 
  mutate(zone = ordered(zone, 
                        levels = c("High", "Mid", "Low")), 
         short_term = str_to_title(short_term), 
         type = ordered(type, 
                        levels = c("Fixed Effects", 
                                   "Younger than 400ka", 
                                   "Older than 4ka",
                                   "Random Effects"))) %>% 
  ggplot(aes(beta_coef, short_term,
             colour = type)) +
  geom_vline(xintercept = 0) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.6) +
  geom_point(position = position_dodge(width = 0.5), 
             size = 2) +
  facet_wrap(~ zone, 
             ncol = 1) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6)) +
  scale_color_brewer(name = NULL,
                     type = "qual", 
                     palette = 6) +
  labs(y = NULL, 
       x = "Beta Coefficient") +
  theme(legend.position = "bottom", 
        panel.border = element_rect(fill = NA,
                                    colour = "grey60"), 
        panel.grid.major.x = element_line(colour = "grey60", 
                                          linetype = "dotted"))

# save plot
ggsave(plot_comparison, filename = here("figures", 
                                        "supplemental",
                                        "sampling_intensity_test.png"), 
       width = image_width, height = image_height*1.5, units = image_units, 
       bg = "white", device = ragg::agg_png)

























dat_trends %>% 
  ggplot(aes(temp_change, climatic_debt)) +
  geom_hline(yintercept = c(-1, 1)) +
  geom_smooth()

# fit overall model 
mod1 <- lm(climatic_debt ~ temp_change, data = dat_trends)

# create equal grid for inference (accounting for differential sampling in bins)
new_data <- tibble(temp_change = seq(-3, 3, by = 0.1), 
                   bin = 0)  

# estimate over this equal grid  
dat_trends_pred <- predict(mod1,
                           newdata = new_data) %>%
  as_tibble() %>% 
  add_column(new_data) %>% 
  rename(predicted_debt = value)


dat_trends_pred %>% 
  ggplot(aes(temp_change, predicted_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_line(data = tibble(predicted_debt = seq(-3 ,
                                               3,
                                               length.out = 20),
                          temp_change = seq(-3, 3,
                                            length.out = 20)),
            colour = "grey40") +
  geom_hline(yintercept = 0, 
             colour = "grey40") +
  geom_line(colour = colour_coral, 
            lwd = 1, 
            alpha = 0.8) +
  labs(y = "Climatic Lag [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.5, 2.5)) +
  theme(axis.ticks = element_blank())

dat_trends_pred %>% 
  filter(between(predicted_debt, -1, 1)) %>% 
  summary()

dat_mod <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(climatic_debt ~ temp_change,
                             data = df)
                      })
  )



dat_mod %>%
  mutate(new_data = if_else(short_term == "warming", 
                            list(tibble(temp_change = seq(0, 3, by = 0.1), 
                                        bin = 0)), 
                            list(tibble(temp_change = seq(-3, 0, by = 0.1), 
                                        bin = 0))), 
         predicted_debt = map2(.x = lm_mod,
                               .y = new_data, 
                               .f = ~ predict(.x,
                                              newdata = .y,
                                              interval =  "confidence"))) %>% 
  select(-c(data, lm_mod)) %>% 
  mutate(predicted_debt = map(predicted_debt, as_tibble)) %>% 
  unnest(cols = c(predicted_debt, new_data)) %>% 
  rename(predicted_debt = fit) %>% 
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
  drop_na(short_term) %>% 
  filter(bin != 4)


dat_trends %>% 
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
  geom_smooth(method = "lm", 
              colour = "grey20", 
              fill = "grey20", 
              alpha = 0.15) +
  geom_smooth(alpha = 0, colour = "grey60", 
              linetype = "dotted") 
