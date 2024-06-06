library(here)
library(merTools)
library(lme4)
library(tidyverse)



# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# debt data
dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) 


# prepare trend data
dat_trends <- read_rds(here("data",
                            "trend_data.rds"))


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


# same for spatial sampling over time
plot_heatmap <- dat_trends %>% 
  count(bin, zone) %>% 
  mutate(zone = ordered(zone, 
                        levels = c("Low", "Mid", "High"))) %>% 
  ggplot(aes(bin, zone, fill = log1p(n))) +
  geom_tile() +
  scale_fill_gradient(low = colour_lavender, 
                      high = colour_coral, 
                      breaks = c(2, 4, 6), 
                      labels = round(expm1(c(2, 4, 6))), 
                      name = "# Assemblages") +
  scale_x_reverse() +
  labs(x = "Age [ka]", 
       y = "Latitudinal zone")

# save plot
ggsave(plot_heatmap, filename = here("figures", 
                                        "supplemental",
                                        "sampling_heatmap.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)
