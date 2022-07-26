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



# prepare trend data
dat_trends <- dat_debt %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)


# model climatic debt versus temperature ----------------------------------


# latitudinal wise
# split data into latitudinal zones and then apply mixed effect models,
# accounting for sampling

dat_mod <- dat_trends %>% 
  filter(short_term == "warming") %>% 
  group_by(zone) %>%
  nest() %>% 
  mutate(lm_mod = map(
    .x = data,
    .f = function(df) {
      lmer(climatic_debt ~ temp_change + (1 | bin),
           data = df)
    }
  ))


# predict under climate change scenarios
new_data_lat <- dat_mod %>%
  mutate(one_degree = map(data, ~ mutate(.x, temp_change = temp_change + 1)),
         one_degree = map2(
           .x = lm_mod,
           .y = one_degree,
           .f = ~ predict(.x,
                          newdata = .y)), 
         two_degree = map(data, ~ mutate(.x, temp_change = temp_change + 2)),
         two_degree = map2(
           .x = lm_mod,
           .y = two_degree,
           .f = ~ predict(.x,
                          newdata = .y)),
         three_degree = map(data, ~ mutate(.x, temp_change = temp_change + 3)),
         three_degree = map2(
           .x = lm_mod,
           .y = three_degree,
           .f = ~ predict(.x,
                          newdata = .y))) %>%
  unnest(cols = c(data, one_degree, two_degree, three_degree)) %>% 
  select(-c(bin:temp_surface, cti, abs_lat:lm_mod)) %>% 
  pivot_longer(cols = one_degree:three_degree, 
               names_to = "scenario", 
               values_to = "pred_lag") %>% 
  # calculate percentage change in lag
  mutate(across(c(pred_lag, climatic_debt), abs),
         perc_change = ((pred_lag - climatic_debt) / climatic_debt)*10)



# visualize ---------------------------------------------------------------


plot_perc_change <- new_data_lat %>% 
  group_by(zone, scenario) %>% 
  summarise(mean_cl_boot(perc_change)) %>% 
  full_join(tibble(scenario = c("one_degree", 
                                "two_degree", 
                                "three_degree"), 
                   temp_label = c("+1°C", 
                                  "+2°C", 
                                  "+3°C"))) %>% 
  mutate(zone = factor(zone,
                       levels = c("High", "Mid", "Low"))) %>%
  ggplot(aes(zone, y, fill = zone)) +
  geom_hline(yintercept = 0,
             colour = "grey70",
             linetype = "dotted") +
  geom_linerange(aes(ymin = ymin, ymax = ymax), 
             colour = "grey20", 
             alpha = 0.8) + 
  geom_point(shape = 21, 
             size = 4, 
             colour = "grey30", 
             alpha = 0.8) +
  scale_fill_manual(values = c(colour_lavender,
                               colour_brown,
                               colour_green),
                    name = "Latitude") +
  facet_wrap(~ temp_label) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3)), 
         colour = "none") +
  labs(y = "Change in Climatic Lag [%]",
       x = NULL) +
  theme(legend.position = c(0.1, 0.8),
        panel.border = element_rect(fill = NA,
                                    colour = "grey90"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# save plot
ggsave(plot_perc_change, filename = here("figures",
                                         "supplemental",
                                         "percentage_change.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)
