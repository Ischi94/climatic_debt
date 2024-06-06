library(here)
library(tidyverse)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))




# load data ---------------------------------------------------------------

# debt data
dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) 

# world map outlines
world <- map_data("world")

# prepare trend data
dat_trends <- read_rds(here("data",
                            "trend_data.rds"))


# models ------------------------------------------------------------------


# split data into latitudinal zones and then apply fixed effect models
dat_mod <- dat_debt %>% 
  # filter out temperature where temp changes
  # are not associated with a response 
  filter(between(temp_change, -0.3, 0.3)) %>%
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
  rename_with(~paste0(.x, "_new")) %>% 
  ungroup() %>% 
  arrange(short_term_new, zone_new) %>% 
  bind_cols(read_csv(here("data",
                          "beta_coefficient_per_latitude.csv")) %>% 
              arrange(short_term, zone)) %>% 
  write(here("data", 
             "beta_coefficient_per_latitude_above_03.csv"))


