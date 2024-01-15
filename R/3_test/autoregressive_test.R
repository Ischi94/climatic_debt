library(here)
library(tidyverse)
library(patchwork)


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


# latitudinal autoregressive model ----------------------------------------


# split data into latitudinal zones and then apply autoregressive models
dat_mod <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(climatic_debt ~ temp_change + lag(temp_change),
                           data = df)
                      })
  )

# extract
dat_mod %>% 
  mutate(beta_coef = map_dbl(lm_mod, 
                             ~ coef(.x) %>% pluck(2)), 
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 2),  
         ci_high = map_dbl(ci, pluck, 4)) %>% 
  select(zone, short_term, beta_coef, ci_low, ci_high) %>%
  ungroup() %>% 
  arrange(zone, short_term) %>% 
  write_csv(here("data",
                 "beta_coefficient_per_latitude_autoregressive.csv"))


