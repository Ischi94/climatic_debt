library(here)
library(tidyverse)
library(vegan)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                         "spp_by_depth.rds")) %>% 
  # right format for WAPLS
  mutate(rel_abund = rel_abund * 100) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292))



# preprocess data ---------------------------------------------------------


# get the species data in the right format
dat_spec <- dat_spp %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  nest_by(bin, zone, DepthHabitat) %>% 
  ungroup() %>% 
  mutate(data = map(data, 
                    ~ pivot_wider(.x, 
                                  id_cols = core_uniq, 
                                  names_from = species, 
                                  values_from = rel_abund, 
                                  values_fn = mean, 
                                  values_fill = 0) %>% 
                      select(-core_uniq) %>% 
                      colSums() %>% 
                      enframe() %>% 
                      mutate(rel_abund = value/sum(value)) %>% 
                      select(-value) %>% 
                      pivot_wider(names_from = name, values_from = rel_abund, 
                                  values_fill = 0))) %>% 
  unnest(cols = c(data)) %>% 
  replace(is.na(.), 0) %>% 
  group_by(DepthHabitat) %>% 
  nest() %>% 
  ungroup()


# calculate turnover ------------------------------------------------------


# convert to matrix and calculate turnover
dat_turnover <- dat_spec %>%
  mutate(dist_mat = map(data, 
                        ~ .x %>% 
                          select(- c(bin, zone)) %>% 
                          vegdist(method = "chisq") %>% 
                          as.matrix()),
         turnover = map(dist_mat, 
                        ~ c(.x[row(.x) == col(.x) + 1], 0)), 
         data = map(data, 
                    ~ select(.x, bin, zone))) %>% 
  select(-dist_mat) %>% 
  unnest(cols = c(data, turnover))


# add temperature
dat_final <- dat_spp %>%
  # add latitudinal zones
  mutate(
    abs_lat = abs(pal.lat),
    zone = case_when(
      abs_lat >= 60 ~ "High",
      between(abs_lat, 30, 60) ~ "Mid",
      between(abs_lat, 0, 30) ~ "Low"
    )
  ) %>%
  # preprocess temperature data
  group_by(DepthHabitat) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(data = map(data, 
                    ~ .x %>% 
                      group_by(zone, bin) %>%
                      summarise(mean_temp = mean(temp_depth)) %>%
                      mutate(temp_change = mean_temp - lead(mean_temp)) %>%
                      ungroup())) %>% 
  unnest(cols = data) %>% 
  # merge with turnover dataframe
  right_join(dat_turnover) %>% 
  drop_na(temp_change)

