library(here)
library(raster)
library(tidyverse)



# load data ---------------------------------------------------------------

# cleaned and processed species data
dat_spp <- read_rds(here("data", 
                         "cleaned_spp_data.rds"))

# age codes for aogcm models
dat_age_stack <- read_csv(here("data",
                               "gcm_model_codes.csv"))


# average living depth ----------------------------------------------------

# calculate the average living depth per habitat and add it to the original data
dat_spp_depth <- dat_spp %>% 
  group_by(DepthHabitat) %>% 
  summarise(ald = mean(ALD)) %>% 
  drop_na(DepthHabitat) %>% 
  right_join(dat_spp)

# the equivalent AOGCM depth levels from which to extract temp for the ald are
# 1, 4, 6, 8
dat_spp_depth <- tibble(DepthHabitat = c("Surface", 
                                         "Surface.subsurface", 
                                         "Subsurface"), 
                        raster_band = c(4, 6, 8)) %>% 
  right_join(dat_spp_depth) 



# extract from raster -----------------------------------------------------



# load tif files as raster stack at the water surface and each average living
# depth and then extract temperature for each core within each bin
dat_temp <- list.files(here("data", "gcm_annual_mean"), # get all downloaded tif files
                        pattern = ".tif$",
                        full.names = TRUE) %>%
  as_tibble() %>% 
  # filter out those bins not needed
  mutate(raster_name = str_match(value, "mean/\\s*(.*?)\\s*.tif")[,2], 
         id = str_extract(raster_name, "^[^_]+(?=_)")) %>% 
  left_join(dat_age_stack) %>% 
  filter(age_1000ka %in% unique(dat_spp_depth$bin)) %>% 
  select(pwd = value, bin = age_1000ka) %>% 
  # add the focal directory of the raster file to the spp data
  full_join(dat_spp_depth) %>% 
  # only the distinct combinations of core within bin
  distinct(bin, core_uniq, pal.lat, pal.long, raster_band, pwd) %>% 
  # load in surface raster if depth is unknown
  # reduce memory
  group_by(bin, raster_band, pwd) %>% 
  nest() %>% 
  ungroup() %>% 
  # load in surface raster if depth is unknown
  mutate(raster_band = replace_na(raster_band, 4)) %>% 
  # load raster files into tibble as list
  mutate(depth_raster = map2(.x = pwd,
                             .y = raster_band,
                             .f = ~ raster(x = .x, band = .y)),  
         surface_raster = map(.x = pwd, # raster at water surface
                              .f = ~ raster(x = .x, band = 1))) %>% 
  # extract temperature from raster for surface
  mutate(coord = map(data, 
                     ~ cbind(.x$pal.long, .x$pal.lat)), 
         temp_surface = map2(.x = surface_raster,
                             .y = coord,
                             .f = ~ raster::extract(.x, .y))) %>% 
  # same for depth
  mutate(temp_depth = map2(.x = depth_raster,
                           .y = coord,
                           .f = ~ raster::extract(.x, .y))) %>%
  unnest(cols = c(data, temp_surface, temp_depth))  %>% 
  # Infer environment if it's missing and some of the adjacent 9 cells have values
  mutate(coord = map2(pal.long, pal.lat,
                      cbind),  
         temp_surface = if_else(is.na(temp_surface), 
                                # distance to corner cell is 196 km for
                                # 1.25-degree resolution (~111 km/degree)
                                map2_dbl(
                                  .x = surface_raster,
                                  .y = coord,
                                  .f = ~ raster::extract(.x, .y,
                                                         buffer = 200 * 1000,
                                                         fun = mean)), 
                                temp_surface), 
         temp_depth = if_else(is.na(temp_depth),
                              map2_dbl(
                                .x = depth_raster,
                                .y = coord,
                                .f = ~ raster::extract(.x, .y,
                                                       buffer = 200 * 1000,
                                                       fun = mean)), 
                              temp_depth)) %>% 
  select(core_uniq, bin, temp_surface, temp_depth, raster_band)


# clean and save ----------------------------------------------------------


# combine with species data
dat_final <- dat_temp %>%
  distinct() %>%  
  right_join(dat_spp_depth, by = c("bin", "core_uniq", "raster_band")) %>% 
  # remove where temperature is unknown
  drop_na(temp_surface, temp_depth) %>% 
  select(-c(N:ref, raster_band, ald)) %>% 
  # restrict to the last 700 ka, to trim edge effects of range-through cores
  filter(bin <= 700) %>% 
  # retain only those cores with at least 4 successive steps 
  group_by(core_uniq) %>% 
  mutate(n_occ = n()) %>% 
  ungroup() %>% 
  filter(n_occ >= 4) %>% 
  select(-n_occ)



# save as rds
dat_final %>% 
  write_rds(here("data", "spp_by_depth.rds"), 
            compress = "gz")


