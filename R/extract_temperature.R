library(here)
library(raster)
library(tidyverse)



# load data ---------------------------------------------------------------

# cleaned and processed species data
dat_spp <- read_rds(here("data", 
                         "cleaned_spp_data.rds"))



# average living depth ----------------------------------------------------

dat_spp %>% 
  group_by(DepthHabitat) %>% 
  summarise(ald = mean(ALD)) %>% 
  drop_na(DepthHabitat)

# the equivalent AOGCM depth levels from which to extract temp for the ald are
# 1, 4, 6, 8
dpths <- c(1, 4, 6, 8)

raster1 <- stack(here("data", 
                       "gcm_annual_mean", 
                       "teit1_140_ann_temp_ym_dpth.tif"), bands = c(1, 4, 6, 8))

plot(raster1)
