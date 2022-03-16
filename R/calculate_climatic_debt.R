library(raster)
library(here)
library(tidyverse)
library(patchwork)



# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- readRDS(here("data", 
                        "spp-and-sampling-data_list-by-depth.rds"))

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m)




# species temperature index -------------------------------------

# using surface water temperature
dat_spp_surf <- dat_spp %>% 
  pluck("temp_ym_0m") %>% 
  pluck("sp") %>% 
  as_tibble() 

# convert cell numbers to coordinates
rEmpt <- raster(ncols=288, nrows=144, xmn=-180, xmx=180, ymn=-90, ymx=90)
dat_spp_surf[,c('centroid_long', 'centroid_lat')] <- xyFromCell(rEmpt, dat_spp_surf$cell)

# calculate climate debt per cell and per bin
dat_sti_surf <- dat_spp_surf %>% 
  group_by(species) %>%
  summarise(sti = mean(temp_ym)) %>% 
  full_join(dat_spp_surf) %>% 
  arrange(cell) %>% 
  group_by(bin, cell, centroid_long, centroid_lat) %>% 
  summarise(cti = mean(sti), 
            act_temp = mean(temp_ym)) %>% 
  ungroup() %>% 
  select(cell, centroid_long, centroid_lat, bin, act_temp, cti) %>% 
  mutate(temp_lag = act_temp - cti, 
         cell = as.factor(cell))


# save data
write_rds(dat_sti_surf, 
          here("data", 
               "cleaned_debt_raster.rds"))

# visualise over time -----------------------------------------------------


# quick plot
dat_sti_surf %>%
  mutate(dist_equ = abs(centroid_lat)) %>%
  ggplot(aes(bin, temp_lag)) +
  geom_line(aes(group = cell, colour = dist_equ),
            alpha = 0.4) +
  geom_smooth(colour = "coral3") +
  scale_x_reverse() +
  theme_minimal()




