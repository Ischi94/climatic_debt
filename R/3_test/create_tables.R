library(here)
library(tidyverse)
library(flextable)
library(officer)


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                         "spp_by_depth.rds")) %>% 
  # right format for WAPLS
  mutate(rel_abund = rel_abund * 100) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292))

# debt data
dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low"), 
         zone = factor(zone, levels = c("High", "Mid", "Low")))

  
dat_niche <- dat_spp %>% 
  group_by(species) %>% 
  summarise(mean_temp = mean(temp_surface), 
            sd_temp = sd(temp_surface))


# create tables -----------------------------------------------------------


# number of assemblages per zone and north/south
dat_zone_tbl <- dat_debt %>% 
  mutate(north_south = if_else(pal.lat > 0, "North", 
                               "South")) %>% 
  count(zone, north_south) %>% 
  pivot_wider(names_from = north_south, values_from = n) %>% 
  arrange(zone) %>% 
  flextable() %>% 
  theme_vanilla() %>% 
  # set column names
  compose(j = 1, part = "header", 
          value = as_paragraph("Latitudinal Zone")) %>%
  compose(j = 2, part = "header", 
          value = as_paragraph("Global North")) %>% 
  compose(j = 3, part = "header", 
          value = as_paragraph("Global South")) 


# number of species through
dat_spp_tbl <- dat_spp %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low"), 
         zone = factor(zone, levels = c("High", "Mid", "Low"))) %>% 
  count(species, zone) %>% 
  pivot_wider(names_from = zone, values_from = n) %>% 
  rowwise() %>% 
  mutate(sum = sum(High, Mid, Low, na.rm = TRUE), 
         .before = High) %>% 
  flextable() %>% 
  theme_vanilla() %>% 
  # set column names
  compose(j = 1, part = "header", 
          value = as_paragraph("Species")) %>%
  compose(j = 2, part = "header", 
          value = as_paragraph("Total"))



# summary of species temperature niches
dat_niche_tbl <- dat_niche %>%
  select(species, mean_temp, sd_temp) %>% 
  mutate(across(c(mean_temp, sd_temp), 
         round, 1)) %>% 
  arrange(species) %>% 
  flextable() %>% 
  theme_vanilla() %>% 
  # set column names
  compose(j = 1, part = "header", 
          value = as_paragraph("Species")) %>%
  compose(j = 2, part = "header", 
          value = as_paragraph("Mean [°C]")) %>% 
  compose(j = 3, part = "header", 
          value = as_paragraph("SD [°C]")) 


# create word document ----------------------------------------------------

# open docx-file and add flextable
my_doc <- read_docx() %>% 
  body_add_flextable(dat_zone_tbl) %>% 
  body_add_break() %>% 
  body_add_flextable(dat_spp_tbl, pos = "after") %>% 
  body_add_break() %>% 
  body_add_flextable(dat_depth_tbl, pos = "after") %>% 
  body_add_break() %>% 
  body_add_flextable(dat_niche_tbl, pos = "after")

# convert to word file/ add input to empty docx
print(my_doc, target = here("figures",
                            "supplemental", 
                            "supplemental_tables.docx"))
