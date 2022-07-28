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

# depth test
dat_depth <- read_csv(here("data", 
                           "beta_coefficient_per_latitude_by_debth.csv")) %>% 
  mutate((across(beta_coef:ci_high, round, 2)), 
         beta_coef_depth = paste0(beta_coef, " [", ci_low, ", ", ci_high, "]")) %>% 
  select(zone, short_term, beta_coef_depth) %>% 
  full_join(read_csv(here("data", 
                          "beta_coefficient_per_latitude.csv")) %>% 
              mutate((across(beta_coef:ci_high, round, 2)),
                     beta_coef = paste0(beta_coef, " [", ci_low, ", ", ci_high, "]")) %>% 
              select(zone, short_term, beta_coef))
  


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


# summary of trend estimates based on surface versus depth temperature
dat_depth_tbl <- dat_depth %>% 
  mutate(zone = factor(zone, levels = c("High", "Mid", "Low"))) %>% 
  arrange(zone) %>% 
  select(zone, short_term, beta_coef, beta_coef_depth) %>% 
  flextable() %>% 
  theme_vanilla() %>% 
  # set column names
  compose(j = 1, part = "header", 
          value = as_paragraph("Latitudinal Zone")) %>%
  compose(j = 2, part = "header", 
          value = as_paragraph("Temperature Change")) %>% 
  compose(j = 3, part = "header", 
          value = as_paragraph("Surface")) %>% 
  compose(j = 4, part = "header", 
          value = as_paragraph("Preferred Depth")) %>% 
  merge_v(j = c("zone"))


# create word document ----------------------------------------------------

# open docx-file and add flextable
my_doc <- read_docx() %>% 
  body_add_flextable(dat_zone_tbl) %>% 
  body_add_break() %>% 
  body_add_flextable(dat_spp_tbl, pos = "after") %>% 
  body_add_break() %>% 
  body_add_flextable(dat_depth_tbl, pos = "after") 

# convert to word file/ add input to empty docx
print(my_doc, target = here("figures",
                            "supplemental", 
                            "supplemental_tables.docx"))
