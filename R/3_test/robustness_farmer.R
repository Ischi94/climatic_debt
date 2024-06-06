library(here)
library(tidyverse)
library(patchwork)
library(rioja)



# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))




# load data ---------------------------------------------------------------

# North Atlantic temperature data from Farmer et al 2008 for past 12ka
dat_temp <- read_table(here("data",   
                            "farmer_2008_MD99_2251.txt"))

# raw fossil occurrence compilation 
dat_raw <- read_rds(here("data",
                         "raw_foraminifera_data.rds"))

# species age range
dat_ranges <- read_csv(here("data", 
                            "final_species_age_ranges.csv"))

# modern depth ranges: Rebotim et al 2017, table 4
dat_depth <- read_rds(here("data",
                           "cleaned_depth_data.rds"))

# wapls model
mod_wapls <-read_rds(here("data", "wapls_model.rds"))

# preliminary data cleaning -----------------------------------------------


dat_clean <- dat_raw %>% 
  mutate(age = age * 1000, # transform to ka
         species = str_replace(species, # replace type
                               "Globototalia anfracta", 
                               "Globorotalia anfracta")) %>% 
  filter(water.depth >= 100) %>%  # remove depth where no breeding populations are viable
  # fix suffices, dashed and underscores
  mutate(core_uniq = str_replace(coreID, "-","_"), 
         core_uniq = str_replace(core_uniq, "181_1119B", "181_1119"), 
         core_uniq = str_split(core_uniq, "[.]"), 
         core_uniq = map_chr(core_uniq, pluck, 1))  %>%
  # remove cores which contain duplicate sets of entries, with different
  # coordinates for each set
  filter(!core_uniq %in% c("RC8C80", "RC8C81",
                           "RC8C82", "RC8C83",
                           "RC9C98", "RC9C99")) %>% 
  select(-coreID)

# average coordinates for cores with same entries but with a small
# location difference (measurment error)
dat_clean <- dat_clean %>%
  select(core_uniq, pal.lat, pal.long) %>% 
  distinct(core_uniq, pal.lat, pal.long) %>% 
  count(core_uniq, name = "dupl_core") %>% 
  left_join(dat_clean) %>% 
  mutate(pal.long = if_else(dupl_core > 1,
                            mean(pal.long),
                            pal.long),
         pal.lat = if_else(dupl_core > 1,
                           mean(pal.lat),
                           pal.lat))




# calculate range-through data --------------------------------------------

dat_rt <- dat_clean %>% 
  select(core_uniq, latitude, longitude, age) %>% 
  group_by(core_uniq) %>% 
  summarise(fad = max(age), 
            lad = min(age), 
            lat = mean(latitude), 
            long = mean(longitude))



# remove absence data -----------------------------------------------------

# remove records which have an abundance entry but no original abundance
dat_clean_pres <- dat_clean %>% 
  filter(!(is.na(orig.abundance) & abundance > 0)) %>% 
  # remove unclear/ ambigous original abundance codings
  filter(!orig.abundance %in% c("-", "/", 
                                "*", "rw?", 
                                "#0.0", "?1", 
                                "p?")) %>% 
  # modify the abundance where original abundance indicates presence, but
  # entered abundance is 0
  mutate(abundance = if_else(orig.abundance %in% c("r","P", "R", "F", 
                                                   "A", "T", "S", "C", 
                                                   "D", "X", "+", "VR",
                                                   "C/A", "x", "T/R", ">0",
                                                   "(F)", "B", "(R)", "F/C", 
                                                   "R/F", "(X)")
                             & abundance == 0, 
                             1, 
                             abundance)) %>% 
  # remove real absences
  filter(abundance != 0)



# clean on species-level data ------------------------------------------------------

# clean taxonomy
dat_clean_pres <- dat_clean_pres %>% 
  mutate(species = str_replace_all(species, 
                                   "Globigerinella radians",
                                   "Globigerinella siphonifera"), 
         species = str_replace_all(species, 
                                   "Menardella fimbriata",
                                   "Menardella menardii")) %>% 
  # remove T. cavernula and excelsa, which originated <800 ka
  filter((!earliest < 0.8) %>% replace_na(TRUE))


# remove species from reworked sediment, present younger than their LAD
dat_clean_pres_nrew <- dat_ranges %>% 
  # there are 5 species in the dataset not present in the FAD/LAD file
  # because of taxonomic revision
  mutate(Species.name = str_replace_all(Species.name, 
                                        "Globoturborotalita druryi",
                                        "Globigerina druryi"),
         Species.name = str_replace_all(Species.name,
                                        "Dentigloborotalia anfracta",
                                        "Globorotalia anfracta"),
         Species.name = str_replace_all(Species.name,
                                        "Pearsonites broedermanni",
                                        "Igorina broedermanni")) %>% 
  # remove reworked occurrences
  select(species = Species.name, lad = End) %>% 
  full_join(dat_clean_pres) %>% 
  filter(!age < lad) 


# remove records with very imprecise age estimates
dat_clean_pres_nrew <- dat_clean_pres_nrew %>% 
  filter(!(rng.age > 8/1000 & !is.na(rng.age))) %>% 
  filter(!age.model %in% c("Berggren1977","Ericson1968","GTSBlow1969","Raffi2006")) %>% 
  filter(age <= 12) 

# bin and summarize -------------------------------------------------------

# set up bins
brk <- seq(0, 12, by = 1)
bins <- data.frame(t = brk, b = brk + 1, mid = brk + 1/2)

# age 'zero' = 1950, and some observations are more recent (i.e. negative age)
bins$t[1] <- -0.1

dat_clean_binned <- dat_clean_pres_nrew %>%
  mutate(bin = map(age, ~ bins$t < .x & bins$b >= .x),
         bin = map_dbl(bin, ~ bins$mid[.x])) %>% 
  # calculate the relative abundance of each species within each core per bin
  group_by(bin, core_uniq, species) %>% 
  summarise(av_rel_ab = mean(rel.abun)) %>% 
  # standardize to relative abundance within bin within core
  mutate(rel_abund = av_rel_ab/sum(av_rel_ab)) %>% 
  ungroup() %>% 
  select(-av_rel_ab) %>% 
  left_join(dat_clean_pres_nrew %>% 
              distinct(core_uniq, pal.lat, pal.long)) %>% 
  # filter to high latitude
  filter(pal.lat >= 60) 


# bin temperature ---------------------------------------------------------


# bin temp data
dat_temp_bin <- dat_temp %>% 
  mutate(bin = map(age, ~ bins$t < .x & bins$b >= .x),
         bin = map_dbl(bin, ~ bins$mid[.x])) %>% 
  group_by(bin) %>% 
  summarise(mean_temp = mean(sst)) %>% 
  # calculate  change
  mutate(temp_change = mean_temp - lead(mean_temp)) 


# predict community temperature index -------------------------------------


# predict temperature on whole dataset using WAPLS
# bring data in right format
dat_spp_full <- dat_clean_binned %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) 

# predict community temperature index (cti)
dat_pred <- dat_spp_full %>% 
  select(-c(core_uniq, bin)) %>% 
  as.data.frame() %>% 
  predict(mod_wapls, newdata = ., npls = 4,
          nboot = 1000) %>% 
  pluck("fit") %>% 
  .[, 4] %>% 
  as_tibble() %>% 
  bind_cols(dat_spp_full) %>% 
  select(core_uniq, bin, cti = value)

# add estimated temperature from aogcm's
dat_trends <- dat_temp_bin %>% 
  left_join(dat_pred) %>% 
  # calculate offset between cti
  mutate(climatic_debt = mean_temp - cti, 
         short_term = if_else(temp_change >= 0, "warm", "cool")) 

# latitudinal wise
# split data into latitudinal zones and then apply fixed effect models
dat_mod <- dat_trends %>% 
  group_by(short_term) %>%
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
  select(short_term, beta_coef, ci_low, ci_high) %>% 
  write_csv(here("data",
                 "beta_coefficient_per_latitude_high_subset.csv"))
