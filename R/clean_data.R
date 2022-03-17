library(here)
library(tidyverse)


# load data ---------------------------------------------------------------

# raw fossil occurrence compilation 
dat_raw <- read_rds(here("data",
                         "raw_foraminifera_data.rds"))

# species age range
dat_ranges <- read_csv(here("data", 
                            "final_species_age_ranges.csv"))

# modern depth ranges: Rebotim et al 2017, table 4
dat_depth <- read_rds(here("data",
                           "cleaned_depth_data.rds"))



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
  
# average coordinates for cores with same entries but with a 0.1 degree
  # location difference
dat_clean <- dat_clean %>%
  group_by(core_uniq) %>%
  mutate(pal.long = if_else(core_uniq == "NA87_22",
                            mean(pal.long),
                            pal.long),
         pal.lat = if_else(core_uniq == "NA87_22",
                            mean(pal.lat),
                           pal.lat)) %>%
  ungroup()
  


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
  mutate(abundance = if_else(!is.na(orig.abundance) & abundance == 0, 
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
  filter(!earliest < 0.8)


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




# Add modern depth ranges to data -----------------------------------------

dat_clean_pres_nrew %>% 
  mutate(Species = str_replace(species, 
                               " ",
                               ".")) %>% 
  left_join(dat_depth, by = "Species") %>% 
  # remove duplicates
  distinct(.keep_all = TRUE) %>% 
  # clean data
  select()

# cut down file size
irrel <- c('orig.species','abundance','orig.abundance','abun.units',
           'age.st','age.en','AM.type','rel.abun','round.age',
           'site','hole','core','section','sample.top','sample.type',
           'total.IDd','preservation','processing','Comments',
           'Average.Area..µm2.','log.area','Approximate.diameter')
cols2toss <- colnames(occ) %in% c(traits, irrel)
occ <- occ[,!cols2toss]

  


