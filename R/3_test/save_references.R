library(here)
library(tidyverse)


# load data ---------------------------------------------------------------

# raw fossil occurrence compilation with references
dat_raw <- dat_raw %>% 
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
                           "RC9C98", "RC9C99"))


# processed debt data
dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds"))


# get references ----------------------------------------------------------

# get cores used in the manuscript
cores_used <- unique(dat_debt$core_uniq)

# filter out references
dat_raw %>% 
  filter(core_uniq %in% cores_used) %>% 
  distinct(source) %>% 
  # save file
  write_delim(here("data", 
                   "references.txt"))
