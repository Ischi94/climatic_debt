library(here)
library(tidyverse)


# load data ---------------------------------------------------------------

dat_raw <- read_rds(here("data", 
                 "raw_foraminifera_data.rds"))



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
  filter(core_uniq == "NA87_22") %>% 
  summarise(new_long = mean(pal.long), 
            new_lat = mean(pal.lat), 
            core_uniq = "NA87_22") %>% 
  full_join(dat_clean %>% 
              filter(core_uniq == "NA87_22")) %>% 
  mutate(pal.long = new_long, 
         pal.lat = new_lat) %>% 
  select(-c(new_long, new_lat)) %>% 
  bind_rows(dat_clean %>% 
              filter(!core_uniq == "NA87_22"))
  

  


