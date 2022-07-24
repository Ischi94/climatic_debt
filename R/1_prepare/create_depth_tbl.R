library(here)
library(tidyverse)


# load data ---------------------------------------------------------------


# modern depth ranges: Rebotim et al 2017, table 4
dat_depth <- read_table(here("data", 
                             "Rebotim_et_al_depth_ranges.txt"))



# clean depth data --------------------------------------------------------


dat_depth_clean <- dat_depth %>% 
  mutate(MaxAbund = str_sub(MaxAbund, end = -2), 
         MaxAbund = as.numeric(MaxAbund)) %>% 
  # add references
  add_column(ref = "Rebotim et al. 2017") %>% 
  # Schiebel et al. 2004 also confirm that minuta is within top 60m
  mutate(DepthHabitat = if_else(Species == "Globigerinita.minuta",
                                "Surface",
                                DepthHabitat),  
         # specify depth range based on Retailleau et al. 2011
         DepthHabitat = if_else(Species == "Globigerinita.uvula",
                                "Surface",
                                DepthHabitat),
         ALD = if_else(Species == "Globigerinita.uvula",
                       26,
                       ALD), 
         ref = if_else(Species == "Globigerinita.uvula",
                       "Retailleau et al. 2011",
                       ref), 
         DepthHabitat = if_else(Species == "Hastigerinella.digitata",
                                "Subsurface",
                                DepthHabitat),
         ref = if_else(Species == "Hastigerinella.digitat",
                       "Hull et al. 2011",
                       ref), 
         # Schiebel et al. sampled ALD of 60.5 m, but quite variable
         # Schmuker and Schiebel 2002 also got depth of 60-80 m
         # but Schiebel and Hemleben 2017 considered the species shallow
         DepthHabitat = if_else(Species == "Dentigloborotalia.anfracta",
                                "Surface",
                                DepthHabitat),
         ref = if_else(Species == "Dentigloborotalia.anfractat",
                       "Schiebel et al. 2004, Rebotim et al. 2017",
                       ref), 
         # B digitata is subsurface
         DepthHabitat = if_else(Species == "Beella.digitata",
                                "Subsurface",
                                DepthHabitat),
         ALD = if_else(Species == "Beella.digitata",
                       150,
                       ALD), 
         VD = if_else(Species == "Beella.digitata",
                      75,
                      VD), 
         ref = if_else(Species == "Beella.digitata",
                       "Schiebel and Hemleben 2017",
                       ref), 
         # G menardii ranges from 20-60m, ALD of 25 read off fig 6c
         # Schiebel et al. 2004 called menardii an upwelling indicator sp
         # and report densities that give ALD=58.6 and Vd=21.9 for one leg,
         # ALD=46.2 and VD=25.2 for another leg (with higher concentrations)
         DepthHabitat = if_else(Species == "Globorotalia.menardii",
                                "Surface",
                                DepthHabitat),
         ALD = if_else(Species == "Globorotalia.menardii",
                       46.2,
                       ALD), 
         VD = if_else(Species == "Globorotalia.menardii",
                      25.2,
                      VD), 
         MaxAbund = if_else(Species == "Globorotalia.menardii",
                            51,
                            MaxAbund), 
         ref = if_else(Species == "Globorotalia.menardii",
                       "Watkins et al. 1996",
                       ref), 
         # Re-assign crassaformis as mid-depth as a compromise among conflicting studies. 
         # Meilland et al (2019) and Ezard et al (2015) consider the species sub-thermocline,
         # but Rebotim et al (2017) and Schiebel (2017) consider it shallow
         DepthHabitat = if_else(Species == "Globorotalia.crassaformis",
                                "Surface.subsurface",
                                DepthHabitat),
         ref = if_else(Species == "Globorotalia.crassaformis",
                       "Ezard et al. 2015, Rebotim et al. 2017, Schiebel 2017, Meilland et al. 2019",
                       ref)) %>% 
  # add depth data
  # add G trilobus data from Loncaric et al. 2006
  add_row(tibble(Species = "Globigerinoides.trilobus",
                 N = NA,
                 MaxAbund = 208,
                 ALD = 17.4,
                 ALDse = NA,
                 VD = 7.0,
                 VDse = NA,
                 DepthHabitat = "Surface",
                 DepthHabitatVar = NA,
                 ref = "Loncaric et al. 2006")) %>% 
  # add G (H) theyeri data from Schiebel et al. 2004
  add_row(tibble(Species = "Hirsutella.theyeri",
                 N = NA,
                 MaxAbund = 0.2,
                 ALD = 50,
                 ALDse = NA,
                 VD = 15,
                 VDse = NA,
                 DepthHabitat = "Surface",
                 DepthHabitatVar = NA,
                 ref = "Schiebel et al. 2004")) %>% 
  # add depth ranges from Schiebel and Hemleben 2017
  add_row(tibble(Species = "Sphaeroidinella.dehiscens",
                 N = NA,
                 MaxAbund = NA,
                 ALD = 150,
                 ALDse = NA,
                 VD = 75,
                 VDse = NA,
                 DepthHabitat = "Subsurface",
                 DepthHabitatVar = NA,
                 ref = "Schiebel and Hemleben 2017")) %>% 
  add_row(tibble(Species = "Globorotalia.ungulata",
                 N = NA,
                 MaxAbund = NA,
                 ALD = 46.2,
                 ALDse = NA,
                 VD = 25.2,
                 VDse = NA,
                 DepthHabitat = "Surface",
                 DepthHabitatVar = NA,
                 ref = "Schiebel and Hemleben 2017")) %>% 
  # Add depth ranges from Watkins et al. 1996
  add_row(tibble(Species = "Globorotalia.tumida",
                 N = NA,
                 MaxAbund = 30,
                 ALD = 50,
                 ALDse = NA,
                 VD = 15,
                 VDse = NA,
                 DepthHabitat = "Surface",
                 DepthHabitatVar = NA,
                 ref = "Watkins et al. 1996")) %>% 
  add_row(tibble(Species = "Globigerinoides.conglobatus",
                 N = NA,
                 MaxAbund = 3,
                 ALD = 50,
                 ALDse = NA,
                 VD = 15,
                 VDse = NA,
                 DepthHabitat = "Surface",
                 DepthHabitatVar = NA,
                 ref = "Watkins et al. 1996")) %>% 
  add_row(tibble(Species = "Globorotaloides.hexagonus",
                 N = NA,
                 MaxAbund = 2,
                 ALD = 150,
                 ALDse = NA,
                 VD = 75,
                 VDse = NA,
                 DepthHabitat = "Subsurface",
                 DepthHabitatVar = NA,
                 ref = "Watkins et al. 1996")) %>% 
  add_row(tibble(Species = "Globorotaloides.hexagonus",
                 N = NA,
                 MaxAbund = 2,
                 ALD = 150,
                 ALDse = NA,
                 VD = 75,
                 VDse = NA,
                 DepthHabitat = "Subsurface",
                 DepthHabitatVar = NA,
                 ref = "Watkins et al. 1996")) %>% 
  add_row(tibble(Species = "Globigerinella.adamsi",
                 N = NA,
                 MaxAbund = 2,
                 ALD = 125,
                 ALDse = NA,
                 VD = 75,
                 VDse = NA,
                 DepthHabitat = "Surface.subsurface",
                 DepthHabitatVar = NA,
                 ref = "Watkins et al. 1996, Schiebel et al. 2004")) %>% 
  add_row(tibble(Species = "Globoquadrina.conglomerata",
                 N = NA,
                 MaxAbund = 2.4,
                 ALD = 42.5,
                 ALDse = NA,
                 VD = 30.8,
                 VDse = NA,
                 DepthHabitat = "Surface.subsurface",
                 DepthHabitatVar = NA,
                 ref = "Watkins et al. 1996, Schiebel et al. 2004")) %>% 
  # synonymize according to Microtax for congruence with occurrence database
  mutate(Species = str_replace(Species, 
                               "Dentigloborotalia.anfracta", 
                               "Globorotalia.anfracta"), 
         Species = str_replace(Species, 
                               "Globorotalia.scitula", 
                               "Hirsutella.scitula"), 
         Species = str_replace(Species, 
                               "Globorotalia.inflata", 
                               "Globoconella.inflata"), 
         Species = str_replace(Species, 
                               "Globorotalia.hirsuta", 
                               "Hirsutella.hirsuta"), 
         Species = str_replace(Species, 
                               "Berggrenia.pumillio", 
                               "Berggrenia.pumilio"), 
         Species = str_replace(Species, 
                               "Globorotalia.crassaformis", 
                               "Truncorotalia.crassaformis"),
         Species = str_replace(Species, 
                               "Globorotalia.menardii", 
                               "Menardella.menardii"), 
         Species = str_replace(Species, 
                              "Globigerinoides.trilobus", 
                              "Trilobatus.trilobus"), 
         Species = str_replace(Species, 
                               "Globigerinoides.ruber.white", 
                               "Globigerinoides.ruber"),
         Species = str_replace(Species, 
                               "Globigerinoides.tenellus", 
                               "Globoturborotalita.tenella"),
         Species = str_replace(Species, 
                               "Globorotalia.truncatulinoides", 
                               "Truncorotalia.truncatulinoides"),
         )


# save data ---------------------------------------------------------------

dat_depth_clean %>% 
  write_rds(here("data", 
                 "cleaned_depth_data.rds"))
