library(here)
library(tidyverse)
library(patchwork)

# load data ---------------------------------------------------------------

dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low"))

plot1 <- dat_debt %>% 
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) / 
           sd(temp_surface, na.rm = TRUE)) %>% 
  ggplot(aes(climatic_debt, temp_anom)) +
  geom_jitter(aes(colour = zone), 
              width = 0, height = 0.1, 
              alpha = 0.1) +
  geom_smooth(aes(colour = zone),
              method = "lm") +
  theme_minimal() +
  theme(legend.position = "none")


plot2 <- dat_debt %>% 
  group_by(bin, zone) %>% 
  summarise(av_temp = mean(temp_surface)) %>% 
  pivot_wider(names_from = zone, values_from = av_temp) %>% 
  ungroup() %>% 
  mutate(across(high:mid, ~ if_else(.x > lead(.x), 
                                    "warming", 
                                    "cooling"))) %>% 
  full_join(dat_debt %>% select(-zone)) %>% 
  pivot_longer(cols = high:mid, 
               names_to = "zone", 
               values_to = "short_term") %>% 
  drop_na(short_term) %>% 
  mutate(abs_debt = abs(climatic_debt)) %>% 
  group_by(short_term, zone) %>% 
  summarise(mean_cl_boot(abs_debt)) %>% 
  ggplot(aes(short_term, y, 
             ymin = ymin, 
             ymax = ymax, 
             colour = zone)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none")



# world map outlines
world <- map_data("world")

plot3 <- dat_debt %>% 
  distinct(pal.lat, pal.long) %>% 
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  ggplot() +
  geom_point(aes(x = pal.long, y = pal.lat,
                 fill = zone),
             na.rm = TRUE,
             alpha = 0.7, shape = 21,
             colour = "grey25",
             stroke = 0.3,
             size = 2) +
  geom_map(aes(map_id = region), 
           data = world, map = world, fill = "grey40") + 
  scale_x_continuous(name = '', limits = c(-180, 180), expand = c(0,0), 
                     labels = NULL, breaks = seq(-180, 180, by = 45)) +
  scale_y_continuous(name = '', limits = c(-90, 90),   expand = c(0,0), 
                     labels = NULL) +
  scale_fill_discrete(name = "Latitude") +
  theme_minimal() +
  coord_map(projection = "mollweide") +
  theme(legend.position = 'bottom')


plot3 / (plot1 + plot2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(2, 1, 1))
  