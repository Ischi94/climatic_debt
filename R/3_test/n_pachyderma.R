library(here)
library(tidyverse)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                         "spp_by_depth.rds")) %>% 
  # right format for WAPLS
  mutate(rel_abund = rel_abund * 100) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292)) %>% 
  filter(species == "Neogloboquadrina pachyderma") 


# get annotations
lab_surf <- paste0(
  "Surface\n",
  "Optimum = ",
  round(density(dat_spp$temp_surface)$x[which.max(density(dat_spp$temp_surface)$y)], 1), "°C",
  "\n", "Mean = ", round(mean(dat_spp$temp_surface), 1), "°C")

lab_depth <- paste0(
  "Depth\n",
  "Optimum = ",
  round(density(dat_spp$temp_depth)$x[which.max(density(dat_spp$temp_depth)$y)], 1), "°C",
  "\n", "Mean = ", round(mean(dat_spp$temp_depth), 1), "°C")

# niche 
plot_niche <- dat_spp %>% 
  pivot_longer(cols = contains("temp"), 
               values_to = "temp", 
               names_to = "temp_type") %>% 
  ggplot(aes(temp, fill = temp_type, 
             colour = temp_type)) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c(colour_lavender, 
                                colour_coral)) +
  scale_fill_manual(values = c(colour_lavender, 
                                colour_coral)) +
  annotate("text", 
           label = lab_depth, 
           x = 8, 
           y = 0.057, 
           size = 8/.pt, 
           hjust = 0,
           colour = colour_lavender) +
  annotate("text", 
           label = lab_surf, 
           x = 12, 
           y = 0.048, 
           hjust = 0,
           size = 8/.pt, 
           colour = colour_coral) +
  labs(x = "Temperature [°C]", 
       y = "Density", 
       fill = NULL, 
       colour = NULL, 
       title = "Thermal niche of N. pachyderma") +
  theme(legend.position = "none", 
        plot.title = element_text(face  = "bold"))

# save
ggsave(plot_niche, filename = here("figures", 
                                   "supplemental",
                                   "niche_n_pachyderma.png"), 
       width = image_width, height = image_height, 
       units = image_units, 
       bg = "white", device = ragg::agg_png)


