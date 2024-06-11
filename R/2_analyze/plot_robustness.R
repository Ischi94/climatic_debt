library(here)
library(tidyverse)

# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# main results
dat_main <- read_csv(here("data", 
                          "beta_coefficient_per_latitude.csv")) %>% 
  add_column(type = "Main")

# autoregressive model
dat_autor <- read_csv(here("data",
                          "beta_coefficient_per_latitude_autoregressive.csv")) %>% 
  add_column(type = "Autoregressive")

# Niche-depth model
dat_depth <- read_csv(here("data",
                           "beta_coefficient_per_latitude_by_depth.csv")) %>% 
  rename(type = "DepthHabitat")

# hemispheric split
dat_hemis <- read_csv(here("data",
                           "beta_coefficient_per_latitude_hemisph.csv")) 

# proxy data
dat_proxy <- read_csv(here("data",
                           "beta_coefficient_per_latitude_high_subset.csv")) %>% 
  add_column(type = "Proxy Data") %>% 
  add_column(zone = "High")

# sampling test
dat_samp <- read_csv(here("data", 
                 "beta_coefficient_per_latitude_sampling.csv")) %>% 
  filter(type %in% c("Older than 4ka", 
                     "Random Effects", 
                     "Younger than 400ka")) 


# wapls test
dat_wapls <- read_csv(here("data", 
                          "beta_coefficient_per_latitude_wapls.csv")) %>% 
  filter(source_data != "Original") %>% 
  rename(type = source_data)



# merge data --------------------------------------------------------------

dat_full <- dat_main %>% 
  full_join(dat_autor) %>% 
  full_join(dat_depth) %>% 
  full_join(dat_hemis) %>% 
  full_join(dat_proxy) %>% 
  full_join(dat_samp) %>% 
  full_join(dat_wapls)

# visualise
plot_comparison <- dat_full %>%
  filter(zone != "Mid") %>% 
  mutate(zone = paste0(zone, " Latitude"),
         zone = ordered(zone, 
                        levels = c("High Latitude", "Low Latitude")), 
         short_term = str_to_title(short_term), 
         ori_id = if_else(type == "Original", 
                          "Main results", 
                          "Sensitivity analyses")) %>% 
  arrange(desc(ori_id)) %>% 
  ggplot(aes(beta_coef, short_term,
             group = type, 
             colour = ori_id)) +
  geom_vline(xintercept = 0) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.6) +
  geom_point(position = position_dodge(width = 0.5), 
             size = 2) +
  facet_wrap(~ zone, 
             ncol = 1) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6)) +
  scale_color_manual(values = c(colour_coral, 
                                colour_grey)) +
  geom_label(data = tibble(zone = "High Latitude",
                          short_term = "Cooling",
                          type = "Original",
                          ori_id = "Main results",
                            beta_coef = 2),
            position = position_nudge(y = -0.1),
            size = 10/.pt,
            alpha = 0.7, 
            label.size = 0,
            aes(label = ori_id)) +
  geom_curve(data = tibble(zone = "High Latitude",
                           short_term = "Cooling",
                           type = "Original",
                           ori_id = "Main results",
                           beta_coef_min = 0.5, 
                           beta_coef_max = 2),
             aes(x = beta_coef_min, 
                 xend = beta_coef_max, 
                 yend = short_term), 
             arrow = arrow(length = unit(0.05, "inch"),
                           ends = "first"),
             lwd = 0.2, 
             curvature = -0.2) +
  labs(y = NULL, 
       x = "Beta Coefficient", 
       colour = NULL) +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA,
                                    colour = "grey60"), 
        panel.grid.major.x = element_line(colour = "grey60", 
                                          linetype = "dotted"))
plot_comparison

# save plot
ggsave(plot_comparison, filename = here("figures",
                                        "fig2_1_sensitivity_test.png"),
       width = image_width, height = image_height, units = image_units,
       bg = "white", device = ragg::agg_png)

