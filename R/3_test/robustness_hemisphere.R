library(here)
library(tidyverse)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))




# load data ---------------------------------------------------------------

# debt data
dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) 


# prepare trend data
dat_trends <- read_rds(here("data",
                            "trend_data.rds"))


# model climatic debt versus temperature ----------------------------------


# latitudinal wise
# split data into latitudinal zones and then apply fixed effect models
dat_mod <- dat_trends %>%
  mutate(hemisph = if_else(pal.lat>= 0, "North", "South")) %>% 
  group_by(short_term, zone, hemisph) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(climatic_debt ~ temp_change,
                           data = df)
                      })
  )

# summarize the beta coefficient 
dat_beta <- dat_mod %>%
  mutate(beta_coef = map_dbl(lm_mod, 
                             ~ coef(.x) %>% pluck(2)), 
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 2),  
         ci_high = map_dbl(ci, pluck, 4)) %>% 
  select(zone, short_term,hemisph, beta_coef, ci_low, ci_high) %>%
  arrange(zone) %>% 
  ungroup() %>% 
  rename(type = hemisph) %>% 
  drop_na()
 
# save data
dat_beta %>% 
  write_csv(here("data", 
                 "beta_coefficient_per_latitude_hemisph.csv"))

# visualize ---------------------------------------------------------------



# build coefficent plot
plot_comparison <- dat_beta %>%
  # add original 
  full_join( read_csv(here("data",
                           "beta_coefficient_per_latitude.csv")) %>% 
               add_column(type = "Original") %>% 
               ungroup()) %>% 
  mutate(zone = ordered(zone, 
                        levels = c("High", "Mid", "Low")), 
         short_term = str_to_title(short_term), 
         type = ordered(type, 
                        levels = c("Original",
                                   "North",
                                   "South"))) %>% 
  ggplot(aes(beta_coef, short_term,
             colour = type)) +
  geom_vline(xintercept = 0) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), 
                 position = position_dodge(width = 0.5), 
                 alpha = 0.6) +
  geom_point(position = position_dodge(width = 0.5), 
             size = 2) +
  facet_wrap(~ zone, 
             ncol = 1) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6)) +
  scale_color_brewer(name = NULL,
                     type = "qual", 
                     palette = 6) +
  labs(y = NULL, 
       x = "Beta Coefficient") +
  theme(legend.position = "bottom", 
        panel.border = element_rect(fill = NA,
                                    colour = "grey60"), 
        panel.grid.major.x = element_line(colour = "grey60", 
                                          linetype = "dotted"))
# save plot
ggsave(plot_comparison, filename = here("figures", 
                                        "supplemental",
                                        "hemisphere_test.png"), 
       width = image_width, height = image_height*1.5, units = image_units, 
       bg = "white", device = ragg::agg_png)
