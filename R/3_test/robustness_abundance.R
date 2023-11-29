library(here)
library(tidyverse)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                         "spp_by_depth.rds")) %>% 
  # right format for WAPLS
  mutate(rel_abund = rel_abund * 100) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292))


# climatic debt based on weighted averages --------------------------------

dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  # convert temperature in temperature anomaly
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) / 
           sd(temp_surface, na.rm = TRUE)) %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)

# climatic debt based on presence only ------------------------------------


# calculate climate debt per cell and per bin
dat_debt_surf <- dat_spp %>%
  group_by(species) %>%
  summarise(sti = mean(temp_surface)) %>%
  full_join(dat_spp) %>%
  group_by(bin, core_uniq, pal.lat, pal.long) %>%
  summarise(cti = mean(sti),
            act_temp = mean(temp_surface)) %>%
  ungroup() %>%
  mutate(temp_debt = act_temp - cti)



# robustness check  -----------------------------------------------------


# prepare trend data
dat_trends <- dat_debt_surf %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  # convert temperature in temperature anomaly
  mutate(temp_anom = (temp_debt - mean(temp_debt, na.rm = TRUE)) / 
           sd(temp_debt, na.rm = TRUE)) %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)

# fit overall model 
mod1 <- lm(temp_debt ~ temp_change,
           data = dat_trends)


# create equal grid for inference (accounting for differential sampling in bins)
new_data <- tibble(temp_change = seq(-3, 3, by = 0.2),
                   bin = 0)  

# estimate over this equal grid  
dat_trends_pred <- predict(mod1,
                           newdata = new_data,
                           interval = "confidence") %>%
  as_tibble() %>% 
  add_column(new_data) %>% 
  rename(predicted_debt = fit)


# in total
plot_trends_comp <- dat_trends_pred %>% 
  ggplot(aes(temp_change, predicted_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  # add the slope based on relative abundance
  geom_smooth(aes(temp_change, climatic_debt), 
              data = dat_debt, 
              method = "lm", 
              colour = colour_coral,
              fill = colour_coral,
              alpha = 0.8, linewidth = 1) +

  geom_ribbon(aes(temp_change, ymin = lwr, ymax = upr), 
              fill = colour_grey,
              alpha = 0.2) +
  geom_abline(intercept = summary(mod1)$coefficients[1,1], 
              slope = summary(mod1)$coefficients[2,1], 
              colour = colour_grey, 
              lwd = 1) +
  annotate(geom = "text", 
           x = -1.4, y = 0, 
           label = "Relative abundance", 
           angle = 2, 
           colour = "grey20",
           size = 10/.pt) +
  annotate(geom = "text", 
           x = -1, y = -7.8, 
           label = "Occurrence only", 
           colour = "grey20",
           angle = 15, 
           size = 10/.pt) +
  labs(y = "Thermal Deviance [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(xlim = c(-2.5, 2.5))

plot_trends_comp

# save plot
ggsave(plot_trends_comp, filename = here("figures", 
                                   "supplemental",
                                   "debt_trends_comparison.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)



# latitudinal wise  -----------------------------------------------------

# split data into latitudinal zones and then apply mixed effect models,
# accounting for sampling
dat_mod <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(temp_debt ~ temp_change,
                             data = df)
                      })
  )

# predict to get trend lines
new_data_lat <- dat_mod %>%
  mutate(new_data = if_else(short_term == "warming", 
                            list(tibble(temp_change = seq(0, 3, by = 0.2), 
                                        bin = 0)), 
                            list(tibble(temp_change = seq(-3, 0, by = 0.2), 
                                        bin = 0))), 
         predicted_debt = map2(.x = lm_mod,
                               .y = new_data, 
                               .f = ~ predict(.x,
                                              newdata = .y,
                                              interval = "prediction") %>% 
                                 as_tibble())) %>% 
  select(-c(lm_mod)) %>% 
  unnest(cols = c(predicted_debt, new_data)) %>% 
  rename(predicted_debt = fit)

# visualize
plot_trends_lat <- new_data_lat %>% 
  mutate(zone = factor(zone, levels = c("High",
                                        "Mid",
                                        "Low"))) %>% 
  ggplot(aes(temp_change, predicted_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_ribbon(aes(fill = zone, 
                  group = interaction(short_term, zone),
                  ymin = lwr, ymax = upr)) +
  geom_smooth(aes(colour = zone,
                  group = interaction(short_term, zone)),
              lwd = 1, 
              se = FALSE) + 
  scale_color_manual(name = "Latitude", 
                     values = alpha(c(colour_lavender,
                                      colour_brown, 
                                      colour_green), 0.8)) +
  scale_fill_manual(name = "Latitude", 
                    values = alpha(c(colour_lavender,
                                     colour_brown, 
                                     colour_green), 0.1)) +
  scale_y_continuous(breaks = seq(-12, 4, 4)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  labs(y = "Thermal Deviance [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  theme(legend.position = c(0.1, 0.85)) +
  theme(axis.ticks = element_blank())

# save plot
ggsave(plot_trends_lat, filename = here("figures",
                                        "supplemental",
                                        "debt_trends_occurrence_latitude.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)


