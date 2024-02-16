library(here)
library(tidyverse)
library(rioja)


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

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 



# fit WAPLS -------------------------------------


# get the species data in the right format
dat_train_spec <- dat_spp %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) %>% 
  arrange(bin, core_uniq) %>% 
  select(-c(core_uniq, bin)) %>% 
  as.data.frame() 


# same with temperature data
dat_train_temp <- dat_spp %>% 
  distinct(core_uniq, bin, temp_surface, temp_depth) %>% 
  arrange(bin, core_uniq) %>% 
  pull(temp_depth) 

# fit the weighted average partial least squares model (WAPLS)
mod_wapls <- WAPLS(dat_train_spec, dat_train_temp, npls = 10)



# predict community temperature index -------------------------------------


# predict temperature on whole dataset using WAPLS
# bring data in right format
dat_spp_full <- dat_spp %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) 

# predict community temperature index (cti)
dat_pred <- dat_spp_full %>% 
  select(-c(core_uniq, bin)) %>% 
  as.data.frame() %>% 
  predict(mod_wapls, newdata = ., npls = 7,
          nboot = 1000) %>% 
  pluck("fit") %>% 
  .[, 7] %>% 
  as_tibble() %>% 
  bind_cols(dat_spp_full) %>% 
  select(core_uniq, bin, cti = value)

# add estimated temperature from aogcm's
dat_debt <- dat_spp %>% 
  distinct(core_uniq, bin, 
           temp_depth, DepthHabitat, 
           pal.lat, pal.long) %>% 
  left_join(dat_pred) %>% 
  # calculate offset between cti
  mutate(climatic_debt = temp_depth - cti) 



# save data
dat_debt %>% 
  write_rds(here("data", 
                 "cleaned_debt_wapls_depth.rds"))




# prepare trend data
dat_trends <- dat_debt %>% 
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat), 
         zone = case_when(
           abs_lat >= 60 ~ "High",
           between(abs_lat, 30, 60) ~ "Mid", 
           between(abs_lat, 0, 30) ~ "Low")) %>% 
  # convert temperature in temperature anomaly
  mutate(temp_anom = (temp_depth - mean(temp_depth, na.rm = TRUE)) / 
           sd(temp_depth, na.rm = TRUE)) %>% 
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)




# get trends
dat_mod <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(climatic_debt ~ temp_change,
                           data = df)
                      })
  )

new_data_lat <- dat_mod %>%
  mutate(new_data = if_else(short_term == "warming", 
                            list(tibble(temp_change = seq(0, 3, by = 0.1), 
                                        bin = 0)), 
                            list(tibble(temp_change = seq(-3, 0, by = 0.1), 
                                        bin = 0))), 
         predicted_debt = map2(.x = lm_mod,
                               .y = new_data, 
                               .f = ~ predict(.x,
                                              newdata = .y,
                                              interval = "confidence") %>% 
                                 as_tibble)) %>% 
  select(-c(data, lm_mod)) %>% 
  unnest(cols = c(predicted_debt, new_data)) %>% 
  rename(predicted_debt = fit)

# plot temperature anomaly versus climatic debt per latitudinal zone
plot_trends_lat <- new_data_lat %>%
  ggplot(aes(temp_change, predicted_debt)) +
  geom_point(aes(temp_change, climatic_debt, 
                 colour = zone), 
             data = dat_trends, 
             alpha = 0.1, 
             size = 0.6) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_ribbon(aes(fill = zone, 
                  group = interaction(short_term, zone),
                  ymin = lwr, ymax = upr)) +
  geom_line(aes(colour = zone, 
                group = interaction(short_term, zone)), 
            lwd = 1) + 
  scale_color_manual(values = alpha(c(colour_lavender,
                                      colour_green,
                                      colour_brown), 0.6)) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                     colour_green,
                                     colour_brown), 0.15)) +
  scale_y_continuous(breaks = seq(-20, 15, 5)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-20, 15), 
                  xlim = c(-2.5, 2.5)) +
  labs(y = "Thermal Deviance [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  theme(legend.position = "none") +
  theme(axis.ticks = element_blank())


# save plot_final
ggsave(plot_trends_lat, filename = here("figures", 
                                        "supplemental",
                                        "debt_trends_depth.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)




# individual models per depth ---------------------------------------------

# get trends
dat_mod_ind <- dat_trends %>% 
  group_by(short_term, zone, 
           DepthHabitat) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lm(climatic_debt ~ temp_change,
                           data = df)
                      })
  )

new_data_lat_ind <- dat_mod_ind %>%
  mutate(new_data = if_else(short_term == "warming", 
                            list(tibble(temp_change = seq(0, 3, by = 0.1), 
                                        bin = 0)), 
                            list(tibble(temp_change = seq(-3, 0, by = 0.1), 
                                        bin = 0))), 
         predicted_debt = map2(.x = lm_mod,
                               .y = new_data, 
                               .f = ~ predict(.x,
                                              newdata = .y,
                                              interval = "confidence") %>% 
                                 as_tibble)) %>% 
  select(-c(data, lm_mod)) %>% 
  unnest(cols = c(predicted_debt, new_data)) %>% 
  rename(predicted_debt = fit)

# visualise
new_data_lat_ind %>%
  ggplot(aes(temp_change, predicted_debt)) +
  geom_point(aes(temp_change, climatic_debt, 
                 colour = zone), 
             data = dat_trends, 
             alpha = 0.1, 
             size = 0.6) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_ribbon(aes(fill = zone, 
                  group = interaction(short_term, zone),
                  ymin = lwr, ymax = upr)) +
  geom_line(aes(colour = zone, 
                group = interaction(short_term, zone)), 
            lwd = 1) + 
  scale_color_manual(values = alpha(c(colour_lavender,
                                      colour_green,
                                      colour_brown), 0.6)) +
  scale_fill_manual(values = alpha(c(colour_lavender,
                                     colour_green,
                                     colour_brown), 0.15)) +
  scale_y_continuous(breaks = seq(-20, 15, 5)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-20, 15), 
                  xlim = c(-2.5, 2.5)) +
  labs(y = "Thermal Deviance [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  theme(legend.position = "none") +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~ DepthHabitat)


# original (reported) beta coefficients
dat_ori_beta <- read_csv(here("data",
                              "beta_coefficient_per_latitude.csv")) %>% 
  add_column(DepthHabitat = "Original")

# add new estimates
# summarize the beta coefficient and save in csv
dat_beta <- dat_mod_ind %>% 
  mutate(beta_coef = map_dbl(lm_mod,
                           ~ coef(.x) %>% pluck(2)),
       ci = map(lm_mod, confint),
       ci_low = map_dbl(ci, pluck, 2),
       ci_high = map_dbl(ci, pluck, 4)) %>%
  select(DepthHabitat, zone, short_term, beta_coef, ci_low, ci_high) %>%
  ungroup() %>%  
  # add original beta coefficients
  full_join(dat_ori_beta)


# build coefficent plot
plot_comparison <- dat_beta %>%
  mutate(zone = ordered(zone, 
                        levels = c("High", "Mid", "Low")), 
         short_term = str_to_title(short_term), 
         source_data = ordered(DepthHabitat, 
                               levels = c("Original", 
                                          "Surface",
                                          "Surface.subsurface", 
                                          "Subsurface"))) %>% 
  ggplot(aes(beta_coef, short_term,
             colour = source_data)) +
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
                                        "depth_test.png"), 
       width = image_width, height = image_height*1.5, units = image_units, 
       bg = "white", device = ragg::agg_png)


# on a global scale -------------------------------------------------------

dat_trends %>%
  mutate(climatic_debt = climatic_debt/3) %>% 
  ggplot(aes(temp_change, climatic_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  geom_line(data = tibble(climatic_debt = seq(-3, 3,
                                              length.out = 20),
                          temp_change = seq(-3, 3,
                                            length.out = 20)),
            colour = "grey40") +
  geom_hline(yintercept = 0, 
             colour = "grey40") +
  geom_smooth(aes(colour = DepthHabitat, 
                  fill = DepthHabitat), 
              linewidth = 0.6, 
              alpha = 0.2) +
  annotate(geom = "text", 
           x = 1.9, y = 2.25, 
           label = "No response", 
           angle = 28, 
           colour = "grey20",
           size = 10/.pt) +
  annotate(geom = "text", 
           x = -1.4, y = 0.25, 
           label = "Equilibrium", 
           colour = "grey20",
           size = 10/.pt) +
  annotate(geom = "curve",
           x = -0.25, xend = 0.8,
           y = -0.9, yend = -2.6,
           curvature = 0.3,
           arrow = arrow(length = unit(0.05, "inch"),
                         ends = "first"),
           colour = colour_grey, lwd = 0.3) +
  annotate(geom = "text",
           x = 1.27, y = -2.6,
           label = "GAM",
           colour = colour_grey,
           size = 10/.pt) +
  labs(y = "Thermal Deviance [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.2, 2.3)) +
  theme(axis.ticks = element_blank())
