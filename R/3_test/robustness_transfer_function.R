library(here)
library(tidyverse)
library(rioja)
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
  filter(!(core_uniq == "124_767B" & bin == 292))

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 

# original (reported) beta coefficients
dat_ori_beta <- read_csv(here("data",
                              "beta_coefficient_per_latitude.csv")) %>% 
  add_column(source_data = "Original")



# prepare data ------------------------------------------------------------

# bring data in right format
dat_spp_full <- dat_spp %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) %>% 
  select(-c(core_uniq, bin)) %>% 
  as.data.frame()


# select training data based on interquartile range
dat_quant <- dat_mean_temp %>% 
  reframe(temp_quant = quantile(temp_ym_0m)) %>% 
  pull(temp_quant) 


# subset based on selected bins
dat_train <- dat_mean_temp %>% 
  filter(between(temp_ym_0m, dat_quant[2], dat_quant[4])) %>% 
  { filter(.data = dat_spp, bin %in% unique(.$bin)) } %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) %>%
  arrange(bin, core_uniq) %>%
  select(-c(core_uniq, bin)) %>%
  as.data.frame()


# same with temperature data
dat_train_temp <- dat_mean_temp %>% 
  filter(between(temp_ym_0m, dat_quant[2], dat_quant[4])) %>% 
  { filter(.data = dat_spp, bin %in% unique(.$bin)) } %>% 
        distinct(core_uniq, bin, temp_surface) %>% 
        arrange(bin, core_uniq) %>% 
        pull(temp_surface)


# imbrie and kipp factor analysis -----------------------------------------

# fit transfer function
mod_ikfa <- IKFA(dat_train, dat_train_temp)


# cross-validation - identify number of components to use for wapls 
nr_comp <- mod_ikfa %>% 
  crossval() %>%
  performance() %>%
  pluck("crossval") %>%
  .[, 1] %>%
  which.min()

# predict temperature values
dat_beta_ikfa <- predict(mod_ikfa, newdata = dat_spp_full)$fit[, nr_comp] %>%
  as_tibble_col(column_name = "cti") %>% 
  bind_cols(dat_spp %>% 
              pivot_wider(id_cols = c(core_uniq, bin), 
                          names_from = species, 
                          values_from = rel_abund, 
                          values_fn = mean, 
                          values_fill = 0)) %>% 
  select(core_uniq, bin, cti) %>% 
  # add estimated temperature from aogcm`s
  right_join(dat_spp %>%
               distinct(core_uniq, bin,
                        temp_surface,
                        pal.lat, pal.long)) %>%
  # calculate mismatch
  mutate(climatic_debt = temp_surface - cti) %>% 
  # split data into latitudinal zones and then apply fixed effect models
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat),
         zone = case_when(abs_lat >= 60 ~ "High",
                          between(abs_lat, 30, 60) ~ "Mid",
                          between(abs_lat, 0, 30) ~ "Low")) %>%
  # convert temperature in temperature anomaly
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) /
           sd(temp_surface, na.rm = TRUE)) %>%
  # prepare trend data
  mutate(
    short_term = if_else(temp_anom > lead(temp_anom),
                         "warming",
                         "cooling"),
    temp_change = temp_anom - lead(temp_anom,
                                   default = mean(temp_anom))
  ) %>%
  drop_na(short_term) %>% 
  # apply fixed effect models
  group_by(short_term, zone) %>%
  nest() %>%
  mutate(
    lm_mod = map(
      .x = data,
      .f = function(df) {
        lm(climatic_debt ~ temp_change,
           data = df)
      }))  %>% 
  # summarize the beta coefficient and save in csv
  mutate(beta_coef = map_dbl(lm_mod,
                             ~ coef(.x) %>% pluck(2)),
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 2),
         ci_high = map_dbl(ci, pluck, 4)) %>%
  select(zone, short_term, beta_coef, ci_low, ci_high) %>%
  ungroup() %>% 
  add_column(source_data = "IKFA") %>% 
  # add original beta coefficients
  full_join(dat_ori_beta)



# maximum likelihood response curves --------------------------------------

# fit transfer function
mod_mlrc <- dat_train %>% 
  mutate_all(~ .x / 100) %>% 
  MLRC(., dat_train_temp)


# cross-validation - identify number of components to use for wapls 
nr_comp <- mod_mlrc %>% 
  crossval() %>%
  performance() %>%
  pluck("crossval") %>%
  .[, 1] %>%
  which.min()

# predict temperature values
dat_beta_mlrc <- dat_spp_full %>%
  mutate_all(~ .x / 100) %>% 
  predict(mod_mlrc, newdata = .) %>% 
  pluck("fit") %>% 
  as_tibble_col(column_name = "cti") %>% 
  bind_cols(dat_spp %>% 
              pivot_wider(id_cols = c(core_uniq, bin), 
                          names_from = species, 
                          values_from = rel_abund, 
                          values_fn = mean, 
                          values_fill = 0)) %>% 
  select(core_uniq, bin, cti) %>% 
  # add estimated temperature from aogcm`s
  right_join(dat_spp %>%
               distinct(core_uniq, bin,
                        temp_surface,
                        pal.lat, pal.long)) %>%
  # calculate mismatch
  mutate(climatic_debt = temp_surface - cti) %>% 
  # split data into latitudinal zones and then apply fixed effect models
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat),
         zone = case_when(abs_lat >= 60 ~ "High",
                          between(abs_lat, 30, 60) ~ "Mid",
                          between(abs_lat, 0, 30) ~ "Low")) %>%
  # convert temperature in temperature anomaly
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) /
           sd(temp_surface, na.rm = TRUE)) %>%
  # prepare trend data
  mutate(
    short_term = if_else(temp_anom > lead(temp_anom),
                         "warming",
                         "cooling"),
    temp_change = temp_anom - lead(temp_anom,
                                   default = mean(temp_anom))
  ) %>%
  drop_na(short_term) %>% 
  # apply fixed effect models
  group_by(short_term, zone) %>%
  nest() %>%
  mutate(
    lm_mod = map(
      .x = data,
      .f = function(df) {
        lm(climatic_debt ~ temp_change,
           data = df)
      }))  %>% 
  # summarize the beta coefficient and save in csv
  mutate(beta_coef = map_dbl(lm_mod,
                             ~ coef(.x) %>% pluck(2)),
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 2),
         ci_high = map_dbl(ci, pluck, 4)) %>%
  select(zone, short_term, beta_coef, ci_low, ci_high) %>%
  ungroup() %>% 
  add_column(source_data = "MLRC") %>% 
  # add original beta coefficients
  full_join(dat_beta_ikfa)



# modern analogue technique -----------------------------------------------


# fit transfer function
mod_mat <- MAT(dat_train, dat_train_temp, 
               lean = FALSE)


# predict temperature values
dat_beta_full <- predict(mod_mat, newdata = dat_spp_full)$fit[, 1] %>%
  as_tibble_col(column_name = "cti") %>% 
  bind_cols(dat_spp %>% 
              pivot_wider(id_cols = c(core_uniq, bin), 
                          names_from = species, 
                          values_from = rel_abund, 
                          values_fn = mean, 
                          values_fill = 0)) %>% 
  select(core_uniq, bin, cti) %>% 
  # add estimated temperature from aogcm`s
  right_join(dat_spp %>%
               distinct(core_uniq, bin,
                        temp_surface,
                        pal.lat, pal.long)) %>%
  # calculate mismatch
  mutate(climatic_debt = temp_surface - cti) %>% 
  # split data into latitudinal zones and then apply fixed effect models
  # add latitudinal zones
  mutate(abs_lat = abs(pal.lat),
         zone = case_when(abs_lat >= 60 ~ "High",
                          between(abs_lat, 30, 60) ~ "Mid",
                          between(abs_lat, 0, 30) ~ "Low")) %>%
  # convert temperature in temperature anomaly
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) /
           sd(temp_surface, na.rm = TRUE)) %>%
  # prepare trend data
  mutate(
    short_term = if_else(temp_anom > lead(temp_anom),
                         "warming",
                         "cooling"),
    temp_change = temp_anom - lead(temp_anom,
                                   default = mean(temp_anom))
  ) %>%
  drop_na(short_term) %>% 
  # apply fixed effect models
  group_by(short_term, zone) %>%
  nest() %>%
  mutate(
    lm_mod = map(
      .x = data,
      .f = function(df) {
        lm(climatic_debt ~ temp_change,
           data = df)
      }))  %>% 
  # summarize the beta coefficient and save in csv
  mutate(beta_coef = map_dbl(lm_mod,
                             ~ coef(.x) %>% pluck(2)),
         ci = map(lm_mod, confint),
         ci_low = map_dbl(ci, pluck, 2),
         ci_high = map_dbl(ci, pluck, 4)) %>%
  select(zone, short_term, beta_coef, ci_low, ci_high) %>%
  ungroup() %>% 
  add_column(source_data = "MAT") %>% 
  # add original beta coefficients
  full_join(dat_beta_mlrc)




# visualise ---------------------------------------------------------------


# build coefficent plot
plot_comparison <- dat_beta_full %>%
  mutate(zone = ordered(zone, 
                        levels = c("High", "Mid", "Low")), 
         short_term = str_to_title(short_term), 
         source_data = ordered(source_data, 
                               levels = c("Original", 
                                          "IKFA",
                                          "MLRC", 
                                          "MAT"))) %>% 
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
                                        "transfer_function_test.png"), 
       width = image_width, height = image_height*1.5, units = image_units, 
       bg = "white", device = ragg::agg_png)


