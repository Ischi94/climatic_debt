library(here)
library(tidyverse)
library(rioja)
library(lme4)

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


# select training data for WA-PLS
# select bins with average temperature based on quantiles
dat_temp_quant <- dat_mean_temp %>% 
  summarize(temp_quant = quantile(temp_ym_0m)) %>% 
  pull(temp_quant) 

dat_bins <- dat_mean_temp %>% 
  filter(between(temp_ym_0m, dat_temp_quant[2], dat_temp_quant[4]))

dat_train <- dat_spp %>% 
  filter(bin %in% unique(dat_bins$bin))


# get the species data in the right format
dat_train_spec <- dat_train %>% 
  pivot_wider(id_cols = c(core_uniq, bin), 
              names_from = species, 
              values_from = rel_abund, 
              values_fn = mean, 
              values_fill = 0) %>% 
  arrange(bin, core_uniq) %>% 
  select(-c(core_uniq, bin)) %>% 
  as.data.frame() 


# same with temperature data
dat_train_temp <- dat_train %>% 
  distinct(core_uniq, bin, temp_surface, temp_depth) %>% 
  arrange(bin, core_uniq) %>% 
  pull(temp_depth) 

# fit the weighted average partial least squares model (WAPLS)
mod_wapls <- WAPLS(dat_train_spec, dat_train_temp)

# cross-validate model
mod_cv <- crossval(mod_wapls, cv.method = "loo")

# How many components to use based on RMSE via loo
nr_comp <- performance(mod_cv)$crossval[, 1] %>% 
  which.min(.)



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
  predict(mod_wapls, newdata = ., npls = nr_comp,
          nboot = 1000) %>% 
  pluck("fit") %>% 
  .[, nr_comp] %>% 
  as_tibble() %>% 
  bind_cols(dat_spp_full) %>% 
  select(core_uniq, bin, cti = value)

# add estimated temperature from aogcm's
dat_debt <- dat_spp %>% 
  distinct(core_uniq, bin, 
           temp_depth, 
           pal.lat, pal.long) %>% 
  left_join(dat_pred) %>% 
  # calculate offset between cti
  mutate(climatic_debt = temp_depth - cti) 



# save data
write_rds(dat_debt, 
          here("data", 
               "cleaned_debt_wapls_by_depth.rds"))


# robustness check  -----------------------------------------------------


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

# fit overall model 
mod1 <- lmer(climatic_debt ~ temp_change + (1 | bin),
             data = dat_trends)

# summarize beta coefficient
bootMer(mod1, fixef, nsim = 1000)



# latitudinal wise
# split data into latitudinal zones and then apply mixed effect models,
# accounting for sampling
dat_mod <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lmer(climatic_debt ~ temp_change + (1 | bin),
                             data = df)
                      })
  )

# summarize the beta coefficient
dat_mod %>% 
  mutate(fix_eff = map(lm_mod, ~ bootMer(.x, fixef, nsim = 1000)), 
         beta_coef = map(lm_mod, fixef),
         beta_coef = map_dbl(beta_coef, pluck(2)),
         beta_coef_sd = map_dbl(fix_eff, ~ sd(.x$t[, 2])), 
         ci_low = beta_coef - 1.96 * beta_coef_sd, 
         ci_high = beta_coef + 1.96 * beta_coef_sd) %>%
  select(zone, short_term, beta_coef, ci_low, ci_high) %>% 
  write_csv(here("data", 
                 "beta_coefficient_per_latitude_by_debth.csv"))

