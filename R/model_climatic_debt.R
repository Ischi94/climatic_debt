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

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 





# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# model comparison --------------------------------------------------------


# climatic debt as a function of temperature
# calculate short-term changes
dat_debt_trends <- dat_debt %>% 
  # summarize per latitudinal zone
  group_by(bin, zone) %>% 
  summarise(av_temp = mean(temp_surface)) %>% 
  pivot_wider(names_from = zone, values_from = av_temp) %>% 
  ungroup() %>% 
  # calculate lags
  mutate(High_lag = lead(High, default = mean(High)), 
         Low_lag = lead(Low, default = mean(Low)), 
         Mid_lag = lead(Mid, default = mean(Mid))) %>% 
  # calculate short-term changes
  mutate(High_st = (High - High_lag)/8, 
         Low_st = (Low - Low_lag)/8, 
         Mid_st = (Mid - Mid_lag)/8) %>% 
  # calculate first long-term trend
  mutate(High_lt1 = (High_lag - lead(High_lag, default = mean(High_lag)))/8, 
         Low_lt1 = (Low_lag - lead(Low_lag, default = mean(Low_lag)))/8, 
         Mid_lt1 = (Mid_lag - lead(Mid_lag, default = mean(Mid_lag)))/8) %>% 
  # calculate second long-term trend
  mutate(High_lt2 = (High_lag - lead(High_lag, default = mean(High_lag), n = 2))/8, 
         Low_lt2 = (Low_lag - lead(Low_lag, default = mean(Low_lag), n = 2))/8, 
         Mid_lt2 = (Mid_lag - lead(Mid_lag, default = mean(Mid_lag), n = 2))/8) %>% 
  # calculate third long-term trend
  mutate(High_lt3 = (High_lag - lead(High_lag, default = mean(High_lag), n = 3))/8, 
         Low_lt3 = (Low_lag - lead(Low_lag, default = mean(Low_lag), n = 3))/8, 
         Mid_lt3 = (Mid_lag - lead(Mid_lag, default = mean(Mid_lag), n = 3))/8) %>% 
  # bring into right format
  select(-c(High:Mid_lag)) %>% 
  pivot_longer(-bin,
               names_to = c("zone", ".value"),
               names_pattern = "(.+)_(.+)") %>% 
  # merge back with climatic debt data
  full_join(dat_debt %>% 
              select(-c(temp_surface:cti, abs_lat)))



# simple linear models

# use leave-one-out cross-validation 
aic_loo <- function(lm_frml,
                    data = dat_debt_trends,
                    parameters = "climatic_debt") {
  
  # get observations
  nr_rows <- data %>% 
    drop_na(all_of(parameters)) %>% 
    nrow() 
  
  # preallocate vector
  dat_aic <- vector(mode = "double", 
                    length = nr_rows)
  
  for (i in 1:nr_rows) {
    
    dat_aic[i] <- data %>% 
      drop_na(all_of(parameters)) %>% 
      add_column(nr = 1:nr_rows) %>% 
      filter(nr != i) %>% 
      lm(formula = lm_frml, data = .) %>% 
      AIC()
  }
  
  return(dat_aic)
  
}



# null model with intercept only
mod_null <- aic_loo(lm_frml = formula(climatic_debt ~ 1))

# short-term only
mod_st <- aic_loo(lm_frml = formula(climatic_debt ~ st),
                  parameters = "st")

# short-term interacting with 8ka long-term trend
mod_lt1 <- aic_loo(lm_frml = formula(climatic_debt ~ st:lt1),
                   parameters = c("st", "lt1"))

# short-term interacting with 16ka long-term trend
mod_lt2 <- aic_loo(lm_frml = formula(climatic_debt ~ st:lt2),
                   parameters = c("st", "lt2"))

# short-term interacting with 24ka long-term trend
mod_lt3 <- aic_loo(lm_frml = formula(climatic_debt ~ st:lt3),
                   parameters = c("st", "lt3"))


# put in a list to call iteratively 
mod_list <- list(mod_null, mod_st,
                 mod_lt1, mod_lt2,
                 mod_lt3)

# summarise leave-one-out AIC

# create comparison table
dat_mod_comp <- tibble(model = c("Null Model",
                                 "Temperature",
                                 "Temperature\n + 8ka",
                                 "Temperature\n + 16ka",
                                 "Temperature\n + 24ka"),
                       aic = map_dbl(mod_list, mean),
                       aic_sd = map_dbl(mod_list, sd)) %>%
  # include uncertainty based on loo
  mutate(lwr = aic - aic_sd, 
         upr = aic + aic_sd) %>% 
  # calculate delta_aic
  mutate(delta_aic = aic - .[[4, 2]],
         lwr = lwr - .[[4, 2]], 
         upr = upr - .[[4, 2]], 
         radius = (upr - lwr)/2)

# summarize uncertainty for r-squared for best performing lag

# preallocate vector
dat_rsq <- vector(mode = "double", 
                  length = 1000)

# iterate in blocks of 500 (bootstrapping)
for (i in 1:1000) {
  
  dat_rsq[i] <- dat_debt_trends %>% 
    drop_na(st, lt2) %>% 
    slice_sample(n = 500) %>% 
    lm(climatic_debt ~ st:lt2, data = .) %>% 
    summary() %>% 
    pluck("r.squared")
}

# summarize via normal approximation
dat_rsq %>% 
  as_tibble() %>% 
  summarise(mean_cl_normal(value))


# climatic debt through time ----------------------------------------------

# calculate average trend and confidence interval via bootstrapping
boot_debt <- function() {
  
  dat_debt %>% 
    slice_sample(n = 1000) %>% 
    group_by(bin) %>% 
    summarise(mean_debt = mean(climatic_debt)) %>% 
    { spline(.$bin, .$mean_debt, ties = min, 
             xout = unique(dat_debt$bin),
             method = "natural") } %>% 
    pluck("y")
  
}

# replicate 10.000 times and get quantile based uncertainty intervals
dat_debt_boot <- replicate(1e4, boot_debt()) %>% 
  as_tibble() %>%
  add_column(bin = unique(dat_debt$bin)) %>% 
  pivot_longer(cols = - bin, 
               names_to = "iteration", 
               values_to = "climatic_debt") %>% 
  group_by(bin) %>% 
  summarise(lwr = quantile(climatic_debt, probs = 0.25), 
            upr = quantile(climatic_debt, probs = 0.75),
            climatic_debt = mean(climatic_debt)) 




# visualizations ----------------------------------------------------------



# global temperature through time
plot_temp <- dat_mean_temp %>% 
  ggplot(aes(bin, temp_ym_0m)) +
  geom_line(colour = "grey20") +
  labs(x = "Age [ka]", 
       y = "Average Global\nTemperature [°C]") +
  scale_x_reverse() +
  coord_cartesian(xlim = c(700, 0), 
                  ylim = c(12, 15)) +
  scale_y_continuous(breaks = c(12, 13, 14))

# increase alpha when less observations, create dataset to map alpha onto
dat_alpha <- dat_debt %>% 
  group_by(bin, core_uniq) %>%
  count() %>%
  ungroup() 

# climatic debt through time
plot_debt_time <- dat_debt_boot %>% 
  ggplot(aes(bin, climatic_debt)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              colour = "white", 
              fill = colour_grey, 
              alpha = 0.7) +
  geom_line(colour = alpha(colour_coral, 0.7), lwd = 1) +
  geom_line(colour = "white", 
            size = 1.1,
            data = dat_debt_boot %>% 
              filter(between(bin, 420, 480)), 
            alpha = 0.8) +
  annotate("curve", 
           x = 350, xend = 410, 
           y = -6.5, yend = -4.5, curvature = -0.3,
           arrow = arrow(ends = "last", 
                         length = unit(.2,"cm")), 
           colour = "grey70") +
  annotate("label",
           x = 230, y = -6.5,
           label = "sampling artefact",
           colour = "grey70",
           size = 10/.pt, 
           label.size = 0) +
  labs(x = "Age [ka]", 
       y = "Average Global\nClimatic Debt [°C]") +
  scale_y_continuous(breaks = seq(-4, 2, by = 2)) +
  scale_x_reverse() +
  theme()


# model comparison
plot_mod_comp <- dat_mod_comp %>% 
  mutate(model = fct_reorder(model, c(5:1))) %>% 
  ggplot(aes(delta_aic, model)) +
  geom_segment(aes(x = 0, xend = delta_aic, 
                   yend = model), 
               colour = "grey20") +
  geom_point(aes(size = radius), 
             colour = colour_grey) +
  geom_point(size = 1.5, shape = 21, 
             stroke = 0.15, 
             fill = colour_coral,
             alpha = 0.7,
             colour = "grey20") +
  labs(y = NULL, x = expression(paste(Delta, "  AIC"))) +
  theme(legend.position = "none")




# combine and save --------------------------------------------------------


# set out layout for final plot
layout <- "
AAA#
AAAC
AAAC
BBBC
BBB#
"

# patch together
plot_final <- plot_debt_time + plot_temp +  plot_mod_comp +
  plot_annotation(tag_levels = "a") +
  plot_layout(design = layout) & 
  coord_cartesian(clip = "off")


# save plot_final
ggsave(plot_final, filename = here("figures", 
                                   "fig4_comparison.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)
