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


# world map outlines
world <- map_data("world")



# climatic debt through time ----------------------------------------------

dat_debt_time <- dat_debt %>%
  group_by(bin) %>% 
  summarise(mean_cl_boot(climatic_debt)) %>% 
  select(bin, climatic_debt = y)
  
# visualise
plot_debt_time <- dat_debt_time %>% 
  ggplot(aes(bin, climatic_debt)) +
  geom_hline(yintercept = 0) +
  geom_line(colour = "coral", lwd = 1.3) +
  labs(x = "Age [ka]", 
       y = "Average Global\nClimatic Debt [°C]") +
  scale_x_reverse() +
  theme_minimal() +
  theme(panel.grid = element_blank())



# global temperature -----------------------------------------------------------


# global temperature through time
plot_temp <- dat_mean_temp %>% 
  ggplot(aes(bin, temp_ym_0m)) +
  geom_line() +
  labs(x = "Age [ka]", 
       y = "Average Global\nTemperature [°C]") +
  scale_x_reverse() +
  coord_cartesian(xlim = c(700, 0)) +
  theme_minimal() +
  theme(panel.grid = element_blank())




# climatic debt as a function of temperature ------------------------------

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
lm_boot <- function(formula, data = dat_debt_trends) {
  
  dat_model <- data 
    }

# short-term only
nr_rows <- dat_debt_trends %>% 
  drop_na(st) %>% 
  nrow()

# preallocate vector
dat_aic <- vector(mode = "double", 
                  length = nr_rows)

for (i in 1:nr_rows) {
  
  dat_aic[i] <- dat_debt_trends %>% 
    drop_na(st) %>% 
    add_column(nr = 1:nr_rows) %>% 
    filter(nr != i) %>% 
    lm(climatic_debt ~ st, data = .) %>% 
    AIC()
  
}


# null model with intercept only
mod_null <- lm(climatic_debt ~ 1, data = dat_debt_trends)

# short-term only
mod_st <- lm(climatic_debt ~ st, data = dat_debt_trends)

# short-term interacting with 8ka long-term trend
mod_lt1 <- lm(climatic_debt ~ st:lt1, data = dat_debt_trends)

# short-term interacting with 16ka long-term trend
mod_lt2 <- lm(climatic_debt ~ st:lt2, data = dat_debt_trends)

# short-term interacting with 24ka long-term trend
mod_lt3 <- lm(climatic_debt ~ st:lt3, data = dat_debt_trends)


# put in a list to call iteratively 
mod_list <- list(mod_null, mod_st,
                 mod_lt1, mod_lt2,
                 mod_lt3)

# create comparison table
dat_mod_comp <- tibble(model = c("Null Model",
                                 "Short Term",
                                 "Short Term:Long Term 8ka",
                                 "Short Term:Long Term 16ka",
                                 "Short Term:Long Term 24ka"),
                       aic = map_dbl(mod_list, AIC),
                       bic = map_dbl(mod_list, BIC)) %>%
  mutate(delta_aic = aic - .[[4, 2]],
         delta_bic = bic - .[[4, 3]])



  
# visualize model comparison
dat_mod_comp %>% 
  mutate(model = fct_reorder(model, c(5:1))) %>% 
  ggplot(aes(delta_aic, model)) +
  geom_segment(aes(x = 0, xend = delta_aic, 
                   yend = model)) +
  geom_point() +
  labs(y = NULL, x = expression(paste(Delta, "  AIC"))) +
  theme_minimal()




plot_debt_time/plot_temp +
  plot_layout(heights = c(2, 1))
