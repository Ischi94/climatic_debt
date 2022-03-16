library(here)
library(tidyverse)
library(patchwork)

# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- readRDS(here("data", 
                        "spp-and-sampling-data_list-by-depth_2020-11-15.rds"))

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m)




# species temperature index -------------------------------------

# using surface water temperature
dat_spp_surf <- dat_spp %>% 
  pluck("temp_ym_0m") %>% 
  pluck("sp") %>% 
  as_tibble()

dat_sti_surf <- dat_spp_surf %>% 
  group_by(species) %>%
  summarise(sti = mean(temp_ym), 
            sti_sd = sd(temp_ym))


# community temperature index -------------------------------------

dat_temp_lag <- dat_spp_surf %>% 
  full_join(dat_sti_surf %>% 
              select(-sti_sd)) %>% 
  group_by(bin, cell) %>%
  summarise(act_temp = mean(temp_ym), 
            cti = mean(sti)) %>% 
  mutate(temp_lag = act_temp - cti) %>% 
  summarise(mean_cl_normal(temp_lag)) %>% 
  rename(temp_lag = y, 
         lwr_ci = ymin, 
         upr_ci = ymax)


plot1 <- dat_temp_lag %>% 
  ggplot(aes(bin, temp_lag)) +
  geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci), 
              alpha = 0.7) +
  geom_line(colour = "coral2", lwd = 1) +
  geom_hline(yintercept = 0) +
  scale_x_reverse() +
  coord_cartesian(xlim = c(700, 0)) +
  labs(x = "Time (ka)", 
       y = "Lag (°C)") +
  theme_minimal()


plot2 <- dat_mean_temp %>% 
  ggplot(aes(bin, temp_ym_0m)) +
  geom_line(colour = "steelblue") +
  scale_x_reverse() +
  coord_cartesian(xlim = c(700, 0)) +
  labs(x = "Time (ka)", 
       y = "Global\nTemperature (°C)") +
  theme_minimal()

plot1/plot2


# corellation -------------------------------------------------------------


dat_cor <- dat_temp_lag %>% 
  select(bin, temp_lag) %>% 
  full_join(dat_mean_temp) 

dat_cor %>% 
  ggplot(aes(temp_ym_0m, temp_lag)) +
  geom_point() +
  geom_smooth(method = "lm")

lm(dat_cor$temp_ym_0m ~ dat_cor$temp_lag) %>% 
  summary()



# ondra's table -----------------------------------------------------------

plot3 <- dat_spp_surf %>% 
  full_join(dat_sti_surf %>% 
              select(-sti_sd)) %>% 
  arrange(cell) %>% 
  group_by(bin, cell) %>% 
  summarise(cti = mean(sti), 
            act_temp = mean(temp_ym)) %>% 
  ungroup() %>% 
  select(cell, bin, act_temp, cti) %>% 
  mutate(temp_lag = act_temp - cti, 
         cell = as.character(cell)) %>% 
  ggplot(aes(bin, temp_lag)) +
  geom_line(aes(group = cell), 
             show.legend = FALSE, alpha = 0.4) +
  geom_point(aes(group = cell), 
            show.legend = FALSE, alpha = 0.4) +
  geom_smooth(colour = "coral3") +
  scale_x_reverse() +
  theme_minimal()

plot3/plot2


