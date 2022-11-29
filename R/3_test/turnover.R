library(here)
library(tidyverse)
library(vegan)
library(patchwork)


# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# cleaned foraminifera occurrence time series along with temperature
dat_spp <- read_rds(here("data", 
                         "spp_by_depth.rds")) %>% 
  # there is one row with na's
  filter(!(core_uniq == "124_767B" & bin == 292))

# average global temperature in each bin
dat_mean_temp <- read_csv(here("data", 
                               "global-MAT_10-deg-grid_8ka.csv")) %>% 
  select(bin, temp_ym_0m) 

dat_proxy <- read_rds(here("data", 
                           "proxy_temperature.rds")) %>% 
  # nas are encode as -999
  filter(temp != -999)



# preprocess data ---------------------------------------------------------


# get the species data in the right format
dat_spec <- dat_spp %>% 
  nest_by(bin) %>% 
  ungroup() %>% 
  mutate(data = map(data, 
                    ~ pivot_wider(.x, 
                                  id_cols = core_uniq, 
                                  names_from = species, 
                                  values_from = rel_abund, 
                                  values_fn = mean, 
                                  values_fill = 0) %>% 
                      select(-core_uniq) %>% 
                      colSums() %>% 
                      enframe() %>% 
                      mutate(rel_abund = value/sum(value)) %>% 
                      select(-value) %>% 
                      pivot_wider(names_from = name, values_from = rel_abund, 
                                  values_fill = 0))) %>% 
  unnest(cols = c(data)) %>% 
  replace(is.na(.), 0) 
  

# calculate turnover ------------------------------------------------------


# convert to matrix and calculate turnover
dat_mat <- dat_spec %>% 
  column_to_rownames(var = "bin") %>% 
  vegdist(method = "chisq") %>% 
  as.matrix()  
  
# extract diagonal
dat_diag <- dat_mat[row(dat_mat) == col(dat_mat) + 1]
  
# join with dataset
dat_final <- dat_mean_temp %>% 
  filter(bin %in% dat_spec$bin) %>% 
  mutate(temp_change = temp_ym_0m - lead(temp_ym_0m)) %>% 
  add_column(turnover = c(dat_diag, 0)) 

# visualise
plot_turnover <- dat_final %>%
  mutate(temp = if_else(temp_change >= 0, "warm", "cool")) %>% 
  drop_na(temp) %>% 
  ggplot(aes(temp_change, turnover)) +
  geom_hline(yintercept = 0, colour = "grey80", 
             linetype = "dotted") +
  geom_vline(xintercept = 0, colour = "grey80", 
             linetype = "dotted") +
  geom_point() +
  geom_smooth(aes(group = temp), 
              method = "lm", 
              colour = colour_coral, 
              fill = colour_coral, 
              alpha = 0.2) +
  labs(y = "Compositional Turnover", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  theme(legend.position = "none")

# save
ggsave(plot_turnover, filename = here("figures",
                                      "supplemental",
                                      "turnover.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)


# same for proxy data -----------------------------------------------------

# join with dataset
dat_final_proxy <- dat_proxy %>% 
  rename(bin = age) %>% 
  filter(bin %in% dat_spec$bin) %>% 
  mutate(temp_change = temp - lead(temp)) %>% 
  add_column(turnover = dat_diag[1:15]) 

# visualise
plot_turnover_proxy <- dat_final_proxy %>%
  mutate(temp = if_else(temp_change >= 0, "warm", "cool")) %>% 
  drop_na(temp) %>% 
  ggplot(aes(abs(temp_change), turnover)) +
  geom_hline(yintercept = 0, colour = "grey80", 
             linetype = "dotted") +
  geom_point() +
  geom_smooth(method = "lm", 
              colour = colour_coral, 
              fill = colour_coral, 
              alpha = 0.2) +
  labs(y = "Compositional Turnover", 
       x = expression(paste("Absolute " ,Delta, " Temperature [°C]"))) +
  theme(legend.position = "none")

# save
ggsave(plot_turnover_proxy, filename = here("figures",
                                      "supplemental",
                                      "turnover_proxy.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)

