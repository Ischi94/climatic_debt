library(here)
library(tidyverse)

# source plotting configurations ------------------------------------------

source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

dat_debt <- read_rds(here("data",
                          "cleaned_debt_wapls.rds")) %>% 
  # convert temperature in temperature anomaly
  mutate(temp_anom = (temp_surface - mean(temp_surface, na.rm = TRUE)) / 
           sd(temp_surface, na.rm = TRUE)) %>% 
  # prepare trend data
  mutate(short_term = if_else(temp_anom > lead(temp_anom),
                              "warming",
                              "cooling"), 
         temp_change = temp_anom - lead(temp_anom, 
                                        default = mean(temp_anom))) %>% 
  drop_na(short_term)


# temperature estimates from proxy data
dat_proxy <- read_rds(here("data", 
                           "proxy_temperature.rds")) %>% 
  # nas are encode as -999
  filter(temp != -999)



# time average sample uncertainty -----------------------------------------


# add sampling buffer
dat_sampling <- distinct(dat_debt, bin) %>% 
  arrange(bin) %>% 
  mutate(lead_bin = lead(bin)) %>% 
  drop_na(lead_bin) %>% 
  mutate(age = map2(.x = bin, 
                    .y = lead_bin, 
                    .f = ~ seq(.x, .y, by = 1))) %>% 
  unnest(c(age)) %>% 
  full_join(dat_proxy) %>% 
  select(-c(lead_bin, age)) %>% 
  nest_by(bin) %>% 
  ungroup() %>% 
  full_join(dat_debt %>% 
              group_by(bin) %>% 
              summarise(climatic_debt = mean(climatic_debt))) %>% 
  unnest(data)

# create thousand unique datasets and loop through
dat_beta_proxy <- vector("double",
                         length = 1000)
dat_interc_proxy <- vector("double",
                           length = 1000)

for (i in 1:1000) {
  
  set.seed(i)
  
  mod <- dat_sampling %>% 
    group_by(bin) %>% 
    slice_sample(n = 1) %>% 
    ungroup() %>% 
    mutate(temp_change = temp - lead(temp)) %>% 
    lm(climatic_debt ~ temp_change, data = .) 
  
  dat_beta_proxy[i] <- mod %>% 
    coefficients() %>% 
    pluck(2)
  
  dat_interc_proxy[i] <- mod %>% 
    coefficients() %>% 
    pluck(1)
  
}

# get model reported in the manuscript
mod1 <- lm(climatic_debt ~ temp_change, data = dat_debt) 

# visualize
plot_proxy <- tibble(beta_coef = dat_beta_proxy,
                     inter_coef = dat_interc_proxy) %>%
  ggplot(aes(x = temp_change, 
             y = climatic_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey40", 
             linetype = "dotted") +
  geom_line(data = tibble(climatic_debt = seq(-3, 3,
                                              length.out = 20),
                          temp_change = seq(-3, 3,
                                            length.out = 20)),
            colour = "white") +
  geom_hline(yintercept = 0, 
             colour = "grey40", 
             linetype = "dotted") +
  geom_abline(aes(slope = beta_coef, intercept = inter_coef), 
              alpha = 0.05) +
  geom_abline(slope = coef(mod1) %>% pluck(2), 
              intercept = coef(mod1) %>% pluck(2), 
              colour = colour_coral, 
              linewidth = 1) +
  labs(y = "Climatic Mismatch [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  scale_y_continuous(breaks = seq(-2, 2, 2)) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(ylim = c(-3, 3), 
                  xlim = c(-2.5, 2.5)) +
  theme(axis.ticks = element_blank())

# save
ggsave(plot_proxy, filename = here("figures",
                                   "supplemental",
                                   "timeaveraging_proxy.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)

