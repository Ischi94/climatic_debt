library(here)
library(tidyverse)
library(lme4)


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
mod1 <- lmer(temp_debt ~ temp_change + (1 | bin),
             data = dat_trends)

# summarize beta coefficient
bootMer(mod1, fixef, nsim = 1000)


# create equal grid for inference (accounting for differential sampling in bins)
new_data <- tibble(temp_change = seq(-3, 3, by = 0.1))  

# estimate over this equal grid  
dat_trends_pred <- predict(mod1, newdata = new_data,
                           re.form = ~ 0,
                           allow.new.levels = TRUE) %>%
  as_tibble() %>% 
  add_column(new_data) %>% 
  rename(predicted_debt = value)


# in total
plot_trends_comp <- dat_trends_pred %>% 
  ggplot(aes(temp_change, predicted_debt)) +
  geom_vline(xintercept = 0, 
             colour = "grey70", 
             linetype = "dotted") +
  # add the slope based on relative abundance
  geom_abline(intercept = -0.214, slope = 0.558, 
              colour = colour_coral, 
              alpha = 0.8, lwd = 1) +
  geom_line(colour = colour_grey, 
            lwd = 1) +
  annotate(geom = "text", 
           x = -1.4, y = 0, 
           label = "Relative abundance", 
           angle = 5, 
           colour = "grey20",
           size = 10/.pt) +
  annotate(geom = "text", 
           x = -1, y = -7.8, 
           label = "Occurrence only", 
           colour = "grey20",
           angle = 25, 
           size = 10/.pt) +
  labs(y = "Climatic Debt [°C]", 
       x = expression(paste(Delta, "  Temperature [°C]"))) +
  scale_x_continuous(breaks = seq(-2, 2, 2)) +
  coord_cartesian(xlim = c(-2.5, 2.5))


# save plot
ggsave(plot_trends_comp, filename = here("figures", 
                                   "supplemental",
                                   "debt_trends_comparison.png"), 
       width = image_width, height = image_height, units = image_units, 
       bg = "white", device = ragg::agg_png)



# latitudinal wise
# split data into latitudinal zones and then apply mixed effect models,
# accounting for sampling
dat_mod <- dat_trends %>% 
  group_by(short_term, zone) %>%
  nest() %>% 
  mutate(lm_mod = map(.x = data,
                      .f = function(df) {
                        lmer(temp_debt ~ temp_change + (1 | bin),
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
                 "beta_coefficient_per_latitude_occurrence.csv"))

