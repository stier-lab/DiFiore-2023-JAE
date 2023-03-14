#----------------------------------------------------------------------------
## Set up and functions
#----------------------------------------------------------------------------
source(here::here("code", "1_setup.R"))
source(here::here("code", "5_clean-obsdata.R"))
source(here::here("code", "theme.R"))


# Add line for merge conflict resolution.

# Misc. numbers needed later in code
# Tank size

# Full tank is 137 x 76. I'm going to stick with areal measures (rather than volumetric) because LTER works off of an areal basis. Therefore experimental tanks were 68.5 x 76. 

tsize <- 137/2 * 76 /10000

names <- s %>% mutate(id = paste(year, site, sep = "-")) %>%
  ungroup() %>%
  select(id) %>%
  distinct(id)

names <- as.vector(names$id)



plain <- function(x,...) {
  format(x, ..., scientific = FALSE, trim = TRUE, drop0trailing = T)
}

#---------------------------------
## Data 
#---------------------------------

ndraws = 1000 # reduced the number of draws from the size 

# Bring in the posteriors from the post hoc analysis of a and h ~ body size
post.a <- read.csv(here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv")) %>%
  filter(.variable %in% c("alphaa", "beta1a", "beta2a")) %>% 
  sample_draws(1000) %>%
  pivot_wider(names_from = .variable, values_from = .value)

post.h <- read.csv(here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv")) %>%
  filter(.variable %in% c("alphah", "beta1h", "beta2h")) %>% 
  sample_draws(1000)%>%
  pivot_wider(names_from = .variable, values_from = .value)


s # This dataframe is a tibble with nested body size distributions for lobsters and urchins.

set.seed(112615)
r.s <- s %>% # reproducible random draws from the size frequency distribution
  group_by(year, site) %>%
  mutate(urc_mass = purrr::map(urc.mass, sample_n, ndraws, replace = T), 
         lob_mass = purrr::map(lob.mass, sample_n, ndraws, replace = T)) %>% 
  dplyr::select(-urc.mass, -lob.mass)

#----------------------------------------------------------------------------
## Simulations based on experimental data
#----------------------------------------------------------------------------

# This is a function for the allometric functional response:

allometricFR <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, ...){
  
  loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass)
  logh <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass)
  a <- exp(loga) * tsize * 2 # convert units from per arena per 2 days to per m2 per day 
  h <- exp(logh) / 2 # convert units from 2 day trials to 1 day. 
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
}



# Simulate interactions assuming all sources of uncertainty
full_wparameter <- r.s %>%
  purrr::pmap(allometricFR,
              a0. = post.a$alphaa,
              h0. = post.h$alphah,
              beta1a. = post.a$beta1a,
              beta2a. = post.a$beta2a,
              beta1h. = post.h$beta1h,
              beta2h. = post.h$beta2h) %>%
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction,
         estimate = "full_wparameter")

# Summary stats for paper

summary(full_wparameter$prediction)
rethinking::PI(full_wparameter$prediction, 0.95)


full <- r.s %>%
  purrr::pmap(allometricFR,
              a0. = median(post.a$alphaa),
              h0. = median(post.h$alphah),
              beta1a. = median(post.a$beta1a),
              beta2a. = median(post.a$beta2a),
              beta1h. = median(post.h$beta1h),
              beta2h. = median(post.h$beta2h)) %>%
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction,
         estimate = "full")

# Summary stats for paper

summary(full$prediction)
rethinking::PI(full$prediction, 0.95)




bind_rows(full,full_wparameter) %>%
  group_by(estimate) %>%
  tidybayes::mean_qi(prediction) %>%
  ggplot(aes(y = estimate, x = prediction))+
  tidybayes::geom_pointinterval(aes(color = estimate, xmin = .lower, xmax = .upper), show.legend = F )+
  scale_x_log10(labels = plain)+
  theme_bd()





#-------------------------------------------------------------------------
## Estimates from the literature
#-------------------------------------------------------------------------

#-----------------------------        
# Rall
# Rall et al. 2012 incorporate temp as a covariate. The basic regression equations are: 
# log(a) = a0 + beta1*log(Mc) + beta2*log(Mr) + beta3*(T-T0)/k*T*T0
# log(h) = h0 + beta1*log(Mc) + beta2*log(Mr) + beta3*(T-T0)/k*T*T0
# where T is temp in kelvin, T0 = 293.15, mass are mg, and time is seconds, k is the boltzman's constant

bz_ev_k = 8.61733326 * 10^-5

allometricRall <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., beta3a., beta3h., T = 1, ...){
  
  urc_mass_mg = urc_mass*1000
  lob_mass_mg = lob_mass*1000
  temp = (289.25 - 293.15)/(bz_ev_k*289.25*293.15)
  loga <- a0. + beta1a.*log(lob_mass_mg) + beta2a.*log(urc_mass_mg) + beta3a.*temp
  logh <- h0. + beta1h.*log(lob_mass_mg) + beta2h.*log(urc_mass_mg) + beta3h.*temp
  a <- exp(loga) 
  h <- exp(logh) 
  consumption_mg_m2_s = a*urc_density*lob_density*T / (1 + a*h*urc_density)
  consumption_g_m2_day = consumption_mg_m2_s * (60*60*24) / 1000 # convert mg per m2 per second to g per m2 per day.
  consumption_g_m2_day
  
}

rall <- r.s %>%
  purrr::pmap(allometricRall, 
              beta1a. = 0.85, 
              beta2a. = 0.09, 
              beta1h. = -0.76, 
              beta2h. = 0.76, 
              h0. = 10.38, 
              a0. = -21.23, 
              beta3a. = 0.42, 
              beta3h. = -0.30) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction,
         estimate = "rall")

summary(rall$prediction)        
hist(rall$prediction)        

#-----------------------------




#-----------------------------
# Barrios-Oneil
# BO et al. 2019 use a different model structure with a temperature covariate
# basic model structure: 
# log(a) = a_0 + beta1*I_ACstat + beta2*I_F + beta3*log(Mc) + beta4*log(Mr) + beta5*log(T) + beta6*log(Mr)*log(T)
# log(cmax) = C_0 + beta1*I_ACstat + beta2*I_F + beta3*log(Mc) + beta4*log(Mr) + beta5*log(T) + beta6*log(Mc)*log(T)
# I assume that I_F = 0 and I_ACstat = 1, such that we have the regression equation for the active predators preying on static prey.
# All beta values from Fig. 4 in text and Table S5. Intercepts estimated as the population level intercept + adjustment for active predators foraging on static prey (AC_stat) + the random intercept deviation from the mean for crustaceans. 
# Units in the paper at g, days, and m^2. So no need for unit conversions.

allometricBO <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, temp, beta5a., beta6a., beta5h., beta6h., ...){
  
  loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass) + beta5a.*log(temp) + beta6a.*log(urc_mass)*log(temp)
  logC <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass) + beta5h.*log(temp) + beta6h.*log(lob_mass)*log(temp)
  a <- exp(loga)
  C <- exp(logC)
  h <- 1/C
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}

BO <- r.s %>%
  purrr::pmap(allometricBO, 
              beta1a. = 0.58, 
              beta2a. = 0.59, 
              beta1h. = 1.44, 
              beta2h. = 0.27,
              temp = 16.1, 
              beta5a. = 1.84,
              beta6a. = 0.12,
              beta5h. = 1.78, 
              beta6h. = -0.39,
              h0. = -6.07 + 0.55 + 0.48, 
              a0. = -8.08 + -1.07) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "BO")
#-----------------------------------------

#-----------------------------------------
# Uiterwall and Delong
# UD 2020 account for arena size, dimensionality, and temp. The model structures are:
# log(a) = a0 + beta1*T + beta2*T^2 + beta3*Mc + beta4*Mr + beta5*arena
# log(h) = h0 + beta1*T + beta2*T^2 + beta3*Mc + beta4*Mr + beta5*Dim + beta6*arena
# Values for all beta values and intercepts are from Table 1 and Table 2 in the main text of the manuscript.
# I assume that temperature = the mean temperature over our experiment and arena size = our arena size.
# Units in the manuscript are mg, cm^2, and days

allometricUD <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., beta3a., beta3h., beta4a., beta4h., beta5a., beta5h., beta6h., temp, arena_m2, dim, T = 1, ...){
  lob_mass_mg = lob_mass*1000 # convert g to mg
  urc_mass_mg = urc_mass*1000 # convert g to mg
  arena_cm2 = arena_m2 * 10^4 # convert m2 to cm2
  urc_density_cm2 = urc_density / 10^4 # convert ind./ m2 to ind./ cm2
  lob_density_cm2 = lob_density / 10^4 # convert ind./m2 to ind./cm2
  loga <- a0. + beta1a.*log(lob_mass_mg) + beta2a.*log(urc_mass_mg) + beta3a.*temp + beta4a.*temp^2 + beta5a.*log(arena_cm2)
  logh <- h0. + beta1h.*log(lob_mass_mg) + beta2h.*log(urc_mass_mg) + beta3h.*temp + beta4h.*temp^2 + beta5h.*log(arena_cm2) + beta6h.*dim 
  a <- exp(loga) # per cm2 per day
  h <- exp(logh) # days
  consumption_mg_cm2_d = a*urc_density_cm2*lob_density_cm2*T / (1 + a*h*urc_density_cm2)
  consumption_g_cm2_d = consumption_mg_cm2_d / 1000 # convert mg to g
  consumption_g_m2_d = consumption_g_cm2_d * 10^4 # convert g per cm2 per day to g per m2 per day
  consumption_g_m2_d
}



# Scrap to think about how to resample from the parameter estimates given mean and SE.
dist <- rnorm(1000, mean = 0, sd = 500)
hist(dist)
mean <- mean(dist)
se_mean <- sd(dist)/sqrt(1000)

dist.95 <- qnorm(0.975)*se_mean

df.dist <- data.frame(x = dist)


ggplot(df.dist)+
  geom_histogram(aes(x = x))+
  geom_vline(xintercept = c(mean, mean -dist.95, mean+dist.95 ))



# Here we simulate the prediction based on UD where we account for parameter uncertainty. To do this we resample fro the 95% CI's of the parameter estimates assuming a uniform distribution. 

helper <- function(mean, se, ndraws = 1000){
  low = mean - qnorm(0.975)*se
  high = mean + qnorm(0.975)*se
  runif(ndraws, min = low, max = high)
}

UD_wparameter <- r.s %>%
  purrr::pmap(allometricUD, 
              beta1a. = helper(0.05, 0.03), 
              beta2a. = helper(-0.0005, 0.02), 
              beta1h. = helper(-0.25, 0.03), 
              beta2h. = helper(0.34, 0.02), 
              h0. = helper(0.83, 0.85), 
              a0. = helper(-8.45, 1.05), 
              beta3a. = helper(0.1, 0.03), 
              beta3h. = helper(-0.24, 0.04), 
              beta4a. = helper(-0.003, 0.001),
              beta4h. = helper(0.005, 0.001),
              beta5a. = helper(0.98, 0.05), 
              beta5h. = helper(-0.01, 0.03), 
              beta6h. = helper(-0.64, 0.19), 
              temp = 16.1, 
              arena_m2 = tsize, 
              dim = 2) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "UD_wparameter")

summary(UD_wparameter$prediction)
hist(UD_wparameter$prediction)




uncertainty <- bind_rows(full, full_wparameter, UD_wparameter, rall, BO, UD) %>%
  group_by(estimate) %>%
  tidybayes::mean_qi(prediction) %>%
  mutate(estimate = case_when(estimate == "full" ~ "Experimental", 
                              estimate == "full_wparameter" ~ "Exp. w/ uncertainty", 
                              estimate == "UD_wparameter" ~ "UD w/ uncertainty", 
                              estimate == "rall" ~ "Rall", 
                              estimate == "UD" ~ "UD", 
                              estimate == "BO" ~ "BO")) %>%
  ggplot(aes(y = estimate, x = prediction))+
  tidybayes::geom_pointinterval(aes(color = estimate, xmin = .lower, xmax = .upper), show.legend = F )+
  scale_x_log10(labels = plain)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "")+
  theme_bd()

ggsave("figures/fig_uncertainty.png", uncertainty, device = "png", width = 6, height = 2.5)


