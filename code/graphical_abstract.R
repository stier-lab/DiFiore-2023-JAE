#----------------------------------------------------------------------------
## Set up and functions
#----------------------------------------------------------------------------
source(here::here("code", "1_setup.R"))
source(here::here("code", "10a_clean-obsdata.R"))
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

ndraws = 10000

# Bring in the posteriors from the post hoc analysis of a and h ~ body size
# post.a <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv")) %>% sample_draws(10000)
# post.h <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv")) %>% sample_draws(10000)

# # eliminate parameter uncertainty
post.a <- read.csv(here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv")) %>%
  as_tibble() %>% 
  pivot_wider(names_from = .variable, values_from = .value) %>%
  summarize(alpha = median(alphaa), beta1 = median(beta1a), beta2 = median(beta2a))

post.h <- read.csv(here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv")) %>%
  as_tibble() %>% 
  pivot_wider(names_from = .variable, values_from = .value) %>%
  summarize(alpha = median(alphah), beta1 = median(beta1h), beta2 = median(beta2h))

s # This dataframe is a tibble with nested body size distributions for lobsters and urchins.

set.seed(112615)
r.s <- s %>% # reproducible random draws from the size frequency distribution
  group_by(year, site) %>%
  mutate(urc_mass = purrr::map(urc.mass, sample_n, ndraws, replace = T), 
         lob_mass = purrr::map(lob.mass, sample_n, ndraws, replace = T)) %>% 
  dplyr::select(-urc.mass, -lob.mass)

#----------------------------------------------------------------------------
## Simuations based on experimental data
#----------------------------------------------------------------------------

# This is a function for the allometric functional response:

allometricFR <- function(lob_mass, urc_mass, urc_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, ...){
  
  loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass)
  logh <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass)
  a <- exp(loga) * tsize * 2 # convert units from per arena per 2 days to per m2 per day 
  h <- exp(logh) / 2 # convert units from 2 day trials to 1 day. 
  a*urc_density*T / (1 + a*h*urc_density) #***
}


# Simulate interactions assuming all sources of uncertainty
full <- r.s %>%
  purrr::pmap(allometricFR,
              a0. = post.a$alpha,
              h0. = post.h$alpha,
              beta1a. = post.a$beta1,
              beta2a. = post.a$beta2,
              beta1h. = post.h$beta1,
              beta2h. = post.h$beta2) %>%
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction,
         estimate = "full")



full %>% 
  group_by(site, year) %>%
  summarize(mean = mean(prediction)) %>%
  ungroup() %>%
  filter(mean == min(mean) | mean == max(mean))



urc_mass <- s %>% 
  mutate(id = paste(site, year, sep = "_")) %>%
  filter(id %in% c("CARP_2014", "MOHK_2016")) %>%
  unnest(urc.mass) %>% 
  select(-lob.mass) %>%
  mutate(species = "urchin")

lob_mass <- s %>% 
  mutate(id = paste(site, year, sep = "_")) %>%
  filter(id %in% c("CARP_2014", "MOHK_2016")) %>%
  unnest(lob.mass) %>% 
  select(-urc.mass) %>%
  mutate(species = "lobster")

sizes <- 
  bind_rows(urc_mass, lob_mass)

p1 <- sizes %>%
  ggplot(aes(x = mass)) +
  geom_density(aes(fill = species), adjust = 1.5, show.legend = F)+
  scale_fill_manual(values = c("#c9345880", "#3458c980"))+
  scale_x_log10()+
  facet_wrap(~forcats::fct_rev(id))+
  theme_bd()

p2 <- full %>%
  mutate(id = paste(site, year, sep = "_")) %>%
  filter(id %in% c("CARP_2014", "MOHK_2016")) %>%
  ggplot(aes(x = prediction))+
  geom_density(aes(fill = id), adjust = 2, show.legend = F, stat = "density")+
  scale_fill_manual(values = c("#31a354bf", "#a1d99bbf"))+
  scale_x_log10(labels = plain)+
  theme_bd()

graphical <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v", axis = "lr")

ggsave("figures/graphical_abstract.svg", graphical)



