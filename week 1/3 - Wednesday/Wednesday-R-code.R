library(sp)
library(tidyverse)
library(cowplot)

# with the function from the Thompson model
install.packages('~/Downloads/RandomFields_3.3.14.tar.gz', repos = NULL)

devtools::install_github("plthompson/mcomsimr")

library(mcomsimr)

S = 10 # number of species
env_traits(species = S,
           max_r = 5,                         # max growth rate
           min_env = 0, max_env = 1,          # range of environment
           env_niche_breadth = 0.05,       # niche breadth
           optima_spacing = 'even')           # spacing of z_i (optima)


##--------simulate metacommumity dynamics with the Thompson model
# Number of patches
P = 20
# Number of species
S = 20

# set up landscape for experiment (want the same landscape for both treatments)
meta_landscape <- landscape_generate(patches = P)

# dispersal matrix
# set rate of distance decay
d = 0.1
disp_mat <- dispersal_matrix(landscape = meta_landscape,
                             torus = TRUE,
                             kernel_exp = d,
                             plot = TRUE)

# generate the time series of the environmental conditions for each patch (same for each treatment)
# Project idea: How does temporal environmental autocorrelation impact population, community or metacommunity dynamics / patterns?
env_conditions <- env_generate(landscape = meta_landscape,
                               env1Scale = 1, # temporal autocorrelation in the environment (here the environment is temporally uncorrelated)
                               timesteps = 1000)

# density independent component of model
densInd_niche <- env_traits(species = S,
                            max_r = 5,                         # max growth rate
                            min_env = 0, max_env = 1,          # range of environment
                            env_niche_breadth = 0.2,       # niche breadth
                            optima_spacing = 'even')           # spacing of z_i (optima)

# species interaction matrix:
# Project idea: examine the impact of alternate competition scenarios
# e.g., create  matrices for other dynamics in the paper (mixed, competitive dominance, destabilising competition)
equal_interaction_mat <- species_int_mat(species = S,
                                         intra = 1,
                                         min_inter = 1,
                                         max_inter = 1)

stabilising_interaction_mat <- species_int_mat(species = S,
                                               intra = 1,
                                               min_inter = 0,
                                               max_inter = 0.8)

# use simulateMC() function to simulate dynamics
sim_equal_comp <- simulate_MC(patches=P, species=S,
                              dispersal = 0.1,
                              landscape = meta_landscape,
                              disp_mat = disp_mat,
                              env.df = env_conditions,
                              max_r = densInd_niche$max_r,
                              env_niche_breadth = densInd_niche$env_niche_breadth,
                              env_optima = densInd_niche$optima,
                              int_mat = equal_interaction_mat,
                              initialization=100, burn_in=300, timesteps=600)
?simulate_MC

sim_stabil_comp <- simulate_MC(patches=P, species=S,
                               dispersal = 0.1,
                               landscape = meta_landscape,
                               disp_mat = disp_mat,
                               env.df = env_conditions,
                               max_r = densInd_niche$max_r,
                               env_niche_breadth = densInd_niche$env_niche_breadth,
                               env_optima = densInd_niche$optima,
                               int_mat = stabilising_interaction_mat,
                               initialization=100, burn_in=300, timesteps=600)

str(sim_equal_comp)

# extract data
sim_equalC_dat <- sim_equal_comp$dynamics.df %>%
  as_tibble() %>%
  # reduce to last 100 time steps
  filter(time > 499 & time < 601)


sim_stabilC_dat <- sim_stabil_comp$dynamics.df %>%
  as_tibble() %>%
  # reduce to last 100 time steps
  filter(time > 499 & time < 601)

# visual inspection

# environmental conditions:
ggplot() +
  geom_line(data = sim_equalC_dat,
            aes(x = time, y = env, colour = env,
                group = patch)) +
  scale_colour_viridis_c()

# patch level population dynamics: equal comp
ggplot() +
  facet_wrap(~patch) +
  geom_line(data = sim_equalC_dat,
            aes(x = time, y = N, colour = optima,
                group = interaction(species, patch))) +
  scale_colour_viridis_c()

ggplot() +
  facet_wrap(~patch) +
  geom_line(data = sim_stabilC_dat,
            aes(x = time, y = N, colour = optima,
                group = interaction(species, patch))) +
  scale_colour_viridis_c()

# Exercises:
# Simulate dynamics according to the classic metacommunity paradigms (Fig 2 in paper) and plot
