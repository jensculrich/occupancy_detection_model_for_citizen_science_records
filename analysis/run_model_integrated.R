## run_model_integrated.R
### Run occupancy model (model_integrated.stan) using real pollinator occurrence data from GBIF
# jcu; started nov 24, 2022

library(rstan)

## --------------------------------------------------
# input data preparation choices
# be careful that the (era_end - era_start) is evenly divisible by the n_intervals
era_start = 2009 # must define start date of the GBIF dataset
era_end = 2020 # must define start date of the GBIF dataset
n_intervals = 3 # must define number of intervals to break up the era into
n_visits = 4 # must define the number of repeat obs years within each interval
# note, should introduce throw error if..
# (era_end - era_start) / n_intervals has a remainder > 0,
min_records_per_species = 100
grid_size = 25000 # 25km x 25 km 
min_population_size = 100 # min pop density in the grid cell (per km^2)
# for reference, 38people/km^2 is ~100people/mile^2
# 100/km^2 is about 250/mile^sq
min_species_for_community_sampling_event = 3
# min records_for_community_sampling_event sets a minimum threshold, if the number
# of records for the taxonomic group within a site within a year is 

source("./analysis/prep_data_integrated.R")
my_data <- prep_data(era_start = era_start, # must define start date of the GBIF dataset
                     era_end = era_end, # must define start date of the GBIF dataset
                     n_intervals = n_intervals, # must define number of intervals to break up the era into
                     n_visits = n_visits, # must define the number of repeat obs years within each interval
                     # note, should introduce throw error if..
                     # (era_end - era_start) / n_intervals has a remainder > 0,
                     min_records_per_species = min_records_per_species,
                     grid_size = grid_size, # 25km x 25 km 
                     min_population_size = min_population_size, # min pop density in the grid cell (per km^2)
                     # for reference, 38people/km^2 is ~100people/mile^2
                     # 100/km^2 is about 250/mile^sq
                     min_records_for_community_sampling_event = min_records_for_community_sampling_event
)

# data to feed to the model
V <- my_data$V # detection data
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_intervals <- my_data$n_intervals # number of surveys 
n_visits <- my_data$n_visits

interval_names <- as.vector(as.numeric(my_data$intervals))
site_names <- my_data$sites
species_names <- my_data$species

city_names <- my_data$city_name_vector

pop_densities <- my_data$pop_densities
site_areas <- my_data$site_areas

# intervals will cause issues if you try to run on only 1 interval
# since it's no longer sent in as a vector of intervals (can you force a single
# integer to be a vector if you truly want to treat all as a single interval?)
intervals_raw <- as.vector(seq(1, n_intervals, by=1)) 
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

stan_data <- c("V", 
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", "species", "sites",
               "pop_densities", "site_areas")

# Parameters monitored
params <- c("mu_psi_0",
            #"psi_species",
            "sigma_psi_species",
            # "psi_interval",
            "mu_psi_interval",
            "sigma_psi_interval",
            "psi_pop_density",
            "psi_site_area",
            "mu_p_0",
            # "p_species",
            "sigma_p_species",
            # "p_site",
            "sigma_p_site",
            "p_interval"
)


# MCMC settings
n_iterations <- 400
n_thin <- 1
n_burnin <- 200
n_chains <- 3
n_cores <- n_chains

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(mu_psi_0 = runif(1, -1, 1),
       sigma_psi_species = runif(1, 0, 1),
       mu_psi_interval = runif(1, -1, 1),
       sigma_psi_interval = runif(1, 0, 1),
       psi_pop_dens = runif(1, -1, 1),
       psi_site_area = runif(1, -1, 1),
       mu_p_0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       sigma_p_site = runif(1, 0, 1),
       p_interval = runif(1, -1, 1)
       
  )
)

## --------------------------------------------------
### Run model
stan_model <- "./models/model0.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 seed = 1,
                 open_progress = FALSE,
                 cores = n_cores)

print(stan_out, digits = 3)

saveRDS(stan_out, "./model_outputs/stan_out_model0.rds")
stan_out <- readRDS("./model_outputs/stan_out_model0.rds")

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out, pars = c(
  "mu_psi_0",
  #"psi_species",
  "sigma_psi_species",
  # "psi_interval",
  "psi_pop_density",
  "psi_site_area",
  "mu_psi_interval",
  "sigma_psi_interval",
  "mu_p_0",
  # "p_species",
  "sigma_p_species",
  # "p_site",
  "sigma_p_site",
  "p_interval"
))

# pairs plot
pairs(stan_out, pars = c(
  "mu_psi_0",
  #"psi_species",
  "sigma_psi_species",
  # "psi_interval",
  "mu_psi_interval",
  "sigma_psi_interval",
  "mu_p_0",
  # "p_species",
  "sigma_p_species",
  # "p_site",
  "sigma_p_site",
  "p_interval"
))

# should now also write a posterior predictive check into the model
