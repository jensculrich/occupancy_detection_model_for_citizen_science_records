### Run occupancy model using real pollinator occurrence data from GBIF
# jcu; started nov 7, 2022

library(rstan)

source("./analysis/prep_data.R")

my_data <- prep_data(era_start = 2007, # must define start date of the GBIF dataset
                        era_end = 2022, # must define start date of the GBIF dataset
                        n_intervals = 5, # must define number of intervals to break up the era into
                        n_visits = 3 # must define the number of repeat obs years within each interval
                        # note, should introduce throw error if..
                        # (era_end - era_start) / n_intervals has a remainder > 0
                      )

# data to feed to the model
V <- my_data$V # detection data
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_intervals <- my_data$n_intervals # number of surveys 
n_visits <- my_data$n_visits

interval_names <- as.numeric(my_data$intervals)
site_names <- my_data$sites
species_names <- my_data$species

intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

stan_data <- c("V", 
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", "species", "sites")

# Parameters monitored
params <- c("mu_psi_0",
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
       mu_p_0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       sigma_p_site = runif(1, 0, 1),
       p_interval = runif(1, -1, 1)
       
  )
)

#targets <- as.data.frame(cbind(param_names, parameters))

## --------------------------------------------------
### Run model
stan_model <- "./simulation/model_simplest.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores)

print(stan_out_sim, digits = 3)

saveRDS(stan_out_sim, "./simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./simulation/stan_out_sim.rds")