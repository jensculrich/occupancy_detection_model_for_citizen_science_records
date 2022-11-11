### Data simulation for citizen science/museum records of pollinator occurrence
# jcu; started oct 10, 2022

## Data simulation for simplest form of the data
# single season, species considered across all sites (not just range)
# and across all visits (not considering that sites without any detections 
# for any species were probably not visited or that sites with only detections for
# a few species only might be targeted and not fully community wide)

# model fitting using 'model_simplest.stan' should return the parameter inputs

simulate_data <- function(n_species,
                          n_sites,
                          n_intervals,
                          n_visits,
                          mu_psi_0,
                          sigma_psi_species,
                          mu_psi_interval,
                          sigma_psi_interval,
                          psi_pop_dens,
                          psi_site_area,
                          mu_p_0,
                          p_interval, 
                          sigma_p_site,
                          sigma_p_species
                          ){
  
  ## ilogit and logit functions
  ilogit <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) log(x/(1-x))
  
  ## visit
  # will need to come back to this portion once simple part is worked out
  # mu.v.0 = 0 
  # mu.v.interval = -0.5
  
  ## type sym
  # will need to come back to this portion once simple part is worked out
  # type.range = "equal"
  # type.visit = 'visit_miss'
  # prop.visits.same = 1
  
  
  # Interval values: numeric vector (will act as covariate data for psi.interval)
  intervals <- seq(1, n_intervals, by=1)
  intervals <- intervals - 1
  sites <- seq(1, n_sites, by=1)
  species <- seq(1, n_species, by=1)
  
  ## --------------------------------------------------
  ### Generate covariate data
  
  ## --------------------------------------------------
  ### Population density
  
  # create a vector of site population density of length = number of sites
  # the model takes z-score scaled data (with mean of 0) so it's ok to center at 0 here
  pop_density <- rnorm(n_sites, mean = 0, sd = 1)
  
  # create a vector of site area of length = number of sites
  # the model takes z-score scaled data (with mean of 0) so it's ok to center at 0 here
  # this simulates the realism that some sites are e.g. partially on ocean or  
  # partially outside the administrative area from which we are drawing collection data.
  site_area <- rnorm(n_sites, mean = 0, sd = 1)
  
  ## --------------------------------------------------
  ### specify species-specific occupancy probabilities
  
  ## species-specific random intercepts
  psi_species <- rnorm(n=n_species, mean=0, sd=sigma_psi_species)
  # species baseline occupancy is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma.psi.sp
  
  ## effect of interval on occupancy (species-specific random slopes)
  psi_interval <- rnorm(n=n_species, mean=mu_psi_interval, sd=sigma_psi_interval)
  # change in each species occupancy across time is drawn from a distribution defined
  # by a community mean (mu_psi_interval) with 
  # species specific variation defined by sigma_psi_interval
  
  ## --------------------------------------------------
  ### specify species-specific detection probabilities
  
  ## species-specific random intercepts
  p_species  <- rnorm(n=n_species, mean = 0, sd=sigma_p_species)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma.p.sp
  
  ## effect of site and interval on detection (site,interval-specific random slopes)
  p_site <- rnorm(n=n_sites, mean=0, sd=sigma_p_site)
  
  # spatiotemporal variability in detection probability (changing across sites and 
  # occupancy intervals), and helps account for the variation that is inherent in 
  # sample effort across space and time in unstructured historical datasets. 
  # spatiotemporal variability at each site*interval is drawn from a distribution defined
  # by a mean of 0 and with site*interval specific variation of sigma.p.site
  
  ## --------------------------------------------------
  ## Create arrays for psi and p
  psi_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals)) 
  # a psi value for each species, at each site, in each interval 
  
  p_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  
  for(species in 1:n_species) { # for each site
    for(site in 1:n_sites) { # for each interval
      for(interval in 1:n_intervals) { # for each species
        
        psi_matrix[species, site, interval] <- ilogit( # occupancy is equal to
          mu_psi_0 + # a baseline intercept
            psi_species[species] + # a species specific intercept
            psi_interval[species]*intervals[interval] + # a species specific temporal change
            psi_pop_dens*pop_density[site] + # a fixed effect of population density 
            psi_site_area*site_area[site] # a fixed effect of site area
        )
        
        for(visit in 1:n_visits) { # for each visit
          
          p_matrix[species, site, interval, visit] <- ilogit( # detection is equal to 
            mu_p_0 + # a baseline intercept
              p_species[species] + # a species specific intercept
              p_site[site] + # a spatiotemporally specific intercept
              p_interval*intervals[interval] # an overall effect of time on detection
          )
          
        } # for each visit
      } # for each species
    } # for each interval
  } # for each site
  
  # preview the psi and p arrays
  # head(psi_matrix[1:n_species, 1:n_sites,1])
  # head(p_matrix[1:n_species, 1:n_sites,7,1])
  
  # (p_matrix[1,1,1:n_intervals,1]) # if p.interval is >0 these should generally be increasing from low to high
  # (p_matrix[2,1,1:n_intervals,1]) # if p.interval is >0 these should generally be increasing from low to high
  # (p_matrix[1,2,1:n_intervals,1]) # if p.interval is >0 these should generally be increasing from low to high
  
  ## --------------------------------------------------
  ## Generate presence and absence
  
  Z <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        
        # eventually here we will need to specify which(site) to restrict to actual range
        Z[species,site,interval] <- rbinom(n = 1, size = 1, 
                                           prob = psi_matrix[species,site,interval])
        
      }
    }
  }
  
  ## --------------------------------------------------
  ## Generate detection non detection data 
  
  V <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals,
                       n_visits=n_visits))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        for(visit in 1:n_visits){
          
          # eventually here we will need to specify which(site) to restrict to actual range
          V[species,site,interval,visit] <- Z[species,site,interval] * # occupancy state * detection prob
            rbinom(n = 1, size = 1, prob = p_matrix[species,site,interval,visit])
          
        }
      }
    }
  } # end simulate detection data
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    V = V, # detection data
    n_species = n_species, # number of species
    n_sites = n_sites, # number of sites
    n_intervals = n_intervals, # number of surveys 
    n_visits = n_visits, # number of visits
    pop_density = pop_density, # vector of pop densities
    site_area = site_area # vector of site areas
  ))
  
} # end simulate_data function


## --------------------------------------------------
### Variable values for data simulation
## study dimensions
n_species = 20 ## number of species
n_sites = 15 ## number of sites
n_intervals = 4 ## number of occupancy intervals
n_visits = 3 ## number of samples per year

## occupancy
mu_psi_0 = -0.5
sigma_psi_species = 0.5
mu_psi_interval = 0.5
sigma_psi_interval = 0.2
psi_pop_dens = -0.5 # fixed effect of population density on occupancy
psi_site_area = 1 # fixed effect of site area on occupancy

## detection
mu_p_0 = -0.5
p_interval = 0.25 # detection probability increasing with time
sigma_p_site = 0.3
sigma_p_species = 0.3

## --------------------------------------------------
### Simulate data

my_simulated_data <- simulate_data(n_species,
                                   n_sites,
                                   n_intervals,
                                   n_visits,
                                   mu_psi_0,
                                   sigma_psi_species,
                                   mu_psi_interval,
                                   sigma_psi_interval,
                                   psi_pop_dens,
                                   psi_site_area,
                                   mu_p_0,
                                   p_interval, 
                                   sigma_p_site,
                                   sigma_p_species)

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V <- my_simulated_data$V # detection data
n_species <- my_simulated_data$n_species # number of species
n_sites <- my_simulated_data$n_sites # number of sites
n_intervals <- my_simulated_data$n_intervals # number of surveys 
n_visits <- my_simulated_data$n_visits

intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

pop_densities <- my_simulated_data$pop_density
site_areas <- my_simulated_data$site_area

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

parameter_value <- c(mu_psi_0,
                      sigma_psi_species,
                      mu_psi_interval,
                      sigma_psi_interval,
                      psi_pop_dens,
                      psi_site_area,
                      mu_p_0,
                      sigma_p_species,
                      sigma_p_site,
                      p_interval
)

# MCMC settings
n_iterations <- 600
n_thin <- 1
n_burnin <- 300
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

targets <- as.data.frame(cbind(params, parameter_value))

## --------------------------------------------------
### Run model
library(rstan)
stan_model <- "./models/model0.stan"

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
View(targets)

saveRDS(stan_out_sim, "./simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./simulation/stan_out_sim.rds")

library(shinystan)
shinystan::launch_shinystan(stan_out_sim)
