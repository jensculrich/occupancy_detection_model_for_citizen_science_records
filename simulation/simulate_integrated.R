### Data simulation for citizen science/museum records of pollinator occurrence
# jcu; started nov 21, 2022

## Data simulation for integrated model form
# single season, species considered across all sites (not just range)
# and across all visits (not considering that sites without any detections 
# for any species were probably not visited or that sites with only detections for
# a few species only might be targeted and not fully community wide)

# model fitting using 'model_simplest.stan' should return the parameter inputs
simulate_data <- function(n_species,
                          n_sites,
                          n_intervals,
                          n_visits,
                          
                          ## ecological process
                          mu_psi_0,
                          mu_psi_interval,
                          sigma_psi_interval,
                          psi_pop_dens,
                          psi_site_area,
                          
                          ## observation process
                          # citizen science observation process
                          mu_p_citsci_0,
                          p_citsci_species,
                          sigma_p_citsci_species,
                          p_citsci_site,
                          sigma_p_citsci_site,
                          p_citsci_interval,
                          p_citsci_pop_density, 
                          
                          # museum record observation process
                          mu_p_museum_0,
                          p_museum_species,
                          sigma_p_museum_species,
                          p_museum_site,
                          sigma_p_museum_site,
                          p_museum_interval,
                          p_museum_pop_density, 
                          
                          introduce_NAs,
                          sites_missing,
                          intervals_missing,
                          visits_missing
){
  
  ## ilogit and logit functions
  ilogit <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) log(x/(1-x))
  
  
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
  p_citsci_species  <- rnorm(n=n_species, mean = 0, sd=sigma_p_citsci_species)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p_species
  
  ## species-specific random intercepts
  p_museum_species  <- rnorm(n=n_species, mean = 0, sd=sigma_p_museum_species)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p_species
  
  ## effect of site and interval on detection (site,interval-specific random slopes)
  p_citsci_site <- rnorm(n=n_sites, mean=0, sd=sigma_p_citsci_site)
  
  ## effect of site and interval on detection (site,interval-specific random slopes)
  p_museum_site <- rnorm(n=n_sites, mean=0, sd=sigma_p_museum_site)
  
  # spatiotemporal variability in detection probability (changing across sites and 
  # occupancy intervals), and helps account for the variation that is inherent in 
  # sample effort across space and time in unstructured historical datasets. 
  # spatiotemporal variability at each site*interval is drawn from a distribution defined
  # by a mean of 0 and with site*interval specific variation of sigma.p.site
  
  ## --------------------------------------------------
  ## Create arrays for psi and p
  psi_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals)) 
  # a psi value for each species, at each site, in each interval 
  
  p_matrix_citsci <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  p_matrix_museum <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
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
          
          p_matrix_citsci[species, site, interval, visit] <- ilogit( # detection is equal to 
            mu_p_citsci_0 + # a baseline intercept
              p_citsci_species[species] + # a species specific intercept
              p_citsci_site[site] + # a spatiotemporally specific intercept
              p_citsci_interval*intervals[interval] + # an overall effect of time on detection
              p_citsci_pop_density*pop_density[site] # an effect of population density on detection ability
          )
          
          p_matrix_museum[species, site, interval, visit] <- ilogit( # detection is equal to 
            mu_p_museum_0 + # a baseline intercept
            p_museum_species[species] + # a species specific intercept
            p_museum_site[site] + # a spatiotemporally specific intercept
            p_museum_interval*intervals[interval] + # an overall effect of time on detection
            p_museum_pop_density*pop_density[site] # an effect of population density on detection ability
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
  
  V_citsci <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals,
                       n_visits=n_visits))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        for(visit in 1:n_visits){
          
          # eventually here we will need to specify which(site) to restrict to actual range
          V_citsci[species,site,interval,visit] <- Z[species,site,interval] * # occupancy state * detection prob
            rbinom(n = 1, size = 1, prob = p_matrix_citsci[species,site,interval,visit])
          
        }
      }
    }
  } # end simulate detection data
  
  V_museum <- array(NA, dim=c(n_species=n_species,
                              n_sites=n_sites,
                              n_intervals=n_intervals,
                              n_visits=n_visits))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        for(visit in 1:n_visits){
          
          # eventually here we will need to specify which(site) to restrict to actual range
          V_museum[species,site,interval,visit] <- Z[species,site,interval] * # occupancy state * detection prob
            rbinom(n = 1, size = 1, prob = p_matrix_museum[species,site,interval,visit])
          
        }
      }
    }
  } # end simulate detection data
  
  sum(V_citsci == 1)
  sum(V_museum == 1)
  sum(V_citsci == 0)
  sum(V_museum == 0)
  
  # simulate NA data for the model to work around
  # for our real data we will have sites that weren't sampled during some visits in some intervals
  if(introduce_NAs == TRUE){
    
    # choose random sites that didn't get visited (for all species) by musuem collecting visits
    site_missed = sample.int(n_sites, sites_missing)
    interval_missed = sample.int(n_intervals, intervals_missing)
    visit_missed = sample.int(n_visits, visits_missing)
    
    # we will make an array that holds values of 1 if sampling occurred
    # or 0 if sampling did not occur at the site*interval*visit.
    V_museum_NA <- V_museum 
    V_museum_NA[1:n_species, site_missed, interval_missed, visit_missed] <- NA
    
    # replace all other values with 1 (was sampled)
    V_museum_NA <- replace(V_museum_NA, V_museum_NA==0, 1)
    # and now replace all NAs with 0, which will act as an indicator for the likelihood function
    # to skip over this sample by contracting the total possible number of observations that could have occurred
    V_museum_NA[is.na(V_museum_NA)] <- 0
    
    # now we want to replace the real data with 0's where sampling did not occur
    # so that we are saying that a a species was not observed at the visit to the site in the interval
    # but the model will remove this from contributing to the probability density by removing
    # the max number of sightings that could have occurred for 
    # each visit in the site*interval with a 0 in V_museum_NA
    V_museum[1:n_species, site_missed, interval_missed, visit_missed] <- 0
    
  }
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    V_citsci = V_citsci, # detection data from citizen science records
    V_museum = V_museum, # detection data from museum records
    V_museum_NA = V_museum_NA, # array indicating whether sampling occurred in a site*interval*visit
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
n_sites = 20 ## number of sites
n_intervals = 3 ## number of occupancy intervals
n_visits = 6 ## number of samples per year

## occupancy
mu_psi_0 = -0.5
sigma_psi_species = 0.5
mu_psi_interval = 0.5
sigma_psi_interval = 0.2
psi_pop_dens = -0.5 # fixed effect of population density on occupancy
psi_site_area = 1 # fixed effect of site area on occupancy

## detection
# citizen science observation process
mu_p_citsci_0 = -0.25
p_citsci_species = 0
sigma_p_citsci_species = 0.5
p_citsci_site = 0
sigma_p_citsci_site = 0.3
p_citsci_interval = 1
p_citsci_pop_density = 1 

# museum record observation process
mu_p_museum_0 = -1
p_museum_species = 0
sigma_p_museum_species = 0.5
p_museum_site = 0
sigma_p_museum_site = 0.3
p_museum_interval = 0
p_museum_pop_density = 0 

## --------------------------------------------------
### Simulate data
set.seed(1)
my_simulated_data <- simulate_data(n_species,
                                   n_sites,
                                   n_intervals,
                                   n_visits,
                                   
                                   # ecological process
                                   mu_psi_0,
                                   mu_psi_interval,
                                   sigma_psi_interval,
                                   psi_pop_dens,
                                   psi_site_area,
                                  
                                   # citizen science observation process
                                   mu_p_citsci_0,
                                   p_citsci_species,
                                   sigma_p_citsci_species,
                                   p_citsci_site,
                                   sigma_p_citsci_site,
                                   p_citsci_interval,
                                   p_citsci_pop_density, 
                                   
                                   # museum record observation process
                                   mu_p_museum_0,
                                   p_museum_species,
                                   sigma_p_museum_species,
                                   p_museum_site,
                                   sigma_p_museum_site,
                                   p_museum_interval,
                                   p_museum_pop_density, 
                                   
                                   # introduce NAs (missed visits)?
                                   introduce_NAs = TRUE,
                                   sites_missing = 0.25*n_sites, 
                                   intervals_missing = 2,
                                   visits_missing = 2)

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V_citsci <- my_simulated_data$V_citsci # detection data
V_museum <- my_simulated_data$V_museum # detection data
V_museum_NA <- my_simulated_data$V_museum_NA # indicator of whether sampling occurred
n_species <- my_simulated_data$n_species # number of species
n_sites <- my_simulated_data$n_sites # number of sites
n_intervals <- my_simulated_data$n_intervals # number of surveys 
n_visits <- my_simulated_data$n_visits

#View(as.data.frame(V_citsci[1:10,1:10,,]))
#View(as.data.frame(V_museum[1:10,1:10,,]))
sum(my_simulated_data$V_citsci == 1)
sum(my_simulated_data$V_museum == 1)

intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

pop_densities <- my_simulated_data$pop_density
site_areas <- my_simulated_data$site_area

stan_data <- c("V_citsci", "V_museum", "V_museum_NA",
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", "species", "sites",
               "pop_densities", "site_areas") 

# Parameters monitored
params <- c("mu_psi_0",
            "mu_psi_interval",
            "sigma_psi_interval",
            "psi_pop_density",
            "psi_site_area",
            
            "mu_p_citsci_0",
            "sigma_p_citsci_species",
            "sigma_p_citsci_site",
            "p_citsci_interval",
            "p_citsci_pop_density", 
            
            "mu_p_museum_0",
            "sigma_p_museum_species",
            "sigma_p_museum_site",
            "p_museum_interval",
            "p_museum_pop_density"
)

parameter_value <- c(mu_psi_0,
                     mu_psi_interval,
                     sigma_psi_interval,
                     psi_pop_dens,
                     psi_site_area,
                     
                     mu_p_citsci_0,
                     sigma_p_citsci_species,
                     sigma_p_citsci_site,
                     p_citsci_interval,
                     p_citsci_pop_density,
                     
                     mu_p_museum_0,
                     sigma_p_museum_species,
                     sigma_p_museum_site,
                     p_museum_interval,
                     p_museum_pop_density
)

# MCMC settings
n_iterations <- 800
n_thin <- 1
n_burnin <- 400
n_chains <- 3
n_cores <- n_chains

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(mu_psi_0 = runif(1, -1, 1),
       mu_psi_interval = runif(1, -1, 1),
       sigma_psi_interval = runif(1, 0, 1),
       psi_pop_density = runif(1, -1, 1),
       psi_site_area = runif(1, -1, 1),
       
       mu_p_citsci_0 = runif(1, -1, 1),
       sigma_p_citsci_species = runif(1, 0, 1),
       sigma_p_citsci_site = runif(1, 0, 1),
       p_citsci_interval = runif(1, -1, 1),
       p_citsci_pop_density = runif(1, -1, 1),
       
       mu_p_museum_0 = runif(1, -1, 1),
       sigma_p_museum_species = runif(1, 0, 1),
       sigma_p_museum_site = runif(1, 0, 1),
       p_museum_interval = runif(1, -1, 1),
       p_museum_pop_density = runif(1, -1, 1)
       
  )
)

targets <- as.data.frame(cbind(params, parameter_value))

## --------------------------------------------------
### Run model
library(rstan)
stan_model <- "./models/model_integrated.stan"

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

saveRDS(stan_out_sim, "./simulation/stan_out_sim_integrated.rds")
stan_out_sim <- readRDS("./simulation/stan_out_sim_integrated.rds")

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out_sim, pars = c(
  "mu_psi_0",
  "mu_p_citsci_0",
  "mu_p_museum_0"
))

# pairs plot
pairs(stan_out, pars = c(
  "mu_psi_0",
  "mu_p_citsci_0",
  "mu_p_museum_0"
))

# should now also write a posterior predictive check into the model


