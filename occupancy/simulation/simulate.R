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
                          sigma_psi_species,
                          sigma_psi_site,
                          mu_psi_open_developed,
                          sigma_psi_open_developed,
                          mu_psi_herb_shrub,
                          sigma_psi_herb_shrub,
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
                          
                          sites_missing,
                          intervals_missing,
                          visits_missing,
                          
                          mu_sites_in_range,
                          sigma_sites_in_range
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
  
  # create a vector of site area of length = number of sites
  # the model takes z-score scaled data (with mean of 0) so it's ok to center at 0 here
  # this simulates the realism that some sites are e.g. partially on ocean or  
  # partially outside the administrative area from which we are drawing collection data.
  site_area <- rnorm(n_sites, mean = 0, sd = 1)
  
  # Realistically the plant cover is negatively correlated with both imp surface and pop density
  # while pop density and imp surface are positively correlated.
  # I will simluate data with these realistic correlations
  #  A correlates with B corMat[1,2], and A with C with corMat[1,3], and B with C with corMat[2,3]
  mu <- c(0, 0, 0)
  stddev <- c(1, 1, 1)
  corMat <- matrix(c(1, 0, -0.25,
                     0, 1, -0.1,
                     -0.25, -0.1, 1),
                   ncol = 3)
  covMat <- stddev %*% t(stddev) * corMat
  correlated_data <- MASS::mvrnorm(n = n_sites, mu = mu, Sigma = covMat, empirical = FALSE)
  
  pop_density <- correlated_data[,1]
  open_developed <- correlated_data[,2]
  herb_shrub <- correlated_data[,3]
  
  ## --------------------------------------------------
  ### specify species-specific occupancy probabilities
  
  ## species-specific random intercepts
  psi_species <- rnorm(n=n_species, mean=0, sd=sigma_psi_species)
  # species baseline occupancy is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma.psi.sp
  
  ## site-specific random intercepts
  psi_site <- rnorm(n=n_sites, mean=0, sd=sigma_psi_site)
  # species baseline occupancy is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma.psi.sp
  
  ## effect of interval on occupancy (species-specific random slopes)
  psi_open_developed <- rnorm(n=n_species, mean=mu_psi_open_developed, sd=sigma_psi_open_developed)
  # change in each species occupancy across time is drawn from a distribution defined
  # by a community mean (mu_psi_open_developed) with 
  # species specific variation defined by sigma_psi_open_developed
  
  ## effect of pop density on occupancy (species-specific random slopes)
  psi_herb_shrub <- rnorm(n=n_species, mean=mu_psi_herb_shrub, sd=sigma_psi_herb_shrub)
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
  logit_psi_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals)) 
  # a psi value for each species, at each site, in each interval 
  
  logit_p_matrix_citsci <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  logit_p_matrix_museum <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  
  for(species in 1:n_species) { # for each site
    for(site in 1:n_sites) { # for each interval
      for(interval in 1:n_intervals) { # for each species
        
        logit_psi_matrix[species, site, interval] <- # occupancy is equal to
          mu_psi_0 + # a baseline intercept
            psi_species[species] + # a species specific intercept
            psi_site[site] + # a site specific intercept
            psi_open_developed[species]*open_developed[interval] + # a species specific temporal change
            psi_herb_shrub[species]*herb_shrub[site] + # a fixed effect of population density 
            psi_site_area*site_area[site] # a fixed effect of site area
        
        for(visit in 1:n_visits) { # for each visit
          
          logit_p_matrix_citsci[species, site, interval, visit] <- # detection is equal to 
            mu_p_citsci_0 + # a baseline intercept
              p_citsci_species[species] + # a species specific intercept
              p_citsci_site[site] + # a spatiotemporally specific intercept
              p_citsci_interval*intervals[interval] + # an overall effect of time on detection
              p_citsci_pop_density*pop_density[site] # an effect of population density on detection ability
          
          logit_p_matrix_museum[species, site, interval, visit] <- # detection is equal to 
            mu_p_museum_0 + # a baseline intercept
            p_museum_species[species] + # a species specific intercept
            p_museum_site[site] + # a spatiotemporally specific intercept
            p_museum_interval*intervals[interval] + # an overall effect of time on detection
            p_museum_pop_density*pop_density[site] # an effect of population density on detection ability
          
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
  
  # range intersection is the only thing considered to be affecting 
  # whether or not a site has potential to be sampled by citizen science efforts
  
  ## --------------------------------------------------
  # Generate species ranges
  ranges <- array(dim = c(n_species, n_sites, n_intervals, n_visits), NA)
  
  sites_in_range_beta1 = 2
  sites_in_range_beta2 = 2
  prob_site_in_range <-  rbeta(n_species, sites_in_range_beta1, sites_in_range_beta2)
  
  for(i in 1:n_species){
    
    ranges[i,,1:n_intervals,1:n_visits] <- rbinom(1, n=n_sites, prob=prob_site_in_range)
    
  }
  
  ## --------------------------------------------------
  ## Generate presence and absence
  
  Z <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        
        # if site is in the species's range, then determine occupancy with some prob psi
        # else occupancy state is 0
        if(ranges[species,site,1,1] == 1) {
          
          Z[species,site,interval] <- rbinom(n = 1, size = 1, 
                                             prob = ilogit(logit_psi_matrix[species,site,interval])
                                             )
        } else{
          
          Z[species,site,interval] <- 0
          
        }
        
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
            rbinom(n = 1, size = 1, 
                   prob = ilogit(logit_p_matrix_citsci[species,site,interval,visit]))
          
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
            rbinom(n = 1, size = 1, 
                   prob = ilogit(logit_p_matrix_museum[species,site,interval,visit]))
          
        }
      }
    }
  } # end simulate detection data
  
  ## --------------------------------------------------
  ## Generate NA indicators (sampling could not have occurred either because
  # a species range does not overlap with the site or no community samples were drawn)
  
  # simulate NA data for the model to work around
  # for our real data we will have some detections that cannot occur because
  # 1) the site is not in the range of the species, or 
  # 2) community sampling didnèt occur at some sites for some visits in some intervals
  
  ## --------------------------------------------------
  # citizen science NAs
  
  # 1 indicates a site was in range and thus the species could be detected there
  V_citsci_NA <- ranges
  # should NEVER have a V_citsci detection outside of the range
  # check <- which(V_citsci>V_citsci_NA)
  
  ## --------------------------------------------------
  # museum NAs
  
  # choose random sites that didn't get visited (for all species) by museum collecting visits
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
  
  # now multiply by whether a site was in a range or not,
  # 1 means the species at the site in the time was BOTH..
  # a target of a community sample AND
  # the site is in the species's range
  V_museum_NA <- V_museum_NA*ranges
  
  # now we want to replace the detection data with 0's where sampling did not occur
  # so that we are saying that a a species occurs at a site*interval..
  # was not observed at the visit to the site in the interval..
  # BUT the model will remove this from contributing to the probability density by removing
  # the max number of sightings that could have occurred for 
  # each visit in the site*interval with a 0 in V_museum_NA
  V_museum[1:n_species, site_missed, interval_missed, visit_missed] <- 0
  
  # should NEVER have a V_museum detection outside of the range and community sampling events
  # check <- which(V_museum>V_museum_NA)  
  
  # sum(V_citsci == 1)
  # sum(V_museum == 1)
  # sum(V_citsci == 0)
  # sum(V_museum == 0)
  
  # sum(V_citsci_NA == 1) # where 1 = sampled plus in range and 0 = unsampled and/or not in range
  # sum(V_museum_NA == 1) # where 1 = sampled plus in range and 0 = unsampled and/or not in range
  # sum(V_citsci_NA == 0)
  # sum(V_museum_NA == 0)
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    V_citsci = V_citsci, # detection data from citizen science records
    V_museum = V_museum, # detection data from museum records
    ranges = ranges, # array indicating whether sampling occurred in a site*interval*visit
    V_museum_NA = V_museum_NA, # array indicating whether sampling occurred in a site*interval*visit
    n_species = n_species, # number of species
    n_sites = n_sites, # number of sites
    n_intervals = n_intervals, # number of surveys 
    n_visits = n_visits, # number of visits
    pop_density = pop_density, # vector of pop densities
    open_developed = open_developed, # vector of impervious surface covers
    herb_shrub = herb_shrub, # vector of perennial plant covers
    site_area = site_area # vector of site areas
  ))
  
} # end simulate_data function


## --------------------------------------------------
### Variable values for data simulation
## study dimensions
n_species = 20 ## number of species
n_sites = 20 ## number of sites
n_intervals = 3 ## number of occupancy intervals
n_visits = 3 ## number of samples per year

## occupancy
mu_psi_0 = -0.5
sigma_psi_species = 0.5
sigma_psi_site = 0.5
mu_psi_open_developed = 0.5
sigma_psi_open_developed = 0.2
mu_psi_herb_shrub = -0.5 # random effect of population density on occupancy
sigma_psi_herb_shrub = 0.2
psi_site_area = 1 # fixed effect of site area on occupancy

## detection
# citizen science observation process
mu_p_citsci_0 = -1
p_citsci_species = 0
sigma_p_citsci_species = 0.5
p_citsci_site = 0
sigma_p_citsci_site = 0.3
p_citsci_interval = 1
p_citsci_pop_density = 1 

# museum record observation process
mu_p_museum_0 = -0.5
p_museum_species = 0
sigma_p_museum_species = 0.5
p_museum_site = 0
sigma_p_museum_site = 0.3
p_museum_interval = 0
p_museum_pop_density = 0 

# introduce NAs (missed visits)?
sites_missing = 0.5*n_sites 
intervals_missing = 2
visits_missing = 2

sites_in_range_beta1 = 2
sites_in_range_beta2 = 2

## --------------------------------------------------
### Simulate data
set.seed(1)
my_simulated_data <- simulate_data(n_species,
                                   n_sites,
                                   n_intervals,
                                   n_visits,
                                   
                                   # ecological process
                                   mu_psi_0,
                                   sigma_psi_species,
                                   sigma_psi_site,
                                   mu_psi_open_developed,
                                   sigma_psi_open_developed,
                                   mu_psi_herb_shrub,
                                   sigma_psi_herb_shrub,
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
                                   sites_missing, 
                                   intervals_missing,
                                   visits_missing,
                                   
                                   sites_in_range_beta1,
                                   sites_in_range_beta2)

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V_citsci <- my_simulated_data$V_citsci # detection data
V_museum <- my_simulated_data$V_museum # detection data
ranges <- my_simulated_data$ranges # indicator of whether sampling occurred
V_museum_NA <- my_simulated_data$V_museum_NA # indicator of whether sampling occurred
n_species <- my_simulated_data$n_species # number of species
n_sites <- my_simulated_data$n_sites # number of sites
n_intervals <- my_simulated_data$n_intervals # number of surveys 
n_visits <- my_simulated_data$n_visits

#View(as.data.frame(V_citsci[1:10,1:10,,]))
#View(as.data.frame(V_museum[1:10,1:10,,]))
sum(my_simulated_data$V_citsci == 1)
sum(my_simulated_data$V_museum == 1)

check_citsci <- which(V_citsci>V_citsci_NA)
check_museum <- which(V_museum>V_museum_NA)
 
intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

pop_densities <- my_simulated_data$pop_density
site_areas <- my_simulated_data$site_area
open_developed <- my_simulated_data$open_developed
herb_shrub <- my_simulated_data$herb_shrub

stan_data <- c("V_citsci", "V_museum", 
               "ranges", "V_museum_NA",
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", "species", "sites",
               "pop_densities", "site_areas", "open_developed", "herb_shrub") 

# Parameters monitored
params <- c("mu_psi_0",
            "sigma_psi_species",
            "sigma_psi_site",
            "mu_psi_open_developed",
            "sigma_psi_open_developed",
            "mu_psi_herb_shrub",
            "sigma_psi_herb_shrub",
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
            "p_museum_pop_density",
            
            "fit_citsci",
            "fit_new_citsci",
            "fit_museum",
            "fit_new_museum"
)

parameter_value <- c(mu_psi_0,
                     sigma_psi_species,
                     sigma_psi_site,
                     mu_psi_open_developed,
                     sigma_psi_open_developed,
                     mu_psi_herb_shrub,
                     sigma_psi_herb_shrub,
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
                     p_museum_pop_density,
                     
                     NA,
                     NA,
                     NA,
                     NA
)

# MCMC settings
n_iterations <- 400
n_thin <- 1
n_burnin <- 200
n_chains <- 4
n_cores <- parallel::detectCores()

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(mu_psi_0 = runif(1, -1, 1),
       sigma_psi_species = runif(1, 0, 1),
       sigma_psi_site = runif(1, 0, 1),
       mu_psi_open_developed = runif(1, -1, 1),
       sigma_psi_open_developed = runif(1, 0, 1),
       mu_psi_herb_shrub = runif(1, -1, 1),
       sigma_psi_herb_shrub = runif(1, 0, 1),
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
stan_model <- "./occupancy/models/model.stan"

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

saveRDS(stan_out_sim, "./occupancy/simulation/stan_out_sim.rds")
stan_out_sim <- readRDS("./simulation/stan_out_sim_integrated_ranges.rds")

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


