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
                          omega,
                          phi,
                          
                          mu_lambda_0,
                          
                          ## observation process
                          # citizen science observation process
                          mu_p_citsci_0,
                          
                          # museum record observation process
                          mu_p_museum_0,
                          
                          sites_in_range_beta1,
                          sites_in_range_beta2
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
  #pop_density <- rnorm(n_sites, mean = 0, sd = 1)
  
  # create a vector of site area of length = number of sites
  # the model takes z-score scaled data (with mean of 0) so it's ok to center at 0 here
  # this simulates the realism that some sites are e.g. partially on ocean or  
  # partially outside the administrative area from which we are drawing collection data.
  #site_area <- rnorm(n_sites, mean = 0, sd = 1)
  
  ## --------------------------------------------------
  ### specify species-specific abundance and detection rates
  #mu_uv <- c(0, 0) # centered means since we already have grand intercepts
  
  #Sigma <- matrix(NA, nrow = 2, ncol = 2)
  #Sigma[1,1] = (sigma_u)^2
  #Sigma[2,2] = (sigma_v)^2
  #Sigma[1,2] = sigma_u * sigma_v * rho_uv
  #Sigma[2,1] = sigma_u * sigma_v * rho_uv
  
  #uv <- MASS::mvrnorm(n = n_species, mu = mu_uv, Sigma = Sigma, empirical = FALSE)
  #uv
  #(rho_uv_simmed <- cor(uv[,1], uv[,2])) # correlation of site occupancy and detection
  
  ## --------------------------------------------------
  ## Create arrays for psi and p
  lambda_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals)) 
  # a psi value for each species, at each site, in each interval 
  
  p_matrix_citsci <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  p_matrix_museum <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  
  for(species in 1:n_species) { # for each site
    for(site in 1:n_sites) { # for each interval
      for(interval in 1:n_intervals) { # for each species
        
        lambda_matrix[species, site, interval] <- exp( # occupancy is equal to
          mu_lambda_0 # a baseline intercept
          # uv[species,1] # a species intercept (correlated with museum detection rates)
            #lambda_site[site] + # a site specific intercept
            #lambda_interval[species]*intervals[interval] + # a species specific temporal change
            #lambda_pop_dens[species]*pop_density[site] + # a fixed effect of population density 
            #lambda_site_area*site_area[site] # a fixed effect of site area
        )
        
        for(visit in 1:n_visits) { # for each visit
          
          p_matrix_citsci[species, site, interval, visit] <-  # detection is equal to 
            mu_p_citsci_0 #+ # a baseline intercept
              #p_citsci_species[species] + # a species specific intercept
              #p_citsci_site[site] + # a spatiotemporally specific intercept
              #p_citsci_interval*intervals[interval] + # an overall effect of time on detection
              #p_citsci_pop_density*pop_density[site] # an effect of population density on detection ability
          
          
          p_matrix_museum[species, site, interval, visit] <- # detection is equal to 
            mu_p_museum_0 #+ # a baseline intercept
            # uv[species,2] # a species intercept (correlated with species abundance)            
            #p_museum_site[site] + # a spatiotemporally specific intercept
            #p_museum_interval*intervals[interval] + # an overall effect of time on detection
            #p_museum_pop_density*pop_density[site] # an effect of population density on detection ability
          
          
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
  
  prob_site_in_range <-  rbeta(n_species, sites_in_range_beta1, sites_in_range_beta2)
  
  for(i in 1:n_species){
    
    ranges[i,,1:n_intervals,1:n_visits] <- rbinom(1, n=n_sites, prob=prob_site_in_range)
    
  }
  
  ranges[1:n_species, 1:n_sites,1,1]
  ranges[1:n_species, 1:n_sites,1,2]
  
  # Creating the Sequence
  # beta_distr = seq(0,1, by=0.1)
  
  # Case 3
  #plot(beta_distr, dbeta(gfg, 5,2), xlab = "X",
  #     ylab = "Beta Density", type = "l",
  #     col = "Red")
  
  ## --------------------------------------------------
  # Generate species*site*interval suitability
  
  suitability <- array(dim = c(n_species, n_sites, n_intervals, n_visits), NA)
  
  for(species in 1:n_species){
    for(site in 1:n_sites){
      for(interval in 1:n_intervals){
        
        # rep across visits so the site is open to a species or closed to a species
        # across all visits 1:n_visits
        suitability[species,site,interval,1:n_sites] <- rep(rbinom(1, 1, prob=omega), n_sites)
        
      }
    }
  }
  
  suitability[1:n_species, 1:n_sites,1,1]
  suitability[1:n_species, 1:n_sites,1,2]
  
  ## --------------------------------------------------
  ## Generate abundance data given means lambda and dispersion phi
  
  N_matrix <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        
        # if site is in the species's range, then determine occupancy with some prob psi
        # else occupancy state is 0
        if(ranges[species,site,1,1] == 1 && suitability[species,site,interval,1]) {
          
          N_matrix[species,site,interval] <- rpois(n = 1, 
                                                   mu = lambda_matrix[species,site,interval], 
                                                   size = phi)

        } else{
          
          N_matrix[species,site,interval] <- 0
          
       }
        
      }
    }
  }
  
  # preview N_matrix
  N_matrix[1:n_species,1:n_sites,1]
  
  ## --------------------------------------------------
  ## Generate abundance data given means lambda and dispersion phi
  
  Z_matrix <- N_matrix
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        
        if(Z_matrix[species,site,interval] > 0) {
          
          # if one or more is present, then transform into a binary response of presence
          Z_matrix[species,site,interval] <- 1
       
        }
        
      }
    }
  }
  
  # preview Z_matrix
  Z_matrix[1:n_species,1:n_sites,1]
  N_matrix[1:n_species,1:n_sites,1]
  
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
          
          V_citsci[species,site,interval,visit] <-  rbinom(
            n = 1, 
            size = N_matrix[species,site,interval], 
            prob = ilogit(p_matrix_citsci[species,site,interval,visit]))
          
        }
      }
    }
  } # end simulate detection data
  
  # preview V_citsci
  V_citsci[1:n_species,1:n_sites,1,1]
  
  V_museum <- array(NA, dim=c(n_species=n_species,
                              n_sites=n_sites,
                              n_intervals=n_intervals,
                              n_visits=n_visits))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        for(visit in 1:n_visits){
          
          V_museum[species,site,interval,visit] <- 
            # Z[species,site,interval] * # occupancy state * detection prob
            # N_matrix[species,site,interval]
            rbinom(n = 1, 
                   size = Z_matrix[species,site,interval], 
                   prob = ilogit(p_matrix_museum[species,site,interval,visit]))
          
        }
      }
    }
  } # end simulate detection data
  
  # preview V_museum
  V_museum[1:n_species,1:n_sites,1,1]
  
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
  V_citsci_NA[1:n_species,1:n_sites,1,1]
  # should NEVER have a V_citsci detection outside of the range
  # check <- which(V_citsci>V_citsci_NA)
  
  ## --------------------------------------------------
  # museum NAs
  
  # choose random sites that didn't get visited (for all species) by museum collecting visits
  #site_missed = sample.int(n_sites, sites_missing)
  #interval_missed = sample.int(n_intervals, intervals_missing)
  #visit_missed = sample.int(n_visits, visits_missing)
  
  # we will make an array that holds values of 1 if sampling occurred
  # or 0 if sampling did not occur at the site*interval*visit.
  #V_museum_NA <- V_museum 
  #V_museum_NA[1:n_species, site_missed, interval_missed, visit_missed] <- NA
  
  # replace all other values with 1 (was sampled)
  #V_museum_NA <- replace(V_museum_NA, V_museum_NA==0, 1)
  # and now replace all NAs with 0, which will act as an indicator for the likelihood function
  # to skip over this sample by contracting the total possible number of observations that could have occurred
  #V_museum_NA[is.na(V_museum_NA)] <- 0
  
  # now multiply by whether a site was in a range or not,
  # 1 means the species at the site in the time was BOTH..
  # a target of a community sample AND
  # the site is in the species's range
  #V_museum_NA <- V_museum_NA*ranges
  
  # now we want to replace the detection data with 0's where sampling did not occur
  # so that we are saying that a a species occurs at a site*interval..
  # was not observed at the visit to the site in the interval..
  # BUT the model will remove this from contributing to the probability density by removing
  # the max number of sightings that could have occurred for 
  # each visit in the site*interval with a 0 in V_museum_NA
  #V_museum[1:n_species, site_missed, interval_missed, visit_missed] <- 0
  
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
  # Max search range for each species*site*interval
  K <- array(dim = c(n_species, n_sites, n_intervals))
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_intervals){
        
        K[i,j,k] <- 5*(max(V_citsci[i,j,k,])+5)
        
      }
    }
  }
  
  # preview V_museum
  K[1:n_species,1:n_sites,1]
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    V_citsci = V_citsci, # detection data from citizen science records
    V_museum = V_museum, # detection data from museum records
    V_citsci_NA = V_citsci_NA, # array indicating whether sampling occurred in a site*interval*visit
    #V_museum_NA = V_museum_NA, # array indicating whether sampling occurred in a site*interval*visit
    n_species = n_species, # number of species
    n_sites = n_sites, # number of sites
    n_intervals = n_intervals, # number of surveys 
    n_visits = n_visits, # number of visits
    #pop_density = pop_density, # vector of pop densities
    #site_area = site_area # vector of site areas
    K = K # upper limit search area
  ))
  
} # end simulate_data function


## --------------------------------------------------
### Variable values for data simulation
## study dimensions
n_species = 12 ## number of species
n_sites = 5 ## number of sites
n_intervals = 3 ## number of occupancy intervals
n_visits = 6 ## number of samples per year

## ecological process
omega = 0.8
phi = 1.5

mu_lambda_0 = 3



## detection
# citizen science observation process
mu_p_citsci_0 = 0
#p_citsci_species = 0
#sigma_p_citsci_species = 0.5
#p_citsci_site = 0
#sigma_p_citsci_site = 0.3
#p_citsci_interval = 1
#p_citsci_pop_density = 1 

# museum record observation process
mu_p_museum_0 = -1

# introduce NAs (missed visits)?
#sites_missing = 0.5*n_sites 
#intervals_missing = 2
#visits_missing = 4

sites_in_range_beta1 = 5
sites_in_range_beta2 = 2

## --------------------------------------------------
### Simulate data
set.seed(1)
my_simulated_data <- simulate_data(n_species,
                                   n_sites,
                                   n_intervals,
                                   n_visits,
                                   
                                   # ecological process
                                   omega, 
                                   phi,
                                   
                                   mu_lambda_0,
                                  
                                   # citizen science observation process
                                   mu_p_citsci_0,
                                   
                                   # museum record observation process
                                   mu_p_museum_0,
                                   
                                   
                                   # range dynamics
                                   sites_in_range_beta1,
                                   sites_in_range_beta2)

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V_citsci <- my_simulated_data$V_citsci # detection data
V_museum <- my_simulated_data$V_museum # detection data
V_citsci_NA <- my_simulated_data$V_citsci_NA # indicator of whether sampling occurred
#V_museum_NA <- my_simulated_data$V_museum_NA # indicator of whether sampling occurred
n_species <- my_simulated_data$n_species # number of species
n_sites <- my_simulated_data$n_sites # number of sites
n_intervals <- my_simulated_data$n_intervals # number of surveys 
n_visits <- my_simulated_data$n_visits
K <- my_simulated_data$K
rho_uv_simmed <- my_simulated_data$rho_uv_simmed
sigma_uv <- cbind(sigma_u, sigma_v)

sum(my_simulated_data$V_citsci >= 1)
sum(my_simulated_data$V_citsci)
#sum(my_simulated_data$V_museum == 1)
 
intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

#pop_densities <- my_simulated_data$pop_density
#site_areas <- my_simulated_data$site_area

stan_data <- c("V_citsci", 
               "V_museum",
               "V_citsci_NA", 
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "K",
               "intervals", "species", "sites"
               ) 

# Parameters monitored
params <- c("omega",
            "phi",
            
            "mu_lambda_0",
            
            "mu_p_citsci_0",
            
            "mu_p_museum_0"
            
)

parameter_value <- c(omega,
                     phi,
                     
                     mu_lambda_0,
                     
                     mu_p_citsci_0,
                     
                     mu_p_museum_0
                     
)

# MCMC settings
n_iterations <- 600
n_thin <- 2
n_burnin <- 300
n_chains <- 3
n_cores <- 4

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(phi = runif(1, 0, 1),
       
       mu_lambda_0 = runif(1, -1, 1),
       
       mu_p_citsci_0 = runif(1, -1, 1),
       
       mu_p_museum_0 = runif(1, -1, 1),
       
       omega = runif(1, 0, 1)
       
  )
)

targets <- as.data.frame(cbind(params, parameter_value))

## --------------------------------------------------
### Run model
library(rstan)
stan_model <- "./models/model_abundance_copy_2.stan"

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

saveRDS(stan_out_sim, "./model_outputs/stan_out_sim_abundance.rds")
stan_out_sim <- readRDS("./simulation/stan_out_sim_abundance.rds")

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out_sim, pars = c(
  "mu_lambda_0",
  "mu_p_citsci_0",
  "mu_p_museum_0",
  "phi",
  "omega"
))

# pairs plot
pairs(stan_out, pars = c(
  "mu_lambda_0",
  "mu_p_citsci_0",
  "phi"
))

# should now also write a posterior predictive check into the model


