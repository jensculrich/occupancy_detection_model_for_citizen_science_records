### Data simulation for citizen science/museum records of pollinator occurrence
# jcu; started nov 21, 2022

# If run on AWS server using pre-built AMI from: 
# YOU WILL NEED TO RUN THIS REINSTALL FIRST
# With the editions of rstan and stanheaders on the AMI, 
# the OS will be unable to compile the stan file.
#remove.packages(c("StanHeaders", "rstan"))
#install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# install.packages("rstan")
# install.packages("MASS") # for simulation of correlated observed data
library(rstan)
library(MASS)

## Data simulation for integrated model form
# single season, species considered across all sites (not just range)
# and across all visits (not considering that sites without any detections 
# for any species were probably not visited or that sites with only detections for
# a few species only might be targeted and not fully community wide)

# model fitting should return the parameter inputs ("targets")

simulate_data <- function(
    
  ## Study design
  n_species,
  n_sites,
  n_intervals,
  n_visits,
  
  ## Ecological process
  # dispersion and occupancy
  gamma_0,
  gamma_1,
  phi,
  
  # abundance
  mu_eta_0,
  sigma_eta_species,
  sigma_eta_site,
  mu_eta_impervious_surface,
  sigma_eta_impervious_surface,
  mu_eta_perennial_plants,
  sigma_eta_perennial_plants,
  eta_site_area,
  
  ## Detection process
  # citizen science observation process
  mu_p_citsci_0,
  sigma_p_citsci_species,
  sigma_p_citsci_site,
  mu_p_citsci_interval,
  sigma_p_citsci_interval,
  p_citsci_pop_density,
  
  # museum record observation process
  mu_p_museum_0,
  sigma_p_museum_species,
  sigma_p_museum_site,
  mu_p_museum_interval,
  sigma_p_museum_interval,
  p_museum_pop_density,
  
  ## Range dynamics
  sites_in_range_beta1,
  sites_in_range_beta2,
  
  ## Community surveys
  sites_missing,
  intervals_missing,
  visits_missing
  
){
  
  ## ilogit and logit functions
  inv_logit <- function(x) exp(x)/(1+exp(x))
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
  corMat <- matrix(c(1, -.3, -.2,
                     -.3, 1, .4,
                     -.2, .4, 1),
                   ncol = 3)
  covMat <- stddev %*% t(stddev) * corMat
  correlated_data <- MASS::mvrnorm(n = n_sites, mu = mu, Sigma = covMat, empirical = FALSE)
  
  perennial_plant_cover <- correlated_data[,1]
  impervious_surfaces <- correlated_data[,2]
  pop_density <- correlated_data[,3]
  
  ## --------------------------------------------------
  ### specify species-specific abundance rates
  
  ## species-specific random intercepts
  eta_species <- rnorm(n=n_species, mean=0, sd=sigma_eta_species)
  # species baseline abundance is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_eta_sp
  
  ## species-specific random intercepts
  eta_site <- rnorm(n=n_sites, mean=0, sd=sigma_eta_site)
  # species baseline abundance is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_eta_sp
  
  ## species-specific random intercepts
  eta_impervious_surface <- rnorm(n=n_species, mean=mu_eta_impervious_surface, sd=sigma_eta_impervious_surface)
  # species baseline abundance is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_eta
  
  ## species-specific random intercepts
  eta_perennial_plant_cover <- rnorm(n=n_species, mean=mu_eta_perennial_plants, sd=sigma_eta_perennial_plants)
  # species baseline abundance is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_eta
  
  ## --------------------------------------------------
  ### specify species-specific detection probabilities
  
  ## species-specific random intercepts
  p_citsci_species  <- rnorm(n=n_species, mean = 0, sd=sigma_p_citsci_species)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p_species
  
  ## site-specific random intercepts
  p_citsci_site  <- rnorm(n=n_sites, mean = 0, sd=sigma_p_citsci_site)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p
  
  ## site-specific random slopes
  p_citsci_interval  <- rnorm(n=n_sites, mean = mu_p_citsci_interval, sd=sigma_p_citsci_interval)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p
  
  ## species-specific random intercepts
  p_museum_species  <- rnorm(n=n_species, mean = 0, sd=sigma_p_museum_species)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p_species
  
  ## site-specific random intercepts
  p_museum_site  <- rnorm(n=n_sites, mean = 0, sd=sigma_p_museum_site)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p
  
  ## site-specific random slopes
  p_museum_interval  <- rnorm(n=n_sites, mean = mu_p_museum_interval, sd=sigma_p_museum_interval)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p
  
  ## --------------------------------------------------
  ## Create arrays for abundance rate (eta) and detection rate (p)
  log_eta_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals)) 
  # a psi value for each species, at each site, in each interval 
  
  p_matrix_citsci <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  p_matrix_museum <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  
  for(species in 1:n_species) { # for each site
    for(site in 1:n_sites) { # for each interval
      for(interval in 1:n_intervals) { # for each species
        
        log_eta_matrix[species, site, interval] <- # log abundance rate is equal to
          mu_eta_0 + # a baseline intercept
          eta_species[species] + # a species-specific intercept  
          eta_site[site] + # a site specific intercept
          eta_impervious_surface[species]*impervious_surfaces[site] + # a species specific temporal change
          #eta_perennial_plant_cover[species]*perennial_plant_cover[site] + # a fixed effect of population density 
          eta_site_area*site_area[site] # a fixed effect of site area
        
        for(visit in 1:n_visits) { # for each visit
          
          p_matrix_citsci[species, site, interval, visit] <-  # detection is equal to 
            mu_p_citsci_0 + # a baseline intercept
              p_citsci_species[species] + # a species-specific intercept
              p_citsci_site[site] #+ # a spatiotemporally specific intercept
              #p_citsci_interval[site]*intervals[interval] + # an overall effect of time on detection
              #p_citsci_pop_density*pop_density[site] # an effect of population density on detection ability
          
          
          p_matrix_museum[species, site, interval, visit] <- # detection is equal to 
            mu_p_museum_0 + # a baseline intercept
            p_museum_species[species] + # a species-specific intercept            
            p_museum_site[site] #+ # a spatiotemporally specific intercept
            #p_museum_interval[site]*intervals[interval] + # an overall effect of time on detection
            #p_museum_pop_density*pop_density[site] # an effect of population density on detection ability
          
          
        } # for each visit
      } # for each species
    } # for each interval
  } # for each site
  
  # preview the psi and p arrays
  # head(eta_matrix[1:n_species, 1:n_sites,1])
  head(p_matrix_citsci[1:n_species, 1:n_sites,1,1])
  
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
  
  #ranges[1:n_species, 1:n_sites,1,1]
  #ranges[1:n_species, 1:n_sites,1,2]
  
  # Creating the Sequence
  # beta_distr = seq(0,1, by=0.1)
  
  # Case 3
  #plot(beta_distr, dbeta(gfg, 5,2), xlab = "X",
  #     ylab = "Beta Density", type = "l",
  #     col = "Red")
  
  ## --------------------------------------------------
  # Generate species*site*interval suitability
  
  omega <- array(dim = c(n_species, n_sites, n_intervals), NA)
  
  for(species in 1:n_species){
    for(site in 1:n_sites){
      for(interval in 1:n_intervals){
        
        # if the site is in range then it has some chance to be occupied
        if(ranges[species,site,1,1] == 1){
          
          # rep across visits so the site is open to a species or closed to a species
          # across all visits 1:n_visits
          omega[species,site,interval] <- 
            rbinom(1, 1, inv_logit(gamma_0 + gamma_1 * log_eta_matrix[species, site, interval]))
                 
        } else{ # the site is not in range so must not be occupied
          # these one count for or against omega estimation since the range indicator will
          # tell STAN to skip over these species*site combinations
          
          omega[species,site,interval] <- 0
        
        } # end if/else
      }
    }
  }
  
  # 1 are occupied sites, 0 are unoccupied sites
  #log_eta_matrix[1:n_species, 1:n_sites,1]
  #omega[1:n_species, 1:n_sites,2]
  
  ## --------------------------------------------------
  ## Generate abundance data given means eta and dispersion phi
  
  N_matrix <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals))
  
  for(species in 1:n_species){
    for(site in 1:n_sites){
      for(interval in 1:n_intervals){
        
        # if site is in the species's range, 
        # and is also suitable
        # then sample abundance with some prob eta from a trunctated (>0) count distribution
        # else abundance state is 0
        if(omega[species,site,interval] == 1){
          
          # eta for species, site, interval
          T <- exp(log_eta_matrix[species,site,interval])
          N_matrix[species,site,interval] <- rnbinom(n=1, mu=T, size=phi) 
          Y0 <- N_matrix[species,site,interval][N_matrix[species,site,interval]>0] 
          r <- (1 - length(Y0))
          # to create a truncated count distr we reject any 0 values and continue sampling
          while(r>0){
            N_matrix[species,site,interval] <- rnbinom(n=1, mu=T, size=phi) 
            Y0 <- c(Y0,
                    N_matrix[species,site,interval][N_matrix[species,site,interval]>0])
            r <- (1 - length(Y0))
          }
          

        } else{
          
          N_matrix[species,site,interval] <- 0
          
       }
        
      }
    }
  }
  
  # preview N_matrix
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
            prob = inv_logit(p_matrix_citsci[species,site,interval,visit]))
          
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
          
          # if site is in the species's range
          # then determine abundance (including 0 abundance) with some prob eta
          # else abundance state is 0
          if(ranges[species,site,1,1] == 1) {
          
          V_museum[species,site,interval,visit] <- 
            # Z[species,site,interval] * # occupancy state * detection prob
            # N_matrix[species,site,interval]
            rbinom(n = 1, 
                   size = omega[species,site,interval], 
                   prob = inv_logit(p_matrix_museum[species,site,interval,visit]))
          
          } else{
            
            # if not in the range we won't detect it.
            # This will be handled by treating as an NA using the V_museum_NA indicator generated below
            V_museum[species,site,interval, visit] <- 0
            
          }
          
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
  # museum NAs
  
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
  check <- which(V_museum>V_museum_NA)  
  
  # sum(V_citsci == 1)
  # sum(V_museum == 1)
  # sum(V_citsci == 0)
  # sum(V_museum == 0)
  
  # sum(V_museum_NA == 1) # where 1 = sampled plus in range  
  # sum(V_museum_NA == 0) # and 0 = unsampled and/or not in range
  
  ## --------------------------------------------------
  # Max search range for each species*site*interval
  K <- array(dim = c(n_species, n_sites, n_intervals))
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_intervals){
        
        K[i,j,k] <- 3*(max(V_citsci[i,j,k,])+5)
        
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
    ranges = ranges, # array indicating whether sampling occurred in a site*interval*visit
    V_museum_NA = V_museum_NA, # array indicating whether sampling occurred in a site*interval*visit
    n_species = n_species, # number of species
    n_sites = n_sites, # number of sites
    n_intervals = n_intervals, # number of surveys 
    n_visits = n_visits, # number of visits
    pop_density = pop_density, # vector of pop densities
    impervious_surfaces = impervious_surfaces, # vector of impervious surface covers
    perennial_plant_cover = perennial_plant_cover, # vector of perennial plant covers
    site_area = site_area, # vector of site areas
    K = K # upper limit search area
  ))
  
} # end simulate_data function


## --------------------------------------------------
### Variable values for data simulation
## study dimensions
n_species = 18 ## number of species
n_sites = 18 ## number of sites
n_intervals = 3 ## number of occupancy intervals
n_visits = 5 ## number of samples per year

## ecological process
#omega = 0.8
gamma_0 = 0.25
gamma_1 = 0.75
phi = 2

# abundance
mu_eta_0 = 1.5
sigma_eta_species = 0.75
sigma_eta_site = 0.75
mu_eta_impervious_surface = -1
sigma_eta_impervious_surface = 0.5
mu_eta_perennial_plants = 1
sigma_eta_perennial_plants = 0.5
eta_site_area = 1

## detection
# citizen science observation process
mu_p_citsci_0 = 0
sigma_p_citsci_species = 0.5
sigma_p_citsci_site = 0.5
mu_p_citsci_interval = 2
sigma_p_citsci_interval = 0.5 
p_citsci_pop_density = 1

# museum record observation process
mu_p_museum_0 = -1
sigma_p_museum_species = 0.5
sigma_p_museum_site = 0.5
mu_p_museum_interval = 0
sigma_p_museum_interval = 0.5 
p_museum_pop_density = 0.75

# introduce NAs (visits that did not survey entire community)?
sites_missing = 0.5*n_sites 
intervals_missing = 2
visits_missing = 2

sites_in_range_beta1 = 5
sites_in_range_beta2 = 2

## --------------------------------------------------
### Simulate data
set.seed(1)
my_simulated_data <- simulate_data(
  
        ## Study design
        n_species,
        n_sites,
        n_intervals,
        n_visits,
        
        ## Ecological process
        # dispersion and occupancy
        gamma_0,
        gamma_1,
        phi,
        
        # abundance
        mu_eta_0,
        sigma_eta_species,
        sigma_eta_site,
        mu_eta_impervious_surface,
        sigma_eta_impervious_surface,
        mu_eta_perennial_plants,
        sigma_eta_perennial_plants,
        eta_site_area,
        
        ## Detection process
        # citizen science observation process
        mu_p_citsci_0,
        sigma_p_citsci_species,
        sigma_p_citsci_site,
        mu_p_citsci_interval,
        sigma_p_citsci_interval,
        p_citsci_pop_density,
        
        # museum record observation process
        mu_p_museum_0,
        sigma_p_museum_species,
        sigma_p_museum_site,
        mu_p_museum_interval,
        sigma_p_museum_interval,
        p_museum_pop_density,
        
        ## Range dynamics
        sites_in_range_beta1,
        sites_in_range_beta2,
        
        ## Community surveys
        sites_missing,
        intervals_missing,
        visits_missing
)

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
K <- my_simulated_data$K


sum(my_simulated_data$V_citsci >= 1)
sum(my_simulated_data$V_citsci)
sum(my_simulated_data$V_museum == 1)
 
intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

pop_densities <- my_simulated_data$pop_density
site_areas <- my_simulated_data$site_area
impervious_surfaces <- my_simulated_data$impervious_surfaces
perennial_plant_cover <- my_simulated_data$perennial_plant_cover

stan_data <- c("V_citsci", 
               "V_museum",
               "V_museum_NA",
               "ranges", 
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "K",
               "intervals", "species", "sites",
               "pop_densities", "site_areas",
               "impervious_surfaces", "perennial_plant_cover") 

# Parameters monitored
params <- c(
            "gamma_0",
            "gamma_1",
            "phi",
            
            "mu_eta_0",
            "sigma_eta_species",
            "sigma_eta_site",
            "mu_eta_impervious_surface",
            "sigma_eta_impervious_surface",
            #"mu_eta_perennial_plants",
            #"sigma_eta_perennial_plants",
            "eta_site_area",
            
            "mu_p_citsci_0",
            "sigma_p_citsci_species",
            "sigma_p_citsci_site",
            #"mu_p_citsci_interval",
            #"sigma_p_citsci_interval",
            #"p_citsci_pop_density",
            
            "mu_p_museum_0",
            "sigma_p_museum_species",
            "sigma_p_museum_site"#,
            #"mu_p_museum_interval",
            #"sigma_p_museum_interval",
            #"p_museum_pop_density",
            
            # posterior predictive check
            #"fit",
            #"fit_new",
            #"fit_occupancy",
            #"fit_occupancy_new"
)

parameter_value <- c(
                     gamma_0,
                     gamma_1,
                     phi,
                     
                     mu_eta_0,
                     sigma_eta_species,
                     sigma_eta_site,
                     mu_eta_impervious_surface,
                     sigma_eta_impervious_surface,
                     #mu_eta_perennial_plants,
                     #sigma_eta_perennial_plants,
                     eta_site_area,
                     
                     mu_p_citsci_0,
                     sigma_p_citsci_species,
                     sigma_p_citsci_site,
                     #mu_p_citsci_interval,
                     #sigma_p_citsci_interval,
                     #p_citsci_pop_density,
                     
                     mu_p_museum_0,
                     sigma_p_museum_species,
                     sigma_p_museum_site#,
                     #mu_p_museum_interval,
                     #sigma_p_museum_interval,
                     #p_museum_pop_density,
                     
                     #NA, # posterior predictive check
                     #NA, # posterior predictive check
                     #NA, # posterior predictive check
                     #NA # posterior predictive check
)

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 3
n_cores <- 4

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
# inits here are uniform from -1 to 1 for all slope effects and 
# 0 to 1 (positive values only) for all variance terms
inits <- lapply(1:n_chains, function(i)
  
  list(#omega = runif(1, 0, 1),
       gamma_0 = runif(1, -0.25, 0.25),
       gamma_1 = runif(1, 0, 1),
       phi = runif(1, 0, 1),
       
       mu_eta_0 = runif(1, -1, 1),
       sigma_eta_species = runif(1, 0, 1),
       sigma_eta_site = runif(1, 0, 1),
       mu_eta_impervious_surface = runif(1, -1, 1),
       sigma_eta_impervious_surface = runif(1, 0, 1),
       #mu_eta_perennial_plants = runif(1, -1, 1),
       #sigma_eta_perennial_plants = runif(1, 0, 1),
       eta_site_area = runif(1, -1, 1),
       
       mu_p_citsci_0 = runif(1, -1, 1),
       sigma_p_citsci_species = runif(1, 0, 1),
       sigma_p_citsci_site = runif(1, 0, 1),
       #mu_p_citsci_interval = runif(1, -1, 1),
       #sigma_p_citsci_interval = runif(1, 0, 1),
       #p_citsci_pop_density = runif(1, -1, 1),
       
       mu_p_museum_0 = runif(1, -1, 1),
       sigma_p_museum_species = runif(1, 0, 1),
       sigma_p_museum_site = runif(1, 0, 1)#,
       #mu_p_museum_interval = runif(1, -1, 1),
       #sigma_p_museum_interval = runif(1, 0, 1),
       #p_museum_pop_density = runif(1, -1, 1)
       
  )
)

targets <- as.data.frame(cbind(params, parameter_value))

## --------------------------------------------------
### Run model

stan_model <- "./abundance-occupancy/models/model_abundance.stan"

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

saveRDS(stan_out_sim, "./abundance-occupancy/model_outputs/stan_out_sim_abundance.rds")
stan_out_sim <- readRDS("./abundance-occupancy/model_outputs/stan_out_sim_abundance.rds")

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out_sim, pars = c(
  "mu_eta_0",
  "eta_site_area",
  "mu_p_citsci_0",
  "mu_p_museum_0",
  "phi",
  "gamma_0",
  "gamma_1"
))

# pairs plot
pairs(stan_out_sim, pars = c(
  "mu_eta_0",
  "eta_site_area",
  "mu_p_citsci_0",
  "mu_p_museum_0",
  "phi",
  "gamma_0",
  "gamma_1"
))

# Posterior predictive check
list_of_draws <- as.data.frame(stan_out_sim)

## --------------------------------------------------
# Abundance

# Evaluation of fit
plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 250000),
     xlim = c(0, 250000))
abline(0, 1, lwd = 2, col = "black")

# Should be close to 1. 
# If the mean of actual data is greater (value > 1)
# then the model underpredicts the real variation in counts.
# If the mean of actual data is less (value < 1)
# then the model overpredicts the variation in counts.
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)

# Should be close to 50% or 0.5
# similarly the actual data should be further away from  
# the expected value of the count about half of the time,
# (versus a count generated using the abundance rate and detection rate)
mean(list_of_draws$fit_new > list_of_draws$fit)

## --------------------------------------------------
# Occupancy

# Evaluation of fit
plot(list_of_draws$fit_occupancy, list_of_draws$fit_occupancy_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 500),
     xlim = c(0, 500))
abline(0, 1, lwd = 2, col = "black")

# Should be close to 1. 
# If the mean of actual data is greater (value > 1)
# then the model underpredicts the real variation in counts.
# If the mean of actual data is less (value < 1)
# then the model overpredicts the variation in counts.
mean(list_of_draws$fit_occupancy) / mean(list_of_draws$fit_occupancy_new)

# Should be close to 50% or 0.5
# similarly the actual data should be further away from  
# the expected value of the count about half of the time,
# (versus a count generated using the abundance rate and detection rate)
mean(list_of_draws$fit_occupancy_new > list_of_draws$fit_occupancy)

