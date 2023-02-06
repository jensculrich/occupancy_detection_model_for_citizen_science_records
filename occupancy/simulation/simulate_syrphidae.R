### Data simulation for citizen science/museum records of pollinator occurrence
# jcu; started nov 21, 2022

## Data simulation for integrated model form
# single season, species considered across all sites (not just range)
# and across all visits (not considering that sites without any detections 
# for any species were probably not visited or that sites with only detections for
# a few species only might be targeted and not fully community wide)

# model fitting using 'model_simplest.stan' should return the parameter inputs
simulate_data <- function(n_genera,
                          n_species_per_genera,
                          n_species = n_genera*n_species_per_genera,
                          n_ecoregion_one,
                          n_ecoregion_three_per_one,
                          n_ecoregion_three,
                          n_sites_per_ecoregion_three,
                          n_sites,
                          n_intervals,
                          n_visits,
                          
                          # ecological process
                          mu_psi_0,
                          sigma_psi_species,
                          sigma_psi_genus,
                          sigma_psi_site,
                          sigma_psi_ecoregion_three,
                          sigma_psi_ecoregion_one,
                          mu_psi_open_developed,
                          sigma_psi_open_developed,
                          mu_psi_herb_shrub_forest,
                          sigma_psi_herb_shrub_forest,
                          psi_site_area,
                          
                          # citizen science observation process
                          mu_p_citsci_0,
                          sigma_p_citsci_species,
                          sigma_p_citsci_site,
                          sigma_p_citsci_ecoregion_three,
                          p_citsci_interval,
                          p_citsci_pop_density, 
                          
                          # museum record observation process
                          mu_p_museum_0,
                          sigma_p_museum_species,
                          sigma_p_museum_site,
                          sigma_p_museum_ecoregion_three,
                          p_museum_total_records,
                          
                          # introduce NAs (missed visits)?
                          sites_missing, 
                          intervals_missing,
                          visits_missing,
                          
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
  
  # create a vector of site area of length = number of sites
  # the model takes z-score scaled data (with mean of 0) so it's ok to center at 0 here
  # this simulates the realism that some sites are e.g. partially on ocean or  
  # partially outside the administrative area from which we are drawing collection data.
  site_area <- rnorm(n_sites, mean = 0, sd = 1)
  
  # Realistically the plant cover is negatively correlated with both imp surface and pop density
  # while pop density and imp surface are positively correlated.
  # I will simluate data with these realistic correlations
  #  A correlates with B corMat[1,2], and A with C with corMat[1,3], and B with C with corMat[2,3]
  mu <- c(0, 0, 0, 0)
  stddev <- c(1, 1, 1, 1)
  corMat <- matrix(c(1, 0, -0.25, .75,
                     0, 1, -0.1, 0,
                     -0.25, -0.1, 1, -.5,
                     .75, 0, -.5, 1),
                   ncol = 4)
  covMat <- stddev %*% t(stddev) * corMat
  correlated_data <- MASS::mvrnorm(n = n_sites, mu = mu, Sigma = covMat, empirical = FALSE)
  
  pop_density <- correlated_data[,1]
  open_developed <- correlated_data[,2]
  herb_shrub_forest <- correlated_data[,3]
  developed_med_high <- correlated_data[,4]
  
  
  # create a matrix of total_records of length = number of sites and row = number of intervals
  # the model takes z-score scaled data (with mean of 0) so it's ok to center at 0 here
  
  # correlate records through time by site
  mu2 <- c(0, 0, 0)
  stddev2 <- c(1, 1, 1)
  corMat2 <- matrix(c(1, 0.8, 0.8,
                     0.8, 1, 0.8,
                     0.8, 0.8, 1),
                   ncol = 3)
  covMat2 <- stddev2 %*% t(stddev2) * corMat2

  # this represents the mean number of records (not unique to species), for 
  # any visit in the interval in which 1 or more records were collected
  total_records_museum <- MASS::mvrnorm(n = n_sites, mu = mu2, Sigma = covMat2, empirical = FALSE)
  
  ## --------------------------------------------------
  ### specify species-specific occupancy probabilities
  
  ## ecoregion3-specific random intercepts
  ## site-specific random intercepts
  genus = rep(1:n_genera, each = n_species_per_genera)
  genera_intercepts <- rep(rnorm(n=n_genera, mean=0, sd=sigma_psi_genus),
                                    each=n_species_per_genera)
  
  ## species-specific random intercepts
  species_intercepts <- rnorm(n=n_species, mean=0, sd=sigma_psi_species)
  # species baseline occupancy is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma.psi.sp
  
  genus_lookup <- rep(1:n_genera, each = n_species_per_genera)
  
  psi_species_nested <- vector(length=n_species)
  
  for(i in 1:n_species){
    
    psi_species_nested[i] <- genera_intercepts[i] + species_intercepts[i]
    
  }
  
  ## ecoregion1-specific random intercepts
  ## site-specific random intercepts
  ecoregion_one = rep(1:n_ecoregion_one, each = n_ecoregion_three_per_one*n_sites_per_ecoregion_three)
  
  # site baseline success is drawn from a normal distribution with mean 0 and 
  # site specific variation defined by sigma_alpha_site
  ecoregion_one_intercepts <- rep(rnorm(n=n_ecoregion_one, mean=0, sd=sigma_psi_ecoregion_one),
                         each=n_ecoregion_three_per_one*n_sites_per_ecoregion_three)
  
  ## ecoregion3-specific random intercepts
  ## site-specific random intercepts
  ecoregion_three = rep(1:n_ecoregion_three, each = n_sites_per_ecoregion_three)
  
  # site baseline success is drawn from a normal distribution with mean 0 and 
  # site specific variation defined by sigma_alpha_site
  ecoregion_three_intercepts <- rep(rnorm(n=n_ecoregion_three, mean=0, sd=sigma_psi_ecoregion_three),
                                  each=n_sites_per_ecoregion_three)
  
  
  ecoregion_one_lookup <- rep(1:n_ecoregion_one, each=n_ecoregion_three_per_one)
  
  ## site-specific random intercepts
  
  # site baseline success is drawn from a normal distribution with mean 0 and 
  # site specific variation defined by sigma_alpha_site
  site_intercepts <- rep(rnorm(n=n_sites, mean=0, sd=sigma_psi_site))
  
  ecoregion_three_lookup <- rep(1:n_ecoregion_three, each=n_sites_per_ecoregion_three)
  
  psi_site_nested <- vector(length=n_sites)
  
  for(i in 1:n_sites){
    
    psi_site_nested[i] <- ecoregion_one_intercepts[i] + 
      ecoregion_three_intercepts[i] + site_intercepts[i]
    
  }
  
  #View(cbind(sites, ecoregion_three, ecoregion_one, 
       #      site_intercepts, ecoregion_three_intercepts, ecoregion_one_intercepts,
       #      psi_site_nested))
  
  ## effect of pop density on occupancy (species-specific random slopes)
  psi_herb_shrub_forest <- rnorm(n=n_species, mean=mu_psi_herb_shrub_forest, sd=sigma_psi_herb_shrub_forest)
  # change in each species occupancy across time is drawn from a distribution defined
  # by a community mean (mu_psi_interval) with 
  # species specific variation defined by sigma_psi_interval
  
  ## effect of pop density on occupancy (species-specific random slopes)
  psi_open_developed <- rnorm(n=n_species, mean=mu_psi_open_developed, sd=sigma_psi_open_developed)
  # change in each species occupancy across time is drawn from a distribution defined
  # by a community mean (mu_psi_interval) with 
  # species specific variation defined by sigma_psi_interval
  
  ## --------------------------------------------------
  ### specify species-specific detection probabilities
  
  ## species-specific random intercepts
  p_citsci_species  <- rnorm(n=n_species, mean = 0, sd=sigma_p_citsci_species)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p_species
  
  ## effect of site and interval on detection (site,interval-specific random slopes)
  ecoregion_three_intercepts_p_citsci <- rep(rnorm(n=n_ecoregion_three, mean=0, sd=sigma_p_citsci_ecoregion_three),
                                             each=n_sites_per_ecoregion_three)
  
  site_intercepts_p_citsci <- rep(rnorm(n=n_sites, mean=0, sd=sigma_p_citsci_site))
  
  p_citsci_site_nested <- vector(length=n_sites)
  
  for(i in 1:n_sites){
    
    p_citsci_site_nested[i] <- ecoregion_three_intercepts_p_citsci[i] + site_intercepts_p_citsci[i]
    
  }
  
  ## species-specific random intercepts
  p_museum_species  <- rnorm(n=n_species, mean = 0, sd=sigma_p_museum_species)
  # species baseline detection probability is drawn from a normal distribution with mean 0 and 
  # species specific variation defined by sigma_p_species
  
  ## effect of site and interval on detection (site,interval-specific random slopes)
  ## effect of site and interval on detection (site,interval-specific random slopes)
  ecoregion_three_intercepts_p_museum <- rep(rnorm(n=n_ecoregion_three, mean=0, sd=sigma_p_museum_ecoregion_three),
                                             each=n_sites_per_ecoregion_three)
  
  site_intercepts_p_museum <- rep(rnorm(n=n_sites, mean=0, sd=sigma_p_museum_site))
  
  p_museum_site_nested <- vector(length=n_sites)
  
  for(i in 1:n_sites){
    
    p_museum_site_nested[i] <- ecoregion_three_intercepts_p_museum[i] + site_intercepts_p_museum[i]
    
  }
  
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
            psi_species_nested[species] + # a species specific intercept
            psi_site_nested[site] + # a site specific intercept
            psi_herb_shrub_forest[species]*herb_shrub_forest[site] + # a species specific effect
            psi_open_developed[species]*open_developed[site] + # a species specific effect
            psi_site_area*site_area[site] # a fixed effect of site area
        
        for(visit in 1:n_visits) { # for each visit
          
          logit_p_matrix_citsci[species, site, interval, visit] <- # detection is equal to 
            mu_p_citsci_0 + # a baseline intercept
              p_citsci_species[species] + # a species specific intercept
              p_citsci_site_nested[site] + # a spatiotemporally specific intercept
              p_citsci_interval*intervals[interval] + # an overall effect of time on detection
              p_citsci_pop_density*pop_density[site] # an effect of population density on detection ability
          
          logit_p_matrix_museum[species, site, interval, visit] <- # detection is equal to 
            mu_p_museum_0 + # a baseline intercept
              p_museum_species[species] + # a species specific intercept
              p_museum_site_nested[site] + # a spatiotemporally specific intercept
              p_museum_total_records*total_records_museum[site,interval]
          
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
    n_species = n_species, # number of species,
    n_genera = n_genera,
    genus = genus,
    genus_lookup = genus_lookup,
    n_ecoregion_three = n_ecoregion_three,
    n_ecoregion_one = n_ecoregion_one,
    ecoregion_three = ecoregion_three,
    ecoregion_one = ecoregion_one,
    ecoregion_three_lookup = ecoregion_three_lookup,
    ecoregion_one_lookup = ecoregion_one_lookup,
    n_sites = n_sites, # number of sites
    n_intervals = n_intervals, # number of surveys 
    n_visits = n_visits, # number of visits
    pop_density = pop_density, # vector of pop densities
    open_developed = open_developed, # vector of impervious surface covers
    herb_shrub_forest = herb_shrub_forest, # vector of perennial plant covers
    site_area = site_area, # vector of site areas
    total_records_museum = total_records_museum # museum records per interval
  ))
  
} # end simulate_data function


## --------------------------------------------------
### Variable values for data simulation
## study dimensions
n_genera = 10 ## number of genera
n_species_per_genera = 4 ## number of species
n_species = n_genera*n_species_per_genera
n_ecoregion_one = 7
n_ecoregion_three_per_one = 7 # ecoregion3 per ecoregion1
n_ecoregion_three = n_ecoregion_one*n_ecoregion_three_per_one
n_sites_per_ecoregion_three = 8
n_sites = n_sites_per_ecoregion_three*n_ecoregion_three ## number of sites
n_intervals = 3 ## number of occupancy intervals
n_visits = 3 ## number of samples per year

## occupancy
mu_psi_0 = -0.5
sigma_psi_species = 0.5
sigma_psi_genus = 0.25
sigma_psi_site = 0.5 # variation across level2
sigma_psi_ecoregion_three = 0.5 # variation across level3
sigma_psi_ecoregion_one = 0.25 # variation across level4
mu_psi_open_developed = -0.25
sigma_psi_open_developed = 0.5
mu_psi_herb_shrub_forest = 0.5 
sigma_psi_herb_shrub_forest = 0.5
psi_site_area = 0.75 # fixed effect of site area on occupancy

## detection
# citizen science observation process
mu_p_citsci_0 = -1
sigma_p_citsci_species = 0.5
sigma_p_citsci_ecoregion_three = 0.3 # variation across level3
sigma_p_citsci_site = 0.3
p_citsci_interval = 1
p_citsci_pop_density = 1 

# museum record observation process
mu_p_museum_0 = -1
sigma_p_museum_species = 0.5
sigma_p_museum_ecoregion_three = 0.3 # variation across level3
sigma_p_museum_site = 0.3
p_museum_total_records = 1

# introduce NAs (missed visits)?
sites_missing = 0.5*n_sites 
intervals_missing = 2
visits_missing = 1

sites_in_range_beta1 = 2
sites_in_range_beta2 = 2

## --------------------------------------------------
### Simulate data
set.seed(3)
my_simulated_data <- simulate_data(n_genera,
                                   n_species_per_genera,
                                   n_species = n_genera*n_species_per_genera,
                                   n_ecoregion_one,
                                   n_ecoregion_three_per_one,
                                   n_ecoregion_three,
                                   n_sites_per_ecoregion_three,
                                   n_sites,
                                   n_intervals,
                                   n_visits,
                                   
                                   # ecological process
                                   mu_psi_0,
                                   sigma_psi_species,
                                   sigma_psi_genus,
                                   sigma_psi_site,
                                   sigma_psi_ecoregion_three,
                                   sigma_psi_ecoregion_one,
                                   mu_psi_open_developed,
                                   sigma_psi_open_developed,
                                   mu_psi_herb_shrub_forest,
                                   sigma_psi_herb_shrub_forest,
                                   psi_site_area,
                                  
                                   # citizen science observation process
                                   mu_p_citsci_0,
                                   sigma_p_citsci_species,
                                   sigma_p_citsci_site,
                                   sigma_p_citsci_ecoregion_three,
                                   p_citsci_interval,
                                   p_citsci_pop_density, 
                                   
                                   # museum record observation process
                                   mu_p_museum_0,
                                   sigma_p_museum_species,
                                   sigma_p_museum_site,
                                   sigma_p_museum_ecoregion_three,
                                   p_museum_total_records,
                                   
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
n_ecoregion_three <- my_simulated_data$n_ecoregion_three
n_ecoregion_one <- my_simulated_data$n_ecoregion_one
n_intervals <- my_simulated_data$n_intervals # number of surveys 
n_visits <- my_simulated_data$n_visits

#View(as.data.frame(V_citsci[1:10,1:10,,]))
#View(as.data.frame(V_museum[1:10,1:10,,]))
sum(my_simulated_data$V_citsci == 1)
sum(my_simulated_data$V_museum == 1)

check_citsci <- which(V_citsci>ranges)
check_museum <- which(V_museum>V_museum_NA)
 
intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
ecoregion_three <- my_simulated_data$ecoregion_three
ecoregion_one <- my_simulated_data$ecoregion_one
ecoregion_three_lookup <- my_simulated_data$ecoregion_three_lookup
ecoregion_one_lookup <- my_simulated_data$ecoregion_one_lookup
species <- seq(1, n_species, by=1)
n_genera <- my_simulated_data$n_genera
genus_lookup <- my_simulated_data$genus_lookup

pop_densities <- my_simulated_data$pop_density
site_areas <- my_simulated_data$site_area
open_developed <- my_simulated_data$open_developed
herb_shrub_forest <- my_simulated_data$herb_shrub_forest
museum_total_records <- my_simulated_data$total_records_museum

stan_data <- c("V_citsci", "V_museum", 
               "ranges", "V_museum_NA",
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", 
               "species", "n_genera", "genus_lookup",
               "sites", "n_ecoregion_three", "n_ecoregion_one",
               "ecoregion_three_lookup", "ecoregion_one_lookup",
               "pop_densities", "site_areas", "open_developed", 
               "herb_shrub_forest", "museum_total_records") 

# Parameters monitored
params <- c("mu_psi_0",
            "sigma_psi_species",
            "sigma_psi_genus",
            "sigma_psi_site",
            "sigma_psi_ecoregion_three",
            "sigma_psi_ecoregion_one",
            "mu_psi_open_developed",
            "sigma_psi_open_developed",
            "mu_psi_herb_shrub_forest",
            "sigma_psi_herb_shrub_forest",
            "psi_site_area",
            
            "mu_p_citsci_0",
            "sigma_p_citsci_species",
            "sigma_p_citsci_site",
            "sigma_p_citsci_ecoregion_three",
            "p_citsci_interval",
            "p_citsci_pop_density", 
            
            "mu_p_museum_0",
            "sigma_p_museum_species",
            "sigma_p_museum_site",
            "sigma_p_museum_ecoregion_three",
            "p_museum_total_records"#,
            
            #"T_rep_citsci",
            #"T_obs_citsci",
            #"P_species_citsci",
            
            #"T_rep_museum",
            #"T_obs_museum",
            #"P_species_museum"
)

parameter_value <- c(mu_psi_0,
                     sigma_psi_species,
                     sigma_psi_genus,
                     sigma_psi_site,
                     sigma_psi_ecoregion_three,
                     sigma_psi_ecoregion_one,
                     mu_psi_open_developed,
                     sigma_psi_open_developed,
                     mu_psi_herb_shrub_forest,
                     sigma_psi_herb_shrub_forest,
                     psi_site_area,
                     
                     mu_p_citsci_0,
                     sigma_p_citsci_species,
                     sigma_p_citsci_site,
                     sigma_p_citsci_ecoregion_three,
                     p_citsci_interval,
                     p_citsci_pop_density,
                     
                     mu_p_museum_0,
                     sigma_p_museum_species,
                     sigma_p_museum_site,
                     sigma_p_museum_ecoregion_three,
                     p_museum_total_records#,
                     
                     #NA,
                     #NA,
                     #NA,
                     #NA,
                     #NA,
                     #NA
)

# MCMC settings
n_iterations <- 800
n_thin <- 1
n_burnin <- 400
n_chains <- 4
n_cores <- 4
delta = 0.9

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
set.seed(2)
inits <- lapply(1:n_chains, function(i)
  
  list(mu_psi_0 = runif(1, -1, 1),
       sigma_psi_species = runif(1, 0, 1),
       sigma_psi_site = runif(1, 0, 1),
       sigma_psi_ecoregion_three = runif(1, 0, 1),
       sigma_psi_ecoregion_one = runif(1, 0, 1),
       mu_psi_open_developed = runif(1, -1, 1),
       sigma_psi_open_developed = runif(1, 0, 1),
       mu_psi_herb_shrub_forest = runif(1, -1, 1),
       sigma_psi_herb_shrub_forest = runif(1, 0, 1),
       psi_site_area = runif(1, -1, 1),
       
       mu_p_citsci_0 = runif(1, -1, 0),
       sigma_p_citsci_species = runif(1, 0, 1),
       sigma_p_citsci_site = runif(1, 0, 1),
       sigma_p_citsci_ecoregion_three = runif(1, 0, 1),
       p_citsci_interval = runif(1, -1, 1),
       p_citsci_pop_density = runif(1, -1, 1),
       
       mu_p_museum_0 = runif(1, -1, 0),
       sigma_p_museum_species = runif(1, 0, 1),
       sigma_p_museum_site = runif(1, 0, 1),
       sigma_p_museum_ecoregion_three = runif(1, 0, 1),
       p_museum_total_records = runif(1, -1, 1)
       
  )
)

targets <- as.data.frame(cbind(params, parameter_value))

## --------------------------------------------------
### Run model
library(rstan)
stan_model <- "./occupancy/models/model_syrphidae.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     open_progress = FALSE,
                     control=list(adapt_delta=delta),
                     cores = n_cores)

print(stan_out_sim, digits = 3)
View(targets)

saveRDS(stan_out_sim, "./occupancy/simulation/stan_out_sim_syrphidae.rds")
stan_out_sim <- readRDS("./occupancy/simulation/stan_out_sim_syrphidae.rds")

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out_sim, pars = c(
  "mu_psi_0",
  "mu_psi_herb_shrub_forest",
  "mu_psi_open_developed",
  "mu_p_citsci_0",
  "p_citsci_interval",
  "p_citsci_pop_density",
  "mu_p_museum_0",
  "p_museum_total_records"
))

# traceplot
traceplot(stan_out_sim, pars = c(
  "sigma_psi_species",
  "sigma_psi_genus",
  "sigma_psi_site",
  "sigma_psi_ecoregion_three",
  "sigma_psi_ecoregion_one",
  "sigma_psi_open_developed",
  "sigma_psi_herb_shrub_forest",
  "sigma_p_citsci_site",
  "sigma_p_citsci_ecoregion_three",
  "sigma_p_museum_site",
  "sigma_p_museum_ecoregion_three",
  "sigma_p_citsci_species",
  "sigma_p_museum_species"
))

# pairs plot
pairs(stan_out_sim, pars = c(
  "mu_psi_0",
  "mu_p_citsci_0",
  "mu_p_museum_0",
  
  "sigma_psi_open_developed",
  "sigma_psi_herb_shrub_forest"
))

## --------------------------------------------------
### PPC

# print rep and obs
print(stan_out_sim, digits = 3, pars=
        c("T_rep_citsci", "T_obs_citsci", "P_species_citsci"))
print(stan_out_sim, digits = 3, pars=
        c("T_rep_museum", "T_obs_museum", "P_species_museum"))

# as data frame
list_of_draws <- as.data.frame(stan_out_sim)
list_of_draws <- list_of_draws[(n_burnin+1):n_iterations,]

# Citizen Science
# P-values
m <- n_iterations - n_burnin
P_average_citsci = vector(length = n_species)

for(i in 1:n_species){
  P_average_citsci[i] = sum(list_of_draws[,88+i])/m
}

print(P_average_citsci)

# Evaluation of fit
plot(list_of_draws[,54], list_of_draws[,19], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 100),
     xlim = c(0, 100))

abline(0, 1, lwd = 2, col = "black")

# Museum
# P-values
m <- n_iterations - n_burnin
P_average_museum = vector(length = n_species)

for(i in 1:n_species){
  P_average_museum[i] = sum(list_of_draws[,193+i])/m
}

print(P_average_museum)

# Evaluation of fit
plot(list_of_draws[,124], list_of_draws[,159], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 100),
     xlim = c(0, 100))

abline(0, 1, lwd = 2, col = "black")

