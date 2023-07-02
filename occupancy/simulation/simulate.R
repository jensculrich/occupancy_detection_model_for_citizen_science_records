### Data simulation for citizen science/museum records of pollinator occurrence
# jcu; started nov 21, 2022

## Data simulation for integrated model form
# single season, species considered across all sites (not just range)
# and across all visits (not considering that sites without any detections 
# for any species were probably not visited or that sites with only detections for
# a few species only might be targeted and not fully community wide)

library(rstan)

# choose a taxon group
taxon = "syrphidae"
taxon = "bombus"

# then you will enter the sim data function
# then you will enter the parameter values (targets)
# then run the simulation
# then prepare the simmed data for stan
# then run and save the stan model
# then inspect the results in comparison with the targets
# as well as model run diagnostics

# model fitting using 'model_simplest.stan' should return the parameter inputs
simulate_data <- function(taxon,
                          n_genera,
                          n_species_per_genera,
                          n_species,
                          n_level_four,
                          n_level_three_per_one,
                          n_level_three,
                          n_sites_per_level_three,
                          n_sites,
                          n_intervals,
                          n_visits,
                          
                          # ecological process
                          mu_psi_0,
                          sigma_psi_species,
                          sigma_psi_genus,
                          sigma_psi_site,
                          sigma_psi_level_three,
                          sigma_psi_level_four,
                          sigma_psi_income_level_three,
                          sigma_psi_income_level_four,
                          mu_psi_income,
                          sigma_psi_income,
                          mu_psi_natural_habitat,
                          sigma_psi_natural_habitat,
                          delta0,
                          delta1,
                          psi_site_area,
                          
                          # citizen science observation process
                          mu_p_cs_0,
                          sigma_p_cs_species,
                          sigma_p_cs_site,
                          sigma_p_cs_level_three,
                          sigma_p_cs_level_four,
                          p_cs_interval,
                          p_cs_pop_density, 
                          
                          # museum record observation process
                          mu_p_rc_0,
                          sigma_p_rc_species,
                          sigma_p_rc_site,
                          sigma_p_rc_level_three,
                          sigma_p_rc_level_four,
                          p_rc_total_records,
                          
                          # correlation (detection)
                          rho,
                          
                          # introduce NAs (missed visits)?
                          omega_community,
                          omega_species,
                          ignore_community_misses,
                          #sites_missing, 
                          #intervals_missing,
                          #visits_missing,
                          
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
  
  # create a vector of site area
  site_area <- rnorm(n_sites, mean = 0, sd = 1)
  
  # Realistically natural habitat area is positively correlated with income
  # and negatively with population density
  # I will simluate data with these realistic correlations
  mu <- c(0, 0, 0)
  stddev <- c(1, 1, 1)
  corMat <- matrix(c(1, 0, -0.25, 
                     0, 1, 0.15,
                     -0.25, 0.15, 1),
                   ncol = 3)
  covMat <- stddev %*% t(stddev) * corMat
  correlated_data <- MASS::mvrnorm(n = n_sites, mu = mu, Sigma = covMat, empirical = FALSE)
  
  pop_density <- correlated_data[,1]
  income <- correlated_data[,2]
  natural_habitat <- correlated_data[,3]
  
  plot(correlated_data[,1], correlated_data[,2])
  cor(correlated_data[,1], correlated_data[,2])
  
  plot(correlated_data[,1], correlated_data[,3])
  cor(correlated_data[,1], correlated_data[,3])
  
  plot(correlated_data[,2], correlated_data[,3])
  cor(correlated_data[,2], correlated_data[,3])
  
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
  total_records_rc <- MASS::mvrnorm(n = n_sites, mu = mu2, Sigma = covMat2, empirical = FALSE)
  
  ## --------------------------------------------------
  ### Introduce correlations between species occupancy and species detection rates 
  
  species_mu <- c(mu_p_cs_0, mu_p_rc_0)

  species_covMat <- matrix(c(sigma_p_cs_species, rho, 
                     rho, sigma_p_rc_species),
                   ncol = 2)
  
  correlated_data2 <- MASS::mvrnorm(n = n_species, 
                                   mu = species_mu, 
                                   Sigma = species_covMat, 
                                   empirical = FALSE)
  
  plot(correlated_data2[,1], correlated_data2[,2])
  cor(correlated_data2[,1], correlated_data2[,2])
  
  ## --------------------------------------------------
  ### Ecological process random effects
  
  ## --------------------------------------------------
  ### specify species-specific occupancy probabilities
  
  if(taxon == "syrphidae"){
    
    ## genus-specific random intercepts
    #genus = rep(1:n_genera, each = n_species_per_genera)
    #genera_intercepts <- rep(rnorm(n=n_genera, mean=mu_psi_0, sd=sigma_psi_genus),
    #                         each=n_species_per_genera)
    
    ## species-specific random intercepts
    #species_intercepts <- rnorm(n=n_species, mean=0, sd=sigma_psi_species)
    # species baseline occupancy is drawn from a normal distribution with mean 0 and 
    # species specific variation defined by sigma.psi.sp
    
    #genus_lookup <- rep(1:n_genera, each = n_species_per_genera)
    
    #psi_species <- vector(length=n_species)
    
    #for(i in 1:n_species){
      
     # psi_species[i] <- genera_intercepts[i] + species_intercepts[i]
      
    #}
    
    ## species-specific random intercepts
    psi_species <- rnorm(n=n_species, mean=mu_psi_0, sd=sigma_psi_species)
    
    genera_intercepts <- NULL
    genus_lookup <- NULL
    
    ## species-specific random intercepts
    species_intercepts <- psi_species
    # species baseline occupancy is drawn from a normal distribution with mean 0 and 
    # species specific variation defined by sigma.psi.sp
    
  } else { # else taxon is bombus (only one genus so there's no among genus variation)
    
    ## species-specific random intercepts
    psi_species <- rnorm(n=n_species, mean=mu_psi_0, sd=sigma_psi_species)
    
    genera_intercepts <- NULL
    genus_lookup <- NULL
    
    ## species-specific random intercepts
    species_intercepts <- psi_species
    # species baseline occupancy is drawn from a normal distribution with mean 0 and 
    # species specific variation defined by sigma.psi.sp
  
  }
  
  
  ## --------------------------------------------------
  ### specify spatially-specific, spatially-nested occupancy probabilities
  
  ## define spatial clusters
  ecoregion_one = rep(1:n_level_four, each = n_level_three_per_one*n_sites_per_level_three)
  level_three = rep(1:n_level_three, each = n_sites_per_level_three)
  ## provide lookup references for STAN
  ecoregion_one_lookup <- rep(1:n_level_four, each=n_level_three_per_one)
  level_three_lookup <- rep(1:n_level_three, each=n_sites_per_level_three)
  
  ## ecoregion1-specific random intercepts
  ecoregion_one_intercepts <- rep(rnorm(n=n_level_four, mean=0, sd=sigma_psi_level_four),
                         each=n_level_three_per_one*n_sites_per_level_three)
  
  ## ecoregion3-specific random intercepts
  level_three_intercepts <- rep(rnorm(n=n_level_three, mean=0, sd=sigma_psi_level_three),
                                  each=n_sites_per_level_three)
  
  ## site-specific random intercepts
  site_intercepts <- rnorm(n=n_sites, mean=0, sd=sigma_psi_site)
  
  # combine intercept adjustments for each site
  psi_site_nested <- vector(length=n_sites)
  
  for(i in 1:n_sites){
    
    psi_site_nested[i] <-
      ecoregion_one_intercepts[i] + 
      level_three_intercepts[i] + 
      site_intercepts[i]
    
  }
  
  #View(as.data.frame(cbind(site_intercepts, level_three_intercepts, ecoregion_one_intercepts, psi_site_nested)))
  
  ## --------------------------------------------------
  ### random slope effects on occupancy
  
  nativity = NULL
  
  if(taxon == "syrphidae"){
    
    nativity = rbinom(n_species, 1, 0.75)
    
  }
  
  ## effect of natural habitat on occupancy (species-specific random slopes)
  if(taxon == "syrphidae"){
    
    mu_psi_natural_habitat = vector(length = n_species)
    
    sigma_psi_natural_habitat = vector(length = n_species)
    
    for(i in 1:n_species){
      mu_psi_natural_habitat[i] = delta0 + delta1*nativity[i]
    }
    
    for(i in 1:n_species){
      sigma_psi_natural_habitat[i] = gamma0 + gamma1*nativity[i]
    }
    
    psi_natural_habitat = vector(length = n_species)
    
    for(i in 1:n_species){
      psi_natural_habitat[i] <- rnorm(n=1, mean=mu_psi_natural_habitat[i], sd=sigma_psi_natural_habitat[i])
    }
    
  } else {
    psi_natural_habitat <- rnorm(n=n_species, mean=mu_psi_natural_habitat, sd=sigma_psi_natural_habitat)
  }
  
  
  
  ## luxury effect on occupancy (species-specific random slopes)
  psi_income <- rnorm(n=n_species, mean=mu_psi_income, sd=sigma_psi_income)
  
  ## --------------------------------------------------
  ### Observation processes random effects
  
  ## Citsci
  
  ## --------------------------------------------------
  ### specify species-specific detection probabilities
  
  ## species-specific random intercepts
  if(taxon == "bombus"){
    p_cs_species <- correlated_data2[,1]
  } else {
    p_cs_species  <- rnorm(n=n_species, mean = mu_p_cs_0, sd=sigma_p_cs_species)
  }
  
  ## --------------------------------------------------
  ### specify spatially-specific, spatially-nested detection probabilities
  
  ## effect of site on detection 
  ecoregion_one_intercepts_p_cs <- rep(rnorm(n=n_level_four, mean=0, sd=sigma_p_cs_level_four),
                                           each=n_level_three_per_one*n_sites_per_level_three)
  
  ## effect of site on detection 
  level_three_intercepts_p_cs <- rep(rnorm(n=n_level_three, mean=0, sd=sigma_p_cs_level_three),
                                             each=n_sites_per_level_three)
  
  site_intercepts_p_cs <- rnorm(n=n_sites, mean=0, sd=sigma_p_cs_site)
  
  p_cs_site_nested <- vector(length=n_sites)
  
  for(i in 1:n_sites){
    
    p_cs_site_nested[i] <- 
      ecoregion_one_intercepts_p_cs[i] +
      level_three_intercepts_p_cs[i] + 
      site_intercepts_p_cs[i]
    
  }
  
  # View(as.data.frame(cbind(site_intercepts_p_cs, level_three_intercepts_p_cs, ecoregion_one_intercepts_p_cs, p_cs_site_nested)))
  
  ## Museum
  
  ## species-specific random intercepts
  if(taxon == "bombus"){
    p_rc_species <- correlated_data2[,2]
  } else {
    p_rc_species  <- rnorm(n=n_species, mean = mu_p_rc_0, sd=sigma_p_rc_species)
  }
  
  ## --------------------------------------------------
  ### specify spatially-specific, spatially-nested detection probabilities
  
  ## effect of site on detection 
  ecoregion_one_intercepts_p_rc <- rep(rnorm(n=n_level_four, mean=0, sd=sigma_p_rc_level_four),
                                           each=n_level_three_per_one*n_sites_per_level_three)
  
  ## effect of site on detection 
  level_three_intercepts_p_rc <- rep(rnorm(n=n_level_three, mean=0, sd=sigma_p_rc_level_three),
                                             each=n_sites_per_level_three)
  
  site_intercepts_p_rc <- rnorm(n=n_sites, mean=0, sd=sigma_p_rc_site)
  
  p_rc_site_nested <- vector(length=n_sites)
  
  for(i in 1:n_sites){
    
    p_rc_site_nested[i] <- 
      ecoregion_one_intercepts_p_rc[i] +
      level_three_intercepts_p_rc[i] + 
      site_intercepts_p_rc[i]
    
  }
  
  # View(as.data.frame(cbind(site_intercepts_p_rc, level_three_intercepts_p_rc, ecoregion_one_intercepts_p_rc, p_rc_site_nested)))
  
  
  # spatiotemporal variability in detection probability (changing across sites and 
  # occupancy intervals), and helps account for the variation that is inherent in 
  # sample effort across space and time in unstructured historical datasets. 
  # spatiotemporal variability at each site*interval is drawn from a distribution defined
  # by a mean of 0 and with site*interval specific variation of sigma.p.site
  
  ## --------------------------------------------------
  ## Create arrays for psi and p
  logit_psi_matrix <- array(NA, dim =c(n_species, n_sites, n_intervals)) 
  # a psi value for each species, at each site, in each interval 
  
  logit_p_matrix_cs <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  logit_p_matrix_rc <- array(NA, dim =c(n_species, n_sites, n_intervals, n_visits))
  # a p value for each species, at each site, in each interval, AND in each visit
  
  for(species in 1:n_species) { # for each site
    for(site in 1:n_sites) { # for each interval
      for(interval in 1:n_intervals) { # for each species
        
        logit_psi_matrix[species, site, interval] <- # occupancy is equal to
            psi_species[species] + # a species specific intercept
            psi_site_nested[site] + # a site specific intercept
            psi_natural_habitat[species]*natural_habitat[site] + # a species specific effect
            psi_income[species]*income[site] + # a species specific effect
            psi_site_area*site_area[site] # a fixed effect of site area
        
        for(visit in 1:n_visits) { # for each visit (but sim constant rates across visits)
          
          logit_p_matrix_cs[species, site, interval, visit] <- # detection is equal to 
              p_cs_species[species] + # a species specific intercept
              p_cs_site_nested[site] + # a spatiotemporally specific intercept # includes global intercept
              p_cs_interval*(intervals[interval]^2) + # an overall effect of time on detection
              p_cs_pop_density*pop_density[site] # an effect of population density on detection ability
          
          logit_p_matrix_rc[species, site, interval, visit] <- # detection is equal to 
              p_rc_species[species] + # a species specific intercept
              p_rc_site_nested[site] + # a spatiotemporally specific intercept # includes global intercept
              p_rc_total_records*total_records_rc[site,interval]
          
        } # for each visit
      } # for each species
    } # for each interval
  } # for each site
  
  # preview the psi and p arrays
  head(logit_psi_matrix[1:n_species, 1:n_sites,1])
  # head(p_matrix[1:n_species, 1:n_sites,7,1])
  
  #(logit_p_matrix_cs[1,1,1:n_intervals,1]) # if p.interval is >0 these should generally be increasing from low to high
  #(logit_p_matrix_cs[2,1,1:n_intervals,1]) # if p.interval is >0 these should generally be increasing from low to high
  #(logit_p_matrix_cs[1,2,1:n_intervals,1]) # if p.interval is >0 these should generally be increasing from low to high
  
  
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
  
  V_cs <- array(NA, dim=c(n_species=n_species,
                       n_sites=n_sites,
                       n_intervals=n_intervals,
                       n_visits=n_visits))
  
  for(interval in 1:n_intervals){
    for(site in 1:n_sites){
      for(species in 1:n_species){
        for(visit in 1:n_visits){
          
          # detection as a default only occurs at sites within the range of a species
          # those outside of the range will be 'skipped' by the likelihood function
          V_cs[species,site,interval,visit] <- Z[species,site,interval] * # occupancy state * detection prob
            rbinom(n = 1, size = 1, 
                   prob = ilogit(logit_p_matrix_cs[species,site,interval,visit]))
          
        }
      }
    }
  } # end simulate detection data
  
  ## --------------------------------------------------
  # citizen science NAs
  
  # 1 indicates a site was in range and thus the species could be detected there
  V_cs_NA <- ranges
  # should NEVER have a V_cs detection outside of the range
  #check <- which(V_cs>V_cs_NA)
  
  ## --------------------------------------------------
  # Museum detections
  
  # instead let's set up a hurdle for whether a sampling event will occur
  # if the hurdle is passed, then the community is surveyed
  #omega_community <- 0.3 # rate of "passing" and doing a community survey
  #omega_species <- 0.05 # rate of "passing" and doing a species specific survey when a comm sample doesnt happen
  
  V_rc <- array(NA, dim=c(n_species=n_species,
                              n_sites=n_sites,
                              n_intervals=n_intervals,
                              n_visits=n_visits))
  
  community_sampled <- array(0, dim=c(n_species=n_species,
                                       n_sites=n_sites,
                                       n_intervals=n_intervals,
                                       n_visits=n_visits))
  
  non_comm_sample <- array(0, dim=c(n_species=n_species,
                                       n_sites=n_sites,
                                       n_intervals=n_intervals,
                                       n_visits=n_visits))
  
  any_sampled <- array(0, dim=c(n_species=n_species,
                                      n_sites=n_sites,
                                      n_intervals=n_intervals,
                                      n_visits=n_visits))
  
  for(site in 1:n_sites){
    for(interval in 1:n_intervals){
      for(visit in 1:n_visits){
        
        # determine if the community was sampled at site during time
        community_sampled[1:n_species,site,interval,visit] <- 
          rbinom(n = 1, size = 1, prob = omega_community)
          
        for(species in 1:n_species){
          
          # if a community sample didn't occur (<1), did a species specific sampling event occur?
          if(community_sampled[species,site,interval,visit] < 1){
            non_comm_sample[species,site,interval,visit] <- 
              rbinom(n = 1, size = 1, prob = omega_species)*
              ranges[species,site,1,1] # if the site is in range
          } 
          
          # if either comm or species sampling occurred then the species was sampled (but maybe not detected)
          if(community_sampled[species,site,interval,visit] > 0 || 
             non_comm_sample[species,site,interval,visit] > 0){
            
            any_sampled[species,site,interval,visit] <- 1
            
          }
          
          # Determine species detections
          V_rc[species,site,interval,visit] <- 
            any_sampled[species,site,interval,visit] * # whether a community survey or species survey occurred 
            ranges[species,site,1,1] * # whether a site is actually in the species range 
            Z[species,site,interval] * # occupancy state 
            rbinom(n = 1, size = 1,
                   prob = ilogit(logit_p_matrix_rc[species,site,interval,visit]))  # detection prob
          
        }
      }
    }
  } # end simulate detection data
  
  # check should always be empty (should not be getting detections where a community survey did not occur)
  #check <- which(V_rc>any_sampled)
  #check <- which(V_rc>ranges)
  
  ## --------------------------------------------------
  # museum NAs
  
  ## --------------------------------------------------
  ## Generate NA indicators (sampling could not have occurred either because
  # a species range does not overlap with the site or no community samples were drawn)
  
  # simulate NA data for the model to work around
  # for our real data we will have some detections that cannot occur because
  # 1) the site is not in the range of the species, or 
  # 2) community sampling didnèt occur at some sites for some visits in some intervals
  
  V_rc_NA <- any_sampled*ranges
  # check <- which(V_rc>V_rc_NA)
  
  # Do we ever have community sampling events where no species were detected?
  # which community_sampled[1:n_species,site,interval,visit] > 0 BUT...
  # V_rc[1:n_species,site,interval,visit] == 0
  non_comm_sample[1:n_species,1:10,1,1]
  community_sampled[1:n_species,1:10,1,1]
  ranges[1:n_species,1:10,1,1]
  any_sampled[1:n_species,1:10,1,1]
  V_rc[1:n_species,1:10,1,1]
  
  counter = 0
  counter2 = 0
  
  for(site in 1:n_sites){
    for(interval in 1:n_intervals){
      for(visit in 1:n_visits){
        
        # if the community was sampled, 
        if(community_sampled[1,site,interval,visit] > 0 &&
           # but no species was detected
           sum(V_rc[1:n_species,site,interval,visit]) == 0){
          
          # increase our counter by 1
          counter <- counter + 1
        }
        
        if(community_sampled[1, site,interval,visit] > 0){
          counter2 <- counter2 + 1
        }
        
      }
    }
  }
  
  # given our detection rates, how frequently does a comm sample occur but no species are detected? 
  (ratio_fully_missed_community_samples <- counter/counter2)
  
  # Now let's try to remove community_sampled for all species (make 0)
  # if no one was detected, and see if it breaks our param estimates
  # This is the form that we get our data in, so we want to see how much it influences the 
  # paramter estimation process if it occurs with ratio_fully_missed_community_samples frequency
  V_rc_NA_new <- V_rc_NA
  for(site in 1:n_sites){
    for(interval in 1:n_intervals){
      for(visit in 1:n_visits){
        
        # if the community was sampled, 
        if(community_sampled[1,site,interval,visit] > 0 &&
           # but no species was detected
           sum(V_rc[1:n_species,site,interval,visit]) == 0){
          
          # replace V_rc_NA_new with 0's for all species
          V_rc_NA_new[1:n_species,site,interval,visit] <- 0
        }
        
      }
    }
  }
  
  #check <- which(V_rc>V_rc_NA_new)
  
  # if we want to ignore cases where community was surveyed but no species detected
  # i.e. treat them like they never happened, we just replace those comm surveys with 0's
  # this is more similar to what happens naturally in the data collection process because we have no info 
  # on whether a community survey occurred, we just get information on the outcome of them occurring.
  # In truth I suspect it's rare that people go out looking for a community at a grid cell size, year long scale
  # and never detect ANY species from the target community. 
  if(ignore_community_misses == TRUE){
    # use the na structure where we drop 1's from columns where no comm survey would be inferred
    # by our data collection and preparation process.
    V_rc_NA <- V_rc_NA_new
  }
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    V_cs = V_cs, # detection data from citizen science records
    V_rc = V_rc, # detection data from museum records
    ranges = ranges, # array indicating whether sampling occurred in a site*interval*visit
    V_rc_NA = V_rc_NA, # array indicating whether sampling occurred in a site*interval*visit
    n_species = n_species, # number of species,
    n_genera = n_genera,
    genus_lookup = genus_lookup,
    n_level_three = n_level_three,
    n_level_four = n_level_four,
    level_three = level_three,
    ecoregion_one = ecoregion_one,
    level_three_lookup = level_three_lookup,
    ecoregion_one_lookup = ecoregion_one_lookup,
    n_sites = n_sites, # number of sites
    n_intervals = n_intervals, # number of surveys 
    n_visits = n_visits, # number of visits
    pop_density = pop_density, # vector of pop densities
    income = income, # vector of income levels
    natural_habitat = natural_habitat, # vector of nat habitat cover
    site_area = site_area, # vector of site areas
    total_records_rc = total_records_rc, # museum records per interval
    species_intercepts = species_intercepts,
    genera_intercepts = genera_intercepts,
    nativity = nativity
  ))
  
} # end simulate_data function


## --------------------------------------------------
### Variable values for data simulation

if(taxon == "syrphidae"){
  
  ## study dimensions
  n_genera = 1 # us 1 unless you want to introduce generic variation
  n_species_per_genera = 100 ## number of species
  n_species = n_genera*n_species_per_genera
  n_level_four = 10
  n_level_three_per_one = 7 # ecoregion3 per ecoregion1
  n_level_three = n_level_four*n_level_three_per_one
  n_sites_per_level_three = 7
  n_sites = n_sites_per_level_three*n_level_three ## number of sites
  
  ## currently the civariance matrix for total records is set up for 3 intervals
  n_intervals = 3 ## number of occupancy intervals # if not three will have to change columns in cormatrix2
  n_visits = 3 ## number of samples per year
  
  ## occupancy
  mu_psi_0 = 0
  sigma_psi_species = 1.5
  sigma_psi_genus = 0  # this won't be used if taxon is "bombus", keep at 0 unless want to introduce
  sigma_psi_site = 1.5 # variation across level2
  sigma_psi_level_three = 1 # variation across level3
  sigma_psi_level_four = 1 # variation across level4
  mu_psi_income = 0 # make sure to specify as 0 if using a model without income 
  sigma_psi_income = 0 # make sure to specify as 0 if using a model without income 
  mu_psi_natural_habitat = 0.75 
  sigma_psi_natural_habitat = 0.75
  delta0 = -0.5
  delta1 = 1
  gamma0 = 0.75
  gamma1 = 0.15
  psi_site_area = 0.25 # fixed effect of site area on occupancy
  
  ## detection
  # citizen science observation process
  mu_p_cs_0 = -3
  sigma_p_cs_species = 1.5
  sigma_p_cs_site = 1.5
  sigma_p_cs_level_three = 1 # variation across level3
  sigma_p_cs_level_four = 0.75 
  p_cs_interval = 0.5
  p_cs_pop_density = 0.5 
  
  # museum record observation process
  mu_p_rc_0 = 0
  sigma_p_rc_species = 0.75
  sigma_p_rc_site = 1
  sigma_p_rc_level_three = 0.75 
  sigma_p_rc_level_four = 0.5 
  p_rc_total_records = 0.25
  
  # correlations
  rho = 0.5
  
  # introduce NAs (missed visits)?
  omega_community = 0.05
  omega_species = 0.005
  ignore_community_misses = TRUE
  
  #sites_missing = 0.5*n_sites 
  #intervals_missing = 2
  #visits_missing = 1
  
  sites_in_range_beta1 = 2
  sites_in_range_beta2 = 2
  
} else {
  
  ## study dimensions
  n_genera = 1 # us 1 unless you want to introduce generic variation
  n_species_per_genera = 50 ## number of species
  n_species = n_genera*n_species_per_genera
  n_level_four = 10
  n_level_three_per_one = 7 # ecoregion3 per ecoregion1
  n_level_three = n_level_four*n_level_three_per_one
  n_sites_per_level_three = 7
  n_sites = n_sites_per_level_three*n_level_three ## number of sites
  
  ## currently the civariance matrix for total records is set up for 3 intervals
  n_intervals = 3 ## number of occupancy intervals # if not three will have to change columns in cormatrix2
  n_visits = 3 ## number of samples per year
  
  ## occupancy
  mu_psi_0 = 0
  sigma_psi_species = 1.5
  sigma_psi_genus = 0  # this won't be used if taxon is "bombus", keep at 0 unless want to introduce
  sigma_psi_site = 1 # variation across level2
  sigma_psi_level_three = 1 # variation across level3
  sigma_psi_level_four = 0.75 # variation across level4
  mu_psi_income = 0.25 # make sure to specify as 0 if using a model without income 
  sigma_psi_income = 0.5 # make sure to specify as 0 if using a model without income 
  mu_psi_natural_habitat = 0.75 
  sigma_psi_natural_habitat = 1
  delta0 = 0
  delta1 = 0
  gamma0 = 0
  gamma1 = 0
  psi_site_area = 0.25 # fixed effect of site area on occupancy
  
  ## detection
  # citizen science observation process
  mu_p_cs_0 = -2.7
  sigma_p_cs_species = 1.25
  sigma_p_cs_site = 1.25
  sigma_p_cs_level_three = 0.75 # variation across level3
  sigma_p_cs_level_four = 0.35 
  p_cs_interval = 0.55
  p_cs_pop_density = 0.5 
  
  # museum record observation process
  mu_p_rc_0 = 0
  sigma_p_rc_species = 0.75
  sigma_p_rc_site = 1
  sigma_p_rc_level_three = 0.75 
  sigma_p_rc_level_four = 0.25 
  p_rc_total_records = 0.25
  
  # correlations
  rho = 0.75
  
  # introduce NAs (missed visits)?
  omega_community = 0.3
  omega_species = 0.025
  ignore_community_misses = TRUE
  
  #sites_missing = 0.5*n_sites 
  #intervals_missing = 2
  #visits_missing = 1
  
  sites_in_range_beta1 = 2
  sites_in_range_beta2 = 2
}


## --------------------------------------------------
### Simulate data
set.seed(1)
my_simulated_data <- simulate_data(taxon,
                                   n_genera,
                                   n_species_per_genera,
                                   n_species,
                                   n_level_four,
                                   n_level_three_per_one,
                                   n_level_three,
                                   n_sites_per_level_three,
                                   n_sites,
                                   n_intervals,
                                   n_visits,
                                   
                                   # ecological process
                                   mu_psi_0,
                                   sigma_psi_species,
                                   sigma_psi_genus,
                                   sigma_psi_site,
                                   sigma_psi_level_three,
                                   sigma_psi_level_four,
                                   sigma_psi_income_level_three,
                                   sigma_psi_income_level_four,
                                   mu_psi_income,
                                   sigma_psi_income,
                                   mu_psi_natural_habitat,
                                   sigma_psi_natural_habitat,
                                   delta0,
                                   delta1,
                                   psi_site_area,
                                   
                                   # citizen science observation process
                                   mu_p_cs_0,
                                   sigma_p_cs_species,
                                   sigma_p_cs_site,
                                   sigma_p_cs_level_three,
                                   sigma_p_cs_level_four,
                                   p_cs_interval,
                                   p_cs_pop_density, 
                                   
                                   # museum record observation process
                                   mu_p_rc_0,
                                   sigma_p_rc_species,
                                   sigma_p_rc_site,
                                   sigma_p_rc_level_three,
                                   sigma_p_rc_level_four,
                                   p_rc_total_records,
                                   
                                   # correlation (detection)
                                   rho,
                                   
                                   # introduce NAs (missed visits)?
                                   omega_community,
                                   omega_species,
                                   ignore_community_misses,
                                   #sites_missing, 
                                   #intervals_missing,
                                   #visits_missing,
                                   
                                   sites_in_range_beta1,
                                   sites_in_range_beta2)

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V_cs <- my_simulated_data$V_cs # detection data
V_rc <- my_simulated_data$V_rc # detection data
ranges <- my_simulated_data$ranges # indicator of whether sampling occurred
V_rc_NA <- my_simulated_data$V_rc_NA # indicator of whether sampling occurred
n_species <- my_simulated_data$n_species # number of species
n_sites <- my_simulated_data$n_sites # number of sites
n_level_three <- my_simulated_data$n_level_three
n_level_four <- my_simulated_data$n_level_four
n_intervals <- my_simulated_data$n_intervals # number of surveys 
n_visits <- my_simulated_data$n_visits

#View(as.data.frame(V_cs[1:10,1:10,,]))
#View(as.data.frame(V_rc[1:10,1:10,,]))
sum(my_simulated_data$V_cs == 1)
sum(my_simulated_data$V_rc == 1)

check_cs <- which(V_cs>ranges)
check_rc <- which(V_rc>V_rc_NA)
 
intervals_raw <- seq(1, n_intervals, by=1)
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
level_three <- my_simulated_data$level_three
level_four <- my_simulated_data$ecoregion_one
level_three_lookup <- my_simulated_data$level_three_lookup
level_four_lookup <- my_simulated_data$ecoregion_one_lookup
species <- seq(1, n_species, by=1)
genus_lookup <- my_simulated_data$genus_lookup

pop_densities <- my_simulated_data$pop_density
site_areas <- my_simulated_data$site_area
avg_income <- my_simulated_data$income
natural_habitat <- my_simulated_data$natural_habitat
rc_total_records <- my_simulated_data$total_records_rc
nativity = my_simulated_data$nativity

species_intercepts <- my_simulated_data$species_intercepts # see how close we get with the estimates
mean(species_intercepts) # should be really close to mu_psi_0
sd(species_intercepts) # should be really close to sigma_psi_species

# the models for each group are slightly different and tracking different params
if(taxon == "syrphidae"){
  
  # data for model
  stan_data <- c("V_cs", "V_rc", 
                 "ranges", 
                 "n_species", "n_sites", "n_intervals", "n_visits", 
                 "intervals", "species", "sites",
                 "n_level_three", 
                 "level_three_lookup", 
                 "n_level_four",
                 "level_four_lookup",
                 "pop_densities", "site_areas", 
                 "nativity",
                 "natural_habitat" 
  ) 
  
  # Parameters monitored
  params <- c(
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "delta0",
    "delta1",
    "gamma0",
    "gamma1",
    "psi_site_area",
    
    "mu_p_cs_0",
    "sigma_p_cs_species",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "sigma_p_cs_level_four",
    "p_cs_interval",
    "p_cs_pop_density", 
    
    "psi_species",
    "psi_natural_habitat",
    
    "psi_level_three", # track city/fine-ecoregion effects
    "psi_level_four", # track broad eco effects
    
    "T_rep_cs",
    "T_obs_cs",
    "P_species_cs",
    
    "mu_psi_natural_habitat_native",
    "mu_psi_natural_habitat_nonnative",
    "mu_psi_natural_habitat_all_species"
  )
  
  parameter_value <- c(mu_psi_0,
                       sigma_psi_species,
                       sigma_psi_site,
                       sigma_psi_level_three,
                       sigma_psi_level_four,
                       delta0,
                       delta1,
                       gamma0,
                       gamma1,
                       psi_site_area,
                       
                       mu_p_cs_0,
                       sigma_p_cs_species,
                       sigma_p_cs_site,
                       sigma_p_cs_level_three,
                       sigma_p_cs_level_four,
                       p_cs_interval,
                       p_cs_pop_density, 
                       
                       NA,
                       NA,
                       
                       NA, # track city/fine-ecoregion effects
                       NA, # track broad eco effects
                       
                       NA,
                       NA,
                       NA,
                       
                       NA,
                       NA,
                       NA)
  
  
  # MCMC settings
  n_iterations <- 1000
  n_thin <- 1
  n_burnin <- 300
  n_chains <- 4
  n_cores <- 4
  #n_cores <- parallel::detectCores()
  delta = 0.97
  
  ## Initial values
  # given the number of parameters, the chains need some decent initial values
  # otherwise sometimes they have a hard time starting to sample
  set.seed(1)
  inits <- lapply(1:n_chains, function(i)
    
    list(
      mu_psi_0 = runif(1, -0.2, 0.2),
      sigma_psi_site = runif(1, 1, 2),
      sigma_psi_level_three = runif(1, 0, 1),
      sigma_psi_level_four = runif(1, 0, 1),
      delta0 = runif(1, -0.5, 0.5),
      delta1 = runif(1, 0, 0.5),
      gamma0 = runif(1, 0.5, 0.75), # must be a positive value!
      gamma1 = runif(1, 0, 0.1), # gamma0+gamma1 inits must be >0!
      psi_site_area = runif(1, -0.5, 0.5),
      
      mu_p_cs_0 = runif(1, -3.25, -2.75),
      sigma_p_cs_site = runif(1, 0, 1),
      sigma_p_cs_level_three = runif(1, 0, 0.2),
      sigma_p_cs_level_four = runif(1, 0, 0.2),
      p_cs_interval = runif(1, 0.5, 0.6),
      p_cs_pop_density = runif(1, 0.4, 0.6)
      
    )
  )
  
  #View(as.data.frame(params))
  #View(as.data.frame(parameter_value))
  
  #params_2 <- params[-1]
  #View(as.data.frame(params_2))
  
  targets <- as.data.frame(cbind(params, parameter_value))
  
} else { # bombus
  
  stan_data <- c("V_cs", "V_rc", 
                 "ranges", "V_rc_NA",
                 "n_species", "n_sites", "n_intervals", "n_visits", 
                 "intervals", "species", "sites",
                 "n_level_three", 
                 "level_three_lookup", 
                 "n_level_four",
                 "level_four_lookup",
                 "pop_densities", "site_areas", "avg_income", 
                 "natural_habitat", "rc_total_records") 
  
  # Parameters monitored
  params <- c("sigma_species_detection",
              "rho",
              
              "mu_psi_0",
              "sigma_psi_species",
              "sigma_psi_site",
              "sigma_psi_level_three",
              "sigma_psi_level_four",
              "mu_psi_income",
              "sigma_psi_income",
              "mu_psi_natural_habitat",
              "sigma_psi_natural_habitat",
              "psi_site_area",
              
              "mu_p_cs_0",
              "sigma_p_cs_site",
              "sigma_p_cs_level_three",
              "sigma_p_cs_level_four",
              "p_cs_interval",
              "p_cs_pop_density", 
              
              "mu_p_rc_0",
              "sigma_p_rc_site",
              "sigma_p_rc_level_three",
              "sigma_p_rc_level_four",
              "p_rc_total_records",
              
              "psi_species",
              "psi_income",
              "psi_natural_habitat",
              
              "psi_site",
              "psi_level_four",
              "psi_level_three", # track city or eco3 effects
              
              "T_rep_cs",
              "T_obs_cs",
              "P_species_cs",
              
              "T_rep_rc",
              "T_obs_rc",
              "P_species_rc"
  )
  
  parameter_value <- c(sigma_p_cs_species,
                       sigma_p_rc_species,
                       rho,
                       
                       mu_psi_0,
                       sigma_psi_species,
                       sigma_psi_site,
                       sigma_psi_level_three,
                       sigma_psi_level_four,
                       mu_psi_income,
                       sigma_psi_income,
                       mu_psi_natural_habitat,
                       sigma_psi_natural_habitat,
                       psi_site_area,
                       
                       mu_p_cs_0,
                       sigma_p_cs_site,
                       sigma_p_cs_level_three,
                       sigma_p_cs_level_four,
                       p_cs_interval,
                       p_cs_pop_density, 
                       
                       mu_p_rc_0,
                       sigma_p_rc_site,
                       sigma_p_rc_level_three,
                       sigma_p_rc_level_four,
                       p_rc_total_records,
                       
                       NA,
                       NA,
                       NA,
                       
                       NA,
                       NA,
                       NA,
                       
                       NA,
                       NA,
                       NA,
                       
                       NA,
                       NA,
                       NA
  )
  
  # MCMC settings
  n_iterations <- 2000
  n_thin <- 1
  n_burnin <- 500
  n_chains <- 4
  n_cores <- parallel::detectCores()
  delta = 0.95
  
  ## Initial values
  # given the number of parameters, the chains need some decent initial values
  # otherwise sometimes they have a hard time starting to sample
  set.seed(1)
  inits <- lapply(1:n_chains, function(i)
    
    list(
      rho = runif(1, 0, 1),
      
      mu_psi_0 = runif(1, -0.5, 0.5),
      sigma_psi_species = runif(1, 0, 1),
      sigma_psi_site = runif(1, 0, 1),
      sigma_psi_level_three = runif(1, 0, 1),
      sigma_psi_level_four = runif(1, 0, 1),
      mu_psi_income = runif(1, -1, 1),
      sigma_psi_income = runif(1, 0, 1),
      mu_psi_natural_habitat = runif(1, -1, 1),
      sigma_psi_natural_habitat = runif(1, 0, 1),
      psi_site_area = runif(1, 0, 0.25),
      
      mu_p_cs_0 = runif(1, -3, -2.75),
      sigma_p_cs_site = runif(1, 0, 0.5),
      sigma_p_cs_level_three = runif(1, 0, 0.5),
      sigma_p_cs_level_four = runif(1, 0, 0.5),
      p_cs_interval = runif(1, 0.5, 0.6),
      p_cs_pop_density = runif(1, 0.4, 0.6),
      
      # start musuem values close to zero
      mu_p_rc_0 = runif(1, -0.5, 0.5),
      sigma_p_rc_site = runif(1, 0, 0.25),
      sigma_p_rc_level_three = runif(1, 0, 0.25),
      sigma_p_rc_level_four = runif(1, 0, 0.25),
      p_rc_total_records = runif(1, -0.5, 0.5)  
      
    )
  )
  
  params_2 <- params[-1]
  #View(as.data.frame(params_2))
  
  targets <- as.data.frame(cbind(c("sigma_p_cs_species", "sigma_p_rc_species", 
                                   params_2), parameter_value))
  
}

#View(as.data.frame(params))
#View(as.data.frame(parameter_value))

View(targets)

## --------------------------------------------------
### Run model

stan_model <-  paste0("./occupancy/models/model_", taxon, ".stan")

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

saveRDS(stan_out_sim,  paste0("./occupancy/simulation/", taxon, "_stan_out_sim.rds"))

stan_out_sim <- readRDS(paste0("./occupancy/simulation/", taxon, "_stan_out_sim.rds"))

# print results
if(taxon == "syrphidae"){
  print(stan_out_sim, digits = 3, pars = c(
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "delta0",
    "delta1",
    "gamma0",
    "gamma1",
    "psi_site_area",
    
    "mu_p_cs_0",
    "sigma_p_cs_species",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "sigma_p_cs_level_four",
    "p_cs_interval",
    "p_cs_pop_density"
  ))
} else {
  print(stan_out_sim, digits = 3, pars = c(
    "rho", 
    "sigma_species_detection[1]",
    "sigma_species_detection[2]", 
    "mu_psi_0",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "mu_psi_income",
    "sigma_psi_income",
    "mu_psi_natural_habitat",
    "sigma_psi_natural_habitat",
    "psi_site_area",
    
    "mu_p_cs_0",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "sigma_p_cs_level_four",
    "p_cs_interval",
    "p_cs_pop_density", 
    
    "mu_p_rc_0",
    "sigma_p_rc_site",
    "sigma_p_rc_level_three",
    "sigma_p_rc_level_four",
    "p_rc_total_records"
  ))
}

View(targets)

# print some specific parameter if desired
print(stan_out_sim, digits = 3, pars=
        c("species_intercepts"))

View(as.data.frame(fit_summary <- rstan::summary(stan_out_sim)))

## --------------------------------------------------
### Simple diagnostic plots

# traceplots
if(taxon == "syrphidae"){
  traceplot(stan_out_sim, pars = c( # occupancy
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "delta0",
    "delta1",
    "gamma0",
    "gamma1",
  ))
  traceplot(stan_out_sim, pars = c(
    "mu_p_cs_0",
    "p_cs_interval",
    "p_cs_pop_density",
    "sigma_p_cs_species",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "sigma_p_cs_level_four"
  ))
} else{
  traceplot(stan_out_sim, pars = c( # occupancy
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "mu_psi_natural_habitat",
    "sigma_psi_natural_habitat",
    "mu_psi_income",
    "sigma_psi_income"
  ))
  traceplot(stan_out_sim, pars = c( # detection
    "mu_p_cs_0",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "sigma_p_cs_level_four",
    "p_cs_interval",
    "p_cs_pop_density",
    "mu_p_rc_0",
    "sigma_p_rc_site",
    "sigma_p_rc_level_three",
    "sigma_p_rc_level_four",
    "p_rc_total_records"
  ))
  traceplot(stan_out_sim, pars = c( # species detection
    "rho",
    "sigma_species_detection"
  ))
}

# pairs plot (for a sample of parameters)
if(taxon == "syrphidae"){
  pairs(stan_out_sim, pars = c(
    "mu_psi_0",
    "mu_p_cs_0",
    
    "sigma_psi_species",
    "sigma_psi_site",
    "delta1",
    "sigma_p_cs_species",
    "sigma_p_cs_site"
  ))
} else {
  pairs(stan_out_sim, pars = c(
    "mu_psi_0",
    "mu_p_cs_0",
    "mu_p_rc_0",
    "mu_psi_income",
    "mu_psi_natural_habitat",
    
    "sigma_psi_site",
    "sigma_p_cs_site",
    "sigma_p_rc_site"
  ))
}

## --------------------------------------------------
### Plot parameter estimates and targets



## --------------------------------------------------
### PPC

# print rep and obs
print(stan_out_sim, digits = 3, pars=
        c("T_rep_cs", "T_obs_cs", "P_species_cs"))
print(stan_out_sim, digits = 3, pars=
        c("T_rep_rc", "T_obs_rc", "P_species_rc"))

# as data frame
list_of_draws <- as.data.frame(stan_out_sim)
list_of_draws <- list_of_draws[(n_burnin+1):n_iterations,]

# Citizen Science
# P-values
m <- n_iterations - n_burnin
P_average_cs = vector(length = n_species)

for(i in 1:n_species){
  P_average_cs[i] = sum(list_of_draws[,88+i])/m
}

print(P_average_cs)

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
P_average_rc = vector(length = n_species)

for(i in 1:n_species){
  P_average_rc[i] = sum(list_of_draws[,193+i])/m
}

print(P_average_rc)

# Evaluation of fit
plot(list_of_draws[,124], list_of_draws[,159], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 100),
     xlim = c(0, 100))

abline(0, 1, lwd = 2, col = "black")

## --------------------------------------------------
### Check how many species FAIL the P-value (above .95 or below 0.05)