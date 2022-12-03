// multi-species occupancy model for GBIF occurrence data
// jcu, started nov 21, 2022.
// builds on model0 by introducing integrated model structure where
// citizen science data and gbif data may have their own observation processes
// and also allows for missing (NA) data


data {
  
  int<lower=1> n_species;  // observed species
  int<lower=1> species[n_species]; // vector of species
  
  int<lower=1> n_sites;  // sites within region
  int<lower=1> sites[n_sites];  // vector of sites
  
  int<lower=1> n_intervals;  // intervals during which sites are visited
  
  real intervals[n_intervals]; // vector of intervals (used as covariate data for 
                                // species specific effect of occupancy interval (time) on occupancy)
                                // needs to begin with intervals[1] = 0, i.e., 
                                // there is no temporal addition in the first interval
  
  int<lower=1> n_visits; // visits within intervals
  
  int<lower=0> V_citsci[n_species, n_sites, n_intervals, n_visits];  // visits l when species i was detected at site j on interval k
  int<lower=0> V_museum[n_species, n_sites, n_intervals, n_visits];  // visits l when species i was detected at site j on interval k
  
  int<lower=0> V_citsci_NA[n_species, n_sites, n_intervals, n_visits];  // indicator where 1 == sampled, 0 == missing data
  int<lower=0> V_museum_NA[n_species, n_sites, n_intervals, n_visits];  // indicator where 1 == sampled, 0 == missing data
  
  vector[n_sites] pop_densities; // population density of each site
  vector[n_sites] site_areas; // spatial area extent of each site
  
} // end data


parameters {
  
  // OCCUPANCY
  real mu_psi_0; // global intercept for occupancy
  
  // species specific intercept allows some species to occur at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] psi_species; // species specific intercept for occupancy
  real<lower=0> sigma_psi_species; // variance in species intercepts
  
  // site specific intercept allows some sites to be occupied at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all sites.
  vector[n_sites] psi_site; // site specific intercept for occupancy
  real<lower=0> sigma_psi_site; // variance in site intercepts
  
  // random slope for species specific temporal effects on occupancy
  vector[n_species] psi_interval; // vector of species specific slope estimates
  real mu_psi_interval; // community mean of species specific slopes
  real<lower=0> sigma_psi_interval; // variance in species slopes
  
  // random slope for species specific population density effects on occupancy
  vector[n_species] psi_pop_density; // vector of species specific slope estimates
  real mu_psi_pop_density; // community mean of species specific slopes
  real<lower=0> sigma_psi_pop_density; // variance in species slopes
  
  // effect of site are on occupancy
  real psi_site_area;
  
  // DETECTION
  
  // citizen science observation process
  real mu_p_citsci_0; // global detection intercept for citizen science records
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] p_citsci_species; // species specific intercept for detection
  real<lower=0> sigma_p_citsci_species; // variance in species intercepts
  
  // random slope for site specific temporal effects on occupancy
  vector[n_sites] p_citsci_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_citsci_site; // variance in site slopes
  
  real p_citsci_interval; // fixed temporal effect on detection probability
  real p_citsci_pop_density; // fixed effect of population on detection probability
  
  // museum records observation process
  real mu_p_museum_0; // global detection intercept for citizen science records
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] p_museum_species; // species specific intercept for detection
  real<lower=0> sigma_p_museum_species; // variance in species intercepts
  
  // random slope for site specific temporal effects on occupancy
  vector[n_sites] p_museum_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_museum_site; // variance in site slopes
  
  real p_museum_interval; // fixed temporal effect on detection probability
  real p_museum_pop_density; // fixed effect of population on detection probability
  
} // end parameters


transformed parameters {
  
  real psi[n_species, n_sites, n_intervals];  // odds of occurrence
  real p_citsci[n_species, n_sites, n_intervals]; // odds of detection by cit science
  real p_museum[n_species, n_sites, n_intervals]; // odds of detection by museum
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals  
          
          psi[i,j,k] = inv_logit( // the inverse of the log odds of occurrence is equal to..
            mu_psi_0 + // a baseline intercept
            psi_species[species[i]] + // a species specific intercept
            psi_site[sites[j]] + // a site specific intercept
            psi_interval[species[i]]*intervals[k] + // a species specific temporal effect
            psi_pop_density[species[i]]*pop_densities[j] + // an effect of pop density on occurrence
            psi_site_area*site_areas[j] // an effect of spatial area of the site on occurrence
            ); // end psi[i,j,k]
            
      } // end loop across all intervals
    } // end loop across all sites
  }  // end loop across all species
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
          p_citsci[i,j,k] = inv_logit( // the inverse of the log odds of detection is equal to..
            mu_p_citsci_0 + // a baseline intercept
            p_citsci_species[species[i]] + // a species specific intercept
            p_citsci_site[sites[j]] + // a spatially specific intercept
            p_citsci_interval*intervals[k] + // an overall effect of time on detection
            p_citsci_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ); // end p_citsci[i,j,k]
           
          p_museum[i,j,k] = inv_logit( // the inverse of the log odds of detection is equal to..
            mu_p_museum_0 + // a baseline intercept
            p_museum_species[species[i]] + // a species specific intercept
            p_museum_site[sites[j]] + // a spatially specific intercept
            p_museum_interval*intervals[k] + // an overall effect of time on detection
            p_museum_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ); // end p_museum[i,j,k]
           
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
             
  
} // end transformed parameters


model {
  
  // PRIORS
  
  // Occupancy (Ecological Process)
  mu_psi_0 ~ cauchy(0, 2.5); // global intercept for occupancy rate
  
  psi_species ~ normal(0, sigma_psi_species); 
  // occupancy intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_psi_species ~ cauchy(0, 1); //informative prior
  
  psi_site ~ normal(0, sigma_psi_site); 
  // occupancy intercept for each site drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_psi_site ~ cauchy(0, 0.5); // informative prior
  
  // the change in time prior is set to be informative 
  // I suspect more of the variation in time comes from detection rate changes 
  // and so I'd like the model to assume that (given that I don't currently have 
  // a huge amount of detection data to let the model make that inference itself)
  psi_interval ~ normal(mu_psi_interval, sigma_psi_interval);
  // occupancy slope (temporal effect on occupancy) for each species drawn from the 
  // community distribution (variance defined by sigma), centered at mu_psi_interval. 
  // centering on mu (rather than 0) allows us to estimate the average effect of
  // the management on abundance across all species.
  mu_psi_interval ~ cauchy(0, 0.25); // community mean
  sigma_psi_interval ~ cauchy(0, 0.1); // community variance
  
  psi_pop_density ~ normal(mu_psi_pop_density, sigma_psi_pop_density);
  // occupancy slope (population density effect on occupancy) for each species drawn from the 
  // community distribution (variance defined by sigma), centered at mu_psi_interval. 
  // centering on mu (rather than 0) allows us to estimate the average effect of
  // the management on abundance across all species.
  mu_psi_pop_density ~ cauchy(0, 2.5); // community mean
  mu_psi_pop_density ~ cauchy(0, 1); // community variance
  
  psi_site_area ~ cauchy(0, 2.5); // effect of site area on occupancy
  
  // Detection (Observation Process)
  
  // citizen science records
  
  mu_p_citsci_0 ~ cauchy(0, 2.5); // global intercept for detection

  p_citsci_species ~ normal(0, sigma_p_citsci_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_citsci_species ~ cauchy(0, 1);
  
  // should redefine p_site so that it is spatially AND temporally heterogenous 
  p_citsci_site ~ normal(0, sigma_p_citsci_site);
  // detection intercept for each site drawn from the spatially heterogenous
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_citsci_site ~ cauchy(0, 1); // spatial variance
  
  p_citsci_interval ~ cauchy(0, 2.5); // temporal effect on detection probability
  
  p_citsci_pop_density ~ cauchy(0, 2.5); // population effect on detection probability
  
  // museum records
  
  mu_p_museum_0 ~ cauchy(0, 2.5); // global intercept for detection
  
  p_museum_species ~ normal(0, sigma_p_museum_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_museum_species ~ cauchy(0, 1);
  
  // should redefine p_site so that it is spatially AND temporally heterogenous 
  p_museum_site ~ normal(0, sigma_p_museum_site);
  // detection intercept for each site drawn from the spatially heterogenous
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_museum_site ~ cauchy(0, 1); // spatial variance
  
  p_museum_interval ~ cauchy(0, 2.5); // temporal effect on detection probability
  
  p_museum_pop_density ~ cauchy(0, 2.5); // population effect on detection probability

  
  // LIKELIHOOD
  
  // Stan can sample the mean and sd of parameters by summing out the
  // parameter (marginalizing) across likelihood statements
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
          
          // if species is detected at the specific site*interval at least once
          // by citizen science efforts OR museum records
          // then the species occurs there. lp_observed calculates
          // the probability density that species occurs given psi, plus the 
          // probability density that we did/did not observe it on each visit l in 1:nvisit
          if(sum(V_citsci[i, j, k, 1:n_visits]) > 0 || sum(V_museum[i, j, k, 1:n_visits]) > 0) {
            
             // lp_observed:
             target += log(psi[i,j,k]) +
                      binomial_lpmf(sum(V_citsci[i,j,k,1:n_visits]) | sum(V_citsci_NA[i,j,k,1:n_visits]), p_citsci[i,j,k]) + 
                      // sum(V_museum_NA[i,j,k,1:n_visits]) below tells us how many sampling 
                      // events actually occurred for museum records
                      binomial_lpmf(sum(V_museum[i,j,k,1:n_visits]) | sum(V_museum_NA[i,j,k,1:n_visits]), p_museum[i,j,k]);
                          
          // else the species was never detected at the site*interval
          // lp_unobserved sums the probability density of:
          // 1) species occupies the site*interval but was not detected on each visit, and
          // 2) the species does not occupy the site*interval
          } else {
            
            // lp_unobserved
            target += log_sum_exp(log(psi[i,j,k]) +
                    binomial_lpmf(0 | sum(V_citsci_NA[i,j,k,1:n_visits]), p_citsci[i,j,k]) +
                    // sum(V_museum_NA[i,j,k,1:n_visits]) below tells us how many sampling 
                    // events actually occurred for museum records
                    binomial_lpmf(0 | sum(V_museum_NA[i,j,k,1:n_visits]), p_museum[i,j,k]),
                    log1m(psi[i,j,k])); 
            
          } // end if/else
          
          
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
  
} // end model

