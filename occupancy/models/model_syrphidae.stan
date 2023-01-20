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
  
  int<lower=0> ranges[n_species, n_sites, n_intervals, n_visits];  // NA indicator where 1 == site is in range, 0 == not in range
  int<lower=0> V_museum_NA[n_species, n_sites, n_intervals, n_visits];  // indicator where 1 == sampled, 0 == missing data
  
  vector[n_sites] site_areas; // (scaled) spatial area extent of each site
  vector[n_sites] pop_densities; // (scaled) population density of each site
  // vector[n_sites] open_developed; // (scaled) developed open surface cover of each site
  // vector[n_sites] developed_med_high; // (scaled) high develop cover of each site
  vector[n_sites] developed_med_high; // (scaled) UNdeveloped open surface cover of each site
  
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
  
  // random slope for species specific population density effects on occupancy
  vector[n_species] psi_developed_med_high; // vector of species specific slope estimates
  real mu_psi_developed_med_high; // community mean of species specific slopes
  real<lower=0> sigma_psi_developed_med_high; // variance in species slopes
  
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
  
  real logit_psi[n_species, n_sites, n_intervals];  // odds of occurrence
  real logit_p_citsci[n_species, n_sites, n_intervals]; // odds of detection by cit science
  real logit_p_museum[n_species, n_sites, n_intervals]; // odds of detection by museum
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals  
          
          logit_psi[i,j,k] = // the inverse of the log odds of occurrence is equal to..
            mu_psi_0 + // a baseline intercept
            psi_species[species[i]] + // a species specific intercept
            psi_site[sites[j]] + // a site specific intercept
            //psi_open_developed[species[i]]*open_developed[k] + // an effect
            //psi_developed_med_high[species[i]]*developed_med_high[j] + // an effect of pop density on occurrence
            psi_developed_med_high[species[i]]*developed_med_high[j] + // an effect 
            psi_site_area*site_areas[j] // an effect of spatial area of the site on occurrence
            ; // end psi[i,j,k]
            
      } // end loop across all intervals
    } // end loop across all sites
  }  // end loop across all species
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
          logit_p_citsci[i,j,k] = // the inverse of the log odds of detection is equal to..
            mu_p_citsci_0 + // a baseline intercept
            p_citsci_species[species[i]] + // a species specific intercept
            p_citsci_site[sites[j]] + // a spatially specific intercept
            p_citsci_interval*intervals[k] + // an overall effect of time on detection
            p_citsci_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_citsci[i,j,k]
           
          logit_p_museum[i,j,k] = // the inverse of the log odds of detection is equal to..
            mu_p_museum_0 + // a baseline intercept
            p_museum_species[species[i]] + // a species specific intercept
            p_museum_site[sites[j]] + // a spatially specific intercept
            p_museum_interval*intervals[k] + // an overall effect of time on detection
            p_museum_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_museum[i,j,k]
           
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
             
  
} // end transformed parameters


model {
  
  // PRIORS
  
  // Occupancy (Ecological Process)
  mu_psi_0 ~ normal(0, 2); // global intercept for occupancy rate
  
  psi_species ~ normal(0, sigma_psi_species); 
  // occupancy intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_psi_species ~ normal(0, 0.5); //informative prior
  
  psi_site ~ normal(0, sigma_psi_site); 
  // occupancy intercept for each site drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_psi_site ~ normal(0, 0.5); // informative prior
  
  // the change in time prior is set to be informative 
  // I suspect more of the variation in time comes from detection rate changes 
  // and so I'd like the model to assume that (given that I don't currently have 
  // a huge amount of detection data to let the model make that inference itself)
  //psi_open_developed ~ normal(mu_psi_open_developed, sigma_psi_open_developed);
  // occupancy slope (temporal effect on occupancy) for each species drawn from the 
  // community distribution (variance defined by sigma), centered at mu_psi_interval. 
  // centering on mu (rather than 0) allows us to estimate the average effect of
  // the management on abundance across all species.
  //mu_psi_open_developed ~ normal(0, 2); // community mean
  //sigma_psi_open_developed ~ normal(0, 1); // community variance
  
  //psi_developed_med_high ~ normal(mu_psi_developed_med_high, sigma_psi_developed_med_high);
  // occupancy slope (population density effect on occupancy) for each species drawn from the 
  // community distribution (variance defined by sigma), centered at mu_psi_interval. 
  // centering on mu (rather than 0) allows us to estimate the average effect of
  // the management on abundance across all species.
  //mu_psi_developed_med_high ~ normal(0, 2); // community mean
  //sigma_psi_developed_med_high ~ normal(0, 1); // community variance
  
  psi_developed_med_high ~ normal(mu_psi_developed_med_high, sigma_psi_developed_med_high);
  // occupancy slope (population density effect on occupancy) for each species drawn from the 
  // community distribution (variance defined by sigma), centered at mu_psi_interval. 
  // centering on mu (rather than 0) allows us to estimate the average effect of
  // the management on abundance across all species.
  mu_psi_developed_med_high ~ normal(0, 1); // community mean
  sigma_psi_developed_med_high ~ normal(0, 0.5); // community variance
  
  psi_site_area ~ normal(0, 0.5); // effect of site area on occupancy
  
  // Detection (Observation Process)
  
  // citizen science records
  
  mu_p_citsci_0 ~ normal(0, 2); // global intercept for detection

  p_citsci_species ~ normal(0, sigma_p_citsci_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_citsci_species ~ cauchy(0, 0.5);
  
  // should redefine p_site so that it is spatially AND temporally heterogenous 
  p_citsci_site ~ normal(0, sigma_p_citsci_site);
  // detection intercept for each site drawn from the spatially heterogenous
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_citsci_site ~ normal(0, 0.5); // spatial variance
  
  p_citsci_interval ~ normal(0, 2); // temporal effect on detection probability
  
  p_citsci_pop_density ~ normal(0, 2); // population effect on detection probability
  
  // museum records
  
  mu_p_museum_0 ~ normal(0, 2); // global intercept for detection
  
  p_museum_species ~ normal(0, sigma_p_museum_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_museum_species ~ normal(0, 0.5);
  
  // should redefine p_site so that it is spatially AND temporally heterogenous 
  p_museum_site ~ normal(0, sigma_p_museum_site);
  // detection intercept for each site drawn from the spatially heterogenous
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_museum_site ~ normal(0, 0.5); // spatial variance
  
  p_museum_interval ~ normal(0, 2); // temporal effect on detection probability
  
  p_museum_pop_density ~ normal(0, 2); // population effect on detection probability

  
  // LIKELIHOOD
  
  // Stan can sample the mean and sd of parameters by summing out the
  // parameter (marginalizing) across likelihood statements
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
        // If the site is in the range of a species, then evaluate lp, otherwise do not (treat as NA).
        if(sum(ranges[i,j,k]) > 0){ // The sum of the NA vector will be == 0 if site is not in range
        
          // if species is detected at the specific site*interval at least once
          // by citizen science efforts OR museum records
          // then the species occurs there. lp_observed calculates
          // the probability density that species occurs given psi, plus the 
          // probability density that we did/did not observe it on each visit l in 1:nvisit
          if(sum(V_citsci[i, j, k, 1:n_visits]) > 0 || sum(V_museum[i, j, k, 1:n_visits]) > 0) {
            
             // lp_observed:
             target += log_inv_logit(logit_psi[i,j,k]) +
                      binomial_logit_lpmf(sum(V_citsci[i,j,k,1:n_visits]) | n_visits, logit_p_citsci[i,j,k]) + 
                      // sum(V_museum_NA[i,j,k,1:n_visits]) below tells us how many sampling 
                      // events actually occurred for museum records
                      binomial_logit_lpmf(sum(V_museum[i,j,k,1:n_visits]) | sum(V_museum_NA[i,j,k,1:n_visits]), logit_p_museum[i,j,k]);
                          
          // else the species was never detected at the site*interval
          // lp_unobserved sums the probability density of:
          // 1) species occupies the site*interval but was not detected on each visit, and
          // 2) the species does not occupy the site*interval
          } else {
            
            // lp_unobserved
            target += log_sum_exp(log_inv_logit(logit_psi[i,j,k]) +
                    binomial_logit_lpmf(0 | 
                      n_visits, logit_p_citsci[i,j,k]) +
                    // sum(V_museum_NA[i,j,k,1:n_visits]) below tells us how many sampling 
                    // events actually occurred for museum records
                    binomial_logit_lpmf(0 | 
                      sum(V_museum_NA[i,j,k,1:n_visits]), logit_p_museum[i,j,k]),
                    
                    log1m_inv_logit(logit_psi[i,j,k])); 
            
          } // end if/else ever observed
        
        } // end if/ in range
          
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
  
} // end model

generated quantities{
  
  int Z[n_species, n_sites, n_intervals];
  
  int z_rep[n_species, n_sites, n_intervals];
  int y_rep_citsci[n_species, n_sites, n_intervals, n_visits]; // repd detections
  int y_rep_museum[n_species, n_sites, n_intervals, n_visits]; // repd detections
  
  real eval_citsci[n_species,n_sites,n_intervals,n_visits]; // expected values
  real eval_museum[n_species,n_sites,n_intervals,n_visits]; // expected values
  
  real T_rep_citsci[n_species]; // Tukey-Freeman distance from eval (species bin)
  real T_obs_citsci[n_species]; // Tukey-Freeman distance from eval (species bin)
  real T_rep_museum[n_species]; // Tukey-Freeman distance from eval (species bin)
  real T_obs_museum[n_species]; // Tukey-Freeman distance from eval (species bin)
  
  real P_species_citsci[n_species]; // P-value by species
  real P_species_museum[n_species]; // P-value by species
  
  // Initialize T_rep and T_obs and P-values
  for(i in 1:n_species){
    
    T_rep_citsci[i] = 0;
    T_obs_citsci[i] = 0;
    T_rep_museum[i] = 0;
    T_obs_museum[i] = 0;
    
    P_species_citsci[i] = 0;
    P_species_museum[i] = 0;
    
  }
      
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
      
        if(sum(ranges[i,j,k]) > 0){ // The sum of the NA vector will be == 0 if site is not in range
          
          // if occupancy state is certain then the expected occupancy is 1
          if(sum(V_citsci[i, j, k, 1:n_visits]) > 0 || sum(V_museum[i, j, k, 1:n_visits]) > 0) {
          
            Z[i,j,k] = 1;
          
          // else the site could be occupied or not
          } else {
            
            // occupancy but never observed by either dataset
            real ulo = inv_logit(logit_psi[i,j,k]) * 
              ((1 - inv_logit(logit_p_citsci[i,j,k]))^n_visits + 
              (1 - inv_logit(logit_p_museum[i,j,k]))^n_visits);
            // non-occupancy
            real uln = (1 - inv_logit(logit_psi[i,j,k]));
            
            // outcome of occupancy given the likelihood associated with both possibilities
            Z[i,j,k] = bernoulli_rng(ulo / (ulo + uln));
            
          } // end else uncertain occupancy state
        
        // else the site is not in range for the species
        } else {
          
          // and by definition the site is unoccupied
          Z[i,j,k] = 0;
          
        } // end else the site is not in the species range
        
      } // end loop across intervals
    } // end loop across sites
  } // end loop across species
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        for(l in 1:n_visits){
          
          // expected detections
          eval_citsci[i,j,k,l] = Z[i,j,k] * 
            bernoulli_logit_rng(logit_p_citsci[i,j,k]);
          eval_museum[i,j,k,l] = Z[i,j,k] * 
            bernoulli_logit_rng(logit_p_museum[i,j,k]);
          
          // occupancy in replicated data
          // should evaluate to zero if the site is not in range
          z_rep[i,j,k] = bernoulli_logit_rng(logit_psi[i,j,k]); 
          if(sum(ranges[i,j,k]) > 0){ // The sum of the NA vector will be == 0 if site is not in range
            z_rep[i,j,k] = z_rep[i,j,k];
          } else {
            z_rep[i,j,k] = 0;
          }
          
          // detections in replicated data
          y_rep_citsci[i,j,k,l] = z_rep[i,j,k] * bernoulli_logit_rng(logit_p_citsci[i,j,k]);
          y_rep_museum[i,j,k,l] = z_rep[i,j,k] * bernoulli_logit_rng(logit_p_museum[i,j,k]);
          
          // Compute fit statistic (Tukey-Freeman) for replicate data
          // Citizen science records
          // Binned by species
          T_rep_citsci[i] = T_rep_citsci[i] + (sqrt(y_rep_citsci[i,j,k,l]) - 
            sqrt(eval_citsci[i,j,k,l]))^2;
          // Compute fit statistic (Tukey-Freeman) for real data
          // Binned by species
          T_obs_citsci[i] = T_obs_citsci[i] + (sqrt(V_citsci[i,j,k,l]) - 
            sqrt(eval_citsci[i,j,k,l]))^2;
          
          // Compute fit statistic (Tukey-Freeman) for replicate data
          // Binned by species
          // Museum records
          if(V_museum_NA[i,j,k,l] == 0){
            
            T_rep_museum[i] = T_rep_museum[i] + 0;
            // Compute fit statistic (Tukey-Freeman) for real data
            // Binned by species
            T_obs_museum[i] = T_obs_museum[i] + 0;
          
          } else{
            
            T_rep_museum[i] = T_rep_museum[i] + (sqrt(y_rep_museum[i,j,k,l]) - 
              sqrt(eval_museum[i,j,k,l]))^2;
            // Compute fit statistic (Tukey-Freeman) for real data
            // Binned by species
            T_obs_museum[i] = T_obs_museum[i] + (sqrt(V_museum[i,j,k,l]) - 
              sqrt(eval_museum[i,j,k,l]))^2;
            
          }
          
        } // end loop across visits
      } // end loop across intervals
    } // end loop across sites
  } // end loop across species
  
  // bin by species
  for(i in 1:n_species) { // loop across all species
    
    // if the discrepancy is lower for the real data for the species
    // versus the replicated data
    if(T_obs_citsci[i] < T_rep_citsci[i]){
      
      // then increase species P by 1      
      P_species_citsci[i] = P_species_citsci[i] + 1;
      // the ppc will involve averaging P across the number of post-burnin iterations
            
    }
    
    // if the discrepancy is lower for the real data for the species
    // versus the replicated data
    if(T_obs_museum[i] < T_rep_museum[i]){
      
      // then increase species P by 1      
      P_species_museum[i] = P_species_museum[i] + 1;
      // the ppc will involve averaging P across the number of post-burnin iterations
            
    }
    
  }

}
