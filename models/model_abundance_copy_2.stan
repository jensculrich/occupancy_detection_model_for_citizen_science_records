// multi-species integrated abundance-occupancy model for NHC data
// jcu, started nov 21, 2022.

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
  
  int<lower=0> K[n_species, n_sites, n_intervals]; // Upper bound of population size
  
  vector[n_sites] site_areas; // spatial area extent of each site
  
} // end data

transformed data {
  int<lower=0> max_y[n_species, n_sites, n_intervals];

  for (i in 1:n_species) {
    for(j in 1:n_sites){
      for(k in 1:n_intervals){
        
        // Set the floor of the latent state search to be at least as many as the most 
        // that we observed of a species at a site in a time interval (by cit sci records)
        max_y[i,j,k] = max(V_citsci[i,j,k]);
        
        // We only search abundance if it is it must be greater than 1, i.e.,
        // was detected by a museum but not by citizen science
        // or we search in a hypothetical situation where the site is suitable
        // we did not detect and records, but we consider the probabiility that 1:K
        // individuals exist and went undetected.
        // Therefore, replace any max_y == 0 with max_y == 1 
        // to be the floor of the latent state search.
        if(max_y[i,j,k] == 0){
            max_y[i,j,k] = 1;
          } else {
            max_y[i,j,k] = max_y[i,j,k];
          }
      
      } // end loop across intervals
    } // end loop across sites
  } // end loop across species

} // end transformed data

parameters {
  
  // ABUNDANCE
  
  //real<lower=0,upper=1> omega;
  //real<lower=0,upper=1> omega;
  real gamma_0; // occupancy intercept
  real gamma_1; // relationship between abundance and occupancy 
  
  real<lower=0> phi; // abundance overdispersion parameter
  
  real mu_eta_0; // global intercept for occupancy
  
  real eta_site_area; // effect of site are on occupancy
  
  // species specific intercept allows some species to occur at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] eta_species; // species specific intercept for occupancy
  real<lower=0> sigma_eta_species; // variance in species intercepts
  
  // DETECTION
  
  // citizen science observation process
  real mu_p_citsci_0; // global detection intercept for citizen science records
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] p_citsci_species; // species specific intercept for detection
  real<lower=0> sigma_p_citsci_species; // variance in species intercepts
  
  // museum records observation process
  real mu_p_museum_0; // global detection intercept for citizen science records
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] p_museum_species; // species specific intercept for detection
  real<lower=0> sigma_p_museum_species; // variance in species intercepts
  
} // end parameters


transformed parameters {
  
  real log_eta[n_species, n_sites, n_intervals]; // mean of abundance process
  real logit_p_citsci[n_species, n_sites, n_intervals]; // odds of detection by cit science
  real logit_p_museum[n_species, n_sites, n_intervals]; // odds of detection by museum
  
  real omega[n_species, n_sites, n_intervals]; // availability
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals  
          
          log_eta[i,j,k] = // 
            mu_eta_0 + // a baseline intercept
            eta_species[species[i]] + // a species-specific intercept
            //psi_site[sites[j]] + // a site specific intercept
            //psi_interval[species[i]]*intervals[k] + // a species specific temporal effect
            //psi_pop_density[species[i]]*pop_densities[j] + // an effect of pop density on occurrence
            eta_site_area*site_areas[j] // an effect of spatial area of the site on occurrence
            ; // end lambda[i,j,k]
          
          // Smith et al. 2012 Ecology trick for incorporating the abundance-occupancy 
          // relationship into a zero-inflated abundance model  
          // availability is predicted by abundance
          // omega is logit scaled
          omega[i,j,k] = gamma_0 + gamma_1 * log_eta[i,j,k];
            
      } // end loop across all intervals
    } // end loop across all sites
  }  // end loop across all species
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
          logit_p_citsci[i,j,k] =  // logit scaled individual-level detection rate
            mu_p_citsci_0 + // a baseline intercept
            p_citsci_species[species[i]] //+ // a species specific intercept
            //p_citsci_site[sites[j]] + // a spatially specific intercept
            //p_citsci_interval*intervals[k] + // an overall effect of time on detection
            //p_citsci_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_citsci[i,j,k]
           
          logit_p_museum[i,j,k] = // logit scaled species-level detection rate
            mu_p_museum_0 + // a baseline intercept
            p_museum_species[species[i]] //+ // a species specific intercept
            //p_museum_site[sites[j]] + // a spatially specific intercept
            //p_museum_interval*intervals[k] + // an overall effect of time on detection
            //p_museum_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_museum[i,j,k]
           
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
  
} // end transformed parameters


model {
  
  // PRIORS
  
  // Abundance (Ecological Process)
  
  gamma_0 ~ normal(0, 0.5);
  gamma_1 ~ normal(0, 0.5);
  
  phi ~ cauchy(0, 2.5); // abundance overdispersion scale parameter
  
  mu_eta_0 ~ cauchy(0, 2.5); // global intercept for abundance rate
  
  eta_species ~ normal(0, sigma_eta_species); 
  // occupancy intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_eta_species ~ normal(0, 1); // weakly informative prior
  
  eta_site_area ~ cauchy(0, 2.5); // effect of site area on abundance rate
  
  // Detection (Observation Process)
  
  // citizen science records
  
  mu_p_citsci_0 ~ cauchy(0, 2.5); // global intercept for (citizen science) detection
  
  p_citsci_species ~ normal(0, sigma_p_citsci_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_citsci_species ~ cauchy(0, 1);

  // museum records
  mu_p_museum_0 ~ cauchy(0, 2.5); // global intercept for (museum) detection
  
  p_museum_species ~ normal(0, sigma_p_museum_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_museum_species ~ cauchy(0, 1);
  
  // LIKELIHOOD
  
  // Stan can sample the mean and sd of parameters by summing out the
  // parameter (marginalizing) across likelihood statements
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
          
        // If the site is in the range of a species, then evaluate lp, otherwise do not (treat as NA).
        if(sum(ranges[i,j,k]) > 0){ // The sum of the NA vector will be == 0 if site is not in range
        
        // If a species was detected at least once by either data set, Nijk > 0;
        // Evaluate sum probabilty of an abundance generating term (lambda_ijk),
        // an individual-level detection rate by citizen science data collections (p_citsci_ijk),
        // a species-level detection rate by museum data collections (p_museum_ijk),
        // and an occupancy-abundance relationship (omega_ijk)
        if(sum(V_citsci[i,j,k]) > 0 || sum(V_museum[i,j,k]) > 0) {
          
          vector[K[i,j,k] - max_y[i,j,k] + 1] lp; // lp vector of length of possible abundances 
            // (abundance is at least as big as the max observed count but may range up to K)
          
          // for each possible abundance in max observed abundance through K:
          for(abundance in 1:(K[i,j,k] - max_y[i,j,k] + 1)){ 
          
            // lp of abundance given ecological model and observational model
            lp[abundance] = 
              neg_binomial_2_log_lpmf( // generation of abundance given count distribution
                max_y[i,j,k] + abundance - 1 | log_eta[i,j,k], phi) + 
              // with abundance detections vectorized over n visits..
              binomial_logit_lpmf( // individual-level detection, citizen science
                V_citsci[i,j,k] | max_y[i,j,k] + abundance - 1, logit_p_citsci[i,j,k]) +
              binomial_logit_lpmf( // binary, species-level detecion, museums
                // (given the number of community sampling events that occurred)
                sum(V_museum[i,j,k]) | sum(V_museum_NA[i,j,k]), logit_p_museum[i,j,k]) +
              // plus outcome of site being available, given the 
              // abundance-dependent probability of suitability
              bernoulli_logit_lpmf(1 | omega[i,j,k]); 
          
          }
                
          target += log_sum_exp(lp);
        
        } else { // else was never detected and the site may or may not be available
          
          real lp[2];
          
          // probability present at an available site with
          // some unknown latent abundance state; but never observed.
          // In this formulation the latent abundance could include 0,
          // potentially causing underestimates in species-level detection ability?
          for(abundance in 1:(K[i,j,k])){
            
            // outcome of site being unavailable for occupancy, given the 
            // abundance-dependent probability of suitability
            lp[1] = bernoulli_logit_lpmf(0 | omega[i,j,k]) ; // site not available for species in interval
            
            // outcome of site being available for occupancy, given the 
            // abundance-dependent probability of suitability
            lp[2] = bernoulli_logit_lpmf(1 | omega[i,j,k]); // available but not observed
            
            // start at abundance - 1 so that the search has a floor of 0
            lp[2] = lp[2] +
             neg_binomial_2_log_lpmf( // generation of abundance given count distribution
                abundance - 1 | log_eta[i,j,k], phi) +
              binomial_logit_lpmf( // 0 individual-level detections, citizen science
                0 | abundance - 1, logit_p_citsci[i,j,k]) +
              binomial_logit_lpmf( // 0 species-level detecions, museums
                // (given the number of community sampling events that occurred)
                0 | sum(V_museum_NA[i,j,k]), logit_p_museum[i,j,k]);
                
          } // end for abundance in 1:K
          
          // sum lp of both possibilities of availability
          // and all possible abundance states that went unobserved if it's available
          target += log_sum_exp(lp);
        
        } // end else species never detected at site during interval
            
        } // end if in range
          
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
  
} // end model

generated quantities {
  
  // Posterior Predictive Check
  
  int<lower=0> N[n_species,n_sites,n_intervals]; // predicted abundance 

  real eval[n_species,n_sites,n_intervals,n_visits]; // Expected values
  
  int y_new[n_species,n_sites,n_intervals,n_visits]; // new data for counts generated from eval
    
  real E[n_species,n_sites,n_intervals,n_visits]; // squared scaled distance of real data from expected value
  real E_new[n_species,n_sites,n_intervals,n_visits]; // squared scaled distance of new data from expected value
  
  real fit = 0; // sum squared distances of real data across all observation intervals
  real fit_new = 0; // sum squared distances of new data across all observation intervals
 
  // predict abundance given log_eta
  for (i in 1:n_species){ // loop across all species
    for (j in 1:n_sites){ // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
        if(bernoulli_logit_rng(omega[i,j,k]) == 1){ // if the site is suitable
            
            // predict the abundance from the count distribution
            N[i,j,k] = neg_binomial_2_rng(exp(log_eta[i,j,k]), phi);
            
        } else { // else the site is not suitable
          
          // and abundance is an excess zero
          N[i,j,k] = 0;
          
        } // end if/else 
      
      } // loop across all intervals
    } // loop across all sites
  } // loop across all species
    
  // Initialize E and E_new
  for(l in 1:n_visits) {
      E[1,1,1,l] = 0;
      E_new[1,1,1,l] = 0;
  }
  
  for (i in 2:n_species){
    for(j in 2:n_sites){
      for(k in 2:n_intervals){
        
        E[i,j,k] = E[i-1,j-1,k-1];
        E_new[i,j,k] = E_new[i-1,j-1,k-1];
    
      }
    }
  }
  
  for (i in 1:n_species){ // loop across all species
    for (j in 1:n_sites){ // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        for(l in 1:n_visits){ // loop across all visits
         
          // Assess model fit using Chi-squared discrepancy
          // Compute fit statistic E for observed data
          eval[i,j,k,l] = inv_logit(logit_p_citsci[i,j,k]) * neg_binomial_2_rng(exp(log_eta[i,j,k]), phi); // expected value at observation i for visit j 
          // (probabilty across visits is fixed) is = expected detection prob * expected abundance
          // Compute fit statistic E_new for real data (V)
          E[i,j,k,l] = square(V_citsci[i,j,k,l] - eval[i,j,k,l]) / (eval[i,j,k,l] + 0.5);
          // Generate new replicate count data and
          y_new[i,j,k,l] = binomial_rng(N[i,j,k], inv_logit(logit_p_citsci[i,j,k]));
          // Compute fit statistic E_new for replicate data
          E_new[i,j,k,l] = square(y_new[i,j,k,l] - eval[i,j,k,l]) / (eval[i,j,k,l] + 0.5);
      
        } // loop across all visits
    
        fit = fit + sum(E[i,j,k]); // descrepancies for each site*species*interval combo (across 1:l visits)
        fit_new = fit_new + sum(E_new[i,j,k]); // descrepancies for generated data (across 1:l visits)
                                      
      } // loop across all intervals
    } // loop across all sites
  } // loop across all species
  
}
