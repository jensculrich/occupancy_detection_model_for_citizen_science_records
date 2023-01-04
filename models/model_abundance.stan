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
  
  vector[n_sites] site_areas; // (scaled) spatial area extent of each site
  vector[n_sites] pop_densities; // (scaled) population density of each site
  vector[n_sites] impervious_surfaces; // (scaled) impervious surface cover of each site
  vector[n_sites] perennial_plant_cover; // (scaled) perennial plant cover of each site
  
} // end data

transformed data {
  int<lower=0> max_y[n_species, n_sites, n_intervals];

  for (i in 1:n_species) {
    for(j in 1:n_sites){
      for(k in 1:n_intervals){
        
        // Set the floor of the latent state search to be at least as many as the most 
        // that we observed of a species at a site in a time interval (by cit sci records)
        max_y[i,j,k] = max(V_citsci[i,j,k]);
        
        // We only search abundance if it must be greater than 1, i.e.,
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
  
  real gamma_0; // occupancy intercept
  real gamma_1; // relationship between abundance and occupancy 
  
  real<lower=0> phi; // abundance overdispersion parameter
  
  real mu_eta_0; // global intercept for occupancy
  
  // species-specific intercept allows some species abundance at higher rates than others, 
  // but with overall estimates for abundance partially informed by the data pooled across all species.
  vector[n_species] eta_species; // species specific intercept for occupancy
  real<lower=0> sigma_eta_species; // variance in species intercepts
  
  // site-specific intercept allows some site abundance at higher rates than others, 
  // but with overall estimates for abundance partially informed by the data pooled across all species.
  vector[n_sites] eta_site; // species specific intercept for abundance
  real<lower=0> sigma_eta_site; // variance in species intercepts
  
  // random slope for species-specific impervious surface effects on abundance
  vector[n_species] eta_impervious_surface; // vector of species specific slope estimates
  real mu_eta_impervious_surface; // community mean of species specific slopes
  real<lower=0> sigma_eta_impervious_surface; // variance in species slopes
  
  // random slope for species-specific perennial plant cover effects on abundance
  vector[n_species] eta_perennial_plants; // vector of species specific slope estimates
  real mu_eta_perennial_plants; // community mean of species specific slopes
  real<lower=0> sigma_eta_perennial_plants; // variance in species slopes
  
  real eta_site_area; // fixed effect of site area on abundance
  
  // DETECTION
  
  // citizen science observation process
  real mu_p_citsci_0; // global detection intercept for citizen science records
  
  // species-specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for detection partially informed by the data pooled across all species.
  vector[n_species] p_citsci_species; // species specific intercept for detection
  real<lower=0> sigma_p_citsci_species; // variance in species intercepts
  
  // site-specific intercept allows some sites to have higher detection rates than others, 
  // but with overall estimates for detection partially informed by the data pooled across all species.
  vector[n_sites] p_citsci_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_citsci_site; // variance in site intercepts
  
  // random intercept for site-specific temporal effects on detection
  vector[n_sites] p_citsci_interval; // vector of species specific slope estimates
  real mu_p_citsci_interval; // community mean of species specific slopes
  real<lower=0> sigma_p_citsci_interval; // variance in species slopes
  
  real p_citsci_pop_density; // fixed effect of population on detection probability
  
  // museum records observation process
  real mu_p_museum_0; // global detection intercept for citizen science records
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] p_museum_species; // species specific intercept for detection
  real<lower=0> sigma_p_museum_species; // variance in species intercepts
  
  // site-specific intercept allows some sites to have higher detection rates than others, 
  // but with overall estimates for detection partially informed by the data pooled across all species.
  vector[n_sites] p_museum_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_museum_site; // variance in site intercepts
  
  // random intercept for site-specific temporal effects on detection
  vector[n_sites] p_museum_interval; // vector of species specific slope estimates
  real mu_p_museum_interval; // community mean of species specific slopes
  real<lower=0> sigma_p_museum_interval; // variance in species slopes
  
  real p_museum_pop_density; // fixed effect of population on detection probability
  
  
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
            eta_site[sites[j]] + // a site specific intercept
            eta_impervious_surface[species[i]]*impervious_surfaces[j] + // a species specific temporal effect
            eta_perennial_plants[species[i]]*perennial_plant_cover[j] + // an effect of pop density on occurrence
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
            p_citsci_species[species[i]] + // a species specific intercept
            p_citsci_site[sites[j]] + // a spatially specific intercept
            p_citsci_interval[sites[j]]*intervals[k] + // an overall effect of time on detection
            p_citsci_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_citsci[i,j,k]
           
          logit_p_museum[i,j,k] = // logit scaled species-level detection rate
            mu_p_museum_0 + // a baseline intercept
            p_museum_species[species[i]] + // a species specific intercept
            p_museum_site[sites[j]] + // a spatially specific intercept
            p_museum_interval[sites[j]]*intervals[k] + // an overall effect of time on detection
            p_museum_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_museum[i,j,k]
           
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
  
} // end transformed parameters


model {
  
  // PRIORS
  
  // Abundance (Ecological Process)
  
  gamma_0 ~ normal(0, 1); // occupancy-abundance relationship intercept
  gamma_1 ~ normal(0, 1); // effect of expected abundance on occupancy rate
  
  phi ~ normal(0, 2); // abundance overdispersion scale parameter
  
  mu_eta_0 ~ normal(0, 2); // global intercept for abundance rate
  
  eta_species ~ normal(0, sigma_eta_species); 
  // occupancy intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_eta_species ~ normal(0, 1); // weakly informative prior (for strong pooling across species)
  
  eta_site ~ normal(0, sigma_eta_site); 
  // occupancy intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_eta_site ~ normal(0, 1); // weakly informative prior (for strong pooling across species)
  
  // Random slope for species-specfic effect of impervious surfaces on abundance
  eta_impervious_surface ~ normal(mu_eta_impervious_surface, sigma_eta_impervious_surface);
  mu_eta_impervious_surface ~ normal(0, 2); // community mean
  sigma_eta_impervious_surface ~ normal(0, 1); // community variance // weakly informative prior (for strong pooling across species)

  // Random slope for species-specfic effect of perennial plant cover on abundance
  eta_impervious_surface ~ normal(mu_eta_impervious_surface, sigma_eta_impervious_surface);
  mu_eta_impervious_surface ~ normal(0, 2); // community mean
  sigma_eta_impervious_surface ~ normal(0, 1); // community variance // weakly informative prior (for strong pooling across species)
  
  eta_site_area ~ normal(0, 2); // effect of site area on abundance rate
  
  // Detection (Observation Process)
  
  // CITIZEN SCIENCE records
  
  mu_p_citsci_0 ~ normal(0, 2); // global intercept for (citizen science) detection
  
  p_citsci_species ~ normal(0, sigma_p_citsci_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_citsci_species ~ normal(0, 1); // community variance // weakly informative prior (for strong pooling across species)
  
  p_citsci_site ~ normal(0, sigma_p_citsci_site); 
  // detection intercept for each site drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_citsci_site ~ normal(0, 1); // community variance // weakly informative prior (for strong pooling across species)
  
  // Random slope for site-specfic effect of time on detection
  p_citsci_interval ~ normal(mu_p_citsci_interval, sigma_p_citsci_interval);
  mu_p_citsci_interval ~ normal(0, 2); // community mean
  sigma_p_citsci_interval ~ normal(0, 1); // community variance // weakly informative prior (for strong pooling across species)
  
  p_citsci_pop_density ~ normal(0, 2); // effect of population density on detection

  // MUSEUM records
  mu_p_museum_0 ~ normal(0, 2); // global intercept for (museum) detection
  
  p_museum_species ~ normal(0, sigma_p_museum_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_museum_species ~ normal(0, 1); // community variance // weakly informative prior (for strong pooling across species)
  
  p_museum_site ~ normal(0, sigma_p_museum_site); 
  // detection intercept for each site drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_museum_site ~ normal(0, 1); // community variance // weakly informative prior (for strong pooling across species)
  
  // Random slope for site-specfic effect of time on detection
  p_museum_interval ~ normal(mu_p_museum_interval, sigma_p_museum_interval);
  mu_p_museum_interval ~ normal(0, 2); // community mean
  sigma_p_museum_interval ~ normal(0, 1); // community variance
  
  p_museum_pop_density ~ normal(0, 2); // effect of population density on detection

  
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
                max_y[i,j,k] + abundance - 1 | log_eta[i,j,k], phi)  
              - neg_binomial_2_lccdf(0 | // with that count distribution truncated to be greater than zero
                exp(log_eta[i,j,k]), phi)
              // with abundance detections vectorized over n visits..
              + binomial_logit_lpmf( // individual-level detection, citizen science
                V_citsci[i,j,k] | max_y[i,j,k] + abundance - 1, logit_p_citsci[i,j,k]) 
              + binomial_logit_lpmf( // binary, species-level detecion, museums
                // (given the number of community sampling events that occurred)
                sum(V_museum[i,j,k]) | sum(V_museum_NA[i,j,k]), logit_p_museum[i,j,k])
              // plus outcome of site being available, given the 
              // abundance-dependent probability of suitability
              + log(inv_logit(omega[i,j,k])); 
          
          }
                
          target += log_sum_exp(lp);
        
        } else { // else was never detected and the site may or may not be available
          
          real lp[2];
          
          // probability occupying a site with
          // some unknown latent abundance state;
          // but never observed.
          for(abundance in 1:(K[i,j,k])){
            
            // outcome of site unoccupied, given the 
            // abundance-dependent probability of occupancy
            lp[1] = log1m(inv_logit(omega[i,j,k])); // site not occupied by species in interval
            
            // outcome of site occupancy, given the 
            // abundance-dependent probability of occupancy
            lp[2] = log(inv_logit(omega[i,j,k])); // occupied but not observed
            
            // start at abundance = 1 so that the search has a floor of 1
            lp[2] = lp[2] 
             + neg_binomial_2_log_lpmf( // generation of abundance given count distribution
                abundance | log_eta[i,j,k], phi) 
             - neg_binomial_2_lccdf(0 | // with that count distribution truncated to be greater than zero
                exp(log_eta[i,j,k]), phi)
             + binomial_logit_lpmf( // but getting 0 individual-level detections, citizen science
                0 | abundance, logit_p_citsci[i,j,k]) 
             + binomial_logit_lpmf( // and 0 species-level detecions, museums
                // (given the number of community sampling events that occurred)
                0 | sum(V_museum_NA[i,j,k]), logit_p_museum[i,j,k]);
                
          } // end for abundance in 1:K
          
          // sum lp of both possibilities of occupancy
          // AND all possible abundance states that went unobserved if the site is occupied
          target += log_sum_exp(lp);
        
        } // end else species never detected at site during interval
            
        } // end if in range
          
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
  
} // end model

generated quantities {
  
  // Posterior Predictive Check
  
  // Abundance
  int<lower=0> N[n_species,n_sites,n_intervals]; // expected abundance 

  real eval[n_species,n_sites,n_intervals,n_visits]; // expected values
  
  int y_new[n_species,n_sites,n_intervals,n_visits]; // new data for counts generated from eval
    
  real E[n_species,n_sites,n_intervals,n_visits]; // squared scaled distance of real data from expected value
  real E_new[n_species,n_sites,n_intervals,n_visits]; // squared scaled distance of new data from expected value
  
  real fit = 0; // sum squared distances of real data from expected values
  real fit_new = 0; // sum squared distances of new data from expected values
 
  // Occupancy
  int<lower=0,upper=1> psi[n_species,n_sites,n_intervals]; // expected occupancy  

  real eval_binary_detection[n_species,n_sites,n_intervals,n_visits]; // expected values
  
  int y_occupancy_new[n_species,n_sites,n_intervals,n_visits]; // new data for counts generated from eval
    
  real E_occupancy[n_species,n_sites,n_intervals,n_visits]; // squared scaled distance of real data from expected value
  real E_occupancy_new[n_species,n_sites,n_intervals,n_visits]; // squared scaled distance of new data from expected value
  
  real fit_occupancy = 0; // sum squared distances of real data from expected values
  real fit_occupancy_new = 0; // sum squared distances of new data from expected values
 
  // Initialize E and E_new
  for(l in 1:n_visits){
        
    E[1,1,1,l] = 0;
    E_new[1,1,1,l] = 0;
    
    E_occupancy[1,1,1,l] = 0;
    E_occupancy_new[1,1,1,l] = 0;
    
  } 
  
  for (i in 2:n_species){
    for(j in 2:n_sites){
      for(k in 2:n_intervals){
        
        E[i,j,k] = E[i-1,j-1,k-1];
        E_new[i,j,k] = E_new[i-1,j-1,k-1];
        
        E_occupancy[i,j,k] = E_occupancy[i-1,j-1,k-1];
        E_occupancy_new[i,j,k] = E_occupancy_new[i-1,j-1,k-1];
        
      }
    }
  }
  
  // Generare expected values for abundance and occupancy
  // Considering that sites are unoccupied if outside of the species' range, and
  // sites have abundance > 0 only if the site is occupied.
  // Further, will need to settle a discrete abundance state using categorical_rng(softmax(lp))
  for (i in 1:n_species){ // loop across species
    for(j in 1:n_sites){ // loop across sites
      for(k in 1:n_intervals){ // loop across intervals
        
        // if the site is in range and is predicted to be occupied
        if(sum(ranges[i,j,k]) > 0){ 
        
        vector[K[i,j,k] - max_y[i,j,k] + 1] lp; // lp vector of length of possible abundances 
            // (abundance is at least as big as the max observed count but may range up to K)
            // add +1 to also include possibility that max_y is the actual true abundance
          
          // for each possible abundance in max observed abundance through K:
          for(abundance in 1:(K[i,j,k] - max_y[i,j,k] + 1)){ 
            
            // max_y[i,j,k] + abundance - 1 starts the search at max_y
            lp[abundance] = 
              neg_binomial_2_log_lpmf(max_y[i,j,k] + abundance - 1 | log_eta[i,j,k], phi)
              - neg_binomial_2_lccdf(0 | // with that count distribution truncated to be greater than zero
                exp(log_eta[i,j,k]), phi)
              + binomial_logit_lpmf(V_citsci[i,j,k] | // and some observed data with a detection rate p
                max_y[i,j,k] + abundance - 1, logit_p_citsci[i,j,k])
              + binomial_logit_lpmf( // binary, species-level detecion, museums
              // (given the number of community sampling events that occurred)
              sum(V_museum[i,j,k]) | sum(V_museum_NA[i,j,k]), logit_p_museum[i,j,k]);
          
          }
          
          // N is the most likely latent discrete state from the truncated count distribution
          N[i,j,k] = categorical_rng(softmax(lp));
        
        } else { // else the site is not in range
            
          // and the latent discrete abundance is zero
          N[i,j,k] = 0;
            
        } // end if/else site is in range and is occupied
        
        // Expected abundance is abundance*occupancy state (if not occupied, then abundance is 0)
        N[i,j,k] = N[i,j,k] * bernoulli_logit_rng(omega[i,j,k]);
        
        // Occupancy state is 1 only if site is both in range and bernoulli_logit_rng() = 1
        if(N[i,j,k] == 0){
          
          psi[i,j,k] = 0;
        
        } else{
          
          psi[i,j,k] = 1;
          
        }
  
      } // end loop across intervals
    } // end loop acros sites
  } // end loop across species
  
  
  for (i in 1:n_species){ // loop across all species
    for (j in 1:n_sites){ // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        for(l in 1:n_visits){ // loop across all visits
          
          // create a temporary count variable, since the distribution needs to be truncated
          int w[n_species,n_sites,n_intervals];
                      
          // Assess model fit using Chi-squared discrepancy
          // Compute fit statistic E for observed data
          
          // The sum of the NA indicator vector ranges == 0 if site is not in range
          if(sum(ranges[i,j,k]) > 0){ 
            
            // Cit sci count data
            // expected count is... 
            eval[i,j,k,l] =
              N[i,j,k] // value of the latent abundance state with highest probability
              * inv_logit(logit_p_citsci[i,j,k]); // times the estimate for detection rate

            // Compute fit statistic E_new for real data (V_citsci)
            E[i,j,k,l] = square(V_citsci[i,j,k,l] - eval[i,j,k,l]) / (eval[i,j,k,l] + 0.5);
            
            // Generate new replicate count data and

            w[i,j,k] = neg_binomial_2_rng(exp(log_eta[i,j,k]), phi);
            
            // the count is from a truncated distr., and must be greater than 0.
            while(w[i,j,k] == 0){
              w[i,j,k] = neg_binomial_2_rng(exp(log_eta[i,j,k]), phi);
            }
            
            y_new[i,j,k,l] = 
              bernoulli_logit_rng(omega[i,j,k])
              * binomial_rng(w[i,j,k],
                inv_logit(logit_p_citsci[i,j,k]));
            
            // Compute fit statistic E_new for replicate data
            E_new[i,j,k,l] = square(y_new[i,j,k,l] - eval[i,j,k,l]) / (eval[i,j,k,l] + 0.5);
            
            // Musuem detections of occupancy
            // expected detection is...
            eval_binary_detection[i,j,k,l] =
              psi[i,j,k] // value of the latent abundance state with highest probability
              * inv_logit(logit_p_museum[i,j,k]); // times the estimate for detection rate
              
            // Compute fit statistic E_new for real data (V_citsci)
            E_occupancy[i,j,k,l] = square(V_museum[i,j,k,l] - eval_binary_detection[i,j,k,l]) / (eval_binary_detection[i,j,k,l] + 0.5);
            
            y_occupancy_new[i,j,k,l] = 
               binomial_rng(bernoulli_logit_rng(omega[i,j,k]), inv_logit(logit_p_museum[i,j,k]));
               
            // Compute fit statistic E_new for replicate data
            E_occupancy_new[i,j,k,l] = square(y_occupancy_new[i,j,k,l] - eval_binary_detection[i,j,k,l]) / (eval_binary_detection[i,j,k,l] + 0.5);
            
          } else { // end if site is in range and should not be considered
          
          // do not contribute to the sum squared distance from expected value
          // if the site is not in the species range
          E[i,j,k,l] = 0; // for real data
          E_new[i,j,k,l] = 0; // or for new data
          
          E_occupancy[i,j,k,l] = 0; // for real data
          E_occupancy_new[i,j,k,l] = 0; // or for new data
          
          } // end if/else site in range
      
        } // end loop across all visits
    
        fit = fit + sum(E[i,j,k]); // descrepancies for real data
        fit_new = fit_new + sum(E_new[i,j,k]); // descrepancies for generated data
        
        fit_occupancy = fit_occupancy + sum(E_occupancy[i,j,k]); // descrepancies for real data
        fit_occupancy_new = fit_occupancy_new + sum(E_occupancy_new[i,j,k]); // descrepancies for generated data
                                      
      } // loop across all intervals
    } // loop across all sites
  } // loop across all species
  
} // end generated quantities
