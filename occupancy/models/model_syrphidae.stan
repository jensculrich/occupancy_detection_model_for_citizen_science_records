// multi-species occupancy model for GBIF occurrence data
// jcu, started nov 21, 2022.

functions {
  
  // covariance matrix for detection rates from different data sources
  matrix custom_cov_matrix(vector sigma, real rho) {
    matrix[2,2] Sigma;
    Sigma[1,1] = square(sigma[1]); // species variation in community science detection rates
    Sigma[2,2] = square(sigma[2]); // species variation in research collection detection rates
    Sigma[1,2] = sigma[1] * sigma[2] * rho; // correlation between species-specific detection rates
    Sigma[2,1] = Sigma[1,2]; // correlation between species-specific detection rates
    return Sigma;
  }
  
  // mean values for multivariate normal distribution (center species on global intercepts)
  vector mu(real mu_p_cs_0, real mu_p_rc_0){
    vector[2] global_intercepts;
    global_intercepts[1] = mu_p_cs_0; 
    global_intercepts[2] = mu_p_rc_0;
    return global_intercepts;
  }
  
}

data {
  
  int<lower=1> n_species; // number of species
  int<lower=1> species[n_species]; // vector of species identities
  
  int<lower=1> n_sites;  // (number of) sites within region (level-2 clusters)
  int<lower=1, upper=n_sites> sites[n_sites];  // vector of sites identities (level-2 clusters)
  int<lower=1> n_level_three;  // (number of) fine-scale ecoregion areas (level-3 clusters)
  int<lower=1> n_level_four;  // (number of) broad-scale ecoregion areas (level-4 clusters)
  int<lower=1> level_three_lookup[n_sites]; // level-3 cluster look up vector for level-3 cluster
  int<lower=1> level_four_lookup[n_level_three]; // level-4 cluster look up vector for level-4 cluster

  int<lower=1> n_intervals;  // number of intervals during which sites are sampled
  int intervals[n_intervals]; // vector of intervals
  int<lower=1> n_visits; // visits within intervals
  
  int<lower=0> V_cs[n_species, n_sites, n_intervals, n_visits];  // visits l when species i was detected at site j on interval k (by community science process)
  int<lower=0> V_rc[n_species, n_sites, n_intervals, n_visits];  // visits l when species i was detected at site j on interval k (by NHC process)
  
  int<lower=0> ranges[n_species, n_sites, n_intervals, n_visits];  // NA indicator where 1 == site is in range, 0 == not in range
  int<lower=0> V_rc_NA[n_species, n_sites, n_intervals, n_visits];  // indicator where 1 == sampled, 0 == missing data (for NHC process)
  
  vector[n_sites] site_areas; // (scaled) spatial area extent
  vector[n_sites] pop_densities; // (scaled) population density
  vector[n_sites] avg_income; // (scaled) household income 
  vector[n_sites] avg_racial_minority; // (scaled) prop. of racial minority population
  vector[n_sites] natural_habitat; // (scaled) natural greenspace area
  vector[n_sites] open_developed; // (scaled) developed greenspace area
  
  vector[n_species] nativity; // nativity vector, 0 == non-native, 1 == native
  
} // end data


parameters {
  
  // OCCUPANCY //
  
    real mu_psi_0; // global intercept for occupancy
    
    // species-specific random intercepts for occupancy
    vector[n_species] psi_species_raw; // species specific intercept for occupancy
    real<lower=0> sigma_psi_species; // variance in species intercepts// Level-3 spatial random effect
    
    // nested site-specific random intercepts for occupancy
    vector[n_sites] psi_site_raw; // site specific intercept for occupancy
    real<lower=0> sigma_psi_site; // variance in site intercepts
    vector[n_level_three] psi_level_three_raw; // level-three intercept for occupancy
    real<lower=0> sigma_psi_level_three; // variance in level-three intercepts
    vector[n_level_four] psi_level_four_raw; // level-four specific intercept for occupancy
    real<lower=0> sigma_psi_level_four; // variance in level-four intercepts
    
    // species specific random effects of natural greenspace on occupancy (with nativity predictors)
    vector[n_species] psi_natural_habitat; // vector of species specific slope estimates
    real delta0; // baseline effect (mean)
    real delta1; // effect of being native on the expected value of the random effect
    real<lower=0> gamma0; // baseline effect (variance) (negative variance not possible)
    real gamma1; // effect of being native on the expected value of the random effect
    
    // fixed effects for other occupancy predictors
    real mu_psi_open_developed; // effect of developed greenspace area
    real mu_psi_income; // effect of income
    real mu_psi_race; // effect of prop. of racial minorities
    real psi_site_area; // effect of site area
  
  // DETECTION //
  
    // Covararying species-specific detection intercepts
    real<lower=-1,upper=1> rho;  // correlation of (community science and NHC detection)
    vector<lower=0>[2] sigma_species_detection; // variance in species-specific detection rates 
      // (sigma_species_detection[1] == community science and sigma_species_detection[2] == research collections detection)
    vector[2] species_intercepts_detection[n_species];// species-level detection intercepts
      // (species_intercepts_detection[1] == community science and species_intercepts_detection[2] == research collections detection)
    
    // COMMUNITY SCIENCE
    
      real mu_p_cs_0; // global detection intercept for community science records
      
      // nested site-specific random intercepts for occupancy
      vector[n_sites] p_cs_site_raw; // site intercepts
      real<lower=0> sigma_p_cs_site; // variance in site intercepts
      vector[n_level_three] p_cs_level_three_raw; // level-three intercepts 
      real<lower=0> sigma_p_cs_level_three;  // variance in level-three intercepts
      
      real p_cs_interval; // fixed effect of time interval on cs detection probability
      real p_cs_pop_density; // fixed effect of population density on cs detection probability
      real p_cs_income; // fixed effect of income on cs detection probability
      real p_cs_race; // fixed effect of prop. of racial minorities on cs detection probability
  
    // NATURAL HISTORY COLLECTIONS
    
      real mu_p_rc_0; // global detection intercept for research collections records
    
      vector[n_sites] p_rc_site_raw; // site intercepts
      real<lower=0> sigma_p_rc_site; // variance in site intercepts
      vector[n_level_three] p_rc_level_three_raw; // level-three intercepts 
      real<lower=0> sigma_p_rc_level_three; // variance in level-three intercepts
  
} // end parameters


transformed parameters {
  
  // logit scaled expected values
  real logit_psi[n_species, n_sites, n_intervals];  // odds of occurrence
  real logit_p_cs[n_species, n_sites, n_intervals]; // odds of detection by community science
  real logit_p_rc[n_species, n_sites, n_intervals]; // odds of detection by research collections
  
  // non centered parameterizations of random intercepts //
  
  // species intercepts
  vector[n_species] psi_species;
  psi_species = sigma_psi_species * psi_species_raw;
  
  // spatially nested intercepts
  vector[n_sites] psi_site;
  vector[n_level_three] psi_level_three;
  vector[n_level_four] psi_level_four;

  vector[n_sites] p_cs_site;
  vector[n_level_three] p_cs_level_three;

  vector[n_sites] p_rc_site;
  vector[n_level_three] p_rc_level_three;
  
  vector[n_species] mu_psi_natural_habitat; // expected value for species specific slopes
  vector[n_species] sigma_psi_natural_habitat; // expected value for variance among species slopes

  // using for() loops to nest the random intercepts for sites
  
  //
  // compute the varying community science detection intercept at the ecoregion1 level
  // Level-4 (n_level_four level-4 random intercepts)
  psi_level_four = sigma_psi_level_four * psi_level_four_raw;
  // compute the varying community science detection intercept at the ecoregion3 level
  // Level-3 (n_level_three level-3 random intercepts)
  for(i in 1:n_level_three){
    psi_level_three[i] = psi_level_four[level_four_lookup[i]] + 
      sigma_psi_level_three * psi_level_three_raw[i];
  }
  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  for(i in 1:n_sites){
    psi_site[i] = psi_level_three[level_three_lookup[i]] + 
      sigma_psi_site * psi_site_raw[i];
  }
  
  //
  // compute the varying community science detection intercept at the ecoregion3 level
  // Level-3 (n_level_three level-3 random intercepts)
  for(i in 1:n_level_three){
    p_cs_level_three[i] = sigma_p_cs_level_three * p_cs_level_three_raw[i];
  }
  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  for(i in 1:n_sites){
    p_cs_site[i] = p_cs_level_three[level_three_lookup[i]] + 
      sigma_p_cs_site * p_cs_site_raw[i];
  }
  
  //
  // compute the varying community science detection intercept at the ecoregion3 level
  // Level-3 (n_level_three level-3 random intercepts)
  for(i in 1:n_level_three){
    p_rc_level_three[i] = sigma_p_rc_level_three * p_rc_level_three_raw[i];
  }
  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  for(i in 1:n_sites){
    p_rc_site[i] = p_rc_level_three[level_three_lookup[i]] + 
      sigma_p_rc_site * p_rc_site_raw[i];
  }
  
  //
  // hard prior to disallow intercept plus nativity adjustment from being negative
  real<lower=0> gamma0_plus_gamma1;
  gamma0_plus_gamma1 = gamma0 + gamma1;
  
  // model the expected value for the random effect using a linear predictor that includes nativity
  for(i in 1:n_species){
    mu_psi_natural_habitat[i] = delta0 + delta1*nativity[i];
  }
  
  // model the group level variation (allow native and non-native groups to have different amounts of variation among species)
  for(i in 1:n_species){
    sigma_psi_natural_habitat[i] = gamma0 + gamma1*nativity[i];
  }
  
  //
  //
  
  // calculate logit scaled expected values for occurrence and detection
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals  
          
          logit_psi[i,j,k] = // the inverse of the log odds of occurrence is equal to..
            mu_psi_0 + // a global intercept
            psi_species[species[i]] + // a species-specific intercept
            psi_site[sites[j]] + // a spatially nested, site-specific intercept
            psi_natural_habitat[species[i]]*natural_habitat[j] + // a species-specific effect of natural greenspace area
            mu_psi_open_developed*open_developed[j] + // a species-specific effect of developed greenspace area
            mu_psi_income*avg_income[j] + // a species-specific effect of income
            mu_psi_race*avg_racial_minority[j] + // an effect of ethnic composition
            psi_site_area*site_areas[j] // an effect of spatial area of the site on occurrence
            ; // end psi[i,j,k]
            
      } // end loop across all intervals
    } // end loop across all sites
  }  // end loop across all species
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
          logit_p_cs[i,j,k] = // the inverse of the log odds of detection is equal to..
            species_intercepts_detection[species[i],1] + // a correlated species-specific intercept // includes global intercept
            p_cs_site[sites[j]] + // a spatially nested, site-specific intercept
            p_cs_interval*(intervals[k]^2) + // an effect of time on detection
            p_cs_pop_density*pop_densities[j] + // an effect of pop density on detection
            p_cs_income*avg_income[j] + // an effect of income on detection
            p_cs_race*avg_racial_minority[j] // an effect of prop. of racial minorities on detection
           ; // end p_cs[i,j,k]

          logit_p_rc[i,j,k] = // the inverse of the log odds of detection is equal to..
            species_intercepts_detection[species[i],2] + // a species specific intercept // includes global intercept
            p_rc_site[sites[j]] // a spatially nested, site-specific intercept
           ; // end p_rc[i,j,k]
          
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
             
  
} // end transformed parameters


model {
  
  // PRIORS
  
  // correlated species effects for detection
  sigma_species_detection[1] ~ normal(0, 2);
  sigma_species_detection[2] ~ normal(0, 2);
  (rho + 1) / 2 ~ beta(2, 2);
  
  // correlated species-specific detection rates
  // will send the mean (mu), variance and correlation to the covariance matrix
  species_intercepts_detection ~ multi_normal(mu(mu_p_cs_0, mu_p_rc_0), 
    custom_cov_matrix(sigma_species_detection, rho));
  
  // Occupancy (Ecological Process)
  mu_psi_0 ~ normal(0, 1); // global intercept for occupancy rate
  
  psi_species_raw ~ std_normal(); // species intercept effects
  sigma_psi_species ~ normal(0, 0.5); // weakly-informative prior
  
  // level-2 spatial grouping
  psi_site_raw ~ std_normal();
  sigma_psi_site ~ normal(0, 0.5); // weakly-informative prior
  // level-3 spatial grouping
  psi_level_three_raw ~ std_normal();
  sigma_psi_level_three ~ normal(0, 0.5); // weakly-informative prior
  // level-4 spatial grouping
  psi_level_four_raw ~ std_normal();
  sigma_psi_level_four ~ normal(0, 0.5); // weakly-informative prior

  psi_natural_habitat ~ normal(mu_psi_natural_habitat, sigma_psi_natural_habitat);
  // community effect (mu) and variation among species (sigma) is defined as a vector 
  // with intercept delta0 and an effect of nativity (delta1) on community mean
  // and intercept gamma0 and an effect of nativity (gamma1) on variation
  delta0 ~ normal(0, 2); // community mean
  delta1 ~ normal(0, 1); // effect of nativity
  gamma0 ~ normal(0, 1); // community mean variance
  gamma1 ~ normal(0, 0.25); // effect of nativity on variance
  
  mu_psi_open_developed ~ normal(0, 2); // effect of developed greenspace area
  mu_psi_income ~ normal(0, 2); // effect of income
  mu_psi_race ~ normal(0, 2); // effect of prop. of racial minorities
  psi_site_area ~ normal(0, 2); // effect of site area on occupancy
  
  // Detection (Observation Process)
  
  // community science records
  
  mu_p_cs_0 ~ normal(0, 1); // global intercept for detection
  
  // level-2 spatial grouping
  p_cs_site_raw ~ std_normal();
  sigma_p_cs_site ~ normal(0, 0.5); // weakly-informative prior
  // level-3 spatial grouping
  p_cs_level_three_raw ~ std_normal();
  sigma_p_cs_level_three ~ normal(0, 0.5); // weakly-informative prior
  
  p_cs_interval ~ normal(0, 2); // effect if time inteval on detection probability
  p_cs_pop_density ~ normal(0, 2); // effect of popualtion density on detection probability
  p_cs_income ~ normal(0, 2); // effect of income on detection probability
  p_cs_race ~ normal(0, 2); // effect of prop. racial minorities on detection probability
  
  // museum records
  
  mu_p_rc_0 ~ normal(0, 0.5); // global intercept for detection
  
  // level-2 spatial grouping
  p_rc_site_raw ~ std_normal();
  sigma_p_rc_site ~ normal(0, 0.25); // weakly-informative prior
  // level-3 spatial grouping
  p_rc_level_three_raw ~ std_normal();
  sigma_p_rc_level_three ~ normal(0, 0.25); // weakly-informative prior
  
  // LIKELIHOOD
  
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
        // If the site is in the range of a species, then evaluate lp, otherwise do not (treat as NA).
        if(sum(ranges[i,j,k]) > 0){ // The sum of the NA vector will be == 0 if site is not in range
        
          // if species is detected at the specific site*interval at least once
          // by community science efforts OR research collection records
          // then the species occurs there. lp_observed calculates
          // the probability density that species occurs given psi, plus the 
          // probability density that we did/did not observe it on all years within the interval.
          if(sum(V_cs[i, j, k, 1:n_visits]) > 0 || sum(V_rc[i, j, k, 1:n_visits]) > 0) {
            
             // lp_observed:
             target += log_inv_logit(logit_psi[i,j,k]) +
                      binomial_logit_lpmf(sum(V_cs[i,j,k,1:n_visits]) | n_visits, logit_p_cs[i,j,k]) + 
                      binomial_logit_lpmf(sum(V_rc[i,j,k,1:n_visits]) | sum(V_rc_NA[i,j,k,1:n_visits]), logit_p_rc[i,j,k]);

          // else the species was never detected at the site*interval
          // lp_unobserved sums the probability density of:
          // 1) species occupies the site*interval but was not detected on each visit, and
          // 2) the species does not occupy the site*interval
          } else {
            
            // lp_unobserved
              // Stan can sample the mean and sd of parameters by summing out the
              // parameter (marginalizing) across likelihood statements
            target += 
                    // present but never detected
                    log_sum_exp(log_inv_logit(logit_psi[i,j,k]) +
                    binomial_logit_lpmf(0 | n_visits, logit_p_cs[i,j,k]) +
                    binomial_logit_lpmf(0 | sum(V_rc_NA[i,j,k,1:n_visits]), logit_p_rc[i,j,k]),
                    // not present
                    log1m_inv_logit(logit_psi[i,j,k])); 
            
          } // end if/else ever observed
        
        } // end if/ in range
          
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
  
} // end model

generated quantities{
  
  // Effect of nat habitat area for all, native and non-native species is a derived parameter
  // we estimate the expected value for all species, native, or non-native species using our linear predictor.
  // By doing this in each step of the HMC we can get a distribution of outcomes 
  // (propagating our uncertainty in delta0 and delta1 to an uncertainty in the group level effects)
  
  real mu_psi_natural_habitat_native;
  mu_psi_natural_habitat_native = delta0 + delta1*1;
  
  real mu_psi_natural_habitat_nonnative;
  mu_psi_natural_habitat_nonnative = delta0 + delta1*0;
  
  real mu_psi_natural_habitat_all_species;
  mu_psi_natural_habitat_all_species = mean(mu_psi_natural_habitat);
  
  //
  // posterior predictive check (number of detections, binned by species) //
  //
  
  // simulate occurrence of species at each site in each year
  int z_simmed[n_species, n_sites, n_intervals]; // simulate occurrence

  for(i in 1:n_species){
   for(j in 1:n_sites){
     for(k in 1:n_intervals){
          z_simmed[i,j,k] = bernoulli_logit_rng(logit_psi[i,j,k]); 
     } // end loop across intervals
    } // end loop across sites
  } // end loop across species 
  
  // now simulate some detections and see how it compares to the observed number of detections
  int<lower=0> W_species_rep_cs[n_species]; // sum of simulated detections
  int<lower=0> W_species_rep_rc[n_species]; // sum of simulated detections

  // initialize at 0
  for(i in 1:n_species){
    W_species_rep_cs[i] = 0;
    W_species_rep_rc[i] = 0;
  } // end loop across species
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all years

          if(sum(ranges[i,j,k]) > 0){
            
            // detections in replicated data (us z_simmed from above)
            W_species_rep_cs[i] = W_species_rep_cs[i] + 
              (z_simmed[i,j,k] * binomial_rng(n_visits, inv_logit(logit_p_cs[i,j,k])));
            
            // detections in replicated data (us z_simmed from above)
            W_species_rep_rc[i] = W_species_rep_rc[i] + 
              (z_simmed[i,j,k] *  binomial_rng(sum(V_rc_NA[i,j,k,]), inv_logit(logit_p_rc[i,j,k])));
           
          } // end if()
           
      } // end loop across intervals
    } // end loop across sites
  } // end loop across species
  
} // end generated quantities
