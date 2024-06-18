// multi-species occupancy model for BBNA occurrence data
// jcu, started nov 21, 2022.

data {
  
  int<lower=1> n_species; // number of observed species
  int<lower=1> species[n_species]; // vector of species identities
  
  int<lower=1> n_sites;  // number of sites (level-2 clusters)
  int<lower=1, upper=n_sites> sites[n_sites]; // vector of site identities
  int<lower=1> n_level_three;  // (number of) fine-scale ecoregion areas (level-3 clusters)
  int<lower=1> n_level_four;  // (number of) broad-scale ecoregion areas (level-4 clusters)

  int<lower=1> level_three_lookup[n_sites]; // level-3 cluster look up vector for level-3 cluster
  int<lower=1> level_four_lookup[n_level_three]; // level-4 cluster look up vector for level-4 cluster

  int<lower=1> n_intervals;  // intervals during which sites are visited
  
  int intervals[n_intervals]; // vector of intervals (used as covariate data for 
                                // species specific effect of occupancy interval (time) on detection)
  int<lower=1> n_visits; // visits within intervals
  
  int<lower=0> V_cs[n_species, n_sites, n_intervals, n_visits];  // visits l when species i was detected at site j on interval k
  int<lower=0> V_rc[n_species, n_sites, n_intervals, n_visits];  // visits l when species i was detected at site j on interval k
  
  int<lower=0> ranges[n_species, n_sites, n_intervals, n_visits];  // NA indicator where 1 == site is in range, 0 == not in range
  int<lower=0> V_rc_NA[n_species, n_sites, n_intervals, n_visits];  // indicator where 1 == sampled, 0 == missing data
  
  vector[n_sites] site_areas; // (scaled) spatial area extent of each site
  vector[n_sites] pop_densities; // (scaled) population density of each site
  vector[n_sites] avg_income; // (scaled) household income of each site
  vector[n_sites] natural_habitat; // (scaled) undeveloped open surface cover of each site
  vector[n_sites] open_developed; // (scaled) open developed surface cover of each site

} // end data


parameters {
  
  // OCCUPANCY
  real mu_psi_0; // global intercept for occupancy
  
  // species specific intercept allows some species to occur at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] psi_species; // species specific intercept for occupancy
  real<lower=0> sigma_psi_species; // variance in species intercepts
  
  // Spatially nested random effect on occupancy rates
  // Level-2 spatial random effect
  vector[n_sites] psi_site_raw; // site specific intercept for occupancy
  real<lower=0> sigma_psi_site; // variance in site intercepts
  // Level-3 spatial random effect
  vector[n_level_three] psi_level_three_raw; // level-three specific intercept for PL outcome
  real<lower=0> sigma_psi_level_three; // variance in level-three intercepts
  // Level-4 spatial random effect
  vector[n_level_four] psi_level_four_raw; // level-four specific intercept for PL outcome
  real<lower=0> sigma_psi_level_four; // variance in level-four intercepts
  
  // random slope for species specific natural habitat effects on occupancy
  vector[n_species] psi_natural_habitat; // vector of species specific slope estimates
  real mu_psi_natural_habitat; // community mean of species specific slopes
  real<lower=0> sigma_psi_natural_habitat; // variance in species slopes
  
  // random slope for species specific open developed effects on occupancy
  //vector[n_species] psi_open_developed; // vector of species specific slope estimates
  real mu_psi_open_developed; // community mean of species specific slopes
  //real<lower=0> sigma_psi_open_developed; // variance in species slopes
  
  // random slope for species specific household income effects on occupancy
  //vector[n_species] psi_income; // vector of species specific slope estimates
  real mu_psi_income; // community mean of species specific slopes
  //real<lower=0> sigma_psi_income; // variance in species slopes
  
  // fixed effect of site area on occupancy
  real psi_site_area;
  
  // DETECTION
  
  // citizen science observation process
  real mu_p_cs_0; // global detection intercept for citizen science records
  
  vector[n_species] p_cs_species_raw; // species specific intercept for detection
  real<lower=0> sigma_p_cs_species; // species variance on the intercept
  
  // random slope for site specific temporal effects on occupancy
  // level-2 spatial clusters
  vector[n_sites] p_cs_site_raw; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_cs_site; // variance in site slopes
  // level-3 spatial clusters
  vector[n_level_three] p_cs_level_three_raw; // level-three specific intercept for cs detection
  real<lower=0> sigma_p_cs_level_three;  // variance in level-three slopes
  
  real p_cs_interval; // fixed temporal effect on cs detection probability
  real p_cs_pop_density; // fixed effect of population on cs detection probability
  real p_cs_income; // fixed effect of income on cs detection probability
 
  // research collections records observation process
  real mu_p_rc_0; // global detection intercept for research collections records
  
  vector[n_species] p_rc_species_raw; // species specific intercept for detection
  real<lower=0> sigma_p_rc_species; // species variance on the intercept
  
  // random slope for site specific temporal effects on occupancy
  // level-2 spatial clusters
  vector[n_sites] p_rc_site_raw; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_rc_site; // variance in site slopes
  // level-3 spatial clusters
  vector[n_level_three] p_rc_level_three_raw; // level-three specific intercept for rc detection
  real<lower=0> sigma_p_rc_level_three; // variance in level-three slopes
  
} // end parameters


transformed parameters {
  
  real logit_psi[n_species, n_sites, n_intervals];  // odds of occurrence
  real logit_p_cs[n_species, n_sites, n_intervals]; // odds of detection by community science
  real logit_p_rc[n_species, n_sites, n_intervals]; // odds of detection by research collections
  
  // non-centered species random effects
  vector[n_species] p_cs_species;
  vector[n_species] p_rc_species;
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  p_cs_species = mu_p_cs_0 + sigma_p_cs_species * p_cs_species_raw;
  p_rc_species = mu_p_rc_0 + sigma_p_rc_species * p_rc_species_raw;

  // spatially nested intercepts
  vector[n_sites] psi_site;
  vector[n_level_three] psi_level_three;
  vector[n_level_four] psi_level_four;

  vector[n_sites] p_cs_site;
  vector[n_level_three] p_cs_level_three;

  vector[n_sites] p_rc_site;
  vector[n_level_three] p_rc_level_three;

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
  //
  
  // calculate logit scaled expected values for occurrence and detection
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals  
          
          logit_psi[i,j,k] = // the inverse of the log odds of occurrence is equal to..
            psi_species[species[i]] + // a species specific intercept
            psi_site[sites[j]] + // a spatially nested, site-specific intercept
            psi_natural_habitat[species[i]]*natural_habitat[j] + // a species-specific effect of natural habitat area
            mu_psi_income*avg_income[j] + // a species-specific effect of household income
            mu_psi_open_developed*open_developed[j] + // a species-specific effect of open developed land
            psi_site_area*site_areas[j] // an effect of spatial area of the site 
            ; // end psi[i,j,k]
            
      } // end loop across all intervals
    } // end loop across all sites
  }  // end loop across all species
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
          logit_p_cs[i,j,k] = // the inverse of the log odds of detection is equal to..
            p_cs_species[species[i]] + // a species specific intercept // includes global intercept
            p_cs_site[sites[j]] + // a spatially specific intercept 
            p_cs_interval*(intervals[k]^2) + // an overall effect of time on detection
            p_cs_pop_density*pop_densities[j] + // an overall effect of pop density on detection
            p_cs_income*avg_income[j] // an overall effect of income on detection
           ; // end p_cs[i,j,k]
           
          logit_p_rc[i,j,k] = // the inverse of the log odds of detection is equal to..
            p_rc_species[species[i]] + // a species specific intercept // includes global intercept
            p_rc_site[sites[j]] // a spatially specific intercept 
           ; // end p_rc[i,j,k]
           
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
             
  
} // end transformed parameters


model {
  
  // PRIORS
  
  // Occupancy (Ecological Process)
  mu_psi_0 ~ normal(0, 1); // global intercept for occupancy rate
  
  // level-2 spatial grouping
  psi_site_raw ~ std_normal();
  sigma_psi_site ~ normal(0, 0.5); // weakly-informative prior
  // level-3 spatial grouping
  psi_level_three_raw ~ std_normal();
  sigma_psi_level_three ~ normal(0, 0.5); // weakly-informative prior
  // level-4 spatial grouping
  psi_level_four_raw ~ std_normal();
  sigma_psi_level_four ~ normal(0, 0.5); // weakly-informative prior
  
  psi_species ~ normal(mu_psi_0, sigma_psi_species); 
  sigma_psi_species ~ normal(0, 1); // weakly-informative prior
  
  psi_natural_habitat ~ normal(mu_psi_natural_habitat, sigma_psi_natural_habitat);
  mu_psi_natural_habitat ~ normal(0, 2); // community mean
  sigma_psi_natural_habitat ~ normal(0, 1); // community variance
  
  mu_psi_income ~ normal(0, 2); // community mean
  mu_psi_open_developed ~ normal(0, 2); // community mean
  
  psi_site_area ~ normal(0, 2); // effect of site area on occupancy
  
  // Detection (Observation Process)
  
  // community science records
  
  mu_p_cs_0 ~ normal(0, 2); // global intercept for detection
  
  // species random-effects
  p_cs_species_raw ~ std_normal();
  sigma_p_cs_species ~ normal(0, 2);
  
  // level-2 spatial grouping
  p_cs_site_raw ~ std_normal();
  sigma_p_cs_site ~ normal(0, 0.5); // weakly-informative prior
  // level-3 spatial grouping
  p_cs_level_three_raw ~ std_normal();
  sigma_p_cs_level_three ~ normal(0, 0.5); // weakly-informative prior
  
  // a temporal effect on detection probability
  p_cs_interval ~ normal(0, 2); 
  
  // a population effect on detection probability
  p_cs_pop_density ~ normal(0, 2);
  
  // an income effect on detection probability
  p_cs_income ~ normal(0, 2);
  
  // museum records
  
  mu_p_rc_0 ~ normal(0, 1); // global intercept for detection
  
  // species random-effects
  p_rc_species_raw ~ std_normal();
  sigma_p_rc_species ~ normal(0, 2);
  
  // level-2 spatial grouping
  p_rc_site_raw ~ std_normal();
  sigma_p_rc_site ~ normal(0, 0.5); // weakly-informative prior
  // level-3 spatial grouping
  p_rc_level_three_raw ~ std_normal();
  sigma_p_rc_level_three ~ normal(0, 0.5); // weakly-informative prior
  
  // LIKELIHOOD
  
  // Stan can sample the mean and sd of parameters by summing out the
  // parameter (marginalizing) across likelihood statements
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
        // If the site is in the range of a species, then evaluate lp, otherwise do not (treat as NA).
        if(sum(ranges[i,j,k]) > 0){ // The sum of the NA vector will be == 0 if site is not in range
        
          // if species is detected at the specific site*interval at least once
          // by community science efforts OR research collection records
          // then the species occurs there. lp_observed calculates
          // the probability density that species occurs given psi, plus the 
          // probability density that we did/did not observe it on each visit l in 1:nvisit
          if(sum(V_cs[i, j, k, 1:n_visits]) > 0 || sum(V_rc[i, j, k, 1:n_visits]) > 0) {
            
             // lp_observed:
             target += log_inv_logit(logit_psi[i,j,k]) +
                      binomial_logit_lpmf(sum(V_cs[i,j,k,1:n_visits]) | n_visits, logit_p_cs[i,j,k]) + 
                      // sum(V_rc_NA[i,j,k,1:n_visits]) below tells us how many sampling 
                      // events actually occurred for research collection records (values can range between 0 and the number of years per interval)
                      binomial_logit_lpmf(sum(V_rc[i,j,k,1:n_visits]) | sum(V_rc_NA[i,j,k,1:n_visits]), logit_p_rc[i,j,k]);
                          
          // else the species was never detected at the site*interval
          // lp_unobserved sums the probability density of:
          // 1) species occupies the site*interval but was not detected on each visit, and
          // 2) the species does not occupy the site*interval
          } else {
            
            // lp_unobserved
            target += log_sum_exp(log_inv_logit(logit_psi[i,j,k]) +
                    binomial_logit_lpmf(0 | 
                      n_visits, logit_p_cs[i,j,k]) +
                    // sum(V_rc_NA[i,j,k,1:n_visits]) below tells us how many sampling 
                    // events actually occurred for museum records
                    binomial_logit_lpmf(0 | 
                      sum(V_rc_NA[i,j,k,1:n_visits]), logit_p_rc[i,j,k]),
                    
                    log1m_inv_logit(logit_psi[i,j,k])); 
            
          } // end if/else ever observed
        
        } // end if/ in range
          
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
  
} // end model

generated quantities{
  
  // occurrence of species at each site in each year
  int z_simmed[n_species, n_sites, n_intervals]; // simulate occurrence

  for(i in 1:n_species){
   for(j in 1:n_sites){
     for(k in 1:n_intervals){
          z_simmed[i,j,k] = bernoulli_logit_rng(logit_psi[i,j,k]); 
      }    
    }
  }
  
  //
  // posterior predictive check (number of detections, binned by species)
  //
  int<lower=0> W_species_rep_cs[n_species]; // sum of simulated detections
  int<lower=0> W_species_rep_rc[n_species]; // sum of simulated detections

  // initialize at 0
  for(i in 1:n_species){
    W_species_rep_cs[i] = 0;
    W_species_rep_rc[i] = 0;
  }
      
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
              (z_simmed[i,j,k] *  binomial_rng(sum(V_rc_NA[i,j,k,]), inv_logit(logit_p_cs[i,j,k])));
           
          } // end if{}
           
      } // end loop across years
    } // end loop across sites
  } // end loop across species
  
} // end generated quantities
