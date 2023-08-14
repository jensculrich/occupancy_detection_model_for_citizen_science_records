// multi-species occupancy model for BBNA occurrence data
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
  real rc_total_records[n_sites, n_intervals]; // (scaled) number of records
  
} // end data


parameters {
  
  // Covararying Parameters
  real<lower=-1,upper=1> rho;  // correlation of (community science and research collections detection)
  vector<lower=0>[2] sigma_species_detection; // variance in species-specific detection rates 
    // (sigma_species_detection[1] == community science and sigma_species_detection[2] == research collections detection)
  vector[2] species_intercepts_detection[n_species];// species-level detection intercepts
    // (species_intercepts_detection[1] == community science and species_intercepts_detection[2] == research collections detection)

  // OCCUPANCY
  real mu_psi_0; // global intercept for occupancy
  
  // species specific intercept allows some species to occur at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] psi_species; // species specific intercept for occupancy
  real<lower=0> sigma_psi_species; // variance in species intercepts
  
  // Spatially nested random effect on occupancy rates
  // Level-2 spatial random effect
  vector[n_sites] psi_site; // site specific intercept for occupancy
  real<lower=0> sigma_psi_site; // variance in site intercepts
  // Level-3 spatial random effect
  vector[n_level_three] psi_level_three; // level-three specific intercept for PL outcome
  real<lower=0> sigma_psi_level_three; // variance in level-three intercepts
  // Level-4 spatial random effect
  vector[n_level_four] psi_level_four; // level-four specific intercept for PL outcome
  real<lower=0> sigma_psi_level_four; // variance in level-four intercepts
  
  // random slope for species specific natural habitat effects on occupancy
  vector[n_species] psi_natural_habitat; // vector of species specific slope estimates
  real mu_psi_natural_habitat; // community mean of species specific slopes
  real<lower=0> sigma_psi_natural_habitat; // variance in species slopes
  
  // random slope for species specific household income effects on occupancy
  vector[n_species] psi_income; // vector of species specific slope estimates
  real mu_psi_income; // community mean of species specific slopes
  real<lower=0> sigma_psi_income; // variance in species slopes
  
  // fixed effect of site area on occupancy
  real psi_site_area;
  
  // DETECTION
  
  // citizen science observation process
  real mu_p_cs_0; // global detection intercept for citizen science records
  
  // random slope for site specific temporal effects on occupancy
  // level-2 spatial clusters
  vector[n_sites] p_cs_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_cs_site; // variance in site slopes
  // level-3 spatial clusters
  vector[n_level_three] p_cs_level_three; // level-three specific intercept for cs detection
  real<lower=0> sigma_p_cs_level_three;  // variance in level-three slopes
  // level-4 spatial clusters
  vector[n_level_four] p_cs_level_four; // level-four specific intercept for cs detection
  real<lower=0> sigma_p_cs_level_four;  // variance in level-four slopes
  
  real p_cs_interval; // fixed temporal effect on cs detection probability
  real p_cs_pop_density; // fixed effect of population on cs detection probability
  
  // research collections records observation process
  real mu_p_rc_0; // global detection intercept for research collections records
  
  // random slope for site specific temporal effects on occupancy
  // level-2 spatial clusters
  vector[n_sites] p_rc_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_rc_site; // variance in site slopes
  // level-3 spatial clusters
  vector[n_level_three] p_rc_level_three; // level-three specific intercept for rc detection
  real<lower=0> sigma_p_rc_level_three; // variance in level-three slopes
  // level-4 spatial clusters
  vector[n_level_four] p_rc_level_four; // level-four specific intercept for rc detection
  real<lower=0> sigma_p_rc_level_four; // variance in level-four slopes
  
  real p_rc_total_records; // fixed effect of total records on rc detection probability
  
} // end parameters


transformed parameters {
  
  real logit_psi[n_species, n_sites, n_intervals];  // odds of occurrence
  real logit_p_cs[n_species, n_sites, n_intervals]; // odds of detection by community science
  real logit_p_rc[n_species, n_sites, n_intervals]; // odds of detection by research collections
  
  // spatially nested intercepts
  real psi0_site[n_sites];
  real psi0_level_three[n_level_three];

  real p0_cs_site[n_sites];
  real p0_cs_level_three[n_level_three];
  
  real p0_rc_site[n_sites];
  real p0_rc_level_three[n_level_three];

  //
  // compute the varying community science detection intercept at the ecoregion3 level
  // Level-3 (n_level_three level-3 random intercepts)
  for(i in 1:n_level_three){
    psi0_level_three[i] = psi_level_four[level_four_lookup[i]] + 
      psi_level_three[i];
  } 

  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  for(i in 1:n_sites){
    psi0_site[i] = psi0_level_three[level_three_lookup[i]] + 
      psi_site[i];
  }
  
  //
  // compute the varying community science detection intercept at the ecoregion3 level
  // Level-3 (n_level_three level-3 random intercepts)
  for(i in 1:n_level_three){
    p0_cs_level_three[i] = p_cs_level_four[level_four_lookup[i]] + 
      p_cs_level_three[i];
  } 

  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  for(i in 1:n_sites){
    p0_cs_site[i] = p0_cs_level_three[level_three_lookup[i]] + 
      p_cs_site[i];
  }
  
  //
  // compute the varying community science detection intercept at the ecoregion3 level
  // Level-3 (n_level_three level-3 random intercepts)
  for(i in 1:n_level_three){
    p0_rc_level_three[i] = p_rc_level_four[level_four_lookup[i]] + 
      p_rc_level_three[i];
  } 

  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  for(i in 1:n_sites){
    p0_rc_site[i] = p0_rc_level_three[level_three_lookup[i]] + 
      p_rc_site[i];
  } 
  
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals  
          
          logit_psi[i,j,k] = // the inverse of the log odds of occurrence is equal to..
            psi_species[species[i]] + // a species specific intercept
            psi0_site[sites[j]] + // a spatially nested, site-specific intercept
            psi_natural_habitat[species[i]]*natural_habitat[j] + // an effect 
            psi_income[species[i]]*avg_income[j] + // an effect
            psi_site_area*site_areas[j] // an effect of spatial area of the site on occurrence
            ; // end psi[i,j,k]
            
      } // end loop across all intervals
    } // end loop across all sites
  }  // end loop across all species
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
          logit_p_cs[i,j,k] = // the inverse of the log odds of detection is equal to..
            species_intercepts_detection[species[i],1] + // a species specific intercept
            p0_cs_site[sites[j]] + // a spatially specific intercept // includes global intercept
            p_cs_interval*(intervals[k]^2) + // an overall effect of time on detection
            p_cs_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_cs[i,j,k]
           
          logit_p_rc[i,j,k] = // the inverse of the log odds of detection is equal to..
            species_intercepts_detection[species[i],2] + // a species specific intercept
            p0_rc_site[sites[j]] + // a spatially specific intercept // includes global intercept
            p_rc_total_records*rc_total_records[j,k] //records at site in interval
           ; // end p_rc[i,j,k]
           
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
             
  
} // end transformed parameters


model {
  
  // PRIORS
  
  // correlated species effects
  sigma_species_detection[1] ~ normal(0, 2);
  sigma_species_detection[2] ~ normal(0, 2);
  (rho + 1) / 2 ~ beta(2, 2);
  
  species_intercepts_detection ~ multi_normal(mu(mu_p_cs_0, mu_p_rc_0), 
    custom_cov_matrix(sigma_species_detection, rho));
  
  // Occupancy (Ecological Process)
  mu_psi_0 ~ normal(0, 1); // global intercept for occupancy rate
  
  // level-2 spatial grouping
  psi_site  ~ normal(0, sigma_psi_site);
  sigma_psi_site ~ normal(0, 0.5); // weakly-informative prior
  // level-3 spatial grouping
  psi_level_three ~ normal(0, sigma_psi_level_three);
  sigma_psi_level_three ~ normal(0, 0.5); // weakly-informative prior
  // level-5 spatial grouping
  psi_level_four ~ normal(0, sigma_psi_level_four);
  sigma_psi_level_four ~ normal(0, 0.5); // weakly-informative prior
  
  psi_species ~ normal(mu_psi_0, sigma_psi_species); 
  sigma_psi_species ~ normal(0, 1); // weakly-informative prior
  
  psi_natural_habitat ~ normal(mu_psi_natural_habitat, sigma_psi_natural_habitat);
  mu_psi_natural_habitat ~ normal(0, 2); // community mean
  sigma_psi_natural_habitat ~ normal(0, 1); // community variance
  
  psi_income ~ normal(mu_psi_income, sigma_psi_income);
  mu_psi_income ~ normal(0, 2); // community mean
  sigma_psi_income ~ normal(0, 1); // community variance
  
  psi_site_area ~ normal(0, 2); // effect of site area on occupancy
  
  // Detection (Observation Process)
  
  // community science records
  
  mu_p_cs_0 ~ normal(0, 2); // global intercept for detection
  
  // level-2 spatial grouping
  p_cs_site  ~ normal(0, sigma_p_cs_site);
  sigma_p_cs_site ~ normal(0, 0.25); // weakly-informative prior
  // level-3 spatial grouping
  p_cs_level_three ~ normal(0, sigma_p_cs_level_three);
  sigma_p_cs_level_three ~ normal(0, 0.25); // weakly-informative prior
  // level-5 spatial grouping
  p_cs_level_four ~ normal(0, sigma_p_cs_level_four);
  sigma_p_cs_level_four ~ normal(0, 0.25); // weakly-informative prior

  // a temporal effect on detection probability
  p_cs_interval ~ normal(0, 2); 
  
  // a population effect on detection probability
  p_cs_pop_density ~ normal(0, 2);
  
  // museum records
  
  mu_p_rc_0 ~ normal(0, 0.5); // global intercept for detection
  
  // level-2 spatial grouping
  p_rc_site  ~ normal(0, sigma_p_rc_site);
  sigma_p_rc_site ~ normal(0, 0.25); // weakly-informative prior
  // level-3 spatial grouping
  p_rc_level_three ~ normal(0, sigma_p_rc_level_three);
  sigma_p_rc_level_three ~ normal(0, 0.25); // weakly-informative prior
  // level-4 spatial grouping
  p_rc_level_four ~ normal(0, sigma_p_rc_level_four);
  sigma_p_rc_level_four ~ normal(0, 0.25); // weakly-informative prior
  
  // an effect of total records at the site during the interval
  p_rc_total_records ~ normal(0, 0.5);
  
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
          if(sum(V_cs[i, j, k, 1:n_visits]) > 0 || sum(V_rc[i, j, k, 1:n_visits]) > 0) {
            
             // lp_observed:
             target += log_inv_logit(logit_psi[i,j,k]) +
                      binomial_logit_lpmf(sum(V_cs[i,j,k,1:n_visits]) | n_visits, logit_p_cs[i,j,k]) + 
                      // sum(V_rc_NA[i,j,k,1:n_visits]) below tells us how many sampling 
                      // events actually occurred for museum records
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
  
  int Z[n_species, n_sites, n_intervals];
  
  int z_rep[n_species, n_sites, n_intervals];
  int y_rep_cs[n_species, n_sites, n_intervals, n_visits]; // repd detections
  int y_rep_rc[n_species, n_sites, n_intervals, n_visits]; // repd detections
  
  real eval_cs[n_species,n_sites,n_intervals,n_visits]; // expected values
  real eval_rc[n_species,n_sites,n_intervals,n_visits]; // expected values
  
  real T_rep_cs[n_species]; // Freeman-Tukey distance from eval (species bin)
  real T_obs_cs[n_species]; // Freeman-Tukey distance from eval (species bin)
  real T_rep_rc[n_species]; // Freeman-Tukey distance from eval (species bin)
  real T_obs_rc[n_species]; // Freeman-Tukey distance from eval (species bin)
  
  real P_species_cs[n_species]; // P-value by species
  real P_species_rc[n_species]; // P-value by species
  
  // Initialize T_rep and T_obs and P-values
  for(i in 1:n_species){
    
    T_rep_cs[i] = 0;
    T_obs_cs[i] = 0;
    T_rep_rc[i] = 0;
    T_obs_rc[i] = 0;
    
    P_species_cs[i] = 0;
    P_species_rc[i] = 0;
    
  }
      
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
      
        if(sum(ranges[i,j,k]) > 0){ // The sum of the NA vector will be == 0 if site is not in range
          
          // if occupancy state is certain then the expected occupancy is 1
          if(sum(V_cs[i, j, k, 1:n_visits]) > 0 || sum(V_rc[i, j, k, 1:n_visits]) > 0) {
          
            Z[i,j,k] = 1;
          
          // else the site could be occupied or not
          } else {
            
            // occupancy but never observed by either dataset
            real ulo = inv_logit(logit_psi[i,j,k]) * 
              ((1 - inv_logit(logit_p_cs[i,j,k]))^n_visits + 
              (1 - inv_logit(logit_p_rc[i,j,k]))^n_visits);
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
          eval_cs[i,j,k,l] = Z[i,j,k] * 
            bernoulli_logit_rng(logit_p_cs[i,j,k]);
          eval_rc[i,j,k,l] = Z[i,j,k] * 
            bernoulli_logit_rng(logit_p_rc[i,j,k]);
          
          // occupancy in replicated data
          // should evaluate to zero if the site is not in range
          z_rep[i,j,k] = bernoulli_logit_rng(logit_psi[i,j,k]); 
          if(sum(ranges[i,j,k]) > 0){ // The sum of the NA vector will be == 0 if site is not in range
            z_rep[i,j,k] = z_rep[i,j,k];
          } else {
            z_rep[i,j,k] = 0;
          }
          
          // detections in replicated data
          y_rep_cs[i,j,k,l] = z_rep[i,j,k] * bernoulli_logit_rng(logit_p_cs[i,j,k]);
          y_rep_rc[i,j,k,l] = z_rep[i,j,k] * bernoulli_logit_rng(logit_p_rc[i,j,k]);
          
          // Compute fit statistic (Tukey-Freeman) for replicate data
          // Citizen science records
          // Binned by species
          T_rep_cs[i] = T_rep_cs[i] + (sqrt(y_rep_cs[i,j,k,l]) - 
            sqrt(eval_cs[i,j,k,l]))^2;
          // Compute fit statistic (Tukey-Freeman) for real data
          // Binned by species
          T_obs_cs[i] = T_obs_cs[i] + (sqrt(V_cs[i,j,k,l]) - 
            sqrt(eval_cs[i,j,k,l]))^2;
          
          // Compute fit statistic (Tukey-Freeman) for replicate data
          // Binned by species
          // Museum records
          if(V_rc_NA[i,j,k,l] == 0){
            
            T_rep_rc[i] = T_rep_rc[i] + 0;
            // Compute fit statistic (Tukey-Freeman) for real data
            // Binned by species
            T_obs_rc[i] = T_obs_rc[i] + 0;
          
          } else{
            
            T_rep_rc[i] = T_rep_rc[i] + (sqrt(y_rep_rc[i,j,k,l]) - 
              sqrt(eval_rc[i,j,k,l]))^2;
            // Compute fit statistic (Tukey-Freeman) for real data
            // Binned by species
            T_obs_rc[i] = T_obs_rc[i] + (sqrt(V_rc[i,j,k,l]) - 
              sqrt(eval_rc[i,j,k,l]))^2;
            
          }
          
        } // end loop across visits
      } // end loop across intervals
    } // end loop across sites
  } // end loop across species
  
  // bin by species
  for(i in 1:n_species) { // loop across all species
    
    // if the discrepancy is lower for the real data for the species
    // versus the replicated data
    if(T_obs_cs[i] < T_rep_cs[i]){
      
      // then increase species P by 1      
      P_species_cs[i] = P_species_cs[i] + 1;
      // the ppc will involve averaging P across the number of post-burnin iterations
            
    }
    
    // if the discrepancy is lower for the real data for the species
    // versus the replicated data
    if(T_obs_rc[i] < T_rep_rc[i]){
      
      // then increase species P by 1      
      P_species_rc[i] = P_species_rc[i] + 1;
      // the ppc will involve averaging P across the number of post-burnin iterations
            
    }
    
  }

}
