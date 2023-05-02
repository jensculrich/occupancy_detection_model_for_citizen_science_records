// multi-species occupancy model for GBIF occurrence data
// jcu, started nov 21, 2022.

functions {
  // covariance matrix for occupancy and detection rates; and detection rates for different data sources
  matrix custom_cov_matrix(vector sigma, real rho) {
    matrix[2,2] Sigma;
    Sigma[1,1] = square(sigma[1]); // species variation in occupancy rates
    Sigma[2,2] = square(sigma[2]); // species variation in cit sci detection rates
    Sigma[1,2] = sigma[1] * sigma[2] * rho; // correlation between occ and cit sci
    Sigma[2,1] = Sigma[1,2];
    return Sigma;
  }
  
  // mean values for multivariate normal distribution (center species on global intercepts)
  vector mu(real mu_p_citsci_0, real mu_p_museum_0){
    vector[2] global_intercepts;
    global_intercepts[1] = mu_p_citsci_0;
    global_intercepts[2] = mu_p_museum_0;
    return global_intercepts;
  }
  
}

data {
  
  int<lower=1> n_species;  // observed species
  int<lower=1> species[n_species]; // vector of species
  int<lower=1> n_genera;  // (number of) genera (level-3 clusters)
  int<lower=1> genus_lookup[n_species]; // level-3 cluster look up vector for level-2 cluster

  int<lower=1> n_sites;  // (number of) sites within region (level-2 clusters)
  int<lower=1, upper=n_sites> sites[n_sites];  // vector of sites
  int<lower=1> n_ecoregion_three;  // (number of) fine-scale (3) ecoregion areas (level-3 clusters)
  int<lower=1> n_ecoregion_one;  // (number of) broad-scale (1) ecoregion areas (level-4 clusters)
  int<lower=1> ecoregion_three_lookup[n_sites]; // level-3 cluster look up vector for level-2 cluster
  int<lower=1> ecoregion_one_lookup[n_ecoregion_three]; // level-4 cluster look up vector for level-3 cluster
  
  int<lower=1> n_intervals;  // intervals during which sites are visited
  
  int intervals[n_intervals]; // vector of intervals (used as covariate data for 
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
  //vector[n_sites] avg_income; // (scaled) developed open surface cover of each site
  vector[n_sites] herb_shrub_forest; // (scaled) undeveloped open surface cover of each site
  real museum_total_records[n_sites, n_intervals]; // (scaled) number of records
  
} // end data


parameters {
  
  // Covararying Parameters
  real<lower=-1,upper=1> rho;  // correlation of (occupancy, cit sci detection)
  vector<lower=0>[2] sigma_species_detection; // sd of (occupancy, detection)
  vector[2] species_intercepts_detection[n_species];  // species-level (occupancy, citsci detection, museum detection)

  // OCCUPANCY
  real mu_psi_0; // global intercept for occupancy
  
  // species specific intercept allows some species to occur at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] psi_species; // species specific intercept for occupancy
  real<lower=0> sigma_psi_species; // variance in species intercepts// Level-3 spatial random effect
  // Level-3 phylogenetic random effect
  vector[n_genera] psi_genus; // site specific intercept for PL outcome
  real<lower=0> sigma_psi_genus; // variance in site intercepts
  
  // Spatially nested random effect on occupancy rates
  // Level-2 spatial random effect
  // site specific intercept allows some sites to be occupied at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all sites.
  vector[n_sites] psi_site; // site specific intercept for occupancy
  real<lower=0> sigma_psi_site; // variance in site intercepts
  // Level-3 spatial random effect
  // site specific intercept allows some sites to have lower success than others, 
  // but with overall estimates for success partially informed by the data pooled across all sites.
  vector[n_ecoregion_three] psi_ecoregion_three; // site specific intercept for PL outcome
  real<lower=0> sigma_psi_ecoregion_three; // variance in site intercepts
  // Level-4 spatial random effect
  // site specific intercept allows some sites to have lower success than others, 
  // but with overall estimates for success partially informed by the data pooled across all sites.
  vector[n_ecoregion_one] psi_ecoregion_one; // site specific intercept for PL outcome
  real<lower=0> sigma_psi_ecoregion_one; // variance in site intercepts
  
  // random slope for species specific natural habitat effects on occupancy
  vector[n_species] psi_herb_shrub_forest; // vector of species specific slope estimates
  real mu_psi_herb_shrub_forest; // community mean of species specific slopes
  real<lower=0> sigma_psi_herb_shrub_forest; // variance in species slopes
  
  // random slope for species specific open low development effects on occupancy
  //vector[n_species] psi_income; // vector of species specific slope estimates
  //real mu_psi_income; // community mean of species specific slopes
  //real<lower=0> sigma_psi_income; // variance in species slopes
  
  // fixed effect of site area on occupancy
  real psi_site_area;
  
  // DETECTION
  
  // citizen science observation process
  real mu_p_citsci_0; // global detection intercept for citizen science records
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  //vector[n_species] p_citsci_species; // species specific intercept for detection
  //real<lower=0> sigma_p_citsci_species; // variance in species intercepts
  
  // random slope for site specific temporal effects on occupancy
  // level-2 spatial clusters
  vector[n_sites] p_citsci_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_citsci_site; // variance in site slopes
  // level-3 spatial clusters
  vector[n_ecoregion_three] p_citsci_ecoregion_three; // site specific intercept for PL outcome
  real<lower=0> sigma_p_citsci_ecoregion_three;
  // level-4 spatial clusters
  vector[n_ecoregion_one] p_citsci_ecoregion_one; // site specific intercept for PL outcome
  real<lower=0> sigma_p_citsci_ecoregion_one; 
  
  real p_citsci_interval; // fixed temporal effect on detection probability
  real p_citsci_pop_density; // fixed effect of population on detection probability
  
  // museum records observation process
  real mu_p_museum_0; // global detection intercept for citizen science records
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  //vector[n_species] p_museum_species; // species specific intercept for detection
  //real<lower=0> sigma_p_museum_species; // variance in species intercepts
  
  // random slope for site specific temporal effects on occupancy
  // level-2 spatial clusters
  vector[n_sites] p_museum_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_museum_site; // variance in site slopes
  // level-3 spatial clusters
  vector[n_ecoregion_three] p_museum_ecoregion_three; // site specific intercept for PL outcome
  real<lower=0> sigma_p_museum_ecoregion_three; 
  // level-4 spatial clusters
  vector[n_ecoregion_one] p_museum_ecoregion_one; // site specific intercept for PL outcome
  real<lower=0> sigma_p_museum_ecoregion_one;
  
  real p_museum_total_records; // fixed effect of total records on detection probability
  
} // end parameters


transformed parameters {
  
  real logit_psi[n_species, n_sites, n_intervals];  // odds of occurrence
  real logit_p_citsci[n_species, n_sites, n_intervals]; // odds of detection by cit science
  real logit_p_museum[n_species, n_sites, n_intervals]; // odds of detection by museum
  
  // spatially nested intercepts
  vector[n_sites] psi0_site;
  vector[n_ecoregion_three] psi0_ecoregion_three;
  vector[n_ecoregion_one] psi0_ecorgion_one;
  
  vector[n_sites] p0_citsci_site;
  vector[n_ecoregion_three] p0_citsci_ecoregion_three;
  vector[n_ecoregion_one] p0_citsci_ecoregion_one;
  
  vector[n_sites] p0_museum_site;
  vector[n_ecoregion_three] p0_museum_ecoregion_three;
  vector[n_ecoregion_one] p0_museum_ecoregion_one;
  
  // phylogenetically nested intercepts
  real psi0_species[n_species];
  
  // compute the varying intercept at the ecoregion1 level
  // Level-4 (n_ecoregion_one level-4 random intercepts)
  psi0_ecorgion_one = sigma_psi_ecoregion_one*psi_ecoregion_one;

  // compute the varying intercept at the ecoregion3 level
  // Level-3 (n_ecoregion_three level-3 random intercepts, nested in ecoregion1)
  psi0_ecoregion_three = psi0_ecorgion_one[ecoregion_one_lookup] + sigma_psi_ecoregion_three*psi_ecoregion_three;

  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3, nested in ecoregion1)
  psi0_site = psi0_ecoregion_three[ecoregion_three_lookup] + sigma_psi_site*psi_site;
  
  //
  // Nested spatial intercept for detection (including global intercept mu)
  // compute the varying occurrence intercept at the ecoregion1 level
  // Level-4 (n_ecoregion_one level-4 random intercepts) vectorized
  p0_citsci_ecoregion_one = p_citsci_ecoregion_one*sigma_p_citsci_ecoregion_one;
  
  // compute the varying citsci detection intercept at the ecoregion3 level
  // Level-3 (n_ecoregion_three level-3 random intercepts)
  p0_citsci_ecoregion_three = p0_citsci_ecoregion_one[ecoregion_one_lookup] +
      p_citsci_ecoregion_three*sigma_p_citsci_ecoregion_three;

  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  p0_citsci_site = p0_citsci_ecoregion_three[ecoregion_three_lookup] + 
      p_citsci_site*sigma_p_citsci_site;
  
  //
  // Nested spatial intercept for detection (including global intercept mu)
  // compute the varying occurrence intercept at the ecoregion1 level
  // Level-4 (n_ecoregion_one level-4 random intercepts) vectorized
  p0_museum_ecoregion_one = p_museum_ecoregion_one*sigma_p_museum_ecoregion_one;
  
  // compute the varying citsci detection intercept at the ecoregion3 level
  // Level-3 (n_ecoregion_three level-3 random intercepts)
  p0_museum_ecoregion_three = p0_museum_ecoregion_one[ecoregion_one_lookup] +
      p_museum_ecoregion_three*sigma_p_museum_ecoregion_three;

  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  p0_museum_site = p0_museum_ecoregion_three[ecoregion_three_lookup] + 
    p_museum_site*sigma_p_museum_site;
    
  // Phylogenetic clustering for occurrence
  // compute the varying intercept at the level-2 species level
  // by clustering within Level-3 (n_genera level-3 random intercepts)
  for(i in 1:n_species){
    psi0_species[i] = psi_genus[genus_lookup[i]] + psi_species[i];
  }
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals  
          
          logit_psi[i,j,k] = // the inverse of the log odds of occurrence is equal to..
            psi0_species[species[i]] + // a phylogenetically nested, species-specific intercept
            psi0_site[sites[j]] + // a spatially nested, site-specific intercept
            psi_herb_shrub_forest[species[i]]*herb_shrub_forest[j] + // an effect 
            //mu_psi_income*avg_income[j] +
            //psi_income[species[i]]*avg_income[j] + // an effect
            psi_site_area*site_areas[j] // an effect of spatial area of the site on occurrence
            ; // end psi[i,j,k]
            
      } // end loop across all intervals
    } // end loop across all sites
  }  // end loop across all species
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
          logit_p_citsci[i,j,k] = // the inverse of the log odds of detection is equal to..
            species_intercepts_detection[species[i],1] + // a species specific intercept
            p0_citsci_site[sites[j]] + // a spatially specific intercept
            p_citsci_interval*(intervals[k]^2) + // an overall effect of time on detection
            p_citsci_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_citsci[i,j,k]
           
          logit_p_museum[i,j,k] = // the inverse of the log odds of detection is equal to..
            species_intercepts_detection[species[i],2] + // a species specific intercept
            p0_museum_site[sites[j]] + // a spatially specific intercept
            p_museum_total_records*museum_total_records[j,k] //records at site in interval
           ; // end p_museum[i,j,k]
           
      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
             
  
} // end transformed parameters


model {
  
  // PRIORS
  
  // correlated species effects
  sigma_species_detection[1] ~ normal(0, 1.5);
  sigma_species_detection[2] ~ normal(0, 1.5);
  (rho + 1) / 2 ~ beta(2, 2);
  
  species_intercepts_detection ~ multi_normal(mu(mu_p_citsci_0, mu_p_museum_0), 
    custom_cov_matrix(sigma_species_detection, rho));
    
  // Occupancy (Ecological Process)
  mu_psi_0 ~ normal(0, 0.25); // global intercept for occupancy rate
  
  // https://betanalpha.github.io/assets/case_studies/divergences_and_bias.html#3_a_non-centered_eight_schools_implementation
  // level-2 spatial grouping
  psi_site  ~ normal(0, 0.25);
  sigma_psi_site ~ normal(0, 1); // weakly-informative prior
  // level-3 spatial grouping
  psi_ecoregion_three ~ normal(0, 0.25);
  sigma_psi_ecoregion_three ~ normal(0, 1); // weakly-informative prior
  // level-4 spatial grouping
  psi_ecoregion_one ~ normal(0, 0.25);
  sigma_psi_ecoregion_one ~ normal(0, 0.5); // weakly-informative prior
  
  // level-2 phylogenetic grouping
  psi_species ~ normal(0, sigma_psi_species); 
  sigma_psi_species ~ normal(0, 0.5); // weakly-informative prior
  // level-3 phylogenetic grouping
  psi_genus ~ normal(mu_psi_0, sigma_psi_genus); 
  sigma_psi_genus ~ normal(0, 0.5); // weakly-informative prior
  
  psi_herb_shrub_forest ~ normal(mu_psi_herb_shrub_forest, sigma_psi_herb_shrub_forest);
  mu_psi_herb_shrub_forest ~ normal(0, 2); // community mean
  sigma_psi_herb_shrub_forest ~ normal(0.75, 0.1); // community variance
  
  //psi_income ~ normal(mu_psi_income, sigma_psi_income);
  //mu_psi_income ~ normal(0, 2); // community mean
  //plotsigma_psi_income ~ normal(0, 0.1); // community variance
  
  psi_site_area ~ normal(0, 2); // effect of site area on occupancy
  
  // Detection (Observation Process)
  
  // citizen science records
  
  mu_p_citsci_0 ~ normal(0, 1); // global intercept for detection
  
  // level-2 spatial grouping
  p_citsci_site  ~ normal(0, 0.25);
  sigma_p_citsci_site ~ normal(0, 0.25); // weakly-informative prior
  // level-3 spatial grouping
  p_citsci_ecoregion_three ~ normal(0, 0.25);
  sigma_p_citsci_ecoregion_three ~ normal(0, 0.25); // weakly-informative prior
  // level-4 spatial grouping
  p_citsci_ecoregion_one ~ normal(0, 0.25);
  sigma_p_citsci_ecoregion_one ~ normal(0, 0.25); // weakly-informative prior
  
  //p_citsci_species ~ normal(mu_p_citsci_0, sigma_p_citsci_species); 
  //sigma_p_citsci_species ~ normal(0,1);
  
  // a temporal effect on detection probability
  p_citsci_interval ~ normal(0, 2); 
  
  // a population effect on detection probability
  p_citsci_pop_density ~ normal(0, 2);
  
  // museum records
  
  mu_p_museum_0 ~ normal(0, 1); // global intercept for detection
  
  // level-2 spatial grouping
  p_museum_site  ~ normal(0, 0.1);
  sigma_p_museum_site ~ normal(0, 0.1); // weakly-informative prior
  // level-3 spatial grouping
  p_museum_ecoregion_three ~ normal(0, 0.1);
  sigma_p_museum_ecoregion_three ~ normal(0, 0.1); // weakly-informative prior
  // level-4 spatial grouping
  p_museum_ecoregion_one ~ normal(0, 0.1);
  sigma_p_museum_ecoregion_one ~ normal(0, 0.1); // weakly-informative prior
  
  //p_museum_species ~ normal(mu_p_museum_0, sigma_p_museum_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  //sigma_p_museum_species ~ normal(0, 0.25);
  
  // an effect of total records at the site during the interval
  p_museum_total_records ~ normal(0, 0.5);
  
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
  
  real T_rep_citsci[n_species]; // Freeman-Tukey distance from eval (species bin)
  real T_obs_citsci[n_species]; // Freeman-Tukey distance from eval (species bin)
  real T_rep_museum[n_species]; // Freeman-Tukey distance from eval (species bin)
  real T_obs_museum[n_species]; // Freeman-Tukey distance from eval (species bin)
  
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
