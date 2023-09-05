// multi-species occupancy model for GBIF occurrence data
// will include income as a predictor, but will also include a simplified random effects
// structure to see if this helps with non-convergence that occurred when using the original model
// jcu, started nov 21, 2022.

data {
  
  int<lower=1> n_species;  // observed species
  int<lower=1> species[n_species]; // vector of species identities
  
  // ended up not including genus clustering (too many rand effects to estimate)
  //int<lower=1> n_genera;  // (number of) genera (level-3 clusters)
  //int<lower=1> genus_lookup[n_species]; // level-3 cluster look up vector for level-2 cluster

  int<lower=1> n_sites;  // (number of) sites within region (level-2 clusters)
  int<lower=1, upper=n_sites> sites[n_sites];  // vector of sites identities (level-2 clusters)
  int<lower=1> n_level_three;  // (number of) fine-scale (3) ecoregion areas (level-3 clusters)
  int<lower=1> n_level_four;  // (number of) broad-scale (1) ecoregion areas (level-4 clusters)
  int<lower=1> level_three_lookup[n_sites]; // level-3 cluster look up vector for level-2 cluster
  int<lower=1> level_four_lookup[n_level_three]; // level-4 cluster look up vector for level-3 cluster
  
  int<lower=1> n_intervals;  // intervals during which sites are visited
  
  int intervals[n_intervals]; // vector of intervals (used as covariate data for 
                                // species specific effect of occupancy interval (time) on detection)
                                
  int<lower=1> n_visits; // visits within intervals
  
  int<lower=0> V_cs[n_species, n_sites, n_intervals, n_visits];  // visits l when species i was detected at site j on interval k
  int<lower=0> V_rc[n_species, n_sites, n_intervals, n_visits];  // visits l when species i was detected at site j on interval k
  
  int<lower=0> ranges[n_species, n_sites, n_intervals, n_visits];  // NA indicator where 1 == site is in range, 0 == not in range

  vector[n_sites] site_areas; // (scaled) spatial area extent of each site
  vector[n_sites] pop_densities; // (scaled) population density of each site
  vector[n_sites] natural_habitat; // (scaled) undeveloped open surface cover of each site
  vector[n_sites] avg_income; // (scaled) household income of each site (with income relativized to state average income)
  vector[n_species] nativity; // nativity vector, 0 == non-native, 1 == native
  
} // end data


parameters {
  
  // OCCUPANCY
  
  real mu_psi_0; // global intercept for occupancy
  
  // species-specific intercepts allow some species to occur at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] psi_species; // species specific intercept for occupancy
  real<lower=0> sigma_psi_species; // variance in species intercepts// Level-3 spatial random effect
  // Level-3 phylogenetic random effect
  // for hoverflies, I also considered nesting species specific random effects within genus specific re's 
  //vector[n_genera] psi_genus; // site specific intercept for occupancy
  //real<lower=0> sigma_psi_genus; // variance in site intercepts
  
  // spatially nested random effect on occupancy rates
  // Level-2 spatial random effect
  vector[n_sites] psi_site; // site specific intercept for occupancy
  real<lower=0> sigma_psi_site; // variance in site intercepts
  // Level-3 spatial random effect
  vector[n_level_three] psi_level_three; // level-three intercept for occupancy
  real<lower=0> sigma_psi_level_three; // variance in level-three intercepts
  
  // random slope for species specific natural habitat effects on occupancy
  vector[n_species] psi_natural_habitat; // vector of species specific slope estimates
  //vector[n_species] mu_psi_natural_habitat; // community mean of species specific slopes
  //real<lower=0> sigma_psi_natural_habitat; // variance in species slopes
  real delta0; // baseline effect (mean)
  real delta1; // effect of being native on the expected value of the random effect
  real<lower=0> gamma0; // baseline effect (variance) (negative variance not possible)
  real gamma1; // effect of being native on the expected value of the random effect
  
  // random slope for species-specific income effects on occupancy
  //vector[n_species] psi_income; // vector of species specific slope estimates
  real mu_psi_income; // community mean of species specific slopes
  //real<lower=0> sigma_psi_income; // variance in species slopes
    
  // fixed effect of site area on occupancy
  real psi_site_area;
  
  // DETECTION
  
  // community science observation process
  
  real mu_p_cs_0; // global detection intercept for community science records
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[n_species] p_cs_species; // species specific intercept for detection
  real<lower=0> sigma_p_cs_species; // variance in species intercepts
  
  // random slope for site specific temporal effects on occupancy
  // level-2 spatial clusters
  vector[n_sites] p_cs_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_cs_site; // variance in site slopes
  // level-3 spatial clusters
  vector[n_level_three] p_cs_level_three; // level three intercept for cs detection
  real<lower=0> sigma_p_cs_level_three; // variance in level three slopes
  
  real p_cs_interval; // fixed temporal effect on detection probability
  real p_cs_pop_density; // fixed effect of population on detection probability
  
} // end parameters


transformed parameters {
  
  real logit_psi[n_species, n_sites, n_intervals];  // odds of occurrence
  real logit_p_cs[n_species, n_sites, n_intervals]; // odds of detection by community science

  vector[n_species] mu_psi_natural_habitat; // expected value for species specific slopes
  vector[n_species] sigma_psi_natural_habitat; // expected value for variance among species slopes
  
  // spatially nested intercepts
  real psi0_site[n_sites];

  real p0_cs_site[n_sites];

  // phylogenetically nested intercepts
  //real psi0_species[n_species];
  
  // hard prior to disallow intercept plus nativity adjustment from being negative
  real<lower=0> gamma0_plus_gamma1;
  gamma0_plus_gamma1 = gamma0 + gamma1;

  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  for(i in 1:n_sites){
    psi0_site[i] = psi_level_three[level_three_lookup[i]] + 
      psi_site[i];
  }

  // compute varying intercept at the site level
  // Level-2 (n_sites level-2 random intercepts, nested in ecoregion3)
  for(i in 1:n_sites){
    p0_cs_site[i] = p_cs_level_three[level_three_lookup[i]] + 
      p_cs_site[i];
  } 
  
  // Phylogenetic clustering for occurrence
  // compute the varying intercept at the level-2 species level
  // by clustering within Level-3 (n_genera level-3 random intercepts)
  //for(i in 1:n_species){
  //  psi0_species[i] = psi_genus[genus_lookup[i]] + psi_species[i];
  //}
  
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
            psi_species[species[i]] + // a phylogenetically nested, species-specific intercept
            psi0_site[sites[j]] + // a spatially nested, site-specific intercept
            psi_natural_habitat[species[i]]*natural_habitat[j] + // a species-specific effect of natural habitat area
            mu_psi_income*avg_income[j] + // a species-specific effect of income
            psi_site_area*site_areas[j] // an effect of spatial area of the site on occurrence
            ; // end psi[i,j,k]
            
      } // end loop across all intervals
    } // end loop across all sites
  }  // end loop across all species
  
  for (i in 1:n_species){   // loop across all species
    for (j in 1:n_sites){    // loop across all sites
      for(k in 1:n_intervals){ // loop across all intervals
        
          logit_p_cs[i,j,k] = // the inverse of the log odds of detection is equal to..
            p_cs_species[species[i]] + // a species-specific intercept
            p0_cs_site[sites[j]] + // a spatially specific intercept
            p_cs_interval*(intervals[k]^2) + // an overall effect of time on detection
            p_cs_pop_density*pop_densities[j] // an overall effect of pop density on detection
           ; // end p_cs[i,j,k]

      } // end loop across all intervals
    } // end loop across all sites
  } // end loop across all species
             
  
} // end transformed parameters


model {
  
  // PRIORS
    
  // Occupancy (Ecological Process)
  mu_psi_0 ~ normal(0, 0.25); // global intercept for occupancy rate
  
  // level-2 spatial grouping
  psi_site  ~ normal(0, sigma_psi_site);
  sigma_psi_site ~ normal(0, 1); // weakly-informative prior
  // level-3 spatial grouping
  psi_level_three ~ normal(0, sigma_psi_level_three);
  sigma_psi_level_three ~ normal(0, 0.5); // weakly-informative prior
  
  // level-2 phylogenetic grouping
  psi_species ~ normal(mu_psi_0, sigma_psi_species); 
  //psi_species ~ normal(0, sigma_psi_species); 
  sigma_psi_species ~ normal(0, 1); // weakly-informative prior
  // level-3 phylogenetic grouping
  //psi_genus ~ normal(mu_psi_0, sigma_psi_genus); 
  //sigma_psi_genus ~ normal(0, 0.05); // weakly-informative prior
  
  psi_natural_habitat ~ normal(mu_psi_natural_habitat, sigma_psi_natural_habitat);
  // community effect (mu) and variation among species (sigma) is defined as a vector 
  // with intercept delta0 and an effect of nativity (delta1) on community mean
  // and intercept gamma0 and an effect of nativity (gamma1) on variation
  delta0 ~ normal(0, 1); // community mean
  delta1 ~ normal(0, 2); // effect of nativity
  gamma0 ~ normal(0, 1); // community mean
  gamma1 ~ normal(0, 0.25); // effect of nativity
  
  mu_psi_income ~ normal(0, 2); // effect of income on occupancy
  psi_site_area ~ normal(0, 2); // effect of site area on occupancy
  
  // Detection (Observation Process)
  
  // community science records
  mu_p_cs_0 ~ normal(-2, 0.25); // global intercept for detection
  // we used an offset intercept given the low detection rates in early years at time = 0.
  
  p_cs_species ~ normal(mu_p_cs_0, sigma_p_cs_species); 
  sigma_p_cs_species ~ normal(0, 1); // weakly-informative prior
  
  // level-2 spatial grouping
  p_cs_site  ~ normal(0, sigma_p_cs_site);
  sigma_p_cs_site ~ normal(0, 0.5); // weakly-informative prior
  // level-3 spatial grouping
  p_cs_level_three ~ normal(0, sigma_p_cs_level_three);
  sigma_p_cs_level_three ~ normal(0, 0.25); // weakly-informative prior
  
  // a temporal effect on detection probability
  p_cs_interval ~ normal(0, 2); 
  
  // a population effect on detection probability
  p_cs_pop_density ~ normal(0, 2);
  
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
          // probability density that we did/did not observe it on each visit l in 1:nvisit.
          // Even though we don't estimate research collection parameters here (fully integrated model),
          // we still use rc detections to anchor and guide the likelihood function.
          // that is, if we see the species using the research collections but not the comm sci, 
          // we still know that the species is present but that our comm sci detection process is missing it.
          if(sum(V_cs[i, j, k, 1:n_visits]) > 0 || sum(V_rc[i, j, k, 1:n_visits]) > 0) {
            
             // lp_observed:
             target += log_inv_logit(logit_psi[i,j,k]) +
                      binomial_logit_lpmf(sum(V_cs[i,j,k,1:n_visits]) | n_visits, logit_p_cs[i,j,k]); 

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
                    binomial_logit_lpmf(0 | n_visits, logit_p_cs[i,j,k]),
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
  // By doing this in each step of the MCMC we can get a distribution of outcomes 
  // (propagating our uncertainty in delta0 and delta1 to an uncertainty in the group level effects)
  
  real mu_psi_natural_habitat_native;
  mu_psi_natural_habitat_native = delta0 + delta1*1;
  
  real mu_psi_natural_habitat_nonnative;
  mu_psi_natural_habitat_nonnative = delta0 + delta1*0;
  
  real mu_psi_natural_habitat_all_species;
  mu_psi_natural_habitat_all_species = mean(mu_psi_natural_habitat);
  
  // posterior predictive check (Freeman-Tukey posterior pred check, binned by species)
  
  // estimate expected values (occurrence is a partially latent variable so when we don't observe the species, we don't actually know the expected values)
  // create replicated data using the paramter estimates and stochastic process defined by the model
  // gather real detecions
  
  // test the rate at which the number of detections in the real data versus the repped data 
  // are closer to the expected values. The FTP is the rate at which the real data are closer.
  
  int Z[n_species, n_sites, n_intervals];
  
  int z_rep[n_species, n_sites, n_intervals];
  int y_rep_cs[n_species, n_sites, n_intervals, n_visits]; // repd detections

  real eval_cs[n_species,n_sites,n_intervals,n_visits]; // expected values

  real T_rep_cs[n_species]; // Freeman-Tukey distance from eval (species bin)
  real T_obs_cs[n_species]; // Freeman-Tukey distance from eval (species bin)

  real P_species_cs[n_species]; // P-value by species

  // Initialize T_rep and T_obs and P-values
  for(i in 1:n_species){
    
    T_rep_cs[i] = 0;
    T_obs_cs[i] = 0;
    
    P_species_cs[i] = 0;

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
              ((1 - inv_logit(logit_p_cs[i,j,k]))^n_visits);
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

          // Compute fit statistic (Tukey-Freeman) for replicate data
          // community science records
          // Binned by species
          T_rep_cs[i] = T_rep_cs[i] + (sqrt(y_rep_cs[i,j,k,l]) - 
            sqrt(eval_cs[i,j,k,l]))^2;
          // Compute fit statistic (Tukey-Freeman) for real data
          // Binned by species
          T_obs_cs[i] = T_obs_cs[i] + (sqrt(V_cs[i,j,k,l]) - 
            sqrt(eval_cs[i,j,k,l]))^2;
          
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
    
  }

}
