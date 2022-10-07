// multi-species occupancy model for GBIF occurrence data

functions{
  
  // if the species is detected at the site*interval at least once..
  real lp_observed(int V, real logit_psi, real logit_p){ 
    
    return log_inv_logit(logit_psi) + 
            // probability density of getting an occurrence 
            // of a species at a site*interval plus..
            binomial_logit_lpmf(V | 1, logit_p); 
            // probability density of then observing or not observing  
            // that specific species at the site*interval per single visit l
  }
  
  // if the species is never detected at the site*interval..
  real lp_unobserved(real logit_psi, real logit_p){ 

    return log_sum_exp(log_inv_logit(logit_psi) +
           // probability density of the species occupying the site*interval
           // but not being detected at the visit l plus..
           log1m_inv_logit(logit_p),
           // probability density of not detecting it on the visit l
           
           // summed with..
           log1m_inv_logit(logit_psi));
           // probability density of the species NOT occupying the site*interval
           // and therefore detection is not possible
  }
  
}

data {
  
  int<lower=1> nsp;  // observed species
  int<lower=1> species[nsp];      // vector of species
  
  int<lower=1> nsite;  // sites within region
  int<lower=1> sites[nsite];      // vector of sites
  
  int<lower=1> ninterval;  // intervals during which sites are visited
  vector[ninterval] intervals; // vector of intervals (used as covariate data for 
                                // fixed effect of occupancy interval (time) on occupancy)
  
  int<lower=1> nvisit; // visits within intervals
  int<lower=0> V[nsp, nsite, ninterval, nvisit];  // visits when species i was detected at site j on interval k
  
  
} // end data


parameters {
  
  // OCCUPANCY
  real mu_psi_0; // global interecept for occupancy
  
  // species specific intercept allows some species to occur at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[nsp] psi_sp; // species specific intercept for occupancy
  real<lower=0> sigma_psi_species; // variance in species intercepts
  
  // random slope for species specific temporal effects on occupancy
  vector[nsp] psi_interval; // vector of species specific slope estimates
  real mu_psi_interval; // community mean of species specific slopes
  real<lower=0> sigma_psi_interval; // variance in species slopes
  
  // DETECTION
  real mu_p_0;
  
  // species specific intercept allows some species to be detected at higher rates than others, 
  // but with overall estimates for occupancy partially informed by the data pooled across all species.
  vector[nsp] p_sp; // species specific intercept for detection
  real<lower=0> sigma_p_species; // variance in species intercepts
  
  // random slope for site specific temporal effects on occupancy
  vector[nsite] p_site; // vector of spatially specific slope estimates
  real<lower=0> sigma_p_site; // variance in site slopes
  
  real p_interval; // fixed temporal effect on detection probability
  
} // end parameters


transformed parameters {
  
  real logit_psi[nsp, nsite, ninterval];  // log odds  of occurrence
  real logit_p[nsp, nsite, ninterval, nvisit];  // log odds of detection
  
  for (i in 1:nsp){   // loop across all species
    for (j in 1:nsite){    // loop across all sites
      for(k in 1:ninterval){
          
          logit_psi[i, j, k] = // log odds  of occurrence is equal to
            mu_psi_0 + // a baseline intercept
            psi_sp[species[i]] + // a species specific intercept
            psi_interval[species[i]]*(k); // a species specific temporal effect
            
      }
    }
  }
  
  for (i in 1:nsp){   // loop across all species
    for (j in 1:nsite){    // loop across all sites
      for(k in 1:ninterval){
        for(l in 1:nvisit){
        
          logit_p[i, j, k, l] = // log odds of detection is equal to
            mu_p_0 + // a baseline intercept
            p_sp[species[i]] + // a species specific intercept
            p_site[sites[j]] + // a spatially specific intercept
            p_interval*intervals[k]; // an overall effect of time on detection
            
        }
      }
    }
  }
  
} // end transformed parameters


model {
  // PRIORS
  
  // Occupancy (Ecological Process)
  mu_psi_0 ~ cauchy(0, 2.5); // global intercept for occupancy rate
  
  psi_sp ~ normal(0, sigma_psi_species); 
  // occupancy intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_psi_species ~ cauchy(0, 2.5);
  
  psi_interval ~ normal(0, sigma_psi_interval);
  // occupancy slope (temporal effect on occupancy) for each species drawn from the 
  // community distribution (variance defined by sigma), centered at mu_psi_interval. 
  // centering on mu (rather than 0) allows us to estimate the average effect of
  // the management on abundance across all species.
  sigma_psi_interval ~ cauchy(0, 2.5); // community variance
  
  // Detection (Observation Process)
  mu_p_0 ~ cauchy(0, 2.5); // global intercept for detection
  
  p_sp ~ normal(0, sigma_p_species); 
  // detection intercept for each species drawn from the community
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_species ~ cauchy(0, 2.5);
  
  // multivariate hierarchical prior
  p_site ~ normal(0, sigma_p_site);
  // detection intercept for each site*interval drawn from the spatiotemporal
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_site ~ cauchy(0, 2.5); // spatiotemporal variance
  
  p_interval ~ cauchy(0, 2.5);
  
  // LIKELIHOOD
  // Stan can sample mean and sd of parameters by summing out the
  // parameter (marginalizing) across likelihood statements
  for(i in 1:nsp) { // loop across all species
    for(j in 1:nsite) { // loop across all sites
      for(k in 1:ninterval){
        for(l in 1:nvisit){
          
          // if species is detected at the specific site*interval at least once
          // lp_observed calculates the probability density that occurs given logit_psi plus
          // the probability density that we did/did not observe it on each visit l in 1:nvisit
          if(sum(V[i, j, k, 1:nvisit]) > 0){ 
            target += lp_observed(V[i, j, k, l], 
              logit_psi[i, j, k], logit_p[i, j, k, l]);
          
          // else the species was never detected at the site*interval
          // lp_unobserved sums the probability density of:
          // 1) species occupies the site*interval but was not detected on each visit, and
          // 2) the species does not occupy the site*interval
          } else {
            target += lp_unobserved(logit_psi[i, j, k], logit_p[i, j, k, l]);
            
          } // end if/else
          
        }
      }
    }
  }
} // end model
