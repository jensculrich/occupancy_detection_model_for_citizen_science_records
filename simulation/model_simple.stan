// multi-species occupancy model for GBIF occurrence data

functions {
  
  matrix cov_matrix_2d(vector sigma, real rho) {
  matrix[2,2] Sigma;
  Sigma[1,1] = square(sigma[1]);
  Sigma[2,2] = square(sigma[2]);
  Sigma[1,2] = sigma[1] * sigma[2] * rho;
  Sigma[2,1] = Sigma[1,2];
  return Sigma;
  
  
  real lp_observed(int X, int K, real logit_psi, real logit_theta) {
  // if species i is observed at site j:
  return log_inv_logit(logit_psi)
    // where logit_psi[i] = (uv[i, 1] + a1_species_occ * species_cov1[i]) + 
                              // (a1_site_occ * site_cov1[j])
    + binomial_logit_lpmf(X | K, logit_theta); 
    // The log binomial probability mass of x successes (observed)
    // in K trials (survey revisits) given the
    // logit-scaled chance of success logit_theta (detection prob)
    // where logit_theta[i] = uv[i, 2]
}
  
} // end functions

data {
  
  int<lower=1> nsp;  // observed species
  int<lower=1> nsites;  // sites within region
  int<lower=1> ninterval;  // intervals during which sites are visited
  int<lower=1> nvisit; // visits within intervals
  int<lower=0> V[nsp, nsites, ninterval, nvisit];  // visits when species i was detected at site j on interval k
  vector[ninterval] intervals;
  
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
  matrix[nsites, ninterval] p_site; // matrix of spatiotemporally specific slope estimates
  real<lower=-1,upper=1> rho_p_site_interval;  // correlation of (site detection, interval detection) heterogeneity
  real<lower=0> sigma_p_site; // variance in site slopes
  
  real p_interval; // fixed temporal effect on detection probability
  
  corr_matrix[ninterval] Omega;        // prior correlation
  vector<lower=0>[ninterval] tau;      // prior scale
  
} // end parameters

transformed parameters {
  
  real logit_psi[nsp, nsites, ninterval];  // log odds  of occurrence
  real logit_p[nsp, nsites, ninterval, nvisit];  // log odds of detection
  
  for (i in 1:nsp){   // loop across all species
    for (j in 1:nsites){    // loop across all sites
      for(k in 1:ninterval){
          
          logit_psi[i, j, k] = // log odds  of occurrence is equal to
            mu_psi_0 + // a baseline intercept
            psi_sp[i] + // a species specific intercept
            psi_interval[i]*(k); // a species specific temporal effect
            
      }
    }
  }
  
  for (i in 1:nsp){   // loop across all species
    for (j in 1:nsites){    // loop across all sites
      for(k in 1:ninterval){
        for(l in 1:nvisit){
        
          logit_p[i, j, k, l] = // log odds of detection is equal to
            mu_p_0 + // a baseline intercept
            p_sp[i] + // a species specific intercept
            p_site[j, k] + // a spatiotemporally specific intercept
            p_interval*intervals[k]; // an overall effect of time on detection
            
        }
      }
    }
  }
  
} // end transformed parameters

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
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
  (rho_p_site_interval + 1) / 2 ~ beta(2, 2);
  p_site ~ multi_normal(rep_vector(0, 2), cov_matrix_2d(sigma_p_site, rho_p_site_interval));
  // detection intercept for each site*interval drawn from the spatiotemporal
  // distribution (variance defined by sigma), centered at 0. 
  sigma_p_site ~ cauchy(0, 2.5); // spatiotemporal variance
  
  p_interval ~ cauchy(0, 2.5);
  
  // LIKELIHOOD
  // Stan can sample mean and sd of parameters by summing out the
  // parameter (marginalizing) across likelihood statements
  for(i in 1:nsp) { // loop across all species
    for(j in 1:nsites) { // loop across all sites
      for(k in 1:ninterval){
        for(l in 1:nvisit){
          
        if(V[i, j, k] > 0) // if species is detected at the specific site 1 or more times
          target += lp_observed(V[i, j], K, logit_psi[i, j], logit_theta[i, j]);
          
        }
      }
    }
  }
} // end model

