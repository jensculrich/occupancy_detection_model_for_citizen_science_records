## Data simulation for simplest form of the data
# single season, species considered across all sites (not just range)
# and across all visits (not considering that sites without any detections 
# for any species were probably not visited or that sites with only detections for
# a few species only might be targeted and not fully community wide)

## --------------------------------------------------
### Variable values for data simulation
## study dimensions
nsp=20
nsite=18 ## number of sites
ninterval=10 ## number of occupancy intervals
nvisit=6 ## number of samples per year


## detection
mu.p.0 = 0
p.interval = 0.5
sigma.p.site = 0.3
sigma.p.sp = 0.3

## occupancy
mu.psi.0 = 0
sigma.psi.sp = 0.5
mu.psi.interval = 0.5
sigma.psi.interval = 0.2

## visit
mu.v.0 = 0 
mu.v.interval = -0.5

## type sym
type.range = "equal"
type.visit = 'visit_miss'
prop.visits.same = 1

## expit and logit functions
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
### specify species-specific occupancy probabilities

## species-specific random intercepts
psi.sp <- rnorm(n=nsp, mean=0, sd=sigma.psi.sp)
# species baseline occupancy is drawn from a normal distribution with mean 0 and 
# species specific variation defined by sigma.psi.sp

## effect of interval on occupancy (species-specific random slopes)
psi.interval <- rnorm(n=nsp, mean=mu.psi.interval, sd=sigma.psi.interval)
# change in each species occupancy across time is drawn from a distribution defined
# by a community mean (mu.psi.interval) with 
# species specific variation defined by sigma.psi.interval

## --------------------------------------------------
### specify species-specific detection probabilities

## species-specific random intercepts
p.sp  <- rnorm(n=nsp, mean = 0, sd=sigma.p.sp)
# species baseline detection probability is drawn from a normal distribution with mean 0 and 
# species specific variation defined by sigma.p.sp

## effect of site and interval on detection (site,interval-specific random slopes)
p.site <- rnorm(n=nsite, mean=0, sd=sigma.p.site)

# spatiotemporal variability in detection probability (changing across sites and 
# occupancy intervals), and helps account for the variation that is inherent in 
# sample effort across space and time in unstructured historical datasets. 
# spatiotemporal variability at each site*interval is drawn from a distribution defined
# by a mean of 0 and with site*interval specific variation of sigma.p.site

## --------------------------------------------------
## Create arrays for psi and p
psi.mat <- array(NA, dim =c(nsp, nsite, ninterval)) 
# a psi value for each species, at each site, in each interval 

p.mat <- array(NA, dim =c(nsp, nsite, ninterval, nvisit))
# a p value for each species, at each site, in each interval, AND in each visit

for(site in 1:nsite) { # for each site
  for(interval in 1:ninterval) { # for each interval
    for(sp in 1:nsp) { # for each species
      
      psi.mat[sp,site,interval] <- expit( # occupancy is equal to
                                     mu.psi.0 + # a baseline intercept
                                     psi.sp[sp] + # a species specific intercept
                                     psi.interval[sp]*(interval)) # a species specific temporal change
      
      for(visit in 1:nvisit) { # for each visit
        
        p.mat[sp,site,interval,visit] <- expit( # detection is equal to 
                                           mu.p.0 + # a baseline intercept
                                           p.sp[sp] + # a species specific intercept
                                           p.site[site] + # a spatiotemporally specific intercept
                                           p.interval*(interval)) # an overall effect of time on detection
      } # for each visit
    } # for each species
  } # for each interval
} # for each site

# preview the psi and p arrays
head(psi.mat[1:20, 1:18,1])
head(p.mat[1:20, 1:18,1,1])

(p.mat[1,1,1:ninterval,1]) # if p.interval is >0 these should generally be increasing from low to high
(p.mat[2,1,1:ninterval,1]) # if p.interval is >0 these should generally be increasing from low to high
(p.mat[1,2,1:ninterval,1]) # if p.interval is >0 these should generally be increasing from low to high

## --------------------------------------------------
## Generate presence and absence

Z <- array(NA, dim=c(nsp=nsp,
                     nsite=nsite,
                     ninterval=ninterval))

for(interval in 1:ninterval){
  for(site in 1:nsite){
    for(sp in 1:nsp){
    
    # eventually here we will need to specify which(site) to restrict to actual range
    Z[sp,site,interval] <- rbinom(n = 1, size = 1, prob = psi.mat[sp,site,interval])
    
    }
  }
}

## --------------------------------------------------
## Generate detection non detection data 

V <- array(NA, dim=c(nsp=nsp,
                     nsite=nsite,
                     ninterval=ninterval,
                     nvisit = nvisit))

for(sp in 1:nsp){
  for(site in 1:nsite){
    for(interval in 1:ninterval){
      for(v in 1:nvisit){
        
        # eventually here we will need to specify which(site) to restrict to actual range
        V[sp,site,interval,v] <- Z[sp,site,interval] * # occupancy state * detection prob
                              rbinom(n = 1, size = 1, prob = p.mat[sp,site,interval,v])
    
      }
    }
  }
}

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V <- V # detection data
nsp <- nsp # number of species
nsite <- nsite # number of sites
ninterval <- ninterval # number of surveys 
nvisit <- nvisit
# Interval values: numeric vector (will act as covariate data for psi.interval)
intervals <- seq(1, ninterval, by=1)
sites <- seq(1, nsite, by=1)
species <- seq(1, nsp, by=1)

stan_data <- c("V", "nsp", "nsite", "ninterval", "nvisit", 
               "species", "sites", "intervals")

# Parameters monitored
params <- c("mu_psi_0",
            "psi_sp",
            "sigma_psi_species",
            "psi_interval",
            "mu_psi_interval",
            "sigma_psi_interval",
            "mu_p_0",
            "p_sp",
            "sigma_p_species",
            "p_site",
            "sigma_p_site",
            "p_interval"
)


# MCMC settings
n_iterations <- 400
n_thin <- 1
n_burnin <- 200
n_chains <- 3
n_cores <- 3

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  list(mu_psi_0 = runif(1, -1, 1),
       psi_sp = runif(1, -1, 1),
       sigma_psi_species = runif(1, -1, 1),
       psi_interval = runif(1, -1, 1),
       mu_psi_interval = runif(1, -1, 1),
       sigma_psi_interval = runif(1, -1, 1),
       mu_p_0 = runif(1, -1, 1),
       p_sp = runif(1, -1, 1),
       sigma_p_species = runif(1, -1, 1),
       p_site = runif(1, -1, 1),
       sigma_p_site = runif(1, -1, 1),
       p_interval = runif(1, -1, 1)
       
  )
)

#targets <- as.data.frame(cbind(param_names, parameters))

## --------------------------------------------------
### Run model
library(rstan)
stan_model <- "./simulation/model_simple.stan"

## Call Stan from R
stan_out_sim <- stan(stan_model,
                   data = stan_data, 
                   init = inits, 
                   pars = params,
                   chains = n_chains, iter = n_iterations, 
                   warmup = n_burnin, thin = n_thin,
                   seed = 1,
                   open_progress = FALSE,
                   cores = n_cores)

print(stan_out_sim, digits = 3)

saveRDS(stan_out_sim, ".//stan_out_sim.rds")
stan_out_sim <- readRDS(".//stan_out_sim.rds")