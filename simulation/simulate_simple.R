## Data simulation for simplest form of the data
# single season, species considered across all sites (not just range)
# and across all visits (not considering that sites without any detections 
# for any species were probably not visited or that sites with only detections for
# a few species only might be targeted and not fully community wide)


## --------------------------------------------------
### Let's try to reform the data into a series of vectors rather than a big array
# I'm having a more difficult time keeping track of all of the random effects in the array form
# and this format has worked for me before. Are there any drawbacks to having the data in a series of
# vectors? Is it more computationally intensive?
library(tidyverse)

## --------------------------------------------------
### Variable values for data simulation

## study dimensions
n_species=10 ## number of species (number of species included in the studies taxonomic scope)
n_sites=10 ## number of sites (number of cities included in the study)
n_intervals=10 ## number of occupancy intervals (i.e. number of sets of years within which we are estimating occupancy rates)
n_visits=3 ## number of repeat 'surveys' per interval (i.e. number of years within an interval in which data is collected)

# number of unique data observation points 
# (observations of unique species*site*occupancy_interval*visit combinations)
R = n_species*n_sites*n_intervals
n_data_per_species <- R/n_species
n_data_per_site <- R/n_sites
n_data_per_interval <- R/n_intervals

## occupancy
mu_psi_0 = 0
sigma_psi_sp = 0.5
mu_psi_interval = 0.5
sigma_psi_interval = 0.2

## detection
mu_p_0 = 0
sigma_p_species = 0.3
sigma_p_site = 0.3
p_interval = 0.5

## expit and logit functions
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
### Covariates

interval <- rep(1:n_intervals, each = n_data_per_interval)

site <- rep(1:n_sites, each = n_species, 
             times = n_data_per_site/n_sites)

species <- rep(1:n_species, times = 
                 n_data_per_species)

visit <- rep(1:n_visits, R/n_visits)

View(as.data.frame(cbind(interval, site, species)))

## --------------------------------------------------
### Ecological process

## add a species-level random effect for occupancy intercept
# (some species more frequently than others), centered on 0
species_intercepts_psi <- rnorm(n_species, 0, sigma_psi_sp)
(mean(species_intercepts_psi)) # should be near 0 since we are centered on 0
(sd(species_intercepts_psi)) # should be near value of sigma_psi_sp

psi_species <- as.numeric(vector(length=R))
psi_species <- rep(species_intercepts_psi[1:n_species], times = n_data_per_species)

# and a species-level random slope (change in occupancy through time
# allowed to respond uniquely by species)
# centered on community mean and defined variation
species_slopes_psi <- rnorm(n_species, mu_psi_interval, sigma_psi_interval)
(mean(species_slopes_psi)) # should be near mu_alpha1

psi_occupancy_interval <- as.numeric(vector(length=R))
psi_occupancy_interval <- rep(species_slopes_psi[1:n_species], times = n_data_per_species)

test <- cbind(interval, site, species, psi_species, psi_occupancy_interval)

### now generate occupancy probabilities
psi_vector = as.numeric(vector(length=R))

for(i in 1:R) { # for each 
      
  psi_vector[i] <- expit( # occupancy is equal to
          mu_psi_0 + # a baseline intercept
          psi_species[i] + # a species specific intercept
          psi_occupancy_interval[i]*(interval[i])) # a species specific temporal change
      
}

psi_data <- as.data.frame(cbind(
  interval, site, species, psi_species, psi_occupancy_interval, psi_vector))

## --------------------------------------------------
## Generate presence and absence data

# Z[i]is the actual outcome of the occupancy coin flip based on weight captured in psi_vector[i]
Z = as.numeric(vector(length=R))

for(i in 1:R) { # for each 
  
  # eventually here we will need to specify which(site) to restrict to actual range
  Z[i] <- rbinom(n = 1, size = 1, prob = psi_vector[i])
  
}

psi_data_simmed <- as.data.frame(cbind(
  interval, site, species, psi_species, psi_occupancy_interval, psi_vector, Z))

## --------------------------------------------------
### Observation process

## add a species-level random effect for occupancy intercept
# (some species more frequently observed than others), centered on 0
species_intercepts_p <- rnorm(n_species, 0, sigma_p_species)
(mean(species_intercepts_p)) # should be near 0 since we are centered on 0
(sd(species_intercepts_p)) # should be near value of sigma_psi_sp

p_species <- as.numeric(vector(length=R))
p_species <- rep(species_intercepts_p[1:n_species], times = n_data_per_species)

## add a site-level random effect for occupancy intercept
# (some sites more effectively surveyed/searched than others), centered on 0
site_intercepts_p <- rnorm(n_sites, 0, sigma_p_site)
(mean(site_intercepts_p)) # should be near 0 since we are centered on 0
(sd(site_intercepts_p)) # should be near value of sigma_psi_sp

p_site <- as.numeric(vector(length=R))
p_site <- rep(site_intercepts_p[1:n_sites], times = n_data_per_site)

test3 <- cbind(interval, site, species, p_species, p_site)

### now generate detection probabilities
p_vector = as.numeric(vector(length=R))

for(i in 1:R) { # for each visit
  
  p_vector[i] <- expit( # detection is equal to 
      mu_p_0 + # a baseline intercept
      p_species[i] + # a species specific intercept
      p_site[i] + # a spatially specific intercept
      p_interval*(interval[i])) # an overall effect of time on detection
  
} # for each visit

test4 <- cbind(interval, site, species, p_species, p_site, p_vector)

## --------------------------------------------------
## Generate some presence/absence detection data

V = matrix(NA, nrow = length(p_vector), ncol = n_visits) # Array for counts

for(i in 1:R) { # for each 
  for(j in 1:n_visits){
  
  # eventually here we will need to specify which(site) to restrict to actual range
  V[i, j] <- Z[i] * # occupancy state * detection prob
    rbinom(n = 1, size = 1, prob = p_vector[i])
  
  }
}

View(as.data.frame(V))

# species should be detected some of the time when they are there, but never when Z = 0.
# View(cbind(Z, V))

# should capture all above into a function

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
R <- R # number of species*site*survey combinations
V <- V # detection data
n_species <- n_species # number of species
n_sites <- n_sites # number of sites
n_intervals <- n_intervals # number of surveys 
n_visits <- n_visits

interval <- interval
site <- site
species <- species
visit <- visit

stan_data <- c("R", "V", "n_intervals", "n_sites", "n_species", "n_visits", 
               "species", "site", "interval", "visit")

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

## --------------------------------------------------
### In Shirey et al. 2022 they use the data in an array format simulated as below ->

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