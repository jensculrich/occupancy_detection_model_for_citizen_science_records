## run_model_integrated.R
### Run abundance-occupancy model (model_abundance.stan) using real pollinator occurrence data from GBIF
# jcu; started nov 24, 2022

library(rstan)

## --------------------------------------------------
# input data preparation choices
# be careful that the (era_end - era_start) is evenly divisible by the n_intervals
era_start = 2008 # must define start date of the GBIF dataset
era_end = 2022 # must define start date of the GBIF dataset
n_intervals = 5 # must define number of intervals to break up the era into
n_visits = 3 # must define the number of repeat obs years within each interval
# note, should introduce throw error if..
# (era_end - era_start) / n_intervals has a remainder > 0,
min_records_per_species = 15 # filters species with less than this many records (total between both datasets)..
# within the time span defined above
grid_size = 30000 # in metres so, e.g., 25000 = 25km x 25 km 
min_population_size = 300 # min pop density in the grid cell (per km^2)
# for reference, 38people/km^2 is ~100people/mile^2
# 100/km^2 is about 250/mile^sq
min_species_for_community_sampling_event = 2 # community sampling inferred if..
# species depositied in single institution from a site in a single year is >= min_species_for_community_sampling_event
# min records_for_community_sampling_event sets a minimum threshold, if the number
# of records for the taxonomic group within a site within a year is 
min_year_for_species_ranges = 2000 # use all data from after this year to infer species ranges
taxon = "syrphidae" # taxon to analyze, either "syrphidae" or "bombus"
# minimum site area 
# if sites are super tiny, the observation process could likely be very unstable
min_site_area = 0.10
# remove specimens lacking species-level id before calculating summary statistics?
# Note, they will get removed before sending to the model either way, but this turns on/off
# whether they are included in the counts of obs per data set, per species, in museums v cit sci, etc.
remove_unidentified_species = TRUE


source("./abundance-occupancy/data_prep/prep_data.R")

my_data <- prep_data(era_start = era_start, # must define start date of the GBIF dataset
                     era_end = era_end, # must define start date of the GBIF dataset
                     n_intervals = n_intervals, # must define number of intervals to break up the era into
                     n_visits = n_visits, # must define the number of repeat obs years within each interval
                     # note, should introduce throw error if..
                     # (era_end - era_start) / n_intervals has a remainder > 0,
                     min_records_per_species = min_records_per_species,
                     grid_size = grid_size, # 25km x 25 km 
                     min_population_size = min_population_size, # min pop density in the grid cell (per km^2)
                     # for reference, 38people/km^2 is ~100people/mile^2
                     # 100/km^2 is about 250/mile^sq
                     min_records_for_community_sampling_event = min_records_for_community_sampling_event,
                     min_year_for_species_ranges = min_year_for_species_ranges,
                     taxon,
                     min_site_area,
                     remove_unidentified_species
                     
)

# save the data in case you want to make tweaks to the model run
# without redoing the data prep
# saveRDS(my_data, "./abundance-occupancy/analysis/prepped_data_list.rds")
my_data <- readRDS("./analysis/prepped_data_list.rds")

gc()

# data to feed to the model
V_citsci <- my_data$V_citsci # citizen science detection data
V_museum <- my_data$V_museum # museum detection data
V_citsci_NA <- my_data$V_citsci_NA # cit science NA indicator array
V_museum_NA <- my_data$V_museum_NA # museum data NA indicator array
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_intervals <- my_data$n_intervals # number of surveys 
n_visits <- my_data$n_visits

interval_names <- as.vector(as.numeric(my_data$intervals))
site_names <- my_data$sites
species_names <- my_data$species

pop_densities <- my_data$pop_densities
impervious_cover <- my_data$impervious_cover
site_areas <- my_data$site_areas

# check correlation between variables
correlation_matrix <- my_data$correlation_matrix

# intervals will cause issues if you try to run on only 1 interval
# since it's no longer sent in as a vector of intervals (can you force a single
# integer to be a vector if you truly want to treat all as a single interval?)
intervals_raw <- as.vector(seq(1, n_intervals, by=1)) 
intervals <- intervals_raw - 1
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

# Get K values
source("./abundance-occupancy/data_prep/get_K.R")
K <- get_K(V_citsci)

# sum(V_museum_NA) # number of species sampling events by museums
# c <- which(V_museum>V_museum_NA) # this will give you numerical value

stan_data <- c("V_citsci", "V_museum",
               "V_citsci_NA", "V_museum_NA", 
               "n_species", "n_sites", "n_intervals", "n_visits", 
               "intervals", "species", "sites",
               "pop_densities", "site_areas"
              )
               # add correct variables

# Parameters monitored
params <- c("mu_psi_0",
            "psi_species",
            "sigma_psi_species",
            "sigma_psi_site",
            "psi_pop_density",
            "mu_psi_pop_density",
            "sigma_psi_pop_density",
            "psi_interval",
            "mu_psi_interval",
            "sigma_psi_interval",
            "psi_site_area",
            
            "mu_p_citsci_0",
            "p_citsci_species",
            "sigma_p_citsci_species",
            "sigma_p_citsci_site",
            "p_citsci_interval",
            "p_citsci_pop_density", 
            
            "mu_p_museum_0",
            "p_museum_species",
            "sigma_p_museum_species",
            "sigma_p_museum_site",
            "p_museum_interval",
            "p_museum_pop_density"
)


# MCMC settings
n_iterations <- 1000
n_thin <- 2
n_burnin <- 500
n_chains <- 3
n_cores <- n_chains
delta = 0.9

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(mu_psi_0 = runif(1, -1, 1),
       sigma_psi_species = runif(1, 0, 0.5),
       sigma_psi_site = runif(1, 0, 0.5),
       mu_psi_interval = runif(1, -0.1, 0.1),
       sigma_psi_interval = runif(1, 0, 0.1),
       mu_psi_pop_density = runif(1, -1, 1),
       sigma_psi_pop_density = runif(1, 0, 1),
       psi_site_area = runif(1, -1, 1),
       
       mu_p_citsci_0 = runif(1, -1, 1),
       sigma_p_citsci_species = runif(1, 0, 1),
       sigma_p_citsci_site = runif(1, 0, 1),
       p_citsci_interval = runif(1, -1, 1),
       p_citsci_pop_density = runif(1, -1, 1),
       
       mu_p_museum_0 = runif(1, -1, 1),
       sigma_p_museum_species = runif(1, 0, 1),
       sigma_p_museum_site = runif(1, 0, 1),
       p_museum_interval = runif(1, -1, 1),
       p_museum_pop_density = runif(1, -1, 1)
       
  )
)

## --------------------------------------------------
### Run model
stan_model <- "./models/model.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 control=list(adapt_delta=delta),
                 seed = 12,
                 open_progress = FALSE,
                 cores = n_cores)

# print main effects
print(stan_out, digits = 3, pars=
        c("mu_psi_0",
          "sigma_psi_species",
          "sigma_psi_site",
          "mu_psi_pop_density",
          "sigma_psi_pop_density",
          "mu_psi_interval",
          "sigma_psi_interval",
          "psi_site_area",
          
          "mu_p_citsci_0",
          "sigma_p_citsci_species",
          "sigma_p_citsci_site",
          "p_citsci_interval",
          "p_citsci_pop_density", 
          
          "mu_p_museum_0",
          "sigma_p_museum_species",
          "sigma_p_museum_site",
          "p_museum_interval",
          "p_museum_pop_density"))

# print sampled random effects
print(stan_out, digits = 3, pars=
        c("psi_pop_density"))

print(stan_out, digits = 3, pars=
        c("psi_species[8]",
          "psi_species[36]",
          "psi_species[41]",
          "psi_species[12]",
          "psi_species[23]",
          
          "psi_pop_density[8]",
          "psi_pop_density[36]",
          "psi_pop_density[41]",
          "psi_pop_density[12]",
          "psi_pop_density[23]",
          
          "p_citsci_species[8]",
          "p_citsci_species[36]",
          "p_citsci_species[41]",
          "p_citsci_species[12]",
          "p_citsci_species[23]",
          
          "p_museum_species[8]",
          "p_museum_species[36]",
          "p_museum_species[41]",
          "p_museum_species[12]",
          "p_museum_species[23]"))



saveRDS(stan_out, "./model_outputs/stan_out_model_integrated_ranges_200_35km_25records.rds")
# stan_out <- readRDS("./model_outputs/stan_out_model_integrated_ranges_200_35km_25records.rds")

## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out, pars = c(
  "mu_psi_0",
  "mu_psi_interval",
  "mu_psi_pop_density",
  "psi_site_area",
  
  "mu_p_citsci_0",
  "sigma_p_citsci_species",
  "sigma_p_citsci_site",
  "p_citsci_interval",
  "p_citsci_pop_density", 
  
  "mu_p_museum_0",
  "sigma_p_museum_species",
  "sigma_p_museum_site",
  "p_museum_interval",
  "p_museum_pop_density"
))

# pairs plot
pairs(stan_out, pars = c(
  "mu_psi_0",
  "mu_psi_interval",
  "sigma_psi_interval",
  "mu_psi_pop_density",
  "psi_site_area",
  
  "mu_p_citsci_0",
  #"sigma_p_citsci_species",
  #"sigma_p_citsci_site",
  "p_citsci_interval",
  #"p_citsci_pop_density", 
  
  "mu_p_museum_0",
  #"sigma_p_museum_species",
  #"sigma_p_museum_site",
  "p_museum_interval"
  # "p_museum_pop_density"
))

# should now also write a posterior predictive check into the model
