## run_model_integrated.R
### Run occupancy model (model_"taxon".stan), using real pollinator occurrence data from GBIF
# jcu; started nov 24, 2022

## --------------------------------------------------
# input data preparation choices - SYRPHIDAE
# be careful that the (era_end - era_start) is evenly divisible by the n_intervals
era_start = 2011 # must define start date of the GBIF dataset
era_end = 2022 # must define start date of the GBIF dataset
n_intervals = 4 # must define number of intervals to break up the era into
n_visits = 3 # must define the number of repeat obs years within each interval
# note, should introduce throw error if..
# (era_end - era_start) / n_intervals has a remainder > 0,
min_records_per_species = 5 # filters species with less than this many records (total between both datasets)..
# within the time span defined above
grid_size = 15000 # in metres so, e.g., 25000 = 25km x 25 km 
min_population_size = 1000 # min pop density in the grid cell (per km^2)

min_species_for_community_sampling_event = 2 # community sampling inferred if..
# species depositied in single institution from a site in a single year is >= min_species_for_community_sampling_event
# min records_for_community_sampling_event sets a minimum threshold, if the number
# of records for the taxonomic group within a site within a year is 
min_year_for_species_ranges = 2000 # use all data from after this year to infer species ranges
taxon = "syrphidae" # taxon to analyze, either "syrphidae" or "bombus"
min_site_area = 0.25 # if sites are super tiny, the observation process could likely be very unstable

# remove specimens lacking species-level id before calculating summary statistics?
# Note, they will get removed before sending to the model either way, but this turns on/off
# whether they are included in the counts of obs per data set, per species, in museums v cit sci, etc.
remove_unidentified_species = TRUE # default to TRUE
consider_species_occurring_outside_sites = FALSE # default to FALSE # consider species that were detected outside of the sites but not at sites?
min_records_per_species_full = 15 # min rec threshold if above is true
make_range_plot = FALSE # default to FALSE # plot ranges
urban_sites = TRUE # default to TRUE # cut sites to above urban threshold if true - if false will return non urban sites
non_urban_subsample_n = 100 # if urban_sites is true, then how many sites do you want to keep? Keeping all will yield too much site data for computer to handle
infer_detections_at_genus = FALSE # default to FALSE # if true, infer non detections only for species in the same genus as a species detected (as opposed to any in the clade)
generate_temporal_plots = FALSE # default to FALSE

## --------------------------------------------------
# input data preparation choices - BOMBUS
# be careful that the (era_end - era_start) is evenly divisible by the n_intervals
era_start = 2011 # must define start date of the GBIF dataset
era_end = 2022 # must define start date of the GBIF dataset
n_intervals = 4 # must define number of intervals to break up the era into
n_visits = 3 # must define the number of repeat obs years within each interval
# note, should introduce throw error if..
# (era_end - era_start + 1) / n_intervals has a remainder > 0,
min_records_per_species = 10 # filters species with less than this many records (total between both datasets)..
# within the time span defined above (is only from urban sites, should redefine to be from anywhere)
grid_size = 15000 # in metres so, e.g., 25000 = 25km x 25 km 
min_population_size = 1000 # min pop density in the grid cell (per km^2)

min_species_for_community_sampling_event = 2 # community sampling inferred if..
# species depositied in single institution from a site in a single year is >= min_species_for_community_sampling_event
# min records_for_community_sampling_event sets a minimum threshold, if the number
# of records for the taxonomic group within a site within a year is 
min_year_for_species_ranges = 2000 # use all data from after this year to infer species ranges
taxon = "bombus" # taxon to analyze, either "syrphidae" or "bombus"
# minimum site area (proportion of grid_sizeXgrid_size that is in the admin area mask and not open water)
# if sites are super tiny, the observation process could likely be very unstable
min_site_area = 0.25

# remove specimens lacking species-level id before calculating summary statistics?
# Note, they will get removed before sending to the model either way, but this turns on/off
# whether they are included in the counts of obs per data set, per species, in museums v cit sci, etc.
remove_unidentified_species = TRUE # default to TRUE
consider_species_occurring_outside_sites = FALSE # default to FALSE # consider species that were detected outside of the sites but not at sites?
min_records_per_species_full = 15 # min rec threshold if above is true
make_range_plot = FALSE # default to FALSE # plot ranges
urban_sites = TRUE # default to TRUE # cut sites to above urban threshold if true - if false will return non urban sites
non_urban_subsample_n = 600 # if urban_sites is true, then how many sites do you want to keep? Keeping all will yield too much site data for computer to handle
infer_detections_at_genus = FALSE # default to FALSE # if true, infer non detections only for species in the same genus as a species detected (as opposed to any in the clade)
generate_temporal_plots = FALSE # default to FALSE

source("./occupancy/data_prep/prep_data.R")

my_data <- prep_data(era_start = era_start, # must define start date of the GBIF dataset
                     era_end = era_end, # must define start date of the GBIF dataset
                     n_intervals = n_intervals, # must define number of intervals to break up the era into
                     n_visits = n_visits, # must define the number of repeat obs years within each interval
                     min_records_per_species = min_records_per_species,
                     grid_size = grid_size, # 25km x 25 km 
                     min_population_size = min_population_size, # min pop density in the grid cell (per km^2)
                     min_records_for_community_sampling_event = min_records_for_community_sampling_event,
                     min_year_for_species_ranges = min_year_for_species_ranges,
                     taxon,
                     min_site_area,
                     remove_unidentified_species,
                     consider_species_occurring_outside_sites,
                     min_records_per_species_full,
                     make_range_plot,
                     urban_sites,
                     non_urban_subsample_n,
                     infer_detections_at_genus,
                     generate_temporal_plots
                     
)

# save the data in case you want to make tweaks to the model run
# without redoing the data prep
saveRDS(my_data, paste("./occupancy/analysis/prepped_data/", 
                       taxon, "/", grid_size / 1000, 
                       "km_", min_population_size, "minpop_", 
                       min_records_per_species, "minpersp_",
                       n_intervals, "ints_", n_visits, "visits",
                       #"_nonurban",
                       ".rds", sep = ""))

my_data <- readRDS(paste0("./occupancy/analysis/prepped_data/",
                          taxon, "/", grid_size / 1000, "km_",
                          min_population_size, "minpop_",
                          min_records_per_species, "minpersp", "_",
                          n_intervals, "ints_",
                          n_visits, "visits",
                          #"_nonurban",
                          ".rds"))

# best to restart R or offload all of the spatial data packages before running the model
gc()

library(rstan)

# data to feed to the model
V_citsci <- my_data$V_citsci # citizen science detection data
V_museum <- my_data$V_museum # museum detection data
ranges <- my_data$V_citsci_NA # cit science NA indicator array
V_museum_NA <- my_data$V_museum_NA # museum data NA indicator array
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_intervals <- my_data$n_intervals # number of surveys 

n_visits <- my_data$n_visits

interval_names <- as.vector(as.numeric(my_data$intervals))
site_names <- my_data$sites
species_names <- my_data$species

#saveRDS(species_names, "./figures/species_names/bombus_names_30km_nonurban.RDS")
#write.csv(as.data.frame(species_names), "./data/syrphidae_names.csv")

pop_densities <- my_data$pop_densities
avg_income <- my_data$avg_income
open_developed <- my_data$developed_open
herb_shrub <- my_data$herb_shrub_cover
site_areas <- my_data$site_areas
herb_shrub_forest <- my_data$herb_shrub_forest
developed_med_high <- my_data$developed_med_high
museum_total_records <- my_data$museum_total_records

ecoregion_three <- my_data$ecoregion_three_vector
ecoregion_one <- my_data$ecoregion_one_vector
ecoregion_three_lookup <- my_data$ecoregion_three_lookup
ecoregion_one_lookup <- my_data$ecoregion_one_lookup
n_ecoregion_three <- my_data$n_ecoregion_three
n_ecoregion_one <- my_data$n_ecoregion_one

# Other information about the data (not used by the model)
correlation_matrix <- my_data$correlation_matrix

species_counts <- my_data$species_counts
species_detections <- my_data$species_detections
species_counts_full <- my_data$species_counts_full

genus_lookup <- my_data$genus_lookup
genus_lookup <- as.numeric(factor(genus_lookup))
n_genera <- my_data$n_genera
genera <- my_data$genera

raw_pop_density <- my_data$raw_pop_density
raw_avg_income <- my_data$raw_avg_income

total_records <- my_data$total_records_since_2000
total_records_since_study <- my_data$total_records_since_2000
citsci_records <- my_data$citsci_records
citsci_detections <- my_data$citsci_detections
museum_records <- my_data$museum_records
museum_detections <- my_data$museum_detections

total_records_full <- my_data$total_records_since_2000_full
total_records_since_study_full <- my_data$total_records_since_2000_full
citsci_records_full <- my_data$citsci_records_full
museum_records_full <- my_data$museum_records_full

# intervals will cause issues if you try to run on only 1 interval
# since it's no longer sent in as a vector of intervals (can you force a single
# integer to be a vector if you truly want to treat all as a single interval?)
intervals_raw <- as.vector(seq(1, n_intervals, by=1)) 
intervals <- intervals_raw - 1
# intervals <- intervals_raw
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

# sum(V_museum_NA) # number of species sampling events by museums
#c <- which(V_museum>V_museum_NA) # this will give you numerical value

# True ? run model with no covariate data for occurence

if(taxon == "bombus"){
  
  if(urban_sites == TRUE){
   
    stan_data <- c("V_citsci", "V_museum", 
                   "ranges", "V_museum_NA",
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_ecoregion_three", "n_ecoregion_one",
                   "ecoregion_three", "ecoregion_one",
                   "ecoregion_three_lookup", "ecoregion_one_lookup",
                   "pop_densities", "site_areas", "avg_income", 
                   "herb_shrub_forest", "museum_total_records") 
    
    # Parameters monitored
    params <- c("mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_site",
                "sigma_psi_ecoregion_three",
                "sigma_psi_ecoregion_one",
                "mu_psi_income",
                "sigma_psi_income",
                #"sigma_psi_income_ecoregion_three",
                #"sigma_psi_income_ecoregion_one",
                "mu_psi_herb_shrub_forest",
                "sigma_psi_herb_shrub_forest",
                "psi_site_area",
                
                "mu_p_citsci_0",
                "sigma_p_citsci_species",
                "sigma_p_citsci_site",
                "sigma_p_citsci_ecoregion_three",
                "sigma_p_citsci_ecoregion_one",
                "p_citsci_interval",
                "p_citsci_pop_density", 
                
                "mu_p_museum_0",
                "sigma_p_museum_species",
                "sigma_p_museum_site",
                "sigma_p_museum_ecoregion_three",
                "sigma_p_museum_ecoregion_one",
                "p_museum_total_records",
                
                "psi_species",
                "psi_income",
                "psi_herb_shrub_forest",
                
                #"T_rep_citsci",
                #"T_obs_citsci",
                "P_species_citsci",
                
                #"T_rep_museum",
                #"T_obs_museum",
                "P_species_museum"
    )
    
    
    # MCMC settings
    n_iterations <- 1000
    n_thin <- 1
    n_burnin <- 500
    n_chains <- 4
    n_cores <- parallel::detectCores()
    #n_cores <- 4
    delta = 0.9
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      
      list(mu_psi_0 = runif(1, -1, 1),
           sigma_psi_species = runif(1, 0, 0.5),
           sigma_psi_site = runif(1, 0, 0.5),
           sigma_psi_ecoregion_three = runif(1, 0, 0.5),
           sigma_psi_ecoregion_one = runif(1, 0, 0.5),
           mu_psi_income = runif(1, -1, 1),
           sigma_psi_income = runif(1, 0, 1),
           #sigma_psi_income_ecoregion_three = runif(1, 0, 1),
           #sigma_psi_income_ecoregion_one = runif(1, 0, 1),
           mu_psi_herb_shrub_forest = runif(1, -1, 1),
           sigma_psi_herb_shrub_forest = runif(1, 0, 0.5),
           psi_site_area = runif(1, -1, 1),
           
           mu_p_citsci_0 = runif(1, -1, 0),
           sigma_p_citsci_species = runif(1, 0, 0.5),
           sigma_p_citsci_site = runif(1, 0, 0.5),
           sigma_p_citsci_ecoregion_three = runif(1, 0, 0.5),
           sigma_p_citsci_ecoregion_one = runif(1, 0, 0.5),
           p_citsci_interval = runif(1, -1, 1),
           p_citsci_pop_density = runif(1, -1, 1),
           
           mu_p_museum_0 = runif(1, -0.5, 0.5),
           sigma_p_museum_species = runif(1, 0, 0.1),
           sigma_p_museum_site = runif(1, 0, 0.1),
           sigma_p_museum_ecoregion_three = runif(1, 0, 0.1),
           sigma_p_museum_ecoregion_one = runif(1, 0, 0.1),
           p_museum_total_records = runif(1,  -0.5, 0.5)
           
      )
    )
     
  } else {
    
    stan_data <- c("V_citsci", "V_museum", 
                   "ranges", "V_museum_NA",
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_ecoregion_three", "n_ecoregion_one",
                   "ecoregion_three", "ecoregion_one",
                   "ecoregion_three_lookup", "ecoregion_one_lookup",
                   "pop_densities", "site_areas","museum_total_records") 
    
    # Parameters monitored
    params <- c("mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_site",
                "sigma_psi_ecoregion_three",
                "sigma_psi_ecoregion_one",
                "psi_site_area",
                
                "mu_p_citsci_0",
                "sigma_p_citsci_species",
                "sigma_p_citsci_site",
                "sigma_p_citsci_ecoregion_three",
                "sigma_p_citsci_ecoregion_one",
                "p_citsci_interval",
                "p_citsci_pop_density", 
                
                "mu_p_museum_0",
                "sigma_p_museum_species",
                "sigma_p_museum_site",
                "sigma_p_museum_ecoregion_three",
                "sigma_p_museum_ecoregion_one",
                "p_museum_total_records",
                
                "psi_species",
                
                #"T_rep_citsci",
                #"T_obs_citsci",
                "P_species_citsci",
                
                #"T_rep_museum",
                #"T_obs_museum",
                "P_species_museum"
    )
    
    
    # MCMC settings
    n_iterations <- 1200
    n_thin <- 1
    n_burnin <- 600
    n_chains <- 4
    n_cores <- parallel::detectCores()
    delta = 0.9
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      
      list(mu_psi_0 = runif(1, -1, 1),
           sigma_psi_species = runif(1, 0, 0.5),
           sigma_psi_site = runif(1, 0, 0.5),
           sigma_psi_ecoregion_three = runif(1, 0, 0.5),
           sigma_psi_ecoregion_one = runif(1, 0, 0.5),
           psi_site_area = runif(1, -1, 1),
           
           mu_p_citsci_0 = runif(1, -1, 0),
           sigma_p_citsci_species = runif(1, 0, 0.5),
           sigma_p_citsci_site = runif(1, 0, 0.5),
           sigma_p_citsci_ecoregion_three = runif(1, 0, 0.5),
           sigma_p_citsci_ecoregion_one = runif(1, 0, 0.5),
           p_citsci_interval = runif(1, -1, 1),
           p_citsci_pop_density = runif(1, -1, 1),
           
           mu_p_museum_0 = runif(1, -0.5, 0.5),
           sigma_p_museum_species = runif(1, 0, 0.1),
           sigma_p_museum_site = runif(1, 0, 0.1),
           sigma_p_museum_ecoregion_three = runif(1, 0, 0.1),
           sigma_p_museum_ecoregion_one = runif(1, 0, 0.1),
           p_museum_total_records = runif(1,  -0.5, 0.5)
           
      )
    )
    
  }
  
} else { # taxon == syrphidae
  
  if(urban_sites == TRUE){
    
    stan_data <- c("V_citsci", "V_museum", 
                   "ranges", "V_museum_NA",
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_genera", "genus_lookup",
                   "n_ecoregion_three", "n_ecoregion_one",
                   "ecoregion_three", "ecoregion_one",
                   "ecoregion_three_lookup", "ecoregion_one_lookup",
                   "pop_densities", "site_areas", "avg_income", 
                   "herb_shrub_forest", "museum_total_records") 
    
    # Parameters monitored
    params <- c("mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_genus",
                "sigma_psi_site",
                "sigma_psi_ecoregion_three",
                "sigma_psi_ecoregion_one",
                "mu_psi_income",
                "sigma_psi_income",
                "mu_psi_herb_shrub_forest",
                "sigma_psi_herb_shrub_forest",
                "psi_site_area",
                
                "mu_p_citsci_0",
                "sigma_p_citsci_species",
                "sigma_p_citsci_site",
                "sigma_p_citsci_ecoregion_three",
                "sigma_p_citsci_ecoregion_one",
                "p_citsci_interval",
                "p_citsci_pop_density", 
                
                "mu_p_museum_0",
                "sigma_p_museum_species",
                "sigma_p_museum_site",
                "sigma_p_museum_ecoregion_three",
                "sigma_p_museum_ecoregion_one",
                "p_museum_total_records",
                
                "psi_species",
                "psi_income",
                "psi_herb_shrub_forest",
                
                #"T_rep_citsci",
                #"T_obs_citsci",
                "P_species_citsci",
                
                #"T_rep_museum",
                #"T_obs_museum",
                "P_species_museum"
    )
    
    
    # MCMC settings
    n_iterations <- 1600
    n_thin <- 1
    n_burnin <- 600
    n_chains <- 4
    n_cores <- parallel::detectCores()
    delta = 0.9
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      
      list(mu_psi_0 = runif(1, -1, 1),
           sigma_psi_species = runif(1, 0, 0.5),
           sigma_psi_genus = runif(1, 0, 0.5),
           sigma_psi_site = runif(1, 0, 0.5),
           sigma_psi_ecoregion_three = runif(1, 0, 0.5),
           sigma_psi_ecoregion_one = runif(1, 0, 0.5),
           mu_psi_income = runif(1, -1, 1),
           sigma_psi_income = runif(1, 0, 0.5),
           mu_psi_herb_shrub_forest = runif(1, -1, 1),
           sigma_psi_herb_shrub_forest = runif(1, 0, 0.5),
           psi_site_area = runif(1, -1, 1),
           
           mu_p_citsci_0 = runif(1, -1, 0),
           sigma_p_citsci_species = runif(1, 0, 0.5),
           sigma_p_citsci_site = runif(1, 0, 0.5),
           sigma_p_citsci_ecoregion_three = runif(1, 0, 0.5),
           sigma_p_citsci_ecoregion_one = runif(1, 0, 0.5),
           p_citsci_interval = runif(1, -1, 1),
           p_citsci_pop_density = runif(1, -1, 1),
           
           mu_p_museum_0 = runif(1, -0.5, 0.5),
           sigma_p_museum_species = runif(1, 0, 0.1),
           sigma_p_museum_site = runif(1, 0, 0.1),
           sigma_p_museum_ecoregion_three = runif(1, 0, 0.1),
           sigma_p_museum_ecoregion_one = runif(1, 0, 0.1),
           p_museum_total_records = runif(1,  -0.5, 0.5)
           
      )
    )
    
  } else { 
    
    stan_data <- c("V_citsci", "V_museum", 
                 "ranges", "V_museum_NA",
                 "n_species", "n_sites", "n_intervals", "n_visits", 
                 "intervals", "species", "sites",
                 "n_genera", "genus_lookup",
                 "n_ecoregion_three", "n_ecoregion_one",
                 "ecoregion_three", "ecoregion_one",
                 "ecoregion_three_lookup", "ecoregion_one_lookup",
                 "pop_densities", "site_areas", "museum_total_records") 
  
  # Parameters monitored
  params <- c("mu_psi_0",
              "sigma_psi_species",
              "sigma_psi_genus",
              "sigma_psi_site",
              "sigma_psi_ecoregion_three",
              "sigma_psi_ecoregion_one",
              "psi_site_area",
              
              "mu_p_citsci_0",
              "sigma_p_citsci_species",
              "sigma_p_citsci_site",
              "sigma_p_citsci_ecoregion_three",
              "sigma_p_citsci_ecoregion_one",
              "p_citsci_interval",
              "p_citsci_pop_density", 
              
              "mu_p_museum_0",
              "sigma_p_museum_species",
              "sigma_p_museum_site",
              "sigma_p_museum_ecoregion_three",
              "sigma_p_museum_ecoregion_one",
              "p_museum_total_records",
              
              "psi_species",

              #"T_rep_citsci",
              #"T_obs_citsci",
              "P_species_citsci",
              
              #"T_rep_museum",
              #"T_obs_museum",
              "P_species_museum"
  )
  
  
  # MCMC settings
  n_iterations <- 1600
  n_thin <- 1
  n_burnin <- 600
  n_chains <- 4
  n_cores <- parallel::detectCores()
  delta = 0.9
  
  ## Initial values
  # given the number of parameters, the chains need some decent initial values
  # otherwise sometimes they have a hard time starting to sample
  inits <- lapply(1:n_chains, function(i)
    
    list(mu_psi_0 = runif(1, -1, 1),
         sigma_psi_species = runif(1, 0, 0.5),
         sigma_psi_genus = runif(1, 0, 0.5),
         sigma_psi_site = runif(1, 0, 0.5),
         sigma_psi_ecoregion_three = runif(1, 0, 0.5),
         sigma_psi_ecoregion_one = runif(1, 0, 0.5),
         psi_site_area = runif(1, -1, 1),
         
         mu_p_citsci_0 = runif(1, -1, 0),
         sigma_p_citsci_species = runif(1, 0, 0.5),
         sigma_p_citsci_site = runif(1, 0, 0.5),
         sigma_p_citsci_ecoregion_three = runif(1, 0, 0.5),
         p_citsci_interval = runif(1, -1, 1),
         p_citsci_pop_density = runif(1, -1, 1),
         
         mu_p_museum_0 = runif(1, -0.5, 0.5),
         sigma_p_museum_species = runif(1, 0, 0.1),
         sigma_p_museum_site = runif(1, 0, 0.1),
         sigma_p_museum_ecoregion_three = runif(1, 0, 0.1),
         p_museum_total_records = runif(1,  -0.5, 0.5)
         
    )
  )
  
  }
}


## --------------------------------------------------
### Run model

if(urban_sites == TRUE){
  stan_model <- paste0("./occupancy/models/model_", taxon, ".stan")
} else {
  stan_model <- paste0("./occupancy/models/model_", taxon, "_simple.stan")
}

## Call Stan from R
stan_out <- stan(stan_model,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 control=list(adapt_delta=delta),
                 seed = 3,
                 open_progress = FALSE,
                 cores = n_cores)

saveRDS(stan_out, paste0(
  "./occupancy/model_outputs/", taxon, "_", grid_size / 1000,
  "km_", min_population_size, "minpop_", 
  min_records_per_species, "minpersp_",
  n_intervals, "ints_", n_visits, "visits_",
  #"nonurban.RDS"  # use if saving a non-urban model run
  ".RDS"
)
)

stan_out <- readRDS(paste0(
  "./occupancy/model_outputs/", taxon, "_", grid_size / 1000, 
  "km_", min_population_size, "minpop", 
  min_records_per_species, "minpersp",
  n_intervals, "_", n_visits, ".RDS"
)
)

#stan_out <- readRDS("./occupancy/model_outputs/bombus_15km_1000minpop10minpersp4_3wide_priors.RDS")
#stan_out <- readRDS("./occupancy/model_outputs/bombus_15km_1000minpop10minpersp4_3_bbna.RDS")

  
# print main effects
print(stan_out, digits = 3, pars=
        c("mu_psi_0",
          "sigma_psi_site",
          "sigma_psi_ecoregion_three",
          "sigma_psi_ecoregion_one",
          "mu_psi_herb_shrub_forest",
          "sigma_psi_herb_shrub_forest",
          "mu_psi_income",
          "sigma_psi_income",
          "psi_site_area"))


print(stan_out, digits = 3, pars=
        c(
          "mu_p_citsci_0",
          "sigma_p_citsci_species",
          "sigma_p_citsci_site",
          "sigma_p_citsci_ecoregion_three",
          "sigma_p_citsci_ecoregion_one",
          "p_citsci_interval",
          "p_citsci_pop_density", 

          "mu_p_museum_0",
          "sigma_p_museum_species",
          "sigma_p_museum_site",
          "sigma_p_museum_ecoregion_three",
          "sigma_p_museum_ecoregion_one",
          "p_museum_total_records"))

View(as.data.frame(species_names))

# print sampled random effects
print(stan_out, digits = 3, pars=
        c("psi_species"))

print(stan_out, digits = 3, pars=
        c("psi_income"))

print(stan_out, digits = 3, pars=
        c("psi_herb_shrub_forest"))

# print sampled ppc
print(stan_out, digits = 3, pars=
        c("P_species_citsci"))

print(stan_out, digits = 3, pars=
        c("P_species_museum"))

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


## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out, pars = c(
  "mu_psi_0",
  "mu_psi_herb_shrub_forest",
  "mu_psi_income",
  "mu_p_citsci_0",
  "p_citsci_interval",
  "p_citsci_pop_density",
  "mu_p_museum_0",
  "p_museum_total_records"
))

# traceplot
traceplot(stan_out, pars = c(
  "sigma_psi_species",
  #"sigma_psi_genus",
  "sigma_psi_site",
  "sigma_psi_ecoregion_three",
  "sigma_psi_ecoregion_one",
  "sigma_psi_income",
  #"sigma_psi_herb_shrub_forest",
  "sigma_p_citsci_site",
  "sigma_p_citsci_ecoregion_three",
  "sigma_p_citsci_ecoregion_one",
  "sigma_p_museum_site",
  "sigma_p_museum_ecoregion_three",
  "sigma_p_museum_ecoregion_one",
  "sigma_p_citsci_species",
  "sigma_p_museum_species"
))

traceplot(stan_out, pars=
        c("psi_herb_shrub_forest[25]"))

# pairs plot
pairs(stan_out, pars = c(
  "mu_psi_0",
  "sigma_psi_species",
  #"mu_psi_interval",
  #"sigma_psi_interval",
  #"mu_psi_pop_density",
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

x=seq(0,3,1)
y=-4.5+0.4*x^2
plot(x,y, col='violet',type='o',lwd=2,lty=1)
