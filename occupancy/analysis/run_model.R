## run_model_integrated.R
### Run occupancy model (model_"taxon".stan), using real pollinator occurrence data from GBIF
# jcu; started nov 24, 2022

# Would you like to prepare data or run an analysis?
# select and enter the data collection choices for either SYRPHIDAE or BOMBUS below

## --------------------------------------------------
# input data preparation choices - SYRPHIDAE
# be careful that the (era_end - era_start) is evenly divisible by the n_intervals
era_start = 2014 # must define start date of the GBIF dataset
era_end = 2022 # must define start date of the GBIF dataset
n_intervals = 3 # must define number of intervals to break up the era into
n_visits = 3 # must define the number of repeat obs years within each interval
# note, should introduce throw error if..
# (era_end - era_start) / n_intervals has a remainder > 0,
min_records_per_species = 5 # filters species with less than this many records (total between both datasets)..
min_unique_detections = 2
# within the time span defined above
grid_size = 10000 # in metres so, e.g., 25000 = 25km x 25 km 
min_population_size = 1200 # min pop density in the grid cell (per km^2)

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
non_urban_subsample_n = 500 # if urban_sites is true, then how many sites do you want to keep? Keeping all will yield too much site data for computer to handle
infer_detections_at_genus = FALSE # default to FALSE # if true, infer non detections only for species in the same genus as a species detected (as opposed to any in the clade)
generate_temporal_plots = FALSE # default to FALSE

# define level 3 spatial clusters by CBSA metro area (by_city == TRUE)
# or by fine scale ecogeographic region (by_city == FALSE)
by_city = FALSE
remove_city_outliers_5stddev = TRUE

## --------------------------------------------------
# input data preparation choices - BOMBUS
# be careful that the (era_end - era_start) is evenly divisible by the n_intervals
era_start = 2011 # must define start date of the GBIF dataset
era_end = 2022 # must define start date of the GBIF dataset
n_intervals = 4 # must define number of intervals to break up the era into
n_visits = 3 # must define the number of repeat obs years within each interval
# note, should introduce throw error if..
# (era_end - era_start + 1) / n_intervals has a remainder > 0,
min_records_per_species = 5 # filters species with less than this many records (total between both datasets)..
min_unique_detections = 1 # filters species not detected at unique sites in unique years at/below this value
# within the time span defined above (is only from urban sites, should redefine to be from anywhere)
grid_size = 10000 # in metres so, e.g., 25000 = 25km x 25 km 
min_population_size = 1200 # min pop density in the grid cell (per km^2)

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

# define level 3 spatial clusters by CBSA metro area (by_city == TRUE)
# or by fine scale ecogeographic region (by_city == FALSE)
by_city = FALSE

## --------------------------------------------------
# Run the prep_data() function to gather new data (this takes about 5 minutes)
# if you've already gathered and presaved the data, then skip this and
# go to the next section (saveRDS() or readRDS())

source("./occupancy/data_prep/prep_data.R")

my_data <- prep_data(era_start = era_start, # must define start date of the GBIF dataset
                     era_end = era_end, # must define start date of the GBIF dataset
                     n_intervals = n_intervals, # must define number of intervals to break up the era into
                     n_visits = n_visits, # must define the number of repeat obs years within each interval
                     min_records_per_species = min_records_per_species,
                     min_unique_detections = min_unique_detections,
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
                     generate_temporal_plots,
                     by_city,
                     remove_city_outliers_5stddev
)


## --------------------------------------------------
# save data or load pre-saved data

# choose a directory for saving new data or loading old data
dir <- ("/")
if(by_city == TRUE) {
  dir <- "/by_city/"
} 
if(urban_sites == FALSE){
  dir <- "/non_urban/"
}

# save the data in case you want to make tweaks to the model run
# without redoing the data prep

# (default save) urban sites clustered in fine scale ecoregion
saveRDS(my_data, paste("./occupancy/analysis/prepped_data/", 
                       taxon, 
                       dir, 
                       grid_size / 1000, 
                       "km_", min_population_size, "minpop_", 
                       min_unique_detections, "minpersp_",
                       n_intervals, "ints_", n_visits, "visits",
                       ".rds", sep = ""))

# load pre-saved data
my_data <- readRDS(paste0("./occupancy/analysis/prepped_data/",
                          taxon, 
                          dir,
                          grid_size / 1000, "km_",
                          min_population_size, "minpop_",
                          min_unique_detections, "minpersp", "_",
                          n_intervals, "ints_",
                          n_visits, "visits",
                          ".rds"))


## --------------------------------------------------
# Once you have the data you want, you will need to format it
# appropriately to prepare it for STAN to read.

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

pop_densities <- my_data$pop_densities
avg_income <- my_data$avg_income
open_developed <- my_data$developed_open
herb_shrub <- my_data$herb_shrub_cover
site_areas <- my_data$site_areas
herb_shrub_forest <- my_data$herb_shrub_forest
developed_med_high <- my_data$developed_med_high
museum_total_records <- my_data$museum_total_records

level_three_cluster <- my_data$elevel_three_vector
n_level_three <- my_data$n_level_three
level_three_lookup <- my_data$level_three_lookup

ecoregion_one <- my_data$ecoregion_one_vector
n_ecoregion_one <- my_data$n_ecoregion_one
ecoregion_one_lookup <- my_data$ecoregion_one_lookup

level_three_names <- my_data$level_three_names
#View(as.data.frame(level_three_names))
level_three_names_unique <- unique(level_three_names)
#View(as.data.frame(level_three_names_unique))

#saveRDS(species_names, "./figures/species_names/syrphidae_names_50km_nonurban.RDS")
#saveRDS(species_names, "./figures/species_names/bombus_names_10km_urban.RDS")
#write.csv(as.data.frame(species_names), "./data/syrphidae_names.csv")

#saveRDS(as.data.frame(cbind(CBSA_names, CBSA_lookup)), "./figures/species_names/CBSA_names_10km_urban.RDS")

species_counts <- my_data$species_counts # number of species detections in time period at urban grid cells
species_detections <- my_data$species_detections # number of unique species/site/year (in time period at urban grid cells) occasions with detections > 0
species_counts_full <- my_data$species_counts_full # number of records anywhere in time frame (not just urban grid cells)

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
sites <- seq(1, n_sites, by=1)
species <- seq(1, n_species, by=1)

# add nativity for syrphids
if(taxon == "syrphidae"){
  
  library(tidyverse)
  df <- read.csv("./data/occurrence_data/syrphidae_nativity.csv") 
  df <- select(df, species, nativity)
  species_df <- as.data.frame(species_names) %>% 
    rename("species" = "species_names")
  
  species_df <- left_join(species_df, df, by ="species")
  species_df <- species_df %>% dplyr::mutate(nativity = replace_na(nativity, 1))
  
  nativity <- pull(species_df, nativity)

  detach("package:tidyverse", unload = TRUE)

}


rm(my_data)
gc()

## --------------------------------------------------
# Once you have the data formatted for STAN..
# you will also need to enter HMC options for STAN 
# (e.g., which params to monitor, initial values specified, 
# how many chains to run, how many iterations, etc.)

# There are for possible model sets:
# bombus urban, 
# bombus non-urban,
# syrphidae urban,
# or syrphidae non-urban

# if you've selected the data collection requisites desired above then just 
# press enter at the if statement and the HMC options will be loaded into the environment
# you will need to manually edit the options within the if statements if you want 
# to change things such as the number of chains to run or number of iterations

if(taxon == "bombus"){
  
  if(urban_sites == TRUE){
   
    stan_data <- c("V_citsci", "V_museum", 
                   "ranges", "V_museum_NA",
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_level_three", 
                   "level_three_lookup", 
                   "n_ecoregion_one",
                   "ecoregion_one_lookup",
                   "pop_densities", "site_areas", "avg_income", 
                   "herb_shrub_forest", "museum_total_records") 
    
    # Parameters monitored
    params <- c("sigma_species_detection",
                "species_intercepts_detection",
                "rho",
                
                "mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_site",
                "sigma_psi_level_three",
                "sigma_psi_ecoregion_one",
                "mu_psi_income",
                "sigma_psi_income",
                "mu_psi_herb_shrub_forest",
                "sigma_psi_herb_shrub_forest",
                "psi_site_area",
                
                "mu_p_citsci_0",
                "sigma_p_citsci_site",
                "sigma_p_citsci_level_three",
                "sigma_p_citsci_ecoregion_one",
                "p_citsci_interval",
                "p_citsci_pop_density", 
                
                "mu_p_museum_0",
                "sigma_p_museum_site",
                "sigma_p_museum_level_three",
                "sigma_p_museum_ecoregion_one",
                "p_museum_total_records",
                
                "psi_species",
                "psi_income",
                "psi_herb_shrub_forest",
                
                "psi_level_three", # track city or eco3 effects
                "psi_ecoregion_one",
                
                #"T_rep_citsci",
                #"T_obs_citsci",
                "P_species_citsci",
                
                #"T_rep_museum",
                #"T_obs_museum",
                "P_species_museum"
    )
    
    
    # MCMC settings
    n_iterations <- 500
    n_thin <- 1
    n_burnin <- 250
    n_chains <- 4
    #n_cores <- parallel::detectCores()
    n_cores <- 4
    delta = 0.9
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      
      list(
        rho = runif(1, 0, 1),
        
        mu_psi_0 = runif(1, -1, 1),
        sigma_psi_species = runif(1, 0, 1),
        sigma_psi_site = runif(1, 0, 1),
        sigma_psi_level_three = runif(1, 0, 1),
        sigma_psi_ecoregion_one = runif(1, 0, 1),
        mu_psi_income = runif(1, -1, 1),
        sigma_psi_income = runif(1, 0, 1),
        mu_psi_herb_shrub_forest = runif(1, -1, 1),
        sigma_psi_herb_shrub_forest = runif(1, 0, 1),
        psi_site_area = runif(1, -1, 1),
        
        mu_p_citsci_0 = runif(1, -1, 0),
        sigma_p_citsci_site = runif(1, 0, 0.5),
        sigma_p_citsci_level_three = runif(1, 0, 0.5),
        sigma_p_citsci_ecoregion_one = runif(1, 0, 0.5),
        p_citsci_interval = runif(1, 0, 1),
        p_citsci_pop_density = runif(1, -1, 1),
        
        # start musuem values close to zero
        mu_p_museum_0 = runif(1, -0.5, 0.5),
        sigma_p_museum_site = runif(1, 0, 0.25),
        sigma_p_museum_level_three = runif(1, 0, 0.25),
        sigma_p_museum_ecoregion_one = runif(1, 0, 0.25),
        p_museum_total_records = runif(1, -0.5, 0.5)    
      )
    )
      
  } else { # bombus non-urban
    
    stan_data <- c("V_citsci", "V_museum", 
                   "ranges", "V_museum_NA",
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_level_three", 
                   "level_three_lookup", 
                   "n_ecoregion_one",
                   "ecoregion_one_lookup",
                   "ecoregion_three_lookup", "ecoregion_one_lookup",
                   "pop_densities", "site_areas","museum_total_records") 
    
    # Parameters monitored
    params <- c("mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_site",
                "sigma_psi_level_three",
                "sigma_psi_ecoregion_one",
                "psi_site_area",
                
                "mu_p_citsci_0",
                "sigma_p_citsci_species",
                "sigma_p_citsci_site",
                "sigma_p_citsci_level_three",
                "sigma_p_citsci_ecoregion_one",
                "p_citsci_interval",
                "p_citsci_pop_density", 
                
                "mu_p_museum_0",
                "sigma_p_museum_species",
                "sigma_p_museum_site",
                "sigma_p_museum_level_three",
                "sigma_p_museum_ecoregion_one",
                "p_museum_total_records",
                
                "psi_species",
                
                "psi_level_three", # track city effects
                "psi_ecoregion_one",
                
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
           sigma_psi_level_three = runif(1, 0, 0.5),
           sigma_psi_ecoregion_one = runif(1, 0, 0.5),
           psi_site_area = runif(1, -1, 1),
           
           mu_p_citsci_0 = runif(1, -1, 0),
           sigma_p_citsci_species = runif(1, 0, 0.5),
           sigma_p_citsci_site = runif(1, 0, 0.5),
           sigma_p_citsci_level_three = runif(1, 0, 0.5),
           sigma_p_citsci_ecoregion_one = runif(1, 0, 0.5),
           p_citsci_interval = runif(1, -1, 1),
           p_citsci_pop_density = runif(1, -1, 1),
           
           mu_p_museum_0 = runif(1, -0.5, 0.5),
           sigma_p_museum_species = runif(1, 0, 0.1),
           sigma_p_museum_site = runif(1, 0, 0.1),
           sigma_p_museum_level_three = runif(1, 0, 0.1),
           sigma_p_museum_ecoregion_one = runif(1, 0, 0.1),
           p_museum_total_records = runif(1,  -0.5, 0.5)
           
      )
    )
    
  }
  
} else { # taxon == syrphidae, urban
  
  if(urban_sites == TRUE){
    
    stan_data <- c("V_citsci", "V_museum", 
                   "ranges", "V_museum_NA",
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_genera", "genus_lookup",
                   "n_level_three", 
                   "level_three_lookup", 
                   "n_ecoregion_one",
                   "ecoregion_one_lookup",
                   "pop_densities", "site_areas", 
                   "avg_income", 
                   "nativity",
                   "herb_shrub_forest" 
                   #"museum_total_records"
                   ) 
    
    # Parameters monitored
    #params <- c(#"rho",
                #"sigma_species_detection",
      
                #"mu_psi_0",
                #"sigma_psi_species",
                #"sigma_psi_genus",
                #"sigma_psi_site",
                #"sigma_psi_level_three",
                #"sigma_psi_ecoregion_one",
                #"mu_psi_income",
                #"sigma_psi_income",
                #"mu_psi_herb_shrub_forest",
                #"delta0",
                #"delta1",
                #"sigma_psi_herb_shrub_forest",
                #"psi_site_area",
                
                #"mu_p_citsci_0",
                #"sigma_p_citsci_site",
                #"sigma_p_citsci_level_three",
                #"sigma_p_citsci_ecoregion_one",
                #"p_citsci_interval",
                #"p_citsci_pop_density", 
                
                #"mu_p_museum_0",
                #"sigma_p_museum_site",
                #"sigma_p_museum_level_three",
                #"sigma_p_museum_ecoregion_one",
                #"p_museum_total_records",
                
                #"psi_species",
                #"psi_income",
                #"psi_herb_shrub_forest",
                
                #"psi_level_three", # track city/fine-ecoregion effects
                #"psi_ecoregion_one", # track broad eco effects
                
                #"T_rep_citsci",
                #"T_obs_citsci",
                #"P_species_citsci",
                
                #"T_rep_museum",
                #"T_obs_museum",
                #"P_species_museum",
                
                #"mu_psi_nat_habitat_native",
                #"mu_psi_nat_habitat_nonnative"
    #)
    
    # Parameters monitored
    params <- c(#"rho",
                #"sigma_species_intercepts",
                
                "mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_genus",
                "sigma_psi_site",
                "sigma_psi_level_three",
                "sigma_psi_ecoregion_one",
                "mu_psi_herb_shrub_forest",
                "delta0",
                "delta1",
                "sigma_psi_herb_shrub_forest",
                "gamma0",
                "gamma1",
                "psi_site_area",
                
                "mu_p_citsci_0",
                "sigma_p_citsci_species",
                "sigma_p_citsci_site",
                "sigma_p_citsci_level_three",
                "sigma_p_citsci_ecoregion_one",
                "p_citsci_interval",
                "p_citsci_pop_density", 
                
                #"mu_p_museum_0",
                #"sigma_p_museum_site",
                #"sigma_p_museum_level_three",
                #"sigma_p_museum_ecoregion_one",
                #"p_museum_total_records",
                
                "psi_species",
                #"species_intercepts",
                #"psi_income",
                "psi_herb_shrub_forest",
                
                "psi_level_three", # track city/fine-ecoregion effects
                "psi_ecoregion_one", # track broad eco effects
                
                #"T_rep_citsci",
                #"T_obs_citsci",
                "P_species_citsci",
                
                #"T_rep_museum",
                #"T_obs_museum",
                #"P_species_museum",
                
                "mu_psi_nat_habitat_native",
                "mu_psi_nat_habitat_nonnative"
    )
    
    
    # MCMC settings
    n_iterations <- 300
    n_thin <- 1
    n_burnin <- 150
    n_chains <- 4
    n_cores <- 4
    #n_cores <- parallel::detectCores()
    delta = 0.9
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    set.seed(1)
    inits <- lapply(1:n_chains, function(i)
      
      list(
            mu_psi_0 = runif(1, -0.5, 0.5),
            sigma_psi_site = runif(1, 3, 4),
            sigma_psi_level_three = runif(1, 0, 1),
            sigma_psi_ecoregion_one = runif(1, 2, 3),
            #mu_psi_herb_shrub_forest = runif(1, 0, 0.5),
            delta0 = runif(1, -0.5, 0.5),
            delta1 = runif(1, 0, 0.5),
            gamma0 = runif(1, 0.5, 0.75),
            gamma1 = runif(1, 0, 0.1),
            sigma_psi_herb_shrub_forest = runif(1, 0.75, 1),
            psi_site_area = runif(1, -0.5, 0.5),
            
            mu_p_citsci_0 = runif(1, -3.5, -2.5),
            sigma_p_citsci_site = runif(1, 0, 2),
            sigma_p_citsci_level_three = runif(1, 0, 2),
            sigma_p_citsci_ecoregion_one = runif(1, 0, 0.5),
            p_citsci_interval = runif(1, 0.5, 0.6),
            p_citsci_pop_density = runif(1, 0.4, 0.6)
            
            # start musuem values close to zero
            #mu_p_museum_0 = runif(1, -3, -1),
            #sigma_p_museum_site = runif(1, 0, 0.25),
            #sigma_p_museum_level_three = runif(1, 0, 0.25),
            #sigma_p_museum_ecoregion_one = runif(1, 0, 0.25),
            #p_museum_total_records = runif(1, -0.5, 0.5)
           
      )
    )
    
  } else { 
    
    stan_data <- c("V_citsci", "V_museum", 
                 "ranges", "V_museum_NA",
                 "n_species", "n_sites", "n_intervals", "n_visits", 
                 "intervals", "species", "sites",
                 "n_genera", "genus_lookup",
                 "n_level_three", 
                 "level_three_lookup", 
                 "n_ecoregion_one",
                 "ecoregion_one_lookup",
                 "ecoregion_three_lookup", "ecoregion_one_lookup",
                 "pop_densities", "site_areas", "museum_total_records") 
  
  # Parameters monitored
  params <- c("mu_psi_0",
              "sigma_psi_species",
              "sigma_psi_genus",
              "sigma_psi_site",
              "sigma_psi_level_three",
              "sigma_psi_ecoregion_one",
              "psi_site_area",
              
              "mu_p_citsci_0",
              "sigma_p_citsci_species",
              "sigma_p_citsci_site",
              "sigma_p_citsci_level_three",
              "sigma_p_citsci_ecoregion_one",
              "p_citsci_interval",
              "p_citsci_pop_density", 
              
              "mu_p_museum_0",
              "sigma_p_museum_species",
              "sigma_p_museum_site",
              "sigma_p_museum_level_three",
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
  n_iterations <- 600
  n_thin <- 1
  n_burnin <- 300
  n_chains <- 4
  n_cores <- parallel::detectCores()
  delta = 0.9
  
  ## Initial values
  # given the number of parameters, the chains need some decent initial values
  # otherwise sometimes they have a hard time starting to sample
  inits <- lapply(1:n_chains, function(i)
    
    list(mu_psi_0 = runif(1, 0, 1),
         sigma_psi_species = runif(1, 0, 0.5),
         sigma_psi_genus = runif(1, 0, 0.5),
         sigma_psi_site = runif(1, 0, 0.5),
         sigma_psi_level_three = runif(1, 0, 0.5),
         sigma_psi_ecoregion_one = runif(1, 0, 0.5),
         psi_site_area = runif(1, -1, 1),
         
         mu_p_citsci_0 = runif(1, -1, 0),
         sigma_p_citsci_species = runif(1, 0, 0.5),
         sigma_p_citsci_site = runif(1, 0, 0.5),
         sigma_p_citsci_level_three = runif(1, 0, 0.5),
         p_citsci_interval = runif(1, -1, 1),
         p_citsci_pop_density = runif(1, -1, 1),
         
         mu_p_museum_0 = runif(1, -0.5, 0.5),
         sigma_p_museum_species = runif(1, 0, 0.1),
         sigma_p_museum_site = runif(1, 0, 0.1),
         sigma_p_museum_level_three = runif(1, 0, 0.1),
         p_museum_total_records = runif(1,  -0.5, 0.5)
         
    )
  )
  
  }
}


## --------------------------------------------------
### Run model

# load appropriate model file from the directory
if(urban_sites == TRUE){
  if(by_city == FALSE){
    stan_model <- paste0("./occupancy/models/model_", taxon, "_covariance_no_rc.stan")
  } else {
    stan_model <- paste0("./occupancy/models/model_", taxon, "_covariance_by_city.stan")
  }
} else {
  stan_model <- paste0("./occupancy/models/model_", taxon, "_simple.stan")
}

# for bombus 20km and 25km which required narrower prior for mu_psi_0 identifiability
# stan_model <- paste0("./occupancy/models/model_", taxon, "_2.stan")

## Call Stan from R
set.seed(1)
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
  "./occupancy/model_outputs/", taxon, dir, 
  taxon, "_",
  grid_size / 1000,
  "km_", min_population_size, "minpop_", 
  min_records_per_species, "minpersp_",
  n_intervals, "ints_", n_visits, "visits_",
  ".rds"
)
)

stan_out <- readRDS(paste0(
  "./occupancy/model_outputs/", taxon, dir, taxon, "_", grid_size / 1000, 
  "km_", min_population_size, "minpop_", 
  min_records_per_species, "minpersp_",
  n_intervals, "ints_", n_visits, "visits_", 
  ".rds"
)
)

stan_out <- readRDS("./occupancy/model_outputs/syrphidae/old_results/syrphidae_10km_1200minpop_5minpersp_4ints_3visits_.RDS")

# print main effects
print(stan_out, digits = 3, pars=
        c("mu_psi_0",
          "sigma_psi_species",
          "sigma_psi_site",
          "sigma_psi_level_three",
          "sigma_psi_ecoregion_one",
          #"mu_psi_herb_shrub_forest",
          "delta0",
          "delta1",
          "mu_psi_nat_habitat_native",
          "sigma_psi_herb_shrub_forest",
          #"psi_income",
          #mu_psi_income",
          #"sigma_psi_income",
          "psi_site_area"))


print(stan_out, digits = 3, pars=
        c(
          "mu_p_citsci_0",
          #"sigma_p_citsci_species",
          "sigma_p_citsci_site",
          "sigma_p_citsci_level_three",
          "sigma_p_citsci_ecoregion_one",
          "p_citsci_interval",
          "p_citsci_pop_density"

          #"mu_p_museum_0",
          #"sigma_p_museum_species",
          #"sigma_p_museum_site",
          #"sigma_p_museum_ecoregion_three",
          #"sigma_p_museum_ecoregion_one",
          #"p_museum_total_records"
          )
          )

print(stan_out, digits = 3, pars=
        c("rho",
          "sigma_species_detection"))

View(as.data.frame(species_names))

View(as.data.frame(stan_out, pars=
                     c("psi_CBSA")))

# print sampled random effects
print(stan_out, digits = 3, pars=
        c("species_intercepts[1,1]"))

print(stan_out, digits = 3, pars=
        c("psi_income"))

print(stan_out, digits = 3, pars=
        c("mu_psi_nat_habitat_native"))

print(stan_out, digits = 3, pars=
        c("psi_herb_shrub_forest"))

# print sampled ppc
print(stan_out, digits = 3, pars=
        c("P_species_citsci"))

print(stan_out, digits = 3, pars=
        c("P_species_museum"))


## --------------------------------------------------
### Simple diagnostic plots

# traceplot
traceplot(stan_out, pars = c(
  "mu_psi_0",
  "psi_site_area",
  #"mu_psi_herb_shrub_forest",
  "delta0",
  "delta1",
  "gamma0",
  "gamma1",
  #"psi_income",
  #"mu_psi_income",
  "mu_p_citsci_0",
  "p_citsci_interval",
  "p_citsci_pop_density"
  #"mu_p_museum_0",
  #"p_museum_total_records"
))

traceplot(stan_out, pars = c(
  "rho",
  "sigma_species_intercepts"
))

# traceplot
traceplot(stan_out, pars = c(
  #"sigma_psi_species",
  #"sigma_psi_genus",
  "sigma_psi_site",
  "sigma_psi_level_three",
  #"sigma_psi_ecoregion_three",
  "sigma_psi_ecoregion_one",
  #"sigma_psi_income",
  #"sigma_psi_herb_shrub_forest",
  "sigma_p_citsci_site",
  "sigma_p_citsci_level_three",
  "sigma_p_citsci_ecoregion_one"
  #"sigma_p_museum_site",
  #"sigma_p_museum_level_three",
  #"sigma_p_museum_ecoregion_one"#,
  #"sigma_p_citsci_species",
  #"sigma_p_museum_species"
))

traceplot(stan_out, pars=
            c("mu_psi_nat_habitat_native"))

traceplot(stan_out, pars=
            c("mu_psi_nat_habitat_nonnative"))

traceplot(stan_out, pars=
        c("psi_herb_shrub_forest[25]"))

# pairs plot
pairs(stan_out, pars = c(
  "mu_psi_0",
  "sigma_psi_species",
  #"mu_psi_interval",
  #"sigma_psi_interval",
  #"mu_psi_herb_shrub_forest",
  "psi_site_area",
  
  "mu_p_citsci_0",
  #"sigma_p_citsci_species",
  #"sigma_p_citsci_site",
  "p_citsci_interval",
  #"p_citsci_pop_density", 
  
  "mu_p_museum_0",
  #"sigma_p_museum_species",
  "sigma_p_museum_site"
))

x=seq(0,3,1)
y=-4.5+0.4*x^2
plot(x,y, col='violet',type='o',lwd=2,lty=1)


# plot species detections from species counts
# draw plot
ggplot(data = species_counts, aes(x = reorder(species, total_count), y = total_count)) + 
  geom_bar(stat = "identity", show.legend = FALSE) + 
  labs(x = "", 
       y = "Number of detections in urban landscapes (2011-2022)") +
  coord_flip()
