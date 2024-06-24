## run_model_integrated.R
### Run occupancy model (model_"taxon".stan), using real pollinator occurrence data from GBIF
# jcu; started nov 24, 2022

# Would you like to prepare data or run an analysis?
# select and enter the data collection choices for either SYRPHIDAE or BOMBUS below
taxon = "syrphidae" # taxon to analyze, either "syrphidae" or "bombus"
taxon = "bombus" # taxon to analyze, either "syrphidae" or "bombus"

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
min_unique_detections = 2 # species must have more detections at more site/years than this (so 2 == 3 or more unique detections)
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
non_urban_subsample_n = 550 # if urban_sites is true, then how many sites do you want to keep? Keeping all will yield too much site data for computer to handle
infer_detections_at_genus = FALSE # default to FALSE # if true, infer non detections only for species in the same genus as a species detected (as opposed to any in the clade)
generate_temporal_plots = FALSE # default to FALSE

# define level 3 spatial clusters by CBSA metro area (by_city == TRUE)
# or by fine scale ecogeographic region (by_city == FALSE)
by_city = FALSE
remove_city_outliers_5stddev = TRUE
canopy_cover_model = FALSE

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
non_urban_subsample_n = 550 # if urban_sites is true, then how many sites do you want to keep? Keeping all will yield too much site data for computer to handle
infer_detections_at_genus = FALSE # default to FALSE # if true, infer non detections only for species in the same genus as a species detected (as opposed to any in the clade)
generate_temporal_plots = FALSE # default to FALSE

# define level 3 spatial clusters by CBSA metro area (by_city == TRUE)
# or by fine scale ecogeographic region (by_city == FALSE)
by_city = FALSE
# filter out urban landscapes with inordinately high population density (i.e. downtown NYC)
remove_city_outliers_5stddev = TRUE
canopy_cover_model = FALSE

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
                     remove_city_outliers_5stddev,
                     canopy_cover_model
                     
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
                       min_unique_detections, "minUniqueDetections_",
                       n_intervals, "ints_", n_visits, "visits",
                       ".rds", sep = ""))

# load pre-saved data
my_data <- readRDS(paste0("./occupancy/analysis/prepped_data/",
                          taxon, 
                          dir,
                          grid_size / 1000, "km_",
                          min_population_size, "minpop_",
                          min_unique_detections, "minUniqueDetections", "_",
                          n_intervals, "ints_",
                          n_visits, "visits",
                          ".rds"))

# or read in some data manually
#my_data <- readRDS("./occupancy/analysis/prepped_data/bombus/non_urban/40km_1000minpop_1minUniqueDetections_4ints_3visits.rds")

## --------------------------------------------------
# Once you have the data you want, you will need to format it
# appropriately to prepare it for STAN to read.

# best to restart R or offload all of the spatial data packages before running the model
gc()

library(rstan)
# use non-centered spatial random effects model?
use_reparameterized_rand_effects_model = FALSE # use non-centored spatial random effects? default = FALSE
bombus_no_rc = FALSE # do not fit integrated model?  default = FALSE

# data to feed to the model
# cs = community science / "citizen science"
# rc = research collections / "museum"
V_cs <- my_data$V_citsci # citizen science detection data
V_rc <- my_data$V_museum # museum detection data
ranges <- my_data$V_citsci_NA # cit science NA indicator array
V_rc_NA <- my_data$V_museum_NA # museum data NA indicator array
n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_intervals <- my_data$n_intervals # number of surveys 

n_visits <- my_data$n_visits

interval_names <- as.vector(as.numeric(my_data$intervals))
site_names <- my_data$sites
species_names <- my_data$species

pop_densities <- my_data$pop_densities
avg_income <- my_data$avg_income
avg_racial_minority <- my_data$avg_minority
open_developed <- my_data$developed_open
herb_shrub <- my_data$herb_shrub_cover
site_areas <- my_data$site_areas
natural_habitat <- my_data$herb_shrub_forest
developed_med_high <- my_data$developed_med_high
rc_total_records <- my_data$museum_total_records

# level three is either CBSA statistical metro areas or ecoregion 3 (depending on preset options above)
level_three_cluster <- my_data$elevel_three_vector
n_level_three <- my_data$n_level_three
level_three_lookup <- my_data$level_three_lookup

# level four is always ecoregion 1
level_four <- my_data$ecoregion_one_vector
n_level_four <- my_data$n_ecoregion_one
level_four_lookup <- my_data$ecoregion_one_lookup
level_four_lookup_by_site <- my_data$ecoregion_one_lookup_by_grid_cell

level_three_names <- my_data$level_three_names
level_three_names_unique <- unique(level_three_names)
#View(as.data.frame(level_three_lookup))
#View(as.data.frame(level_three_names))
#View(as.data.frame(level_three_names_unique))
level_four_names <- my_data$ecoregion_one_names
#View(as.data.frame(level_four_names))
level_four_names_unique <- unique(level_four_names)
#View(as.data.frame(level_four_lookup))
#View(as.data.frame(level_three_names_unique))

CBSA_names <- my_data$CBSA_names
#View(as.data.frame(unique(CBSA_names)))

#saveRDS(species_names, "./figures/species_names/syrphidae_names_10km_urban.RDS")
#saveRDS(species_names, "./figures/species_names/bombus_names_40km_nonurban.RDS")
#write.csv(as.data.frame(species_names), "./data/syrphidae_names.csv")

#saveRDS(as.data.frame(cbind(CBSA_names, CBSA_lookup)), "./figures/species_names/CBSA_names_10km_urban.RDS")

species_counts <- my_data$species_counts # number of species detections in time period at urban grid cells
species_detections <- my_data$species_detections # number of unique species/site/year (in time period at urban grid cells) occasions with detections > 0
species_counts_full <- my_data$species_counts_full # number of records anywhere in time frame (not just urban grid cells)

#sum(species_counts$total_count)
#sum(species_detections$museum_detections, na.rm = TRUE)

genus_lookup <- my_data$genus_lookup
genus_lookup <- as.numeric(factor(genus_lookup))
n_genera <- my_data$n_genera
genera <- my_data$genera

raw_pop_density <- my_data$raw_pop_density
raw_avg_income <- my_data$raw_avg_income

total_records <- my_data$total_records_since_2000
total_records_since_study <- my_data$total_records_since_2000
cs_records <- my_data$citsci_records
cs_detections <- my_data$citsci_detections
rc_records <- my_data$museum_records
rc_detections <- my_data$museum_detections

total_records_full <- my_data$total_records_since_2000_full
total_records_since_study_full <- my_data$total_records_since_2000_full
cs_records_full <- my_data$citsci_records_full
rc_records_full <- my_data$museum_records_full

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
  library(dplyr) # have to reload this for cluster computing
  df <- read.csv("./data/occurrence_data/syrphidae_nativity.csv") 
  df <- select(df, species, nativity)
  species_df <- as.data.frame(species_names) %>% 
    rename("species" = "species_names")
  
  species_df <- left_join(species_df, df, by ="species")
  #species_df <- species_df %>% dplyr::mutate(nativity = replace_na(nativity, 1))
  species_df$nativity[is.na(species_df$nativity)] <- 1 # cluster won't recognize replace_na()
  
  nativity <- pull(species_df, nativity)
  
  #saveRDS(nativity, "./figures/species_names/syrphidae_nativity_10km_urban.RDS")
  
  detach("package:tidyverse", unload = TRUE)

}


rm(my_data)
gc()

#saveRDS(nativity, "./figures/species_names/syrphidae_nativity_10km_urban.RDS")

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
    
    if(bombus_no_rc == FALSE){
   
    stan_data <- c("V_cs", "V_rc", 
                   "ranges", "V_rc_NA",
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_level_three", 
                   "level_three_lookup", 
                   "n_level_four",
                   "level_four_lookup",
                   "pop_densities", "site_areas", 
                   "avg_income", "avg_racial_minority",
                   "open_developed", "natural_habitat") 
    
    # Parameters monitored
    params <- c("sigma_species_detection",
                "species_intercepts_detection",
                "rho",
                
                "mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_site",
                "sigma_psi_level_three",
                "sigma_psi_level_four",
                "mu_psi_income",
                "mu_psi_race",
                "mu_psi_natural_habitat",
                "mu_psi_open_developed",
                "sigma_psi_natural_habitat",
                "psi_site_area",
                
                "mu_p_cs_0",
                "sigma_p_cs_site",
                "sigma_p_cs_level_three",
                "p_cs_interval",
                "p_cs_pop_density",
                "p_cs_income",
                
                "mu_p_rc_0",
                "sigma_p_rc_site",
                "sigma_p_cs_level_three",
                
                "psi_species",
                "psi_natural_habitat",
                
                "psi_site",
                
                "W_species_rep_cs",
                "W_species_rep_rc"
    )
    
    # MCMC settings
    n_iterations <- 4000
    n_thin <- 1
    n_burnin <- 2000
    n_chains <- 4
    #n_cores <- parallel::detectCores()
    n_cores <- 4
    delta = 0.95
    max_treedepth = 10
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      
      list(
        
        mu_psi_0 = runif(1, -1, 1),
        sigma_psi_species = runif(1, 0, 1),
        sigma_psi_site = runif(1, 1, 2),
        sigma_psi_level_three = runif(1, 0, 1),
        sigma_psi_level_four = runif(1, 0, 1),
        mu_psi_income = runif(1, -1, 1),
        mu_psi_natural_habitat = runif(1, -1, 1),
        sigma_psi_natural_habitat = runif(1, 0, 1),
        psi_site_area = runif(1, -1, 1),
        mu_psi_open_developed = runif(1, -1, 1),
        mu_psi_income = runif(1, -1, 1),
        mu_psi_race = runif(1, -1, 1),
        
        mu_p_cs_0 = runif(1, -1, 0),
        sigma_p_cs_site = runif(1, 0.5, 1),
        sigma_p_cs_level_three = runif(1, 0, 0.5),
        p_cs_interval = runif(1, 0, 1),
        p_cs_pop_density = runif(1, -1, 1),
        
        mu_p_rc_0 = runif(1, -0.5, 0.5),
        sigma_p_rc_site = runif(1, 0, 0.5),
        sigma_p_rc_level_three = runif(1, 0, 0.5)
      )
    )
    
    } else { # do not use fully integrated model
      
      stan_data <- c("V_cs", "V_rc", 
                     "ranges", "V_rc_NA",
                     "n_species", "n_sites", "n_intervals", "n_visits", 
                     "intervals", "species", "sites",
                     "n_level_three", 
                     "level_three_lookup", 
                     "n_level_four",
                     "level_four_lookup",
                     "pop_densities", "site_areas", "avg_income", 
                     "natural_habitat", "rc_total_records") 
      
      # Parameters monitored
      params <- c(
                  
                  "mu_psi_0",
                  "sigma_psi_species",
                  "sigma_psi_site",
                  "sigma_psi_level_three",
                  "sigma_psi_level_four",
                  "mu_psi_income",
                  "sigma_psi_income",
                  "mu_psi_natural_habitat",
                  "sigma_psi_natural_habitat",
                  "psi_site_area",
                  
                  "mu_p_cs_0",
                  "sigma_p_cs_species",
                  "sigma_p_cs_site",
                  "sigma_p_cs_level_three",
                  "sigma_p_cs_level_four",
                  "p_cs_interval",
                  "p_cs_pop_density",
                  
                  "psi_species",
                  "psi_income",
                  "psi_natural_habitat",
                  
                  "psi_site",
                  "psi_level_four",
                  "psi_level_three", # track city or eco3 effects
                  
                  "T_rep_cs",
                  "T_obs_cs",
                  "P_species_cs"
      )
      
      # MCMC settings
      n_iterations <- 400
      n_thin <- 1
      n_burnin <- 200
      n_chains <- 4
      #n_cores <- parallel::detectCores()
      n_cores <- 4
      delta = 0.95
      
      ## Initial values
      # given the number of parameters, the chains need some decent initial values
      # otherwise sometimes they have a hard time starting to sample
      inits <- lapply(1:n_chains, function(i)
        
        list(
          
          mu_psi_0 = runif(1, -1, 1),
          sigma_psi_species = runif(1, 0, 1),
          sigma_psi_site = runif(1, 1, 2),
          sigma_psi_level_three = runif(1, 0, 1),
          sigma_psi_level_four = runif(1, 0, 1),
          mu_psi_income = runif(1, -1, 1),
          sigma_psi_income = runif(1, 0, 1),
          mu_psi_natural_habitat = runif(1, -1, 1),
          sigma_psi_natural_habitat = runif(1, 0, 1),
          psi_site_area = runif(1, -1, 1),
          
          mu_p_cs_0 = runif(1, -1, 0),
          sigma_p_cs_site = runif(1, 0.5, 1),
          sigma_p_cs_level_three = runif(1, 0, 0.5),
          sigma_p_cs_ecoregion_one = runif(1, 0, 0.5),
          p_cs_interval = runif(1, 0, 1),
          p_cs_pop_density = runif(1, -1, 1)

        )
      )
    }
      
  } else { # bombus non-urban
    
    stan_data <- c("V_cs", "V_rc", 
                   "ranges", "V_rc_NA",
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_level_three", 
                   "level_three_lookup", 
                   "n_level_four",
                   "level_four_lookup",
                   "pop_densities", "site_areas", 
                   "rc_total_records") 
    
    # Parameters monitored
    params <- c("sigma_species_detection",
                "species_intercepts_detection",
                "rho",
                
                "mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_site",
                "sigma_psi_level_three",
                "sigma_psi_level_four",
                "psi_site_area",
                
                "mu_p_cs_0",
                "sigma_p_cs_site",
                "sigma_p_cs_level_three",
                "sigma_p_cs_level_four",
                "p_cs_interval",
                "p_cs_pop_density", 
                
                "mu_p_rc_0",
                "sigma_p_rc_site",
                #"sigma_p_rc_level_three",
                #"sigma_p_rc_level_four",
                "p_rc_total_records",
                
                "psi_species",

                "psi_site",
                "psi_level_four",
                "psi_level_three", # track city or eco3 effects
                
                "T_rep_cs",
                "T_obs_cs",
                "P_species_cs",
                
                "T_rep_rc",
                "T_obs_rc",
                "P_species_rc"
    )
    
    # MCMC settings
    n_iterations <- 2000
    n_thin <- 1
    n_burnin <- 500
    n_chains <- 4
    #n_cores <- parallel::detectCores()
    n_cores <- 4
    delta = 0.95
    
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
        sigma_psi_level_four = runif(1, 0, 1),
        psi_site_area = runif(1, -1, 1),
        
        mu_p_cs_0 = runif(1, -1, 0),
        sigma_p_cs_site = runif(1, 0, 0.5),
        sigma_p_cs_level_three = runif(1, 0, 0.5),
        sigma_p_cs_ecoregion_one = runif(1, 0, 0.5),
        p_cs_interval = runif(1, 0, 1),
        p_cs_pop_density = runif(1, -1, 1),
        
        # start musuem values close to zero
        mu_p_rc_0 = runif(1, -0.5, 0.5),
        sigma_p_rc_site = runif(1, 0, 0.25),
        #sigma_p_rc_level_three = runif(1, 0, 0.25),
        #sigma_p_rc_level_four = runif(1, 0, 0.25),
        p_rc_total_records = runif(1, -0.5, 0.5)    
      )
    )
    
  }
  
} else { # taxon == syrphidae, urban
  
  if(urban_sites == TRUE){
    
    stan_data <- c("V_cs", "V_rc", "V_rc_NA",
                   "ranges", 
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_level_three", 
                   "level_three_lookup", 
                   "n_level_four",
                   "level_four_lookup",
                   "pop_densities", "site_areas", 
                   "nativity",
                   "avg_income", "avg_racial_minority",
                   "open_developed", "natural_habitat"
                   ) 
    
    # Parameters monitored
    params <- c("sigma_species_detection",
                "species_intercepts_detection",
                "rho",
                
                "mu_psi_0",
                "sigma_psi_species",
                "sigma_psi_site",
                "sigma_psi_level_three",
                "sigma_psi_level_four",
                "delta0",
                "delta1",
                "gamma0",
                "gamma1",
                "psi_site_area",
                "mu_psi_income",
                "mu_psi_race",
                "mu_psi_open_developed",
                
                "mu_p_cs_0",
                "sigma_p_cs_site",
                "sigma_p_cs_level_three",
                "p_cs_interval",
                "p_cs_pop_density",
                "p_cs_income",
                
                "mu_p_rc_0",
                "sigma_p_rc_site",
                "sigma_p_rc_level_three",
                
                "psi_species",
                "psi_natural_habitat",
                
                "psi_site",
                
                "W_species_rep_cs",
                "W_species_rep_rc",
                
                "mu_psi_natural_habitat_native",
                "mu_psi_natural_habitat_nonnative",
                "mu_psi_natural_habitat_all_species"
    )
    
    # MCMC settings
    n_iterations <- 300
    n_thin <- 1
    n_burnin <- 150
    n_chains <- 4
    n_cores <- 4
    #n_cores <- parallel::detectCores()
    delta = 0.95
    max_treedepth = 12
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    set.seed(1)
    inits <- lapply(1:n_chains, function(i)
      
      list(
            mu_psi_0 = runif(1, 0, 0.5),
            sigma_psi_species = runif(1, 1.5, 2.5),
            sigma_psi_site = runif(1, 1.5, 2.5),
            sigma_psi_level_three = runif(1, 0.75, 1.25),
            sigma_psi_level_four = runif(1, 0.5, 1),
            delta0 = runif(1, -0.5, 0.5),
            delta1 = runif(1, 0, 0.5),
            gamma0 = runif(1, 0.75, 1), # must be a positive value!
            gamma1 = runif(1, 0, 0.1), # gamma0+gamma1 inits must be >0!
            psi_site_area = runif(1, -0.5, 0.5),
            mu_psi_income = runif(1, -0.25, 0.25),
            mu_psi_open_developed = runif(1, -0.25, 0.25),
            
            mu_p_cs_0 = runif(1, -3, -2.5),
            sigma_p_cs_species = runif(1, 0, 1),
            sigma_p_cs_site = runif(1, 1, 1.25),
            sigma_p_cs_level_three = runif(1, 0.75, 1),
            #sigma_p_cs_level_four = runif(1, 0.75, 1),
            p_cs_interval = runif(1, 0.5, 0.6),
            p_cs_pop_density = runif(1, 0.4, 0.6),
            
            mu_p_rc_0 = runif(1, -0.5, 0.5),
            sigma_p_rc_site = runif(1, 0, 0.5),
            sigma_p_rc_level_three = runif(1, 0, 0.5)
           
      )
    )
    
  } else { 
    
    stan_data <- c("V_cs", "V_rc", "V_rc_NA",
                   "ranges", 
                   "n_species", "n_sites", "n_intervals", "n_visits", 
                   "intervals", "species", "sites",
                   "n_genera", "genus_lookup",
                   "n_level_three", 
                   "level_three_lookup", 
                   "n_level_four",
                   "level_four_lookup",
                   "pop_densities", "site_areas"
                   ) 
    
    # Parameters monitored
    params <- c(
      "sigma_species_detection",
      "species_intercepts_detection",
      "rho",
      
      "mu_psi_0",
      "sigma_psi_species",
      "sigma_psi_site",
      "sigma_psi_level_three",
      "sigma_psi_level_four",
      "psi_site_area",
      
      "mu_p_cs_0",
      "sigma_p_cs_site",
      "sigma_p_cs_level_three",
      "p_cs_interval",
      "p_cs_pop_density",
      "p_cs_income",
      
      "mu_p_rc_0",
      "sigma_p_rc_site",
      "sigma_p_rc_level_three",
      
      "psi_species",
      "psi_natural_habitat",
      
      "psi_site",
      "psi_level_four",
      "psi_level_three",
      
      "W_species_rep_cs",
      "W_species_rep_rc"
      
    )
    
    
    # MCMC settings
    n_iterations <- 300
    n_thin <- 1
    n_burnin <- 150
    n_chains <- 4
    n_cores <- 4
    #n_cores <- parallel::detectCores()
    delta = 0.95
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      
      list(
        
        mu_psi_0 = runif(1, -1, 1),
        sigma_psi_species = runif(1, 0, 1),
        sigma_psi_site = runif(1, 1, 2),
        sigma_psi_level_three = runif(1, 0, 1),
        sigma_psi_level_four = runif(1, 0, 1),
        mu_psi_income = runif(1, -1, 1),
        mu_psi_natural_habitat = runif(1, -1, 1),
        sigma_psi_natural_habitat = runif(1, 0, 1),
        psi_site_area = runif(1, -1, 1),
        mu_psi_open_developed = runif(1, -1, 1),
        
        mu_p_cs_0 = runif(1, -1, 0),
        sigma_p_cs_site = runif(1, 0.5, 1),
        sigma_p_cs_level_three = runif(1, 0, 0.5),
        p_cs_interval = runif(1, 0, 1),
        p_cs_pop_density = runif(1, -1, 1),
        
        mu_p_rc_0 = runif(1, -0.5, 0.5),
        sigma_p_rc_site = runif(1, 0, 0.5),
        sigma_p_rc_level_three = runif(1, 0, 0.5)
      )
    )
  
  }
}


## --------------------------------------------------
### Run model

# load appropriate model file from the directory
if(urban_sites == TRUE){
  stan_model <- paste0("./occupancy/models/model_", taxon, ".stan")
} else {
  stan_model <- paste0("./occupancy/models/model_", taxon, "_simple.stan")
}

if(use_reparameterized_rand_effects_model == TRUE){
  level_four_lookup <- level_four_lookup_by_site
  stan_model <- paste0("./occupancy/models/model_", taxon, "_reparameterized_rand_effects.stan")
}

# or manually enter a model name
#stan_model <- paste0("./occupancy/models/model_", taxon, "_no_open_developed.stan")

## Call Stan from R
set.seed(1)
stan_out <- stan(stan_model,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 control=list(adapt_delta=delta, max_treedepth=max_treedepth),
                 seed = 3,
                 open_progress = FALSE,
                 cores = n_cores)

saveRDS(stan_out, paste0(
  "./occupancy/model_outputs/", taxon, dir, 
  taxon, "_",
  grid_size / 1000,
  "km_", min_population_size, "minpop_", 
  min_unique_detections, "minUniqueDetections_",
  n_intervals, "ints_", n_visits, "visits_",
  ".rds"
)
)

stan_out <- readRDS(paste0(
  "./occupancy/model_outputs/", taxon, dir, taxon, "_", grid_size / 1000, 
  "km_", min_population_size, "minpop_", 
  min_unique_detections, "minUniqueDetections_",
  n_intervals, "ints_", n_visits, "visits_", 
  ".rds"
)
)

# read in a model ouput manually
#stan_out <- readRDS("./occupancy/model_outputs/large_files/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits.rds")
#stan_out <- readRDS("./occupancy/model_outputs/large_files/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits.rds")
stan_out <- readRDS("./occupancy/model_outputs/bombus/by_city/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits.RDS")

# print main effects
# print results
# for urban sites
if(taxon == "syrphidae"){
  print(stan_out, digits = 3, pars = c(
    "rho", 
    "sigma_species_detection[1]",
    "sigma_species_detection[2]", 
    
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    #"sigma_psi_level_four",
    "delta0",
    "delta1",
    "gamma0",
    "gamma1",
    "psi_site_area",
    "mu_psi_income",
    "mu_psi_race",
    "mu_psi_open_developed",
    #"sigma_psi_open_developed",
    "mu_psi_natural_habitat_native",
    "mu_psi_natural_habitat_nonnative",
    "mu_psi_natural_habitat_all_species",
    
    "mu_p_cs_0",
    #"sigma_p_cs_species",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    #"sigma_p_cs_level_four",
    "p_cs_interval",
    "p_cs_pop_density",
    "p_cs_income",
    
    "mu_p_rc_0",
    "sigma_p_rc_site",
    "sigma_p_rc_level_three"
  ))
} else {
  print(stan_out, digits = 3, pars = c(
    "rho", 
    "sigma_species_detection[1]",
    "sigma_species_detection[2]", 
    "mu_psi_0",
    "sigma_psi_site",
    #"sigma_psi_level_three",
    #"sigma_psi_level_four",
    #"sigma_psi_income",
    "mu_psi_natural_habitat",
    "sigma_psi_natural_habitat",
    "mu_psi_open_developed",
    "mu_psi_income",
    "mu_psi_race",
    #"sigma_psi_open_developed",
    "psi_site_area",
    
    "mu_p_cs_0",
    #"sigma_p_cs_species",
    "sigma_p_cs_site",
    #"sigma_p_cs_level_three",
    "p_cs_interval",
    "p_cs_pop_density", 
    "p_cs_income",
    
    "mu_p_rc_0",
    #"sigma_p_rc_species",
    "sigma_p_rc_site"
    #"sigma_p_rc_level_three"
  ))
}

# non urban sites
if(taxon == "syrphidae"){
  print(stan_out, digits = 3, pars = c(
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "psi_site_area",
    
    "mu_p_cs_0",
    "sigma_p_cs_species",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "sigma_p_cs_level_four",
    "p_cs_interval",
    "p_cs_pop_density"
  ))
} else {
  print(stan_out, digits = 3, pars = c(
    "rho", 
    "sigma_species_detection[1]",
    "sigma_species_detection[2]", 
    "mu_psi_0",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "psi_site_area",
    
    "mu_p_cs_0",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    #"sigma_p_cs_level_four",
    "p_cs_interval",
    "p_cs_pop_density",
    "p_cs_income",
    
    "mu_p_rc_0",
    "sigma_p_rc_site",
    "sigma_p_rc_level_three"
    #"sigma_p_rc_level_four"
    #"p_rc_natural_habitat"
  ))
}

print(stan_out, digits = 3, pars = c(
  "mu_psi_herb_shrub_forest", 
  "mu_psi_income"
))

# print some specific parameter if desired
print(stan_out, digits = 3, pars=
        c("psi_site"))

View(as.data.frame(rstan::summary(stan_out)))


## --------------------------------------------------
### Simple diagnostic plots

# traceplots
# urban sites
if(taxon == "syrphidae"){
  traceplot(
    stan_out, pars = c( # occupancy
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    #"sigma_psi_level_four",
    "delta0",
    "delta1",
    "gamma0",
    "gamma1",
    "psi_site_area",
    "mu_psi_income",
    "mu_psi_race",
    "mu_psi_open_developed"
  ))
  traceplot(stan_out, pars = c(
    "mu_p_cs_0",
    "p_cs_interval",
    "p_cs_pop_density",
    #"sigma_p_cs_species",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    #"sigma_p_cs_level_four",
    "mu_p_rc_0",
    "sigma_p_rc_site",
    "sigma_p_rc_level_three"
  ))
  traceplot(stan_out, pars=
              c("mu_psi_natural_habitat_native",
                "mu_psi_natural_habitat_nonnative",
                "mu_psi_natural_habitat_all_species"
  ))
} else{
  traceplot(stan_out, pars = c( # occupancy
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "mu_psi_natural_habitat",
    "sigma_psi_natural_habitat",
    "mu_psi_income",
    #"sigma_psi_income",
    "psi_site_area",
    "mu_psi_open_developed"
    #"sigma_psi_open_developed"
  ))
  traceplot(stan_out, pars = c( # detection
    "mu_p_cs_0",
    #"sigma_p_cs_species",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "p_cs_interval",
    "p_cs_pop_density",
    "mu_p_rc_0",
    #"sigma_p_rc_species",
    "sigma_p_rc_site",
    "sigma_p_rc_level_three"
  ))
  traceplot(stan_out, pars = c( # species detection
    "rho",
    "sigma_species_detection"
  ))
}

# non-urban
if(taxon == "syrphidae"){
  traceplot(stan_out, pars = c( # occupancy
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four",
    "psi_site_area"
  ))
  traceplot(stan_out, pars = c(
    "mu_p_cs_0",
    "p_cs_interval",
    "p_cs_pop_density",
    "sigma_p_cs_species",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "sigma_p_cs_level_four"
  ))
} else{
  traceplot(stan_out, pars = c( # occupancy
    "mu_psi_0",
    "sigma_psi_species",
    "sigma_psi_site",
    "sigma_psi_level_three",
    "sigma_psi_level_four"
  ))
  traceplot(stan_out, pars = c( # detection
    "mu_p_cs_0",
    "sigma_p_cs_site",
    "sigma_p_cs_level_three",
    "sigma_p_cs_level_four",
    "p_cs_interval",
    "p_cs_pop_density",
    "mu_p_rc_0",
    "sigma_p_rc_site"
    #"sigma_p_rc_level_three",
    #"sigma_p_rc_level_four",
    #"p_rc_total_records"
  ))
  traceplot(stan_out, pars = c( # species detection
    "rho",
    "sigma_species_detection"
  ))
}

# pairs plot
pairs(stan_out, pars = c(
  "mu_psi_0",
  #"sigma_psi_species",
  "sigma_psi_site",
  #"mu_psi_herb_shrub_forest",
  #"psi_site_area",
  
  "mu_p_cs_0",
  #"sigma_p_citsci_species",
  "sigma_p_cs_site"
  #"p_cs_interval",
  #"p_citsci_pop_density", 
  
  #"mu_p_rc_0"
  #"sigma_p_rc_site",
  #"sigma_p_rc_level_three"
  #"sigma_p_museum_species",
  #"sigma_p_museum_site"
))

sx=seq(0,3,1)
y=-4.5+0.4*x^2
plot(x,y, col='violet',type='o',lwd=2,lty=1)


# plot species detections from species counts
# draw plot
ggplot(data = species_counts, aes(x = reorder(species, total_count), y = total_count)) + 
  geom_bar(stat = "identity", show.legend = FALSE) + 
  labs(x = "", 
       y = "Number of detections in urban landscapes (2011-2022)") +
  coord_flip()

# get an "average" P value
fit_summary <- rstan::summary(stan_out)

fit_summary$summary[1,] # top row (first monitored paramter from the fit)
fit_summary$summary[1,1] # parameter mean
fit_summary$summary[1,3] # parameter sd
fit_summary$summary[1,4] # 2.5 % CI
fit_summary$summary[1,8] # 97.5 % CI

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
# hoverflies
(mean_FTP <- mean(fit_summary$summary[804:944,1]))
# bumble bees cs
mean_FTP <- mean(fit_summary$summary[760:791,1])
# bumble bees rc
mean_FTP <- mean(fit_summary$summary[856:887,1])
