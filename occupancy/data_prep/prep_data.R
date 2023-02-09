### Prepeare data to feed to model_integrated.stan
# jcu; started oct 27, 2022

# will need to assign occupancy intervals and visit numbers
# for occupancy (as opposed to abundance), will need to 
# filter down to one unique occurrence per species*site*interval*year

## The data we will need to prepare to feed to the model:
# V_citsci # an array of presence/absence citizen science detection data
# V_citsci_NA # indicator of whether site is in species range 
# V_museum # an array of presence/absence museum records detection data
# V_museum_NA # indicator of whether a) site is in species range and b) ..
# if a museum community sampling event occurred at the site in the visit time
# n_species # number of species
# n_sites # number of sites
# n_intervals # number of occupancy intervals 
# n_visits # number of visits in each interval

library(tidyverse)

prep_data <- function(era_start, era_end, n_intervals, n_visits, 
                      min_records_per_species,
                      grid_size, min_population_size, 
                      min_records_for_community_sampling_event,
                      min_year_for_species_ranges,
                      taxon,
                      min_site_area,
                      remove_unidentified_species,
                      consider_species_occurring_outside_sites,
                      min_records_per_species_full
) {
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  source("./occupancy/data_prep/get_spatial_data.R")
  
  # retrieve the spatial occurrence record data
  my_spatial_data <- get_spatial_data(
    grid_size, min_population_size, taxon, min_site_area)
  
  # save the data in case you want to make tweaks to the prep data
  # without redoing the raster extractions
  # saveRDS(my_spatial_data, "./occupancy/analysis/prepped_data/spatial_data_list.rds")
  # my_spatial_data <- readRDS("./occupancy/analysis/prepped_data/spatial_data_list.rds")
  gc()
  
  df_id_urban_filtered <- as.data.frame(my_spatial_data$df_id_urban_filtered)
  
  # remove specimens with no species level identification
  if(remove_unidentified_species == TRUE) {
    
    df_id_urban_filtered <- df_id_urban_filtered %>%
      filter(species != "")
    
  }
  
  # spatial covariate data to pass to run model
  urban_grid <- my_spatial_data$urban_grid
  
  raw_pop_density <- urban_grid$pop_density_per_km2
  pop_density_vector <- urban_grid$scaled_pop_den_km2
  site_area_vector <- urban_grid$scaled_site_area
  developed_open_vector <- urban_grid$scaled_developed_open
  herb_shrub_vector <- urban_grid$scaled_herb_shrub_cover
  forest_vector <- urban_grid$scaled_forest
  herb_shrub_forest_vector <- urban_grid$scaled_herb_shrub_forest
  developed_med_high_vector <- urban_grid$scaled_developed_med_high
  ecoregion_three_vector <- urban_grid$ecoregion_three_vector
  ecoregion_one_vector <- urban_grid$ecoregion_one_vector

  site_name_vector <- urban_grid$grid_id
  
  # correlation matrix
  correlation_matrix <- my_spatial_data$correlation_matrix
  
  my_spatial_data
  
  ecoregion_three_lookup <- my_spatial_data$ecoregion_three_lookup
  ecoregion_one_lookup <- my_spatial_data$ecoregion_one_lookup
  n_ecoregion_three <- my_spatial_data$n_ecoregion_three
  n_ecoregion_one <- my_spatial_data$n_ecoregion_one
  
  # other info to pass to output that we may want to keep track of
  # correlation between variables
  
  ## --------------------------------------------------
  # occurrence data ONLY FROM SITES
  
  # total records from time period
  total_records_since_2000 <- nrow(df_id_urban_filtered)
  
  total_records_since_study <- df_id_urban_filtered %>%
    filter(year > era_start) %>%
    nrow()
  
  # citsci and museum records from time period
  citsci_records <- df_id_urban_filtered %>%
    filter(basisOfRecord == "HUMAN_OBSERVATION") %>%
    filter(year >= era_start) %>%
    nrow()
  
  # citsci and museum records from time period
  citsci_detections <- df_id_urban_filtered %>%
    filter(basisOfRecord == "HUMAN_OBSERVATION") %>%
    filter(year >= era_start) %>%
    group_by(species, grid_id, year) %>%
    slice(1) %>%
    nrow()
  
  museum_records <- df_id_urban_filtered %>%
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
    filter(year >= era_start) %>%
    nrow()
  
  museum_detections <- df_id_urban_filtered %>%
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
    filter(year >= era_start) %>%
    group_by(species, grid_id, year) %>%
    slice(1) %>%
    nrow()
  
  # species counts
  species_counts <- df_id_urban_filtered %>%
    filter(year > era_start) %>%
    group_by(species) %>%
    count(name="total_count")
  
  # species counts
  species_counts_citsci <- df_id_urban_filtered %>%
    filter(basisOfRecord == "HUMAN_OBSERVATION") %>%
    filter(year > era_start) %>%
    group_by(species) %>%
    count(name="citsci_count")
  
  # species counts
  species_counts_museum <- df_id_urban_filtered %>%
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
    filter(year > era_start) %>%
    group_by(species) %>%
    count(name="museum_count")
  
  # species counts
  species_counts <- species_counts %>%
    left_join(., species_counts_citsci) %>%
    left_join(., species_counts_museum)
  
  rm(species_counts_citsci, species_counts_museum)
  
  # species detections
  species_detections <- df_id_urban_filtered %>%
    filter(year > era_start) %>%
    group_by(species, grid_id, year) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(species) %>%
    count(name="total_detections")
  
  # species detections
  species_detections_citsci <- df_id_urban_filtered %>%
    filter(basisOfRecord == "HUMAN_OBSERVATION") %>%
    filter(year > era_start) %>%
    group_by(species, grid_id, year) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(species) %>%
    count(name="citsci_detections")
  
  # species detections
  species_detections_museum <- df_id_urban_filtered %>%
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
    filter(year > era_start) %>%
    group_by(species, grid_id, year) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(species) %>%
    count(name="museum_detections")
  
  # species counts
  species_detections <- species_detections %>%
    left_join(., species_detections_citsci) %>%
    left_join(., species_detections_museum)
  
  rm(species_detections_citsci, species_detections_museum)
  
  ## --------------------------------------------------
  # occurrence data FROM ANYWHERE IN CONTINENTAL US
  
  # read either the syrphidae data or the bombus data
  df_full <- read.csv(paste0("./data/occurrence_data/", taxon, "_data_all.csv")) %>%
    filter(species != "") %>% 
    
    # and perform any further initial data filters
    
    # filter out B. impatiens outside of it's recently expanding native range (Looney et al.)
    # (filter out occurrences west of 105 Longitude)
    filter(!(species == "Bombus impatiens" & decimalLongitude < -105)) %>%
    filter(!(species == "Bombus pensylvanicus" & decimalLongitude < -110)) %>%
    
    # filter out records with high location uncertainty (threshold at 10km)
    # assuming na uncertainty (large portion of records) is under threshold
    mutate(coordinateUncertaintyInMeters = replace_na(coordinateUncertaintyInMeters, 0)) %>%
    filter(coordinateUncertaintyInMeters < 10000)
  
  # total records from time period
  total_records_since_2000_full <- nrow(df_full)
  
  total_records_since_study_full <- df_full %>%
    filter(year > era_start) %>%
    nrow()
  
  # citsci and museum records from time period
  citsci_records_full <- df_full %>%
    filter(basisOfRecord == "HUMAN_OBSERVATION") %>%
    nrow()
  
  # citsci and museum records from time period
  #citsci_detections_full <- df_full %>%
    #filter(basisOfRecord == "HUMAN_OBSERVATION") %>%
    #filter(year >= era_start) %>%
    # group_by(species, grid_id, year) %>%
    #slice(1) %>%
    #nrow()
  
  museum_records_full <- df_full %>%
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
    nrow()
  
  #museum_detections <- df_full %>%
    #filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
    #filter(year >= era_start) %>%
    #group_by(species, grid_id, year) %>%
    #slice(1) %>%
    #nrow()
  
  # species counts
  species_counts_full <- df_full %>%
    group_by(species) %>%
    count(name="total_count")
  
  # species counts
  species_counts_citsci_full <- df_full %>%
    filter(basisOfRecord == "HUMAN_OBSERVATION") %>%
    group_by(species) %>%
    count(name="citsci_count")
  
  # species counts
  species_counts_museum_full <- df_full %>%
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>%
    group_by(species) %>%
    count(name="museum_count")
  
  # species counts
  species_counts_full <- species_counts_full %>%
    left_join(., species_counts_citsci_full) %>%
    left_join(., species_counts_museum_full)
  
  ## --------------------------------------------------
  # assign species list based on whether we want occurrences from anywhere in the extent or sites only 
  
  # anywehere in the extent
  if(consider_species_occurring_outside_sites == TRUE){
    
    ## Get unique species 
    # create an alphabetized list of all species encountered across all sites*intervals*visits
    species_list <- df_full %>%
      
      # remove species with total observations (n) < min_records_per_species 
      group_by(species) %>%
      add_tally() %>%
      filter(n >= min_records_per_species_full) %>%
      ungroup() %>%
      
      group_by(species) %>%
      slice(1) %>% # take one row per species (the name of each species)
      ungroup() %>%
      dplyr::select(species) # extract species names column as vector
    
    # get vectors of species, sites, intervals, and visits 
    # these are all species that were observed at least min_records_per_species
    species_vector <- species_list %>%
      pull(species)
    
  } else { # else only from sites
    
    ## Get unique species 
    # create an alphabetized list of all species encountered across all sites*intervals*visits
    species_list <- df_filtered %>%
      group_by(species) %>%
      slice(1) %>% # take one row per species (the name of each species)
      ungroup() %>%
      dplyr::select(species) # extract species names column as vector
    
    # get vectors of species, sites, intervals, and visits 
    # these are all species that were observed at least min_records_per_species
    species_vector <- species_list %>%
      pull(species)
    
  }
  
  ## --------------------------------------------------
  # get genus for phylogenetic clustering
  if(taxon == "syrphidae"){
    
    genus_lookup <- species_list %>%
      mutate(genus = word(species, 1)) %>%
      pull(genus)
    
    genus_vector <- species_list %>%
      mutate(genus = word(species, 1)) %>%
      group_by(genus) %>%
      slice(1) %>%
      ungroup() %>%
      pull(genus) 
    
    n_genera <- length(genus_vector)
    
  } else {
    
    # no generic clustering for bumblebees
    genus_lookup <- NULL
    genus_vector <- NULL
    n_genera <- NULL
    
  }
  
  ## --------------------------------------------------
  # auto generate plot data (hashtag out if not using)
  
  # count for iNat vs Museum
  out <- df_full %>%
    group_by(basisOfRecord) %>% 
    count(year)
  
  out2 <- df_id_urban_filtered %>%
    group_by(basisOfRecord) %>% 
    count(year)
  
  xints <- vector()
  for(i in 1:(n_intervals+1)){
    xints[i] <- (era_start - 0.5) + n_visits*(i - 1)
  }
  
  xints2 <- vector()
  for(i in 1:(n_intervals*n_visits + 1)){
    xints2[i] <- (era_start - 0.5) + (i - 1)
  }
  
  # Chronological record counts split by basis of record (any)
  ggplot(out, aes(x = year, y = n, col = as.factor(basisOfRecord))) + 
    xlim(2000, 2023) + # choose some years for the axes
    geom_line() +
    scale_colour_manual(name = "Basis of records", 
                        labels = c("Citizen science records", "Museum records"),
                        values=c("red","blue")) +
    geom_vline(xintercept=xints, 
               linewidth=3, alpha=0.5) +
    geom_vline(xintercept=xints2, 
      linewidth=1, alpha=0.7, linetype = 'dotted') +
    theme_bw() +
    xlab("Year") +
    ylab("Number of Records \n(continental U.S.)") +
    theme(legend.position = c(0.15, 0.8),
          legend.title = element_text(colour="black", size=14, 
                                      face="bold"),
          legend.text = element_text(colour="black", size=12),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16)
    )
  
  ggsave(paste0("./figures/occurrence_data/", taxon, "_temporal_full.pdf"),
         width = 11, height = 8, units = "in") 
  
  # Chronological record counts split by basis of record (sites only)
  ggplot(out2, aes(x = year, y = n, col = as.factor(basisOfRecord))) + 
    xlim(2000, 2023) + # choose some years for the axes
    geom_line() +
    scale_colour_manual(name = "Basis of Records", 
                        labels = c("Citizen science records", "Museum records"),
                        values=c("red","blue")) +
    geom_vline(xintercept=xints, 
               linewidth=3, alpha=0.5) +
    geom_vline(xintercept=xints2, 
               linewidth=1, alpha=0.7, linetype = 'dotted') +
    theme_bw() +
    xlab("Year") +
    ylab("Number of Records \n(urban sites)") +
    theme(legend.position = c(0.15, 0.8),
          legend.title = element_text(colour="black", size=14, 
                                      face="bold"),
          legend.text = element_text(colour="black", size=12),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16)
    )
  
  ggsave(paste0("./figures/occurrence_data/", taxon, "_temporal_urban_sites.pdf"),
         width = 11, height = 8, units = "in") 
  
  
  ## --------------------------------------------------
  # assign study dimensions
  
  total_years <- era_end - era_start + 1
  remainder <- total_years %% n_intervals
  n_visits <- n_visits
  min_records_per_species <- min_records_per_species
  
  ## --------------------------------------------------
  # assign occupancy visits and intervals and reduce to one
  # sample per species per site per visit
  
  df_filtered <- df_id_urban_filtered %>%
    
    # remove records (if any) missing species level identification
    filter(species != "") %>%
    
    # assign year as - year after era_start
    mutate(occ_year = (year - era_start)) %>% 
    
    # remove data from years that are in the remainder
    # occupancy intervals have to be equal in length for the model to process
    # e.g. we can't have intervals of 3 years and then one interval with only 1 year
    # filter(occ_year %in% (remainder:1)) %>%
    
    # remove years before the start date
    filter(occ_year >= 0) %>%
    # remove years after end date
    filter(occ_year < total_years) %>%
    
    # now assign the years into 1:n_intervals
    mutate(occ_interval = occ_year %/% n_visits) %>%
    
    # add a sampling round (1:n)
    mutate(visit = (occ_year %% n_visits)) %>%
    
    # remove species with total observations (n) < min_records_per_species 
    filter(species %in% species_vector) %>%
    
    # # one unique row per site*species*occ_interval*visit combination
    # group_by(grid_id, species, occ_interval, visit) %>% 
    # slice(1) %>%
    # ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, grid_id, occ_interval, occ_year, visit) %>%
    arrange((occ_year))
  # end pipe
  
  ## --------------------------------------------------
  # filter to citizen science records
  # assign occupancy visits and intervals and reduce to one
  # sample per species per site per visit
  
  df_citsci <- df_id_urban_filtered %>%
    
    # remove records (if any) missing species level identification
    filter(species != "") %>%
    
    # assign year as - year after era_start
    mutate(occ_year = (year - era_start)) %>%
    
    # remove years before the start date
    filter(occ_year >= 0) %>%
    # remove years after end date
    filter(occ_year < total_years) %>%
    
    # now assign the years into 1:n_intervals
    mutate(occ_interval = occ_year %/% n_visits) %>%
    
    # add a sampling round (1:n)
    mutate(visit = (occ_year %% n_visits)) %>%
    
    # remove species with total observations (n) < min_records_per_species
    # remove species with total observations (n) < min_records_per_species 
    filter(species %in% species_vector) %>%
    
    # filter to citizen science data only
    filter(basisOfRecord == "HUMAN_OBSERVATION") %>% 
    
    # one unique row per site*species*occ_interval*visit combination
    # contrast this line with the abundance models!
    group_by(grid_id, species, occ_interval, visit) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, grid_id, occ_interval, occ_year, visit) %>%
    arrange((occ_year))
  # end pipe
  
  ## --------------------------------------------------
  # filter to museum records
  # assign occupancy visits and intervals and reduce to one
  # sample per species per site per visit
  
  df_museum <- df_id_urban_filtered %>%
    
    # remove records (if any) missing species level identification
    filter(species != "") %>%
    
    # assign year as - year after era_start
    mutate(occ_year = (year - era_start)) %>%
    
    # remove years before the start date
    filter(occ_year >= 0) %>%
    # remove years after end date
    filter(occ_year < total_years) %>%
    
    # now assign the years into 1:n_intervals
    mutate(occ_interval = occ_year %/% n_visits) %>%
    
    # add a sampling round (1:n)
    mutate(visit = (occ_year %% n_visits)) %>%
    
    # remove species with total observations (n) < min_records_per_species
    filter(species %in% species_vector) %>%
    
    # filter to citizen science data only
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>% 
    
    # determine whether a community sampling event occurred
    # using collector name might be overly conservative because for example
    # the data includes recordedBy == J. Fulmer *and* recordedBy J. W. Fulmer
    # instead grouping by collections housed in the same institution from the same year
    # within a site
    group_by(institutionCode, year, grid_id) %>%
    mutate(n_species_sampled = n_distinct(species)) %>%
    # filter(n_species_sampled >= min_species_for_community_sampling_event) %>%
    
    # one unique row per site*species*occ_interval*visit combination
    group_by(grid_id, species, occ_interval, visit) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, grid_id, occ_interval, occ_year, visit) %>%
    arrange((occ_year))
  # end pipe
  
  ## --------------------------------------------------
  # museum community samples
  
  community_samples <- df_id_urban_filtered %>%
    
    # remove records (if any) missing species level identification
    filter(species != "") %>%
    
    # assign year as - year after era_start
    mutate(occ_year = (year - era_start)) %>%
    
    # remove years before the start date
    filter(occ_year >= 0) %>%
    # remove years after end date
    filter(occ_year < total_years) %>%
    
    # now assign the years into 1:n_intervals
    mutate(occ_interval = occ_year %/% n_visits) %>%
    
    # add a sampling round (1:n)
    mutate(visit = (occ_year %% n_visits)) %>%
    
    # remove species with total observations (n) < min_records_per_species
    filter(species %in% species_vector) %>%
    
    # filter to citizen science data only
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>% 
    
    # determine whether a community sampling event occurred
    # using collector name might be overly conservative because for example
    # the data includes recordedBy == J. Fulmer *and* recordedBy J. W. Fulmer
    # instead grouping by collections housed in the same institution from the same year
    # within a site
    group_by(institutionCode, year, grid_id) %>%
    mutate(n_species_sampled = n_distinct(species)) %>%
    ungroup() %>%
    dplyr::select(grid_id, occ_interval, occ_year, visit, species, n_species_sampled) %>%
    group_by(grid_id, occ_interval, occ_year, visit) %>%
    # slice max in case there are multiple institutions collecting from a site in a year
    slice_max(n_species_sampled) %>%
    # then take one per year
    slice(1) %>%
    mutate(community_sampled = ifelse(
      n_species_sampled >= min_species_for_community_sampling_event,
      1, 0)) %>%
    dplyr::select(-species, -n_species_sampled) %>%
    mutate(occ_interval = as.character(occ_interval),
           visit = as.character(visit))
  # end pipe
  
  ## --------------------------------------------------
  # records per sampling event
  
  museum_total_records <- df_id_urban_filtered %>%
    
    # remove records (if any) missing species level identification
    filter(species != "") %>%
    
    # assign year as - year after era_start
    mutate(occ_year = (year - era_start)) %>%
    
    # remove years before the start date
    filter(occ_year >= 0) %>%
    # remove years after end date
    filter(occ_year < total_years) %>%
    
    # now assign the years into 1:n_intervals
    mutate(occ_interval = occ_year %/% n_visits) %>%
    
    # add a sampling round (1:n)
    mutate(visit = (occ_year %% n_visits)) %>%
    
    # remove species with total observations (n) < min_records_per_species
    filter(species %in% species_vector) %>%
    
    # filter to citizen science data only
    filter(basisOfRecord == "PRESERVED_SPECIMEN") %>% 
    
    # determine whether a community sampling event occurred
    # using collector name might be overly conservative because for example
    # the data includes recordedBy == J. Fulmer *and* recordedBy J. W. Fulmer
    # instead grouping by collections housed in the same institution from the same year
    # within a site
    group_by(year, grid_id) %>%
    add_tally(name = "records_per_year_per_site") %>%
    slice(1) %>%
    ungroup() %>%
    # now add an average per sampling interval 
    group_by(occ_interval, grid_id) %>%
    mutate(museum_total_records = mean(records_per_year_per_site)) %>%
    slice(1) %>%
    ungroup() %>%
    
    dplyr::select(grid_id, occ_interval, museum_total_records) %>%
    mutate(occ_interval = as.numeric(occ_interval),
           grid_id = as.integer(grid_id),
           museum_total_records = as.numeric(museum_total_records)) %>%
    # scale the variable (before comparing to sites with no museum data - 
    # since they will be passed in the likelihood function)
    mutate(museum_total_records = center_scale(museum_total_records))
  # end pipe
  
  # join with all siteXintervals after generating that table
  
  ## --------------------------------------------------
  # Extract stan data from df
  # I use the list of all species and all sites 
  # rather than species and sites sampled potentially by only citizen science or by only museums
  # to generate these master lists of all possible sites and species
  
  ## Get unique sites 
  # create an alphabetized list of all sites
  site_vector <- site_name_vector
  
  site_list <- as.data.frame(site_vector) %>%
    rename("grid_id" = "site_vector") %>%
    mutate(grid_id = as.character(grid_id))
  
  interval_vector <- as.vector(levels(as.factor(df_filtered$occ_interval)))
  
  visit_vector <- as.vector(levels(as.factor(df_filtered$visit)))
  
  # find study dimensions from the pooled data
  n_species <- length(species_vector)
  
  n_sites <- length(site_vector)
  
  n_intervals <- length(interval_vector)
  
  n_visits <- length(visit_vector)
  
  ## --------------------------------------------------
  ## Now we are ready to create the detection matrices, V_citsci and V_museum
  
  ## --------------------------------------------------
  ## V_citsci 
  V_citsci <- array(data = NA, dim = c(n_species, n_sites, n_intervals, n_visits))
  
  for(i in 1:n_intervals){
    for(j in 1:n_visits){
      
      # iterate this across visits within intervals
      temp <- df_citsci %>%
        # filter to indices for interval and visit
        filter(occ_interval == (i - 1), # intervals start at 0 so need to subtract 1 from i
               visit == (j - 1)) %>%  # visits start at 0 so need to subtract 1 from j
        # now join with all species (so that we include species not captured during 
        # this interval*visit but which might actually be at some sites undetected)
        full_join(species_list, by="species") %>%
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        mutate(grid_id = as.character(grid_id)) %>%
        full_join(site_list, by="grid_id") %>%
        # have to convert back to integer to get correct ordering from low to high
        # MUST MATCH SITE NAMES VECTOR
        mutate(grid_id = as.integer(grid_id)) %>%
        # group by SPECIES
        group_by(species) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(grid_id, row, fill = 0) %>%
        # replace number of unique site captures of the species (if > 1) with 1.
        mutate_at(5:(ncol(.)), ~replace(., . > 1, 1)) %>%
        # if more columns are added these indices above^ might need to change
        # 5:(n_species+5) represent the columns of each site in the matrix
        # just need the matrix of 0's and 1's
        dplyr::select(5:(ncol(.))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(species  %in% levels(as.factor(species_list$species)))
      
      # convert from dataframe to matrix
      temp_matrix <- as.matrix(temp)
      # remove species names
      temp_matrix <- temp_matrix[,-1]
      
      # replace NAs for the interval i and visit j with the matrix
      V_citsci[1:n_species, 1:n_sites,i,j] <- temp_matrix[1:n_species, 1:n_sites]
      
    }
  }
  
  class(V_citsci) <- "numeric"
  
  ## --------------------------------------------------
  ## V_museum 
  V_museum <- array(data = NA, dim = c(n_species, n_sites, n_intervals, n_visits))
  
  for(i in 1:n_intervals){
    for(j in 1:n_visits){
      
      # iterate this across visits within intervals
      temp <- df_museum %>%
        # filter to indices for interval and visit
        filter(occ_interval == (i - 1), # intervals start at 0 so need to subtract 1 from i
               visit == (j - 1)) %>%  # visits start at 0 so need to subtract 1 from j
        # now join with all species (so that we include species not captured during 
        # this interval*visit but which might actually be at some sites undetected)
        full_join(species_list, by="species") %>%
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        mutate(grid_id = as.character(grid_id)) %>%
        full_join(site_list, by="grid_id") %>%
        # have to convert back to integer to get correct ordering from low to high
        # MUST MATCH SITE NAMES VECTOR
        mutate(grid_id = as.integer(grid_id)) %>%
        # group by SPECIES
        group_by(species) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(grid_id, row, fill = 0) %>%
        # replace number of unique site captures of the species (if > 1) with 1.
        mutate_at(5:(ncol(.)), ~replace(., . > 1, 1)) %>%
        # if more columns are added these indices above^ might need to change
        # 5:(n_species+5) represent the columns of each site in the matrix
        # just need the matrix of 0's and 1's
        dplyr::select(5:(ncol(.))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(species  %in% levels(as.factor(species_list$species)))
      
      
      # convert from dataframe to matrix
      temp_matrix <- as.matrix(temp)
      # remove species names
      temp_matrix <- temp_matrix[,-1]
      
      # replace NAs for the interval i and visit j with the matrix
      V_museum[1:n_species, 1:n_sites,i,j] <- temp_matrix[1:n_species, 1:n_sites]
      
    }
  }
  
  class(V_museum) <- "numeric"
  
  ## --------------------------------------------------
  # Generate an indicator array for whether or not a community-wide museum 
  # sampling event occurred at the site in a year
  
  # first make a df of all possible site visits * species visits
  all_species_site_visits <- as.data.frame(cbind(
    rep(species_vector, times=(n_sites*n_intervals*n_visits)),
    as.integer(rep(site_vector, each=n_species, times=n_intervals*n_visits)),
    rep(interval_vector, each=n_species*n_sites, times=(n_visits)),
    rep(visit_vector, each=n_species*n_sites*n_intervals)
    #rep(0:(total_years-1),)
  )) %>%
    rename("species" = "V1",
           "grid_id" = "V2",
           "occ_interval" = "V3",
           "visit" = "V4") %>%
    #"occ_year" = "V5") %>%
    mutate(grid_id = as.integer(grid_id))
  
  ## --------------------------------------------------
  # Get species ranges
  source("./occupancy/data_prep/get_species_ranges.R")
  
  species_ranges <- get_species_ranges(urban_grid,
                                       site_name_vector,
                                       n_sites,
                                       species_vector,
                                       n_species,
                                       min_year_for_species_ranges,
                                       taxon)
  
  gc()
  
  ## --------------------------------------------------
  # Now infer an NA structure for sites outside of range or 
  # sites* visits (museums only) where no community sampling event occurred
  
  ## --------------------------------------------------
  # Infer citizen science missing data
  # Particularly, whether the species can be sampled or not at a site given its range
  # this is the only condition considered to effect whether a species could potentially
  # be detected by a citizen science survey effort (compare with museum records)
  
  df_citsci_visits <- as.data.frame(cbind(rep(species_vector, each=n_sites),
                                          rep(site_vector, times = n_species),
                                          rep(1, times = n_species*n_sites),
                                          unlist(species_ranges))) %>%
    rename("species" = "V1",
           "grid_id" = "V2",
           "community_sample" = "V3",
           "in_range" = "V4") %>%
    mutate(sampled = as.numeric(community_sample)*
             as.numeric(in_range),
           grid_id = as.integer(grid_id)) %>%
    left_join(all_species_site_visits, .) 
  
  V_citsci_NA <- array(data = df_citsci_visits$sampled, dim = c(n_species, n_sites, n_intervals, n_visits))
  # check <- which(V_citsci>V_citsci_NA) # should NEVER have occurrences where the species can't be sampled,
  # thus check should be empty
  
  ## --------------------------------------------------
  # Infer museum record missing data
  # Particularly, whether the species can be sampled or not at a site given its range
  # AND 
  # whether or not a community sampling event has occurred
  # be detected by a citizen science survey effort (compare with museum records)
  
  # want to keep detections for species that were detected in isolation
  # but not infer that other species could be detected if we just have a few records
  
  all_visits_museum_visits_joined <- left_join(all_species_site_visits, community_samples, 
                                               by=c("grid_id", "occ_interval", "visit")) %>%
    # create an indicator if the site visit was a sample or not
    mutate(community_sampled = replace_na(community_sampled, 0),
           occ_interval = as.integer(occ_interval),
           visit = as.integer(visit)) %>%
    dplyr::select(-occ_year)
  
  ranges <- as.data.frame(cbind(rep(species_vector, each=n_sites),
                                rep(site_vector, times = n_species),
                                unlist(species_ranges))) %>%
    rename("species" = "V1",
           "grid_id" = "V2",
           "in_range" = "V3") %>%
    mutate(grid_id = as.integer(grid_id))
  
  all_visits_museum_visits_joined <- left_join(all_visits_museum_visits_joined, ranges) %>%
    mutate(sampled = as.numeric(community_sampled)*
             as.numeric(in_range))
  
  # for < min_species_for_community_sampling_event with a sampled indicator
  # i.e., keep the data for the singletons and doubletons, while not inferring abscense for the rest of the community
  df_museum <- df_museum %>%
    # if the species was sampled than at least that species was sampled (may also be a comm sample)
    mutate(non_comm_sample = 1) 
  
  # sampled for each row if either a community sample or a non comm sample
  all_visits_museum_visits_joined <- left_join(all_visits_museum_visits_joined, df_museum) %>%
    mutate(non_comm_sample = replace_na(non_comm_sample, 0)) %>%
    mutate(any_sampled = ifelse(non_comm_sample == 1, 1, sampled))
  
  # now spread into 4 dimensions
  V_museum_NA <- array(data = all_visits_museum_visits_joined$any_sampled, dim = c(n_species, n_sites, n_intervals, n_visits))
  check <- which(V_museum>V_museum_NA) # this will give you numerical value
  # thus check should be empty
  
  ## --------------------------------------------------
  # Construct records per visit covariate
  
  all_interval_site_visits <- all_species_site_visits %>%
    group_by(occ_interval, grid_id) %>%
    slice(1) %>%
    dplyr::select(occ_interval, grid_id) %>%
    mutate(occ_interval = as.numeric(occ_interval)) %>%
    ungroup()
  
  all_interval_sites_museum_joined <- left_join(all_interval_site_visits, museum_total_records, 
                                               by=c("grid_id", "occ_interval")) %>%
    # create an indicator if the site visit was a sampled or not
    mutate(museum_total_records = replace_na(museum_total_records, 0),
           occ_interval = as.integer(occ_interval)) %>%
    arrange(grid_id)
  
  # now spread into 2 dimensions
  museum_total_records <- matrix(data = all_interval_sites_museum_joined$museum_total_records, 
                        nrow= n_sites, ncol= n_intervals)
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    V_citsci = V_citsci, # citizen science detection data
    V_museum = V_museum, # museum detection data
    V_citsci_NA = V_citsci_NA, # array indicating whether sampling occurred in a site*interval*visit
    V_museum_NA = V_museum_NA, # array indicating whether sampling occurred in a site*interval*visit
    n_species = n_species, # number of species
    n_sites = n_sites, # number of sites
    n_intervals = n_intervals, # number of surveys 
    n_visits = n_visits,
    
    intervals = interval_vector,
    sites = site_vector,
    species = species_vector,
    n_genera = n_genera,
    genera = genus_vector,
    genus_lookup = genus_lookup,
    
    raw_pop_density = raw_pop_density,
    pop_densities = pop_density_vector,
    site_areas = site_area_vector,
    herb_shrub_cover = herb_shrub_vector,
    developed_open = developed_open_vector,
    forest = forest_vector,
    herb_shrub_forest = herb_shrub_forest_vector,
    developed_med_high = developed_med_high_vector,
    museum_total_records = museum_total_records,
    
    ecoregion_three_vector = ecoregion_three_vector,
    ecoregion_one_vector = ecoregion_one_vector,
    ecoregion_three_lookup = ecoregion_three_lookup,
    ecoregion_one_lookup = ecoregion_one_lookup,
    n_ecoregion_three = n_ecoregion_three,
    n_ecoregion_one = n_ecoregion_one,
    
    correlation_matrix = correlation_matrix,
    
    total_records = total_records_since_2000,
    total_records_since_study = total_records_since_2000,
    citsci_records = citsci_records,
    citsci_detections = citsci_detections,
    museum_records = museum_records,
    museum_detections = museum_detections,
    species_counts = species_counts,
    species_detections = species_detections,
    
    total_records_full = total_records_since_2000_full,
    total_records_since_study_full = total_records_since_2000_full,
    citsci_records_full = citsci_records_full,
    museum_records_full = museum_records_full,
    species_counts_full = species_counts_full
    
  ))
  
}
