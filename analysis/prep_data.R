### Prepeare data to feed to model
# jcu; started oct 27, 2022

# will need to assign occupancy intervals and visit numbers
# for occupancy (as opposed to abundance), will need to 
# filter down to one unique occurrence per species*site*interval*year

## The data we will need to prepare to feed to the model:
# V <- V # an array of detection data
# n_species # number of species
# n_sites # number of sites
# n_intervals # number of occupancy intervals 
# n_visits # number of visits in each interval

library(tidyverse)

prep_data <- function(era_start, era_end, n_intervals, n_visits, min_records_per_species,
                      grid_size, min_population_size) {
  
  # spatially explicit occurrence data
  # df <- read.csv("./data/data_urban_occurrences.csv")
  
  source("./data/get_spatial_data.R")
  
  my_spatial_data <- get_spatial_data(grid_size, min_population_size)
  df <- my_spatial_data$df_id_urban_filtered
  
  ## --------------------------------------------------
  # assign study dimensions
  
  total_years <- era_end - era_start + 1
  remainder <- total_years %% n_intervals
  n_visits <- n_visits
  min_records_per_species <- min_records_per_species
  
  ## --------------------------------------------------
  # assign occupancy visits and intervals and reduce to one
  # sample per species per site per visit
  
  df_filtered <- df %>%
    
    # remove records (if any) missing species level identification
    filter(species != "")
  
  # assign year as - year after era_start
  mutate(occ_year = (year - era_start)) %>% # need to -1 so the start year is year 0
    
    
    # remove data from years that are in the remainder
    # occupancy intervals have to be equal in length for the model to process
    # e.g. we can't have intervals of 3 years and then one interval with only 1 year
    # filter(occ_year %in% (remainder:1)) %>%
    
    # remove years before the start date
    filter(occ_year >= 0) %>%
    
    # now assign the years into 1:n_intervals
    mutate(occ_interval = occ_year %/% n_visits) %>%
    
    # add a sampling round (1:n)
    mutate(visit = (occ_year %% n_visits)) %>%
    
    # remove species with total observations (n) < min_records_per_species 
    group_by(species) %>%
    add_tally() %>%
    filter(n >= min_records_per_species) %>%
    ungroup() %>%
    
    # one unique row per visit*site*species combination
    group_by(visit, grid_id, species) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    select(species, grid_id, occ_interval, occ_year, visit) %>%
    arrange((occ_year))
  # end pipe
  
  ## --------------------------------------------------
  # Extract stan data from df
  
  ## Get unique species and sites
  # create an alphabetized list of all species encountered across all sites*intervals*visits
  species_list <- df_filtered %>%
    group_by(species) %>%
    slice(1) %>% # take one row per species (the name of each species)
    ungroup() %>%
    select(species) # extract species names column as vector
  
  # create an alphabetized list of all sites
  site_list <- df_filtered %>%
    group_by(grid_id) %>%
    slice(1) %>% # take one row per species (the name of each species)
    ungroup() %>%
    select(grid_id) # extract species names column as vector
  
  # get vectors of species, sites, intervals, and visits 
  species_vector <- species_list %>%
    pull(species)
  
  site_vector <- as.character(site_list %>%
                                pull(grid_id))
  
  interval_vector <- as.vector(levels(as.factor(df_filtered$occ_interval)))
  
  visit_vector <- as.vector(levels(as.factor(df_filtered$visit)))
  
  # find study dimensions from the pooled data
  n_species <- length(species_vector)
  
  n_sites <- length(site_vector)
  
  n_intervals <- length(interval_vector)
  
  n_visits <- length(visit_vector)
  
  ## --------------------------------------------------
  ## Now we are ready to create the detection matrix, V
  
  V <- array(data = NA, dim = c(n_species, n_sites, n_intervals, n_visits))
  
  for(i in 1:n_intervals){
    for(j in 1:n_visits){
      
      # iterate this across visits within intervals
      temp <- df_filtered %>%
        # filter to indices for interval and visit
        filter(occ_interval == (i - 1), # intervals start at 0 so need to subtract 1 from i
               visit == (j - 1)) %>%  # visits start at 0 so need to subtract 1 from j
        # now join with all species (so that we include species not captured during 
        # this interval*visit but which might actually be at some sites undetected)
        full_join(species_list, by="species") %>%
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        full_join(site_list, by="grid_id") %>%
        # separate_rows(grid_id, sep = ",") %>%
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
        select(levels(as.factor(site_vector))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(species  %in% levels(as.factor(df$species)))
      
      
      # convert from dataframe to matrix
      temp_matrix <- as.matrix(temp)
      # remove species names
      temp_matrix <- temp_matrix[,-1]
      
      # replace NAs for the interval i and visit j with the matrix
      V[1:n_species, 1:n_sites,i,j] <- temp_matrix[1:n_species, 1:n_sites]
      
    }
  }
  
  class(V) <- "numeric"
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    V = V, # detection data
    n_species = n_species, # number of species
    n_sites = n_sites, # number of sites
    n_intervals = n_intervals, # number of surveys 
    n_visits = n_visits,
    
    intervals = interval_vector,
    sites = site_vector,
    species = species_vector
    
  ))
  
}
