### Gather and prep spatial data for pollinator occurrence
# jcu; started oct 27, 2022

# Collect spatial data to describe where the observations occur i.e.
# We will construct a grid of cells across the spatial extent (California)
# Each cell will be a site, where observations of any species can
# be said to have occurred in any of the years of our timeline.
# Choosing the spatial size of the grid cells will be important,
# and also importantly we will plan to limit the grid cells to those 
# that contain 'urban areas'.

# THEN intersect the filtered pollinator occurrence data against the urban grid
# to retain only occurrence data from our sites. This is the data that we will 
# format and then feed to the model (data_urban_occurrences.csv)

# Finally, extract spatial covariate data including:
# population density
# impervious surface cover
# site area

# This file contains a function 'get_spatial_data' which
# can be called from a 'run_model_X_.R' file with specific data collection
# choices of grid size and minimum pop in the cell to be chosen or varied

# The file 'explore_spatial_data.R' walks through this process, with plots to visualize
# how the spatial data is gathered and how data 'collection' choices 
# like grid size and minimum urbanization intensity to include a site
# shape the data that we will analyze with our model

library(tidyverse) # data carpentry
library(tigris) # get state shapefile
library(sf) # spatial data processing
library(raster) # read and format raster data
library(exactextractr) # quick extraction of raster data

get_spatial_data <- function(
  grid_size, # square edge dimensions of a site, in meters
  min_population_size, # minimum pop density of a site to be considered "urban"
  taxon, # prepare data for syrphidae or bombus
  min_site_area, # min land area (in the state admin areas) of a site to be included
  urban_sites,
  non_urban_subsample_n,
  min_records_per_species
){
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  ## To save time and computing power, if study dimensions match preprocessing reqts.
  ## then use the preprocessed data to save time, 
  #if(grid_size == 30000 && min_population_size == 300 && min_site_area == 0.10){
    
    # read the occurrence data for the given taxon
    #df_id_dens <- readRDS(paste0(
      #"./preprocessed_data/",
      #taxon,
      #"_occurrence_records_30km_300minpop_.RDS"))
    
    #grid_pop_dens <- readRDS("./preprocessed_data/site_data_30km_300minpop_.RDS")
    
    ## --------------------------------------------------
    # Extract variables
    
    #scaled_pop_density <- grid_pop_dens %>% 
      #pull(scaled_pop_den_km2)
    
    #scaled_site_area <- grid_pop_dens %>% 
      #pull(scaled_site_area)
    
    #scaled_developed_open <- grid_pop_dens %>% 
      #pull(scaled_developed_open)
    
    #scaled_herb_shrub_cover <- grid_pop_dens %>% 
      #pull(scaled_herb_shrub_cover)
    
    #scaled_forest <- grid_pop_dens %>% 
      #pull(scaled_forest)
    
    #scaled_developed_med_high <- grid_pop_dens %>% 
      #pull(scaled_developed_med_high)
    
    ## --------------------------------------------------
    # Calculate correlations between variables 
    
    #my_variables <- as.data.frame(cbind(scaled_pop_density, 
                                        #scaled_site_area, 
                                        #scaled_developed_open,
                                        #scaled_developed_med_high,
                                        #scaled_herb_shrub_cover, 
                                        #scaled_forest))
    
    #correlation_matrix <- cor(my_variables)
    
  #} else{
    
    ## --------------------------------------------------
    # Envrionmental raster data
    
    # USA_Contiguous_Albers_Equal_Area_Conic
    crs <- 5070
    
    # spatial data - state shapefile
    states <- tigris::states() %>%
      # lower 48 + DC
      filter(REGION != 9) %>%
      filter(!NAME %in% c("Alaska", "Hawaii"))
    
    #str(states)
    #st_crs(states)
    
    states_trans <- states  %>% 
      st_transform(., crs) # USA_Contiguous_Albers_Equal_Area_Conic
    
    ## --------------------------------------------------
    # Overlay the shapefile with a grid of sites of size == 'grid_size' 
    
    # create _km grid - here you can substitute by specifying grid_size above
    grid <- st_make_grid(states_trans, cellsize = c(grid_size, grid_size)) %>% 
      st_sf(grid_id = 1:length(.))
    
    ## --------------------------------------------------
    # Extract mean population density in each grid cell
    
    # pop density raster
    # https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
    # 2015 pop density at 1km resolution
    pop_raster=raster("./data/spatial_data/population_density/gpw_v4_population_density_rev11_2015_30_sec.tif")
    
    pop_raster <- crop(pop_raster, states)
    pop_raster <- mask(pop_raster, states)
    
    # project the grid to the raster
    crs_raster <- sf::st_crs(raster::crs(pop_raster))
    prj1 <- st_transform(grid, crs_raster)
    
    # Extract raster values to list object, and then summarize by the mean value
    r.vals <- exactextractr::exact_extract(pop_raster, prj1, 'mean')
    
    ## --------------------------------------------------
    # filter sites to urban areas using the r.vals for pop density
    
    grid_pop_dens <- cbind(grid, r.vals) %>% 
      rename("pop_density_per_km2" = "r.vals")
    
    # now filter out the non-urban areas (areas below our pop density threshold)
    # or if urban_sites is false, filter to some non-urban areas and subsample n rows
    if(urban_sites == FALSE){
      grid_pop_dens <- grid_pop_dens %>%
        filter(pop_density_per_km2 < min_population_size) %>%
            slice_sample(., n = non_urban_subsample_n, replace = FALSE) %>%
        # need to make grid id in ascending numeric order otherwise chaos 
        # will ensue later if trying to set up the NA indicator array! 
        # when you use slice_sample (or any form of sample) they are randomly ordered
        mutate(grid_id = row_number())
    } else {
      grid_pop_dens <- grid_pop_dens %>%
        filter(pop_density_per_km2 > min_population_size) 
    }
    
    
    # free unused space
    rm(pop_raster, prj1, r.vals, crs_raster)
    gc()
    
    ## --------------------------------------------------
    # Extract median household income from each census block group and then average out across grid cell
    
    # census block-group income data
    # 2020 income data 
    # B19013. Median Household Income in the Past 12 Months (in 2020 Inflation-Adjusted Dollars)
    #https://data2.nhgis.org/main
    # Typically, Block Groups have a population of 600 to 3,000 people.
    income=st_read("./data/spatial_data/socioeconomic_data/US_blck_grp_2020.shp")
    gc()
    income_data <- read.csv("./data/spatial_data/socioeconomic_data/nhgis0001_ds249_20205_blck_grp.csv")
    
    # calculate median household income RELATIVE to 
    # median household income for the state, 
    # i.e., $80,000 might be really high for WV but relatively low for CA.
    income_data <- income_data %>%
      group_by(STATE) %>%
      mutate(mean_state_AMR8E001 = mean(AMR8E001, na.rm=TRUE)) %>%
      ungroup() %>%
      mutate(relative_AMR8E001 = AMR8E001 / mean_state_AMR8E001) %>%
      dplyr::select(GISJOIN, relative_AMR8E001)
    
    income <- left_join(income, income_data, by="GISJOIN")
    rm(income_data)
    gc()
    
    # transform state shapefile to crs
    income_trans <- st_transform(income, crs) # albers equal area
    
    grid_pop_dens <- st_join(grid_pop_dens, income_trans)
    grid_pop_dens <- grid_pop_dens %>%
      group_by(grid_id) %>%
      mutate(avg_income = mean(relative_AMR8E001, na.rm = TRUE)) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::select(grid_id, pop_density_per_km2, avg_income) %>%
      filter(!is.na(avg_income)) %>%
      mutate(scaled_avg_income = center_scale(avg_income))
    
    rm(income, income_trans)
    gc()
    
    
    ## --------------------------------------------------
    # Extract environmental variables from each remaining site
    
    # land cover raster 30m x 30m
    # https://www.mrlc.gov/data/nlcd-2016-land-cover-conus
    # land cover data from 2016 
    land=raster::raster("./data/spatial_data/land_cover/land_cover/nlcd_2016_land_cover_l48_20210604.img")
    
    # make sure that the grid is still projected to the raster
    crs_raster <- sf::st_crs(raster::crs(land))
    prj1 <- st_transform(grid_pop_dens, crs_raster)
    
    # then extract values and cbind with the grid_pop_dens
    # extract raster values to list object
    
    r.vals_land <- exactextractr::exact_extract(land, prj1)
    gc(verbose = FALSE)
    
    # want to reclassify open water as NA
    r.vals_land_NA <- lapply(r.vals_land, function(x) na_if(x$value, 0))
    r.vals_land_NA <- lapply(r.vals_land_NA, function(x) na_if(x,11))
    gc(verbose = FALSE)
    
    # first find out the proportion of rows that are NA's 
    # (this will be our site area)
    r.site_area <- lapply(r.vals_land_NA, 
                          function(x) { 
                            # which rows are not NAs (are not masked or coded as open water)
                            # divided by the number of rows
                            # yields a site area
                            (length(which(!is.na(x))) / length(x))
                          } 
    ) 
    gc(verbose = FALSE)
    
    # now drop NA values so that the below estimates are the proportion of cover
    # of all land cover in the administrative area
    r.vals_land_NA <- lapply(r.vals_land_NA, na.omit)
    gc(verbose = FALSE)
    
    # now pull out site proportion of each type
    # for legend of category number codes see: 
    # https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
    
    r.mean_herb_shrub <- lapply(r.vals_land_NA, 
                                function(x) { 
                                  (length(which(x %in% c(52,71))) / length(x))
                                } 
    ) 
    
    r.mean_dev_open <- lapply(r.vals_land_NA, 
                              function(x) { 
                                (length(which(x %in% c(21))) / length(x))
                              } 
    ) 
    
    r.mean_high_dev <- lapply(r.vals_land_NA, 
                              function(x) { 
                                (length(which(x %in% c(23,24))) / length(x))
                              } 
    ) 
    
    r.mean_forest <- lapply(r.vals_land_NA, 
                            function(x) { 
                              (length(which(x %in% c(41,42,43))) / length(x))
                            } 
    )
    
    r.mean_herb_shrub_forest <- lapply(r.vals_land_NA, 
                                       function(x) { 
                                         (length(which(x %in% c(52,71,41,42,43,90,95))) / length(x))
                                       }
    )                                
                                       
    grid_pop_dens <- cbind(grid_pop_dens,
                           unlist(r.site_area),
                           unlist(r.mean_herb_shrub),
                           unlist(r.mean_dev_open),
                           unlist(r.mean_high_dev),
                           unlist(r.mean_forest),
                           unlist(r.mean_herb_shrub_forest)) %>%
      rename("site_area" = "unlist.r.site_area.",
             "herb_shrub_cover" = "unlist.r.mean_herb_shrub.",
             "developed_open" = "unlist.r.mean_dev_open.",
             "developed_med_high" = "unlist.r.mean_high_dev.",
             "forest" = "unlist.r.mean_forest.",
             "herb_shrub_forest" = "unlist.r.mean_herb_shrub_forest.") %>%
      # remove sites below filter for minimum site area
      filter(site_area > min_site_area) %>%
      # center-scale the variables
      mutate(scaled_pop_den_km2 = center_scale(pop_density_per_km2),
             scaled_site_area = center_scale(site_area),
             scaled_herb_shrub_cover = center_scale(herb_shrub_cover),
             scaled_developed_open = center_scale(developed_open),
             scaled_developed_med_high = center_scale(developed_med_high),
             scaled_forest = center_scale(forest),
             scaled_herb_shrub_forest = center_scale(herb_shrub_forest)
      )
    
    rm(crs_raster, grid, land, prj1, 
       r.mean_dev_open, r.mean_forest, r.mean_herb_shrub, r.mean_high_dev, r.mean_herb_shrub_forest,
       r.site_area, r.vals_land, r.vals_land_NA, states, states_trans)
    gc()
    
    ## --------------------------------------------------
    # Determine spatial cluster of each site
    # (will need to move this below site filter)
    
    ## --------------------------------------------------
    # Ecoregion data for site clustering
    
    ecoregion3 <- sf::read_sf("./data/spatial_data/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp")
    ecoregion1 <- sf::read_sf("./data/spatial_data/na_cec_eco_l1/NA_CEC_ECO_Level1.shp")
    
    ## level 3 cluster (ecoregion3)
    crs_ecoregion3 <- sf::st_crs(raster::crs(ecoregion3))
    prj1 <- st_transform(grid_pop_dens, crs_ecoregion3)
    
    ecoregion_three_names <- st_join(prj1, ecoregion3) %>%
      group_by(grid_id) %>%
      slice(which.max(Shape_Area)) %>% 
      pull(NA_L3CODE) 
    
    ecoregion_three_vector <- as.numeric(as.factor(ecoregion_three_names))
    
    ecoregion_three_lookup <- as.numeric(as.factor(ecoregion_three_names))
    
    n_ecoregion_three <- unique(ecoregion_three_lookup) %>%
      length()
    
    ## level 4 cluster (ecoregion1)
    crs_ecoregion1 <- sf::st_crs(raster::crs(ecoregion1))
    prj1 <- st_transform(grid_pop_dens, crs_ecoregion1)
    
    ecoregion_one_names <- st_join(prj1, ecoregion1) %>%
      group_by(grid_id) %>%
      slice(which.max(Shape_Area)) %>% 
      pull(NA_L1CODE)
    
    ecoregion_one_vector <- as.numeric(as.factor(ecoregion_one_names))
    
    n_ecoregion_one <- unique(ecoregion_one_vector) %>%
      length()
    
    grid_pop_dens <- cbind(grid_pop_dens, ecoregion_three_vector, ecoregion_one_vector)
    
    ecoregion_one_lookup <- grid_pop_dens %>%
      group_by(ecoregion_three_vector, ecoregion_one_vector) %>%
      slice(1) %>%
      ungroup() %>%
      group_by(ecoregion_three_vector) %>%
      slice(1) %>%
      ungroup() %>%
      pull(ecoregion_one_vector)
    
    rm(prj1, ecoregion3, ecoregion1)
    
    ## --------------------------------------------------
    # Occurrence data
    
    # read either the syrphidae data or the bombus data
    if(taxon == "syrphidae"){
      df <- read.csv(paste0("./data/occurrence_data/", taxon, "_data_all.csv"))
    } else {
      df <- read.csv(paste0("./data/occurrence_data/bbna_private/bbna_trimmed.csv"))
    }
    
    ## --------------------------------------------------
    # Prep the data
    
    # if not doing urban sites, we want to filter the occurrence data so that 
    # we still only fit a model for the species that we see in the urban sites
    # there are more species that occur anywhere than in urban sites alone
    if(urban_sites == FALSE){
      species_names <- readRDS(paste0("./figures/species_names/", 
                                      taxon,
                                      "_names_15km_urban.RDS"))
      
      # then filter out species that never occur at our urban sites
      df <- df %>% 
        filter(species %in% species_names)
    }
    
    # make the df into a spatial file
    df$decimalLongitude <- na_if(df$decimalLongitude, '')
    df$decimalLatitude <- na_if(df$decimalLatitude, '')
    
    df <- df %>% 
      filter(!is.na(decimalLongitude)) %>%
      filter(!is.na(decimalLatitude)) %>%
      
      # filter out species that we have so few records that we can't 
      # confidently say which cities they could occur in
      group_by(species) %>%
      add_tally(name="records_per_species") %>%
      filter(records_per_species > min_records_per_species) %>%
      ungroup()
      
    
    (df_sf <- st_as_sf(df,
                       coords = c("decimalLongitude", "decimalLatitude"), 
                       crs = 4326))
    
    # and then transform it to the crs
    df_trans <- st_transform(df_sf, crs = crs(grid_pop_dens))
    
    ## --------------------------------------------------
    # Now tag records with site ID's based on spatial intersection
    
    # which grid square is each point in?
    df_id_dens <- df_trans %>% 
    
      st_join(grid_pop_dens, join = st_intersects) %>% as.data.frame %>%
      # filter out records from outside of the urban grid
      filter(!is.na(grid_id)) %>%
      left_join(., dplyr::select(
        df, id, decimalLatitude, decimalLongitude), by="id") %>% 
      
      # and perform any further initial data filters
      
      filter(decimalLatitude < 50) %>% # remove any points from alaska (or untagged with state name but from alaska)
    
      # filter out records with high location uncertainty (threshold at 10km)
      # assuming na uncertainty (large portion of records) is under threshold
      # here we assume records with no listed uncertainty are precise or at least precise enough to fall within our big sites
      mutate(coordinateUncertaintyInMeters = replace_na(coordinateUncertaintyInMeters, 0)) %>%
      filter(coordinateUncertaintyInMeters < 10000)
    
    # When we include random sites (lots of remote locations)
    # we end up with many sites with no or only one or two detections ever.
    # We will remove sites with 1 or fewer total records
    # remove sites with no records
    if(urban_sites == FALSE){
      
      df_id_dens <- df_id_dens %>%
        # group by site
        group_by(grid_id) %>%
        # calculate total records (rows) per site
        add_tally(name="total_records_per_site") %>%
        # filter out rows where total records = 0
        ungroup() %>%
        filter(total_records_per_site > 1)
      
    }
    
    # manually crop bombus occurrences from outside of core range
    # later, share references for "core range" and show plots for excluded and included record points
    if(taxon == "bombus"){
      df_id_dens <- df_id_dens %>% 
        filter(!(species == "impatiens" & decimalLongitude < -100)) %>%
        filter(!(species == "affinis" & (!(state.prov %in% 
                                             c("Minnesota", "Iowa", "Wisconsin", "Illinois",
                                               "Indiana", "Ohio", "West Virginia", "Virginia")))))
    }
    
    # free unused space
    rm(df, df_sf, df_trans)
    gc()
    
    
    ## --------------------------------------------------
    # Calculate correlations between variables
    
    correlation_matrix <- cor(as.data.frame(cbind(grid_pop_dens$scaled_pop_den_km2,
                                                  grid_pop_dens$scaled_site_area,
                                                  grid_pop_dens$scaled_developed_open,
                                                  grid_pop_dens$scaled_developed_med_high,
                                                  grid_pop_dens$scaled_herb_shrub_cover,
                                                  grid_pop_dens$scaled_forest,
                                                  grid_pop_dens$scaled_herb_shrub_forest,
                                                  grid_pop_dens$scaled_avg_income)))
    
    colnames(correlation_matrix) <- c("scaled_pop_den_km2", 
                                      "scaled_site_area", 
                                      "scaled_developed_open",
                                      "scaled_developed_med_high",
                                      "scaled_herb_shrub_cover",
                                      "scaled_forest",
                                      "scaled_herb_shrub_forest",
                                      "scaled_avg_income")
    
    rownames(correlation_matrix) <- c("scaled_pop_den_km2", 
                                      "scaled_site_area", 
                                      "scaled_developed_open",
                                      "scaled_developed_med_high",
                                      "scaled_herb_shrub_cover",
                                      "scaled_forest",
                                      "scaled_herb_shrub_forest",
                                      "scaled_avg_income")
    
  #} # end else
  
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    df_id_urban_filtered = df_id_dens,
    urban_grid = grid_pop_dens,
    ecoregion_three_lookup = ecoregion_three_lookup,
    ecoregion_one_lookup = ecoregion_one_lookup,
    n_ecoregion_three = n_ecoregion_three,
    n_ecoregion_one = n_ecoregion_one,
    correlation_matrix = correlation_matrix
    
  ))
}

