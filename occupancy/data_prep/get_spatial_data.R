### Gather and prep spatial data for pollinator occurrence
# jcu; started oct 27, 2022

# Collect spatial data to describe where the observations occur i.e.
# We will construct a grid of cells across the spatial extent (Continental U.S.)
# Each cell will be a site, where observations of any species can
# be said to have occurred in any of the years of our timeline.
# Choosing the spatial size of the grid cells will be important,
# and also importantly we will plan to limit the grid cells to those 
# that contain 'urban areas'.

# Then intersect the filtered pollinator occurrence data against the urban grid
# to retain only occurrence data from our sites. This is the data that we will 
# format and then feed to the model

# Finally, extract spatial covariate data including:
# population density
# natural habitat area
# household income
# site area

# The function 'get_spatial_data' which is contained in this file
# can be called from a 'run_model.R' file with specific data collection
# choices of grid size and minimum pop in the cell to be chosen or varied

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
  urban_sites, # are we going to look at only urban areas or any areas
  non_urban_subsample_n, # if not restricting to urban, how many sites to choose (grid of whole country is REALLY computationally unfeasible)
  min_records_per_species, # only include species with total detections greater than this (for purposes of defining range)
  min_unique_detections, # only include species detected at urban sites more times than this (so we can figure out if the species is rare v hard to detect)
  era_start, # starting year for when to start counting number of detections in urban landscapes (to filter by minUniqueDetections)
  by_city, # organize random effects clusters by metropolitan statistical area rather than ecoregion 3
  remove_city_outliers_5stddev, # remove REALLY dense urban landscapes (i.e., downtown NYC)
  canopy_cover_model
){
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  
  center_scale <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) # z-score scale
  }
  
  ## --------------------------------------------------
  # Envrionmental raster data
  
  # USA_Contiguous_Albers_Equal_Area_Conic
  crs <- 5070
  
  # spatial data - state shapefile
  states <- tigris::states() %>%
    # lower 48 + DC
    filter(REGION != 9) %>%
    filter(!NAME %in% c("Alaska", "Hawaii"))
  
  states_trans <- states  %>% 
    st_transform(., crs) # USA_Contiguous_Albers_Equal_Area_Conic
  
  if(canopy_cover_model == TRUE){
    ## --------------------------------------------------
    # Ecoregion data for filtering to eastern US
    
    # read in ecoregion 1 data
    ecoregion1 <- sf::read_sf("./data/spatial_data/na_cec_eco_l1/NA_CEC_ECO_Level1.shp")
    
    ecoregion1_region8 <- filter(ecoregion1, NA_L1CODE == "8") %>%
      st_transform(., crs)
    
    states_trans <- sf::st_intersection(states_trans, ecoregion1_region8)
    
    ## --------------------------------------------------
    # Overlay the shapefile with a grid of sites of size == 'grid_size' 
    
    # create "grid_size" km grid over the area
    grid <- st_make_grid(ecoregion1_region8, cellsize = c(grid_size, grid_size)) %>% 
      st_sf(grid_id = 1:length(.))
  } else {
    ## --------------------------------------------------
    # Overlay the shapefile with a grid of sites of size == 'grid_size' 
    
    # create "grid_size" km grid over the area
    grid <- st_make_grid(states_trans, cellsize = c(grid_size, grid_size)) %>% 
      st_sf(grid_id = 1:length(.))
  }
  
  
  ## --------------------------------------------------
  # Extract mean population density in each grid cell
  
  # pop density raster
  # https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
  # 2015 pop density at 1km resolution
  pop_raster=raster("./data/spatial_data/population_density/gpw_v4_population_density_rev11_2015_30_sec.tif")
  
  # crop and mask pop density to the area of the study
  pop_raster <- crop(pop_raster, states)
  pop_raster <- mask(pop_raster, states)
  
  # project the grid to the raster
  crs_raster <- sf::st_crs(raster::crs(pop_raster))
  prj1 <- st_transform(grid, crs_raster)
  
  # Extract raster values to list object, and then summarize by the mean value
  r.vals <- exactextractr::exact_extract(pop_raster, prj1, 'mean')
  
  ## --------------------------------------------------
  # filter sites to urban areas using the r.vals for pop density
  
  # join population density values to the grid table
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
  # Clip out any remaining grid cells outside of the extent
  
  #grid <-  sf::st_intersection(grid, ecoregion1_region8)
  states_trans <- states_trans %>% st_union()
  
  grid_pop_dens <- sf::st_intersection(grid_pop_dens, states_trans)  
  
  ## --------------------------------------------------
  # get latitude and longitude of each grid cell
  
  site_latlong <- st_centroid(st_transform(grid_pop_dens)) %>%
    st_transform(., 4326) %>%
    dplyr::select(-grid_id, -pop_density_per_km2) %>%
    dplyr::mutate(longitude = sf::st_coordinates(.)[,1],
                  latitude = sf::st_coordinates(.)[,2]) %>%
    as.data.frame()
  
  grid_pop_dens <- cbind(grid_pop_dens, site_latlong[,2:3]) 
    
  ## --------------------------------------------------
  # Extract canopy cover from each remaining site
  
  
  canopy=raster("./data/spatial_data/science_tcc_CONUS_2016_v2021-4/science_tcc_conus_2016_v2021-4.tif")
  gc()
  
  # project the grid to the raster
  crs_raster <- sf::st_crs(raster::crs(canopy))
  prj1 <- st_transform(grid_pop_dens, crs_raster)
  
  # Extract raster values to list object, and then summarize by the mean value
  # r.vals <- exactextractr::exact_extract(canopy, prj1, 'mean')
  r.vals <- exactextractr::exact_extract(canopy, prj1)
  
  r.vals_NA <- lapply(r.vals, function(x) na_if(x$value, 255)) # water values - make NA
  r.vals_NA <- lapply(r.vals_NA, na.omit)
  
  r.mean_tcc <- lapply(r.vals_NA, 
                       function(x) { 
                         mean(x)
                       } 
  ) 
  
  # join population density values to the grid table
  grid_pop_dens <- cbind(grid_pop_dens, unlist(r.mean_tcc)) %>% 
    rename("mean_canopy_cover" = "unlist.r.mean_tcc.")
  
  # view the grid on the polygons
  #ggplot() +
  #geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  #geom_sf(data = grid_pop_dens, aes(fill = mean_canopy_cover), lwd = 0.05) +
  #coord_sf(datum = NA)  +
  #labs(x = "") +
  #labs(y = "") 
  
  gc()
  
  ## --------------------------------------------------
  # Extract impervious surface cover from each remaining site
  
  imp_surface=raster("./data/spatial_data/nlcd_2016_impervious_l48_20210604/nlcd_2016_impervious_l48_20210604.img")
  
  # project the grid to the raster
  crs_raster <- sf::st_crs(raster::crs(imp_surface))
  prj1 <- st_transform(grid_pop_dens, crs_raster)
  
  # Extract raster values to list object, and then summarize by the mean value
  # r.vals <- exactextractr::exact_extract(canopy, prj1, 'mean')
  r.vals <- exactextractr::exact_extract(imp_surface, prj1)
  
  #r.vals_NA <- lapply(r.vals, function(x) na_if(x$value, 255))
  #r.vals_NA <- lapply(r.vals, na.omit)
  
  r.vals_vector <- vector(length = nrow(grid_pop_dens))
  
  for(i in 1:nrow(grid_pop_dens)){
    r.vals_vector[i] <- mean(r.vals[[i]][["value"]])
  }
  
  # join population density values to the grid table
  grid_pop_dens <- cbind(grid_pop_dens, r.vals_vector) %>% 
    rename("mean_imp_surface" = "r.vals_vector")
  
  cor_tcc_imp <- cor(grid_pop_dens$mean_canopy_cover, grid_pop_dens$mean_imp_surface)
  
  # view the grid on the polygons
  #ggplot() +
  #geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  #geom_sf(data = grid_pop_dens, aes(fill = mean_imp_surface), lwd = 0.05) +
  #coord_sf(datum = NA)  +
  #labs(x = "") +
  #labs(y = "") 
  
  ## --------------------------------------------------
  # Extract median household income from each census block group and then average out across grid cell
  # Simultaneously extract race data
  
  # census block-group income data
  # 2020 income data 
  # B19013. Median Household Income in the Past 12 Months (in 2020 Inflation-Adjusted Dollars)
  #https://data2.nhgis.org/main
  # Typically, Block Groups have a population of 600 to 3,000 people.
  # spatial format identifying where block groups are located
  income=st_read("./data/spatial_data/socioeconomic_data/US_blck_grp_2020.shp")
  gc()
  # table format data indicating the income in block groups
  income_data <- read.csv("./data/spatial_data/socioeconomic_data/nhgis0001_ds249_20205_blck_grp.csv")
  
  # calculate median household income RELATIVE to 
  # median household income for the state, 
  # i.e., $80,000 might translate to high purchasing power in Indiana but not necessarily in California
  income_data <- income_data %>%
    group_by(STATE) %>%
    mutate(mean_state_AMR8E001 = mean(AMR8E001, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(relative_AMR8E001 = AMR8E001 / mean_state_AMR8E001) %>%
    dplyr::select(GISJOIN, GEOID, relative_AMR8E001)
  
  income_data <- income_data %>%
    mutate(census_tract = str_sub(GEOID, 6, -2)) 
  
  # table format data indicating the income in block groups
  race_data <- read.csv("./data/spatial_data/DECENNIALDP2020.DP1_2024-06-19T183658/DECENNIALDP2020.DP1-Data.csv")
  
  race_data <- race_data %>%
    rename("census_tract" = "GEO_ID") %>%
    dplyr::select(census_tract, DP1_0078P) %>%
    slice(55:nrow(race_data)) %>% # first rows are states
    mutate(DP1_0078P = as.numeric(DP1_0078P)) %>% # DP1_0078 == Percentage of the population id'd as white, no other race
    filter(!is.na(DP1_0078P))
  
  race_data <- race_data %>%
    mutate(census_tract = str_sub(census_tract, 8))
  
  income_data <- left_join(income_data, race_data, by = c("census_tract"))
  
  # join table data with spatial data
  income <- left_join(income, income_data, by="GISJOIN")
  
  gc()
  
  # transform state shapefile to crs
  income_trans <- st_transform(income, crs) # albers equal area
  income_trans <- income_trans %>%
    dplyr::select(GEOID.y, census_tract, relative_AMR8E001, DP1_0078P)
  
  # join income data with spatial file holding population density data
  grid_pop_dens <- st_join(grid_pop_dens, income_trans)
  
  # we will take the average household income across all block groups 
  # that intersect each grid cell
  grid_pop_dens <- grid_pop_dens %>%
    group_by(grid_id) %>%
    # ignore the NA's which indicate block group areas where no one / very few people live
    # in urban landscapes these are typically airports, military zones, or low density farmland blocks 
    mutate(avg_income = mean(relative_AMR8E001, na.rm = TRUE)) %>%
    mutate(minority = (100 - DP1_0078P)) %>%
    mutate(avg_minority = mean(minority, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(grid_id, pop_density_per_km2, longitude, latitude,
                  mean_canopy_cover, mean_imp_surface,
                  avg_income, avg_minority) %>%
    filter(!is.na(avg_income),
           !is.na(avg_minority)) %>%
    # scale the variable
    mutate(scaled_canopy_cover = center_scale(mean_canopy_cover),
           scaled_imp_surface = center_scale(mean_imp_surface),
           scaled_avg_income = center_scale(avg_income),
           scaled_avg_minority = center_scale(avg_minority))
  
  gc()
  
  # view the grid on the polygons
  #ggplot() +
  #geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  #geom_sf(data = grid_pop_dens, aes(fill = scaled_avg_minority), lwd = 0.05) +
  #coord_sf(datum = NA)  +
  #labs(x = "") +
  #labs(y = "") 
  
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
  
  r.mean_low_dev <- lapply(r.vals_land_NA, 
                           function(x) { 
                             (length(which(x %in% c(22))) / length(x))
                           } 
  )
  
  r.mean_high_dev <- lapply(r.vals_land_NA, 
                            function(x) { 
                              (length(which(x %in% c(23,24))) / length(x))
                            } 
  )
  
  r.mean_low_med_high_dev <- lapply(r.vals_land_NA, 
                                    function(x) { 
                                      (length(which(x %in% c(22,23,24))) / length(x))
                                    } 
  )
  
  r.mean_forest <- lapply(r.vals_land_NA, 
                          function(x) { 
                            (length(which(x %in% c(41,42,43))) / length(x))
                          } 
  )
  
  # this is what we use as "natural habitat area" covariate
  r.mean_herb_shrub_forest <- lapply(r.vals_land_NA, 
                                     function(x) { 
                                       (length(which(x %in% c(52,71,41,42,43,90,95))) / length(x))
                                     }
  )                                
  
  # provide option to remove site outliers that are extremely different in pop dens (i.e. downtown NYC)
  if(remove_city_outliers_5stddev == TRUE){
    
    grid_pop_dens <- cbind(grid_pop_dens,
                           unlist(r.site_area),
                           unlist(r.mean_herb_shrub),
                           unlist(r.mean_dev_open),
                           unlist(r.mean_low_dev),
                           unlist(r.mean_high_dev),
                           unlist(r.mean_low_med_high_dev),
                           unlist(r.mean_forest),
                           unlist(r.mean_herb_shrub_forest)) %>%
      rename("site_area" = "unlist.r.site_area.",
             "herb_shrub_cover" = "unlist.r.mean_herb_shrub.",
             "developed_open" = "unlist.r.mean_dev_open.",
             "developed_low" = "unlist.r.mean_low_dev.",
             "developed_med_high" = "unlist.r.mean_high_dev.",
             "developed_low_med_high" = "unlist.r.mean_low_med_high_dev.",
             "forest" = "unlist.r.mean_forest.",
             "herb_shrub_forest" = "unlist.r.mean_herb_shrub_forest.") %>%
      # remove sites below filter for minimum site area
      filter(site_area > min_site_area) %>%
      # center-scale the variables
      mutate(scaled_pop_den_km2 = center_scale(pop_density_per_km2),
             scaled_site_area = center_scale(site_area),
             scaled_herb_shrub_cover = center_scale(herb_shrub_cover),
             scaled_developed_open = center_scale(developed_open),
             scaled_developed_low = center_scale(developed_low),
             scaled_developed_med_high = center_scale(developed_med_high),
             scaled_developed_low_med_high = center_scale(developed_low_med_high),
             scaled_forest = center_scale(forest),
             scaled_herb_shrub_forest = center_scale(herb_shrub_forest)
      ) %>%
      # remove extreme population outliers (more than 5 sd from the mean)
      filter(scaled_pop_den_km2 < 5) %>%
      filter(scaled_pop_den_km2 > -5) %>%
      # rescale the variables (in case any sites removed)
      mutate(scaled_pop_den_km2 = center_scale(pop_density_per_km2),
             scaled_site_area = center_scale(site_area),
             scaled_herb_shrub_cover = center_scale(herb_shrub_cover),
             scaled_developed_open = center_scale(developed_open),
             scaled_developed_low = center_scale(developed_low),
             scaled_developed_med_high = center_scale(developed_med_high),
             scaled_developed_low_med_high = center_scale(developed_low_med_high),
             scaled_forest = center_scale(forest),
             scaled_herb_shrub_forest = center_scale(herb_shrub_forest)
      )
    
  } else{
    
    grid_pop_dens <- cbind(grid_pop_dens,
                           unlist(r.site_area),
                           unlist(r.mean_herb_shrub),
                           unlist(r.mean_dev_open),
                           unlist(r.mean_low_dev),
                           unlist(r.mean_high_dev),
                           unlist(r.mean_low_med_high_dev),
                           unlist(r.mean_forest),
                           unlist(r.mean_herb_shrub_forest)) %>%
      rename("site_area" = "unlist.r.site_area.",
             "herb_shrub_cover" = "unlist.r.mean_herb_shrub.",
             "developed_open" = "unlist.r.mean_dev_open.",
             "developed_low" = "unlist.r.mean_low_dev.",
             "developed_med_high" = "unlist.r.mean_high_dev.",
             "developed_low_med_high" = "unlist.r.mean_low_med_high_dev.",
             "forest" = "unlist.r.mean_forest.",
             "herb_shrub_forest" = "unlist.r.mean_herb_shrub_forest.") %>%
      # remove sites below filter for minimum site area
      filter(site_area > min_site_area) %>%
      # center-scale the variables
      mutate(scaled_pop_den_km2 = center_scale(pop_density_per_km2),
             scaled_site_area = center_scale(site_area),
             scaled_herb_shrub_cover = center_scale(herb_shrub_cover),
             scaled_developed_open = center_scale(developed_open),
             scaled_developed_low = center_scale(developed_low),
             scaled_developed_med_high = center_scale(developed_med_high),
             scaled_developed_low_med_high = center_scale(developed_low_med_high),
             scaled_forest = center_scale(forest),
             scaled_herb_shrub_forest = center_scale(herb_shrub_forest)
      )
  }                               
  
  # clear the work space
  rm(crs_raster, grid, land, prj1, 
     r.mean_dev_open, r.mean_forest, r.mean_herb_shrub, r.mean_high_dev, r.mean_herb_shrub_forest,
     r.site_area, r.vals_land, r.vals_land_NA, r.mean_low_med_high_dev, r.mean_low_dev,
     states, states_trans)
  gc()
  
  ## --------------------------------------------------
  # Determine spatial cluster of each site
  # (will need to move this below site filter)
  
  ## --------------------------------------------------
  # Join with Metro areas or Ecoregion 3
  
  if(by_city == TRUE){ # if TRUE, cluster by metro area
    
    # Metropolitan statistical areas (CBSA's)
    # https://catalog.data.gov/dataset/tiger-line-shapefile-2019-nation-u-s-current-metropolitan-statistical-area-micropolitan-statist
    # updated 2018, based on 2010 census data
    level_three_cluster <- sf::read_sf("./data/spatial_data/tl_2019_us_cbsa/tl_2019_us_cbsa.shp")
    
    ## level 3 cluster (city)
    crs_level_three <- sf::st_crs(raster::crs(level_three_cluster))
    prj1 <- st_transform(grid_pop_dens, crs_level_three)
    
    level_three_names <- st_join(prj1, CBSA) %>%
      group_by(grid_id) %>%
      slice(which.max(ALAND)) %>% 
      pull(NAME) 
    
    level_three_vector <- as.numeric(as.factor(level_three_names))
    
    level_three_lookup <- as.numeric(as.factor(level_three_names))
    
    n_level_three <- unique(level_three_lookup) %>%
      length()
    
  } else{ # else cluster by ecoregion 3 
    
    level_three_cluster <- sf::read_sf("./data/spatial_data/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp")
    
    ## level 3 cluster (ecoregion 3 fine scale ecoregion)
    crs_level_three <- sf::st_crs(raster::crs(level_three_cluster))
    prj1 <- st_transform(grid_pop_dens, crs_level_three)
    
    level_three_names <- st_join(prj1, level_three_cluster) %>%
      group_by(grid_id) %>%
      slice(which.max(Shape_Area)) %>% 
      pull(NA_L3CODE) 
    
    level_three_vector <- as.numeric(as.factor(level_three_names))
    
    level_three_lookup <- as.numeric(as.factor(level_three_names))
    
    n_level_three <- unique(level_three_vector) %>%
      length()
  }
  
  # join the level three cluster info to the site data
  grid_pop_dens <- cbind(grid_pop_dens, level_three_vector)
  
  ## --------------------------------------------------
  # Ecoregion data for site clustering
  
  # read in ecoregion 1 data
  ecoregion1 <- sf::read_sf("./data/spatial_data/na_cec_eco_l1/NA_CEC_ECO_Level1.shp")
  
  ## level 4 cluster (ecoregion1)
  crs_ecoregion1 <- sf::st_crs(raster::crs(ecoregion1))
  prj1 <- st_transform(grid_pop_dens, crs_ecoregion1)
  
  # names of the ecoregions in case interested in which ones have higher occurrence/detection
  ecoregion_one_names <- st_join(prj1, ecoregion1) %>%
    group_by(grid_id) %>%
    slice(which.max(Shape_Area)) %>% 
    pull(NA_L1CODE)
  
  # as a simple vector
  ecoregion_one_vector <- as.numeric(as.factor(ecoregion_one_names))
  
  # number of ecoregions
  n_ecoregion_one <- unique(ecoregion_one_vector) %>%
    length()
  
  # join ecoregion one info with the sites
  grid_pop_dens <- cbind(grid_pop_dens, 
                         ecoregion_one_vector)
  
  # script for five tiered clusters
  #ecoregion_three_lookup <- grid_pop_dens %>%
  #group_by(CBSA_vector, ecoregion_three_vector) %>%
  #slice(1) %>%
  #ungroup() %>%
  #group_by(CBSA_vector) %>%
  #slice(1) %>%
  #ungroup() %>%
  #pull(ecoregion_three_vector)
  
  #ecoregion_one_lookup <- grid_pop_dens %>%
  #group_by(ecoregion_three_vector, ecoregion_one_vector) %>%
  #slice(1) %>%
  #ungroup() %>%
  #group_by(ecoregion_three_vector) %>%
  #slice(1) %>%
  #ungroup() %>%
  #pull(ecoregion_one_vector)
  
  # just do level 3 clustered directly in eco1 as level 4
  # which ecoregion 1 is each ecoregion 3 located within?
  ecoregion_one_lookup <- grid_pop_dens %>%
    group_by(level_three_vector, ecoregion_one_vector) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(level_three_vector) %>%
    slice(1) %>%
    ungroup() %>%
    pull(ecoregion_one_vector)
  
  # also get info on which eco 1, every grid cell is in (the default model does not use this info)
  # which ecoregion 1 is each grid cell located within?
  ecoregion_one_lookup_by_grid_cell <- grid_pop_dens %>%
    group_by(grid_id, ecoregion_one_vector) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(grid_id) %>%
    slice(1) %>%
    ungroup() %>%
    pull(ecoregion_one_vector)
  
  
  # pull city names for sites (not used for analysis but interesting to track when plotting
  # out the site-specific random effects)
  CBSA <- sf::read_sf("./data/spatial_data/tl_2019_us_cbsa/tl_2019_us_cbsa.shp")
  
  ## level 3 cluster (city)
  crs_CBSA <- sf::st_crs(raster::crs(CBSA))
  prj1 <- st_transform(grid_pop_dens, crs_CBSA)
  
  # which CBSA (metro area) is each grid cell located within?
  CBSA_names <- st_join(prj1, CBSA) %>%
    group_by(grid_id) %>%
    slice(which.max(ALAND)) %>% 
    pull(NAME) 
  
  # clear the workspace
  rm(prj1, 
     level_three_cluster, 
     ecoregion1)
  
  
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
  #df$decimalLongitude <- na_if(df$decimalLongitude, '')
  #df$decimalLatitude <- na_if(df$decimalLatitude, '')
  
  df <- df %>% 
    # filter out any records with no spatial data
    filter(!is.na(decimalLongitude)) %>%
    filter(!is.na(decimalLatitude)) %>%
    
    # filter out species that we have so few records that we can't 
    # confidently say which cities they could occur in
    group_by(species) %>%
    add_tally(name="records_per_species") %>%
    filter(records_per_species > min_records_per_species) %>%
    ungroup()
  
  # convert to a shapefile
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
  
  # figure out how many unique site/year detections for each species
  # we will filter out species that don't meet a minimum threshold of unique occurrence detections
  species_with_enough_detections <- df_id_dens %>%
    filter(year >= era_start) %>%
    group_by(species, grid_id, year) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(species) %>%
    count(name="total_detections") %>%
    filter(total_detections > min_unique_detections) %>%
    filter(species != "") %>%
    pull(species)
  
  # now filter out species that were detected at sites fewer than a minimum number of times
  df_id_dens <- df_id_dens %>% 
    filter(species %in% species_with_enough_detections)
  
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
  
  # manually crop species occurrences from outside of core range
  # later, share references for "core range" and show plots for excluded and included record points
  if(taxon == "bombus"){
    df_id_dens <- df_id_dens %>% 
      filter(!(species == "impatiens" & decimalLongitude < -100)) %>%
      filter(!(species == "perplexus" & decimalLongitude < -100)) %>%
      filter(!(species == "pensylvanicus" & decimalLongitude < -120)) %>%
      filter(!(species == "pensylvanicus" & ((state.prov %in% c("Idaho", "Montana", "North Dakota")))))  %>%
      filter(!(species == "affinis" & (!(state.prov %in% 
                                           c("Minnesota", "Iowa", "Wisconsin", "Illinois",
                                             "Indiana", "Ohio", "West Virginia", "Virginia"))))) %>%
      filter(!(species == "vosnesenskii" & (!(state.prov %in% c("California", "Oregon", 
                                                                "Washington", "Idaho", "Nevada")))))
  }
  
  # free unused space
  rm(df, df_sf, df_trans)
  gc()
  
  
  ## --------------------------------------------------
  # Generate a matrix that will later be used to 
  # calculate correlations between predictor variables
  # this will be looked at later in a separate file to consider colinearity among variables
  
  correlation_matrix <- as.matrix(as.data.frame(cbind(
    grid_pop_dens$scaled_pop_den_km2,
    grid_pop_dens$scaled_site_area,
    grid_pop_dens$scaled_avg_income,
    grid_pop_dens$scaled_avg_minority,
    grid_pop_dens$scaled_canopy_cover,
    grid_pop_dens$scaled_imp_surface,
    grid_pop_dens$scaled_developed_open,
    grid_pop_dens$scaled_developed_low,
    grid_pop_dens$scaled_developed_med_high,
    grid_pop_dens$scaled_developed_low_med_high,
    grid_pop_dens$scaled_herb_shrub_cover,
    grid_pop_dens$scaled_forest,
    grid_pop_dens$scaled_herb_shrub_forest)))

  #} # end else
  
  #output <- as.data.frame(grid_pop_dens) 
  #output <- output[,-32]
  #output <- cbind(output, CBSA_names)
  #write.csv(as.data.frame(output), "./data/spatial_data/site_data.csv")
  
  #write.csv(as.data.frame(cbind(ecoregion_one_names, level_three_names)), "./data/spatial_data/eco3.csv")
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    df_id_urban_filtered = df_id_dens,
    urban_grid = grid_pop_dens,
    ecoregion_one_lookup = ecoregion_one_lookup,
    n_ecoregion_one = n_ecoregion_one,
    ecoregion_one_names = ecoregion_one_names,
    level_three_names = level_three_names,
    level_three_lookup = level_three_lookup,
    ecoregion_one_lookup_by_grid_cell = ecoregion_one_lookup_by_grid_cell,
    n_level_three = n_level_three,
    correlation_matrix = correlation_matrix,
    CBSA_names = CBSA_names
    
  ))
}

