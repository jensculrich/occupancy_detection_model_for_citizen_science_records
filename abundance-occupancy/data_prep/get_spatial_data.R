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

library(tidyverse)
library(tigris) # get state shapefile
library(sf) # spatial data processing
library(raster) # process pop dens raster

get_spatial_data <- function(
  grid_size,
  min_population_size
){
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  ## Use preprocessed data to save time, if study dimensions match preprocessing reqts.
  if(grid_size == 30000 && min_population_size == 300){
    
    df_id_dens <- readRDS("./preprocessed_data/occurrence_records_30km_200minpop_.RDS")
    
    grid_pop_dens <- readRDS("./preprocessed_data/site_data_30km_200minpop_.RDS")
    
    ## --------------------------------------------------
    # Extract variables
    
    scaled_pop_density <- grid_pop_dens %>% 
      pull(scaled_pop_den_km2)
    
    scaled_impervious_cover <- grid_pop_dens %>% 
      pull(scaled_impervious_cover)
    
    scaled_tin <- grid_pop_dens %>% 
      pull(scaled_tin_cover)
    
    scaled_site_area <- grid_pop_dens %>% 
      pull(scaled_site_area)
    
    ## --------------------------------------------------
    # Calculate correlations between variables 
    
    my_variables <- as.data.frame(cbind(scaled_pop_density, 
                                        scaled_impervious_cover, 
                                        scaled_tin,
                                        scaled_site_area))
    
    correlation_matrix <- cor(my_variables)
    
  } else{
    
    ## --------------------------------------------------
    # occurrence data
    df <- read.csv("./data/occurrence_data/syrphidae/data_unfiltered.csv") 
    
    ## --------------------------------------------------
    # Spatial extent and urban areas
    
    # CRS for NAD83 / UTM Zone 10N
    crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"
    
    # spatial data - California state shapefile
    states <- tigris::states() %>%
      filter(NAME %in% c("California", "Oregon", "Washington", "Arizona", "Nevada"))
    
    states_trans <- states  %>% 
      st_transform(., 26910) # NAD83 / UTM Zone 10N
    
    ## --------------------------------------------------
    # Environmental data rasters
    
    # pop density raster
    # https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
    # 2015 pop density at 1km resolution
    pop_raster=raster("./data/spatial_data/population_density/gpw_v4_population_density_rev11_2015_30_sec.tif")
    
    # (aggregated) impervious surface cover raster
    # https://www.mrlc.gov/data/nlcd-2016-percent-developed-imperviousness-conus
    # note this impervious surface cover has already been aggregated from 30mx30m to 300mx300m
    imp <- raster::raster(
      "./data/spatial_data/impervious_surface/300m_aggregate_crop_nlcd_2016_impervious.tif")
    
    # Time integrated NDVI raster
    # 2014 TIN surface cover (Time-Integrated NDVI)
    # https://earthexplorer.usgs.gov/ 
    # > Vegetation Monitoring > Phenology > eMODIS Phenology > 250m res - 2014 - Western NA V2
    tin=raster::raster("./data/spatial_data/time_integrated_ndvi/TIN2014_wUSAeM250m_v2.tif")
    
    ## --------------------------------------------------
    # Overlay the shapefile with a grid of sites of size == 'grid_size' 
    
    # create _km grid - here you can substitute by specifying grid_size above
    grid <- st_make_grid(states_trans, cellsize = c(grid_size, grid_size)) %>% 
      st_sf(grid_id = 1:length(.))
    
    # Determine which grid cell (site) each occurrence record is from 
    df <- df %>%
      # spatial projection of the occurrence data
      st_as_sf(.,
               coords = c("decimalLongitude", "decimalLatitude"), 
               crs = 4326) %>%
      # transformed to projection of the CA shapefile
      st_transform(., crs = crs) %>%
      st_join(grid, join = st_intersects) %>% as.data.frame
    
    ## --------------------------------------------------
    # Extract mean population density in each grid cell
    
    crs_raster <- sf::st_crs(raster::crs(pop_raster))
    
    # Join mean values to the grid cell data
    grid_pop_dens <- cbind(grid, 
                           
                           # Use list apply to calculate mean raster value for each grid cell
                           # ignore na's (values where the grid cell does not overlap with land area)                       
                           unlist(lapply(
                             raster::extract(raster::mask( # extract values from 
                               raster::crop(pop_raster, states), states), # the cropped and masked population density raster
                               st_transform(grid, crs_raster)), # for each grid cell that has been transformed to the crs of the raster
                             
                             FUN=mean, na.rm=TRUE) # and calculate the mean from all values extracted from a grid cell
                           )) 
    
    grid_pop_dens$unlist.lapply.raster..extract.raster..mask.raster..crop.pop_raster..[
      grid_pop_dens$unlist.lapply.raster..extract.raster..mask.raster..crop.pop_raster.. == "NaN"] <- NA
    
    grid_pop_dens <- grid_pop_dens %>%
      # and rename the extracted variable
      rename("pop_density_per_km2" = "unlist.lapply.raster..extract.raster..mask.raster..crop.pop_raster..") %>%
      
      # and scale the population density
      mutate(scaled_pop_den_km2 = center_scale(pop_density_per_km2))
    
    # free unused space
    rm(pop_raster)
    gc(verbose = FALSE)
    
    ## --------------------------------------------------
    # add pop density covariate data to the occurrence data
    # AND simultaneously filter out all the associated occurrence data
    # from the sites that are NOT CONSIDERED URBAN AREAS
    # so that we can make comparisons about insect populations among cities
    
    # Join the pop dens data back with the df that now also has grid (site) names
    # We also filter out non-urban sites 
    # (and sites that are outside of the admin area and thus have NA population density)
    # and add a scaled density variable
    df_id_dens <- left_join(df_id, grid_pop_dens) %>%
      left_join(., dplyr::select(df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") %>%
      # Filter to the sites that are actually 'urban' 
      filter(pop_density_per_km2 > min_population_size)
    
    # free unused space
    rm(df, df_id)
    gc(verbose = FALSE)
    
    ## --------------------------------------------------
    # Now give just one row per grid cell, with columns including
    # the grid cell name, the pop density and scaled pop density, and the geometry
    
    # first project the data again using the lat and long colums
    urban_grid_prepped <- grid_pop_dens %>%
      filter(pop_density_per_km2 > min_population_size) %>%
      # this will hold a geometry for both the grid AND an intersecting records
      # but we specifically need the grid geometry for the site area calculation below
      dplyr::select(grid_id, 
                    pop_density_per_km2, scaled_pop_den_km2)
    
    ## --------------------------------------------------
    # Pull vector data from the urban grid cells
    
    scaled_pop_density <- urban_grid_prepped %>%
      pull(scaled_pop_den_km2)
    
    site_name_vector <- as.character(
      urban_grid_prepped %>%
        group_by(grid_id) %>%
        slice(1) %>% # take one row per site (the name of each site)
        ungroup() %>%
        dplyr::select(grid_id) %>% # extract species names column as vector
        pull(grid_id)
    )
    
    ## --------------------------------------------------
    # Now get impervious surface cover from each grid cell
    
    # project the urban sites to the impervious surface cover raster
    crs_imp <- sf::st_crs(raster::crs(imp))
    prj2 <- st_transform(urban_grid_prepped, crs_imp)
    
    # Extract raster values to list object
    # this takes a while since there are many raster cells with their own values
    # in each grid cell.
    # Use list apply to calculate mean raster value for each grid cell
    r.vals_imp <- raster::extract(
      imp <- mask(imp, prj2), # mask the impervious surface cover to these sites
      prj2) # before extracting from these sites
    
    # values of 127 are values of water or outside the mask
    # for i in 1: the number of sites
    for(i in 1:length(r.vals_imp)){
      # go through and replace 127's with NA's
      r.vals_imp[[i]][r.vals_imp[[i]]==127] <- NA
    }
    
    # now calculate the mean value per grid cell
    r.mean_imp <- lapply(r.vals_imp, FUN=mean, na.rm=TRUE)
    
    impervious_cover <- unlist(r.mean_imp) 
    
    scaled_impervious_cover <- center_scale(impervious_cover)
    
    # free unused space
    rm(imp, prj2, r.vals_imp, r.mean_imp, impervious_cover)
    gc(verbose = FALSE)
    
    ## --------------------------------------------------
    # Add raster extraction for any additional spatial covariates here...
    
    ## --------------------------------------------------
    # Calculate land area of grid cells 
    # some cells might partially be outside of the area where we are getting records from
    # e.g. a cell half in California and half in Mexico or a cell that is along the
    # coastline and only overlaps slightly with land
    # we would expect fewer species to occur in these smaller areas and therefore should
    # account for site area (extent of grid cell intersection w/ shapefile) in our analysis
    
    # THE VECTOR 'scaled_grid_area' is the output that lists scaled site area 
    # in order from lowest site number to highest.
    
    # intersect between administrative area and the grid
    grid_intersect <- st_intersection(CA_nad83, urban_grid_prepped)
    #plot(CA_nad83$geometry, axes = TRUE)
    #plot(urban_grid_prepped$geometry, add = TRUE)
    #plot(grid_intersect$geometry, add = TRUE, col = 'red')
    
    # for each field, get area overlapping with admin area
    scaled_site_area <- grid_intersect %>% 
      # add in areas in m2
      mutate(area = st_area(.) %>% as.numeric()) %>% 
      # if you have multiple administrative areas, say we add another state
      # we'd want to summarize the total area across these admin areas
      as_tibble() %>% 
      group_by(grid_id) %>% 
      summarize(area = sum(area)) %>%
      # create a scaled variable
      mutate(scaled_site_area = center_scale(area)) %>%
      pull(scaled_site_area)
    
    # free unused space
    rm(grid_intersect)
    gc(verbose = FALSE)
    
    ## --------------------------------------------------
    # Calculate correlations between variables 
    my_variables <- as.data.frame(cbind(scaled_pop_density, 
                                        scaled_impervious_cover, 
                                        scaled_site_area))
    
    correlation_matrix <- cor(my_variables)
    
  }
  
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    df_id_urban_filtered = df_id_dens,
    scaled_pop_density = scaled_pop_density,
    scaled_impervious_cover = scaled_impervious_cover,
    scaled_site_area = scaled_site_area,
    site_name_vector = site_name_vector,
    urban_grid = urban_grid_prepped,
    correlation_matrix = correlation_matrix))
  
}