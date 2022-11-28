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
    (x - mean(x)) / sd(x)
  }
  
  ## --------------------------------------------------
  # Read in the spatial data and occurrence data
  
  # CRS for NAD83 / UTM Zone 10N
  crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"
  
  # spatial data - California state shapefile
  CA <- tigris::states() %>%
    filter(NAME == "California")
  
  CA_nad83 <- CA  %>% 
    st_transform(., 26910) # NAD83 / UTM Zone 10N
  
  # spatial data - 'urban areas'
  # downloaded manually to the directory from:
  # https://maps.princeton.edu/catalog/stanford-zd071bk4213
  urban_areas <- st_read('./data/california_urban_areas/zd071bk4213.shp') %>%
    st_transform(crs)
  
  # pop density raster
  # https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
  # 2015 pop density at 1km resolution
  pop_raster=raster("./data/pden2010_block/pden2010_block/gpw_v4_population_density_rev11_2015_30_sec.tif")
  
  # occurrence data
  df <- read.csv("./data/data_unfiltered.csv") 
  
  df_nad83 <- df %>%
    # spatial projection of the occurrence data
    st_as_sf(.,
             coords = c("decimalLongitude", "decimalLatitude"), 
             crs = 4326) %>%
    # transformed to projection of the CA shapefile
    st_transform(., crs = crs)
  
  ## --------------------------------------------------
  # Overlay the shapefile with a grid of sites of size == 'grid_size' 
  
  # create _km grid - here you can substitute by specifying grid_size above
  grid <- st_make_grid(CA_nad83, cellsize = c(grid_size, grid_size)) %>% 
    st_sf(grid_id = 1:length(.))
  
  # Determine which grid cell (site) each occurrence record is from 
  df_id <- df_nad83 %>% st_join(grid, join = st_intersects) %>% as.data.frame
  
  ## --------------------------------------------------
  # Extract mean population density in each grid cell
  
  ## --------------------------------------------------
  # Prep the population density raster
  r2_pop_dens <- crop(pop_raster, CA)
  r3_pop_dens <- mask(r2_pop_dens, CA)
  
  # project the grid to the raster
  crs_raster <- sf::st_crs(raster::crs(r3_pop_dens))
  prj1 <- st_transform(grid, crs_raster)
  
  # Extract raster values to list object
  # this takes a while since there are many raster cells with their own values
  # in each grid cell.
  r.vals <- raster::extract(r3_pop_dens, prj1)
  
  # Use list apply to calculate mean raster value for each grid cell
  # ignore na.s (values where the grid cell does not overlap with land area)
  r.mean <- lapply(r.vals, FUN=mean, na.rm=TRUE)
  
  # Join mean values to the grid cell data
  grid_pop_dens <- cbind(grid, unlist(r.mean)) %>% 
    rename("pop_density_per_km2" = "unlist.r.mean.")
  
  # free unused space
  rm(r3_pop_dens, r2_pop_dens, pop_raster)
  gc(verbose = FALSE)
  
  ## --------------------------------------------------
  # add pop density covariate data to the occurrence data
  # AND simultaneously FILTER THE SITES TO ONLY THOSE THAT ARE CONSIDERED URBAN
  # so that we can make comparisons across urban areas
  
  # Join the pop dens data back with the df that now also has grid (site) names
  # We also filter out non-urban sites 
  # (and sites that are outside of the admin area and thus have NA population density)
  # and add a scaled density variable
  df_id_dens <- left_join(df_id, grid_pop_dens) %>%
    left_join(., dplyr::select(df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") %>%
    # Filter to the sites that are actually 'urban' 
    filter(pop_density_per_km2 > min_population_size) %>%
    # scale the population density
    mutate(scaled_pop_den_km2 = center_scale(pop_density_per_km2))
  
  ## --------------------------------------------------
  # Add city data (e.g. city name if that end up being of interest)
  
  # first project the data again using the lat and long colums
  df_w_dens_sf <- st_as_sf(df_id_dens,
                           coords = c("decimalLongitude", "decimalLatitude"), 
                           crs = 4326) %>%
      # and then transform it to the correct crs again
      st_transform(df_w_dens_sf, crs = crs)
  
  # join with city data
  # first remove the geometry for the individual collection points
  st_geometry(df_w_dens_sf) <- NULL # remove point geometry, coerce to data.frame
  df_w_dens_sf <- rename(df_w_dens_sf, "geometry" = ".")
  df_w_dens_sf <- st_set_geometry(df_w_dens_sf, df_w_dens_sf$geometry) # set geometry, return sf
  
  urban_grid_w_cities <- st_join(df_w_dens_sf, urban_areas, join = st_intersects)
  
  urban_grid_prepped <- urban_grid_w_cities %>% 
    # we will select one row per grid cell (one city admin area)
    # keeping the city with the largest population size for now
    group_by(grid_id) %>% 
    slice(which.max(pop2010)) %>%
    ungroup() %>%
    # this will hold a geometry for both the grid AND an intersecting record
    # but we specifically need the grid geometry for the site area calculation below
    rename("city" = "name",
           "total_pop_of_largest_intersecting_city" = "pop2010") %>%
    dplyr::select(grid_id, pop_density_per_km2, scaled_pop_den_km2, 
                  city, total_pop_of_largest_intersecting_city, geometry) 
  
  ## --------------------------------------------------
  # Pull vector data from the urban grid cells
  
  scaled_pop_density <- urban_grid_prepped %>%
    pull(scaled_pop_den_km2)
  
  city_names <- urban_grid_prepped %>%
    pull(city)
  
  site_name_vector <- as.character(
    urban_grid_prepped %>%
    group_by(grid_id) %>%
    slice(1) %>% # take one row per site (the name of each site)
    ungroup() %>%
    dplyr::select(grid_id) %>% # extract species names column as vector
    pull(grid_id)
    )
  
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
  scaled_grid_area <- grid_intersect %>% 
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
  
  return(list(df_id_urban_filtered = df_id_dens,
              scaled_pop_density = scaled_pop_density,
              scaled_grid_area = scaled_grid_area,
              city_names = city_names,
              site_name_vector = site_name_vector))
  
}