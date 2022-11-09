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

get_spatial_data <- function(
  grid_size,
  min_population_size
){
  
  ## --------------------------------------------------
  # Read in the spatial data and occurrence data
  
  # CRS for NAD83 / UTM Zone 10N
  crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"
  
  # spatial data - California state shapefile
  CA <- tigris::states() %>%
    filter(NAME == "California")
  str(CA)
  st_crs(CA)
  
  # spatial data - 'urban areas'
  # downloaded manually to the directory from:
  # https://maps.princeton.edu/catalog/stanford-zd071bk4213
  urban_areas <- st_read('./data/california_urban_areas/zd071bk4213.shp') %>%
    st_transform(crs)
  
  # occurrence data
  df <- read.csv("./data/data_unfiltered.csv")
  
  ## --------------------------------------------------
  # Prep the data
  
  ## transform
  # transform state shapefile to crs
  CA_trans <- st_transform(CA, 26910) # NAD83 / UTM Zone 10N
  
  # make the df into a spatial file
  (df_sf <- st_as_sf(df,
                     coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = 4326))
  
  # and then transform it to the crs
  df_trans <- st_transform(df_sf, crs = crs)
  
  ## --------------------------------------------------
  # Overlay spatial polygon with grid 
  
  # create _km grid - here you can substitute by specifying grid_size above
  grid <- st_make_grid(CA_trans, cellsize = c(grid_size, grid_size)) %>% 
    st_sf(grid_id = 1:length(.))
  
  # create labels for each grid_id
  # grid_lab <- st_centroid(grid) %>% cbind(st_coordinates(.))
  
  # which grid square is each point in?
  df_id <- df_trans %>% st_join(grid, join = st_intersects) %>% as.data.frame
  
  ## Let's try again with only the grid cells that overlap with urban areas
  urban_grid <- st_join(grid, urban_areas, join = st_intersects)
  
  urban_grid <- urban_grid %>% 
    filter(!is.na(name)) %>% # remove grid cells that don't overlap with urban areas
    # we will select one row per grid cell (one city admin area)
    # keeping the city with the largest population size for now
    group_by(grid_id) %>% 
    slice(which.max(pop2010)) %>%
    # filter out grid cells without significant urban centers
    filter(pop2010 > min_population_size)
  
  # Which grid square is each point in?
  df_id_urban <- df_trans %>% st_join(urban_grid, join = st_intersects) %>% as.data.frame
  
  # Now drop all specimens that DO NOT originate from urban or urban adjacent cells
  # i.e. grid_id is NA
  df_id_urban_filtered <- df_id_urban %>%
    filter(!is.na(grid_id))
  
  return(list(df_id_urban_filtered = df_id_urban_filtered))
  
}