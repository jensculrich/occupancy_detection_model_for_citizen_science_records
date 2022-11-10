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

# This file 'get_spatial_data.R' contains a function 'get_spatial_data' which
# can be called from a 'run_model_X_.R' file with specific data collection
# choices of grid size and minimum pop in the cell to be chosen or varied

# This file walks through this process, with plots to visualize
# how the spatial data is gathered and how data 'collection' choices 
# like grid size and minimum urbanization intensity to include a site
# shape the data that we will analyze with our model

library(tidyverse)
library(tigris) # get state shapefile
library(sf) # spatial data processing
library(raster) # process pop dens raster

## --------------------------------------------------
# Get Data

## global options
# grid size
grid_size <- 25000 # 25km x 25 km
# CRS for NAD83 / UTM Zone 10N
crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"
# minimum population size
# a grid cell must intersect a city with min_population_size or more
# people living there. This filters out grid cell sites that do not
# have a an urban center above some minimum threshold
min_population_size <- 25000 

# spatial data - California state shapefile
CA <- tigris::states() %>%
  filter(NAME == "California")
str(CA)
st_crs(CA)
crs(CA)

# spatial data - 'urban areas' - with population size in 2015
# downloaded manually to the directory from:
# https://maps.princeton.edu/catalog/stanford-zd071bk4213
urban_areas <- st_read('./data/california_urban_areas/zd071bk4213.shp') %>%
  st_transform(crs)

# pop density raster
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
# 2015 pop density at 1km resolution
pop_raster=raster("./data/pden2010_block/pden2010_block/gpw_v4_population_density_rev11_2015_30_sec.tif")
crs(pop_raster)

r2 <- crop(pop_raster, CA)
r3 <- mask(r2, CA)
maxValue(r3)

plot(log(r3+1),
     col=terrain.colors(10),
     alpha=1,
     legend=T,
     main="Log((Population Density/km^2) + 1)")

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

# and let's just plot 100 points to get a picture of the data and shapefile
df_trans_100 <- df_trans %>%
  sample_n(100)

ggplot() +
  geom_sf(data=CA) + 
  # geom_sf(data=urban_areas) +
  geom_sf(data = df_trans_100, 
          aes(fill = species), 
          size = 4, alpha = 0.5, shape = 23) +
  theme(legend.position="none")

# this takes a long long time to plot
r3_df <- as.data.frame(r3, xy = TRUE) %>%
  rename("pop_dens" = "gpw_v4_population_density_rev11_2015_30_sec") %>%
  mutate(logp1_pop_dens = log(pop_dens + 1))
ggplot(r3_df) +
  geom_tile(aes(x=x, y=y, fill=gpw_v4_population_density_rev11_2015_30_sec))
plot(pop_raster)

## --------------------------------------------------
# Overlay spatial polygon with grid 

# create _km grid - here you can substitute by specifying grid_size above
grid <- st_make_grid(CA_trans, cellsize = c(grid_size, grid_size)) %>% 
  st_sf(grid_id = 1:length(.))

# create labels for each grid_id
grid_lab <- st_centroid(grid) %>% cbind(st_coordinates(.))

# view sampled points transposed to the grid on the polygon
ggplot() +
  geom_sf(data = CA_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = sample_n(df_trans, 500), 
          aes(fill = species), 
          size = 4, alpha = 0.5, shape = 23) +
  geom_sf(data = grid, fill = 'transparent', lwd = 0.3) +
  geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="none")

# which grid square is each point in?
df_id <- df_trans %>% st_join(grid, join = st_intersects) %>% as.data.frame

# project to raster
crs_raster <- sf::st_crs(raster::crs(r3))
prj1 <- st_transform(grid, crs_raster)

# Extract raster values to list object
r.vals <- raster::extract(r3, prj1)

# Use list apply to calculate mean for each polygon
r.mean <- lapply(r.vals, FUN=mean)

# Join mean values to polygon data
unlist(r.mean)
grid_pop_dens <- cbind(grid, unlist(r.mean))

## Now let's join the pop dens back with the df
df_id_dens <- left_join(df_id, grid_pop_dens) %>%
  rename("pop_density" = "unlist.r.mean.")



## Let's try again with only the grid cells that overlap with urban areas
urban_grid <- st_join(grid, urban_areas, join = st_intersects)

test <- st_join(df_id_dens, urban_grid, join = st_intersects)

# not using this now
# urban_grid <- urban_grid %>% 
#  filter(!is.na(name)) %>% # remove grid cells that don't overlap with urban areas
  # we will select one row per grid cell (one city admin area)
  # keeping the city with the largest population size for now
#  group_by(grid_id) %>% 
#  slice(which.max(pop2010)) %>%
#  # filter out grid cells without significant urban centers
#  filter(pop2010 > min_population_size)

# create labels for each grid_id
urban_grid_lab <- st_centroid(urban_grid) %>% cbind(st_coordinates(.))

# view sampled points transposed to the grid on the polygon
ggplot() +
  geom_sf(data = CA_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = sample_n(df_trans, 500), 
          aes(fill = species), 
          size = 4, alpha = 0.5, shape = 23) +
  geom_sf(data = urban_grid, fill = 'transparent', lwd = 0.3) +
  geom_text(data = urban_grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="none")

# Which grid square is each point in?
df_id_urban <- df_trans %>% st_join(urban_grid, join = st_intersects) %>% as.data.frame

# Now drop all specimens that DO NOT originate from urban or urban adjacent cells
# i.e. grid_id is NA
df_id_urban_filtered <- df_id_urban %>%
  filter(!is.na(grid_id))

# view sampled urban only occurrence transposed to the grid on the polygon

(df_id_urban_filtered_sf <- st_as_sf(df_id_urban_filtered,
                                     sf_column_name = "geometry", 
                                     crs = crs))

(p <- ggplot() +
    geom_sf(data = CA_trans, fill = 'white', lwd = 0.05) +
    geom_sf(data = sample_n(df_id_urban_filtered_sf, 1000), 
            aes(fill = species), 
            size = 4, alpha = 0.5, shape = 23) +
    geom_sf(data = urban_grid, fill = 'transparent', lwd = 0.3) +
    geom_text(data = urban_grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    theme(legend.position="none"))

# could save data frame of urban occurrences as a .csv
# don't need to keep this file if the prep_data function calls the get_spatial_data function
# write.csv(df_id_urban_filtered, "./data/data_urban_occurrences.csv")

