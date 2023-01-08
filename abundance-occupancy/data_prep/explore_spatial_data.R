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
## Operation Functions
## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

## --------------------------------------------------
# occurrence data
df <- read.csv("./data/occurrence_data//syrphidae/data_unfiltered.csv")

## --------------------------------------------------
# Study design

## global options
# grid size
grid_size <- 35000 # e.g., 25000 = 25km x 25 km sites
# CRS for NAD83 / UTM Zone 10N
crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"
# minimum population size
# a grid cell must intersect a city with min_population_size or more
# people living there. This filters out grid cell sites that do not
# have a an urban center above some minimum threshold

# min_population_size of 38 (/km^2) is ~ 100/mile^2 which is a typical threshold for 
# considering an area to be 'urban'
# let's up the minimum a bit and go with 100 per sq km, which is about 260/sq mile
min_population_size <- 200 

## --------------------------------------------------
# Spatial extent and urban areas

# spatial data - California state shapefile
states <- tigris::states() %>%
  filter(NAME %in% c("California", "Oregon", "Washington", "Arizona", "Nevada"))
str(states)
st_crs(states)
crs(states)

# spatial data - 'urban areas' - with population size in 2015
# downloaded manually to the directory from:
# https://maps.princeton.edu/catalog/stanford-zd071bk4213
#urban_areas <- st_read('./data/california_urban_areas/zd071bk4213.shp') %>%
#  st_transform(crs)

## --------------------------------------------------
# ecological rasters

# human pop density raster
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
# 2015 pop density at 1km resolution
pop_raster=raster("./data/spatial_data/population_density/gpw_v4_population_density_rev11_2015_30_sec.tif")
crs(pop_raster)

# impervious surface raster
# 2016 urban imp surface cover https://www.mrlc.gov/data/nlcd-2016-percent-developed-imperviousness-conus
imp=raster::raster("D:/urban_spatial_data/nlcd_2016_impervious_l48_20210604/nlcd_2016_impervious_l48_20210604.img")
raster::crs(imp)

## --------------------------------------------------
# Prep the data

## transform
# transform state shapefile to crs
states_trans <- st_transform(states, 26910) # NAD83 / UTM Zone 10N

# make the df into a spatial file
(df_sf <- st_as_sf(df,
                   coords = c("decimalLongitude", "decimalLatitude"), 
                   crs = 4326))

# and then transform it to the crs
df_trans <- st_transform(df_sf, crs = crs)

ggplot() +
  geom_sf(data=states)

# and let's just plot 100 points to get a picture of the data and shapefile
df_trans_100 <- df_trans %>%
  sample_n(100)

ggplot() +
  geom_sf(data=states) + 
  # geom_sf(data=urban_areas) +
  geom_sf(data = df_trans_100, 
          aes(fill = species), 
          size = 4, alpha = 0.5, shape = 23) +
  theme(legend.position="none")

## --------------------------------------------------
# Overlay spatial polygon with grid 

# create _km grid - here you can substitute by specifying grid_size above
grid <- st_make_grid(states_trans, cellsize = c(grid_size, grid_size)) %>% 
  st_sf(grid_id = 1:length(.))

# create labels for each grid_id
grid_lab <- st_centroid(grid) %>% cbind(st_coordinates(.))

# view sampled points transposed to the grid on the polygon
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  #geom_sf(data = sample_n(df_trans, 500), 
  #         aes(fill = species), 
  #        size = 4, alpha = 0.5, shape = 23) +
  geom_sf(data = grid, fill = 'transparent', lwd = 0.3) +
  #geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="none")


## --------------------------------------------------
# Extract pop density from each grid cell

r2 <- crop(pop_raster, states)
r3 <- mask(r2, states)
maxValue(r3)
r3

plot(log(r3+1),
     col=rev(terrain.colors(10)),
     alpha=1,
     legend=T,
     main=expression("(log) Population Density/km"^2),
     xlab="Longitude", ylab="Latitude")

my_window <- extent(-125, -112, 32, 42)
plot(my_window, col=NA, xlab="longitude", ylab = "latitude")
plot(log(r3+1), add=T)

# project the grid to the raster
crs_raster <- sf::st_crs(raster::crs(r3))
prj1 <- st_transform(grid, crs_raster)

# Extract raster values to list object
# this takes a while since there are many raster cells with their own values
# in each grid cell.
r.vals <- raster::extract(r3, prj1)

# Use list apply to calculate mean raster value for each grid cell
r.mean <- lapply(r.vals, FUN=mean, na.rm=TRUE)

# this takes a long long time to plot
#r3_df <- as.data.frame(r3, xy = TRUE) %>%
#  rename("pop_dens" = "gpw_v4_population_density_rev11_2015_30_sec") %>%
#  mutate(logp1_pop_dens = log(pop_dens + 1))
#ggplot(r3_df) +
#  geom_tile(aes(x=x, y=y, fill=gpw_v4_population_density_rev11_2015_30_sec))
#plot(pop_raster)

# Join mean values to polygon data
grid_pop_dens <- cbind(grid, unlist(r.mean)) %>% 
  rename("pop_density_per_km2" = "unlist.r.mean.")

# now filter out the non-urban areas (areas below our pop density threshold)
grid_pop_dens <- grid_pop_dens %>%
  filter(pop_density_per_km2 > min_population_size) %>%
  mutate(scaled_pop_den_km2 = center_scale(pop_density_per_km2))


# and gather the site names and data that will be handled by the model
scaled_pop_density <- grid_pop_dens %>%
  pull(scaled_pop_den_km2)

grid_id_names <- as.character(grid_pop_dens %>%
                                pull(grid_id))

## --------------------------------------------------
# Tag records with site ID's based on spatial intersection

# which grid square is each point in?
df_id_dens <- df_trans %>% 
  st_join(grid_pop_dens, join = st_intersects) %>% as.data.frame %>%
  filter(!is.na(grid_id)) %>%
  left_join(., dplyr::select(
    df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") 


## --------------------------------------------------
# plot the spatial data

# Reproject the occurrence data if you want to plot random sample of NHC records
df_w_dens_sf <- st_as_sf(df_id_dens,
                         coords = c("decimalLongitude", "decimalLatitude"), 
                         crs = 4326)

# and then transform it to the crs again
df_w_dens_trans <- st_transform(df_w_dens_sf, crs = crs)

# create labels for each grid_id
urban_grid_lab <- st_centroid(grid_pop_dens) %>% cbind(st_coordinates(.))

# view sampled points transposed to the grid on the polygon
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, lwd = 0.3) +
  #scale_fill_gradient2() +
  geom_sf(data = sample_n(df_w_dens_trans, 100), 
          aes(fill = species), 
          size = 4, alpha = 0.5, shape = 23) +
  #geom_text(data = urban_grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  ggtitle("Random sample of 100 NHC records from urban site areas", 
          subtitle = "(coloured by species)") +
  theme(legend.position = "none")

# view sites only with population density
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = pop_density_per_km2), lwd = 0.3) +
  scale_fill_gradient2(name = expression("Population/km"^2)) +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Population density in urban areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 
  # coord_sf(datum = NA)

# view sites only with scaled population density
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_pop_den_km2), lwd = 0.3) +
  scale_fill_gradient2(name = expression("Scaled population/km"^2)) +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Population density in urban areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 
# coord_sf(datum = NA)

## --------------------------------------------------
# Calculate land area of grid cells 
# some cells might partially be outside of the area where we are getting records from
# e.g. a cell half in California and half in Mexico or a cell that is along the
# coastline and only overlaps slightly with land
# we would expect fewer species to occur in these smaller areas and therefore should
# account for site area (extent of grid cell intersection w/ shapefile) in our analysis

# THE VECTOR 'scaled_grid_area' is the output that lists scaled site area 
# in order from lowest site number to highest.

# intersect - note that sf is intelligent with attribute data!
grid_intersect <- st_intersection(states_trans, grid_pop_dens)
plot(states_trans$geometry, axes = TRUE, 
     xlab = " x pixel", ylab = "y pixel")
plot(grid_pop_dens$., add = TRUE)
plot(grid_intersect$geometry, add = TRUE, col = 'red')
title("Site x Land Area Intersection")

# add in areas in m2
attArea <- grid_intersect %>% 
  mutate(area = st_area(.) %>% as.numeric())

# for each field, get area overlapping with admin area
attArea <- attArea %>% 
  as_tibble() %>% 
  group_by(grid_id) %>% 
  summarize(area = sum(area)) %>%
  mutate(scaled_site_area = center_scale(area)) 

scaled_grid_area <- attArea %>% 
  pull(scaled_site_area)

grid_pop_dens <- grid_pop_dens %>%
  left_join(., attArea, by="grid_id")

# view sites only with scaled data
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_site_area), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled site area") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Within-extent land area of urban sites") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 
# coord_sf(datum = NA)

## --------------------------------------------------
# Calculate land impervious surface cover value of each grid cell

# The ## code below should be used if you want to reaggregate the raster data.

# project the raster to the first raster used (population density - "r3")
crs_imp <- sf::st_crs(raster::crs(imp))
crs_raster <- sf::st_crs(raster::crs(r3))

# this is too costly for my computer
# imp2 <- projectRaster(from=imp, to=r3)

# try reducing extent first
# aggregate at scale x (to a coarser scale to speed the processing and save comp load)
# 1/30meters = x/500 meters -> x = 16.7
# 1/30meters = x/300 meters -> x = 10
imp_agg <- aggregate(imp, 10, fun=mean)

writeRaster(imp_agg, 
            "D:/urban_spatial_data/nlcd_2016_impervious_l48_20210604/300m_aggregate_nlcd_2016_impervious_35km_CONUS.tif",
            overwrite=TRUE)

# crop to smaller bounding box to make more manageable 
extent(imp_agg)
new_extent <- extent(-2400000, -1100000, 950000, 3200000)
imp_agg_cropped <- crop(x = imp_agg, y = new_extent)

writeRaster(imp_agg_cropped, 
            "D:/urban_spatial_data/nlcd_2016_impervious_l48_20210604/300m_aggregate_crop_nlcd_2016_impervious_35km_CONUS.tif",
            overwrite=TRUE)

# now let's try this again with the smaller file
imp2 <- projectRaster(from=imp_agg_cropped, to=r3)

# project the grid to the raster
crs_raster <- sf::st_crs(raster::crs(imp2))
prj1 <- st_transform(grid, crs_raster)

plot(imp2)
plot(prj1, colour = NA, add = TRUE)


# imp2[imp2 == 127] <- NA


# project the grid to the raster

# originally was: crop, aggregate, mask
crs_imp <- sf::st_crs(raster::crs(imp))
prj2 <- st_transform(grid_pop_dens, crs_imp)
extent(prj2)

# try.. crop extent to western us, mask to sites, crop extent, aggregate
# these are all taking a long time.. pick back up later.
imp <- mask(imp, prj2)
imp2 <- crop(imp, prj2)

# aggregate at scale x (to a coarser scale to speed the processing and save comp load)
# 1/30meters = x/500 meters -> x = 16.7
# 1/30meters = x/300 meters -> x = 10
imp_agg <- aggregate(imp, 10, fun=mean)

writeRaster(imp_agg, "D:/urban_spatial_data/nlcd_2016_impervious_l48_20210604/300m_aggregate_nlcd_2016_impervious_35km_2.tif")

imp_agg <- raster("D:/urban_spatial_data/nlcd_2016_impervious_l48_20210604/300m_aggregate_nlcd_2016_impervious_35km.tif")
imp_agg <- mask(imp_agg, prj2)

# Extract raster values to list object
# this takes a while since there are many raster cells with their own values
# in each grid cell.
# Use list apply to calculate mean raster value for each grid cell

r.vals_imp <- raster::extract(imp_agg, prj2)

for(i in 1:length(r.vals_imp)){
  r.vals_imp[[i]][r.vals_imp[[i]]==127] <- NA
}

r.mean_imp <- lapply(r.vals_imp, FUN=mean, na.rm=TRUE)

grid_pop_dens <- cbind(grid_pop_dens, unlist(r.mean_imp)) %>%
  rename("impervious_cover" = "unlist.r.mean_imp.") %>%
  mutate(scaled_impervious_cover = center_scale(impervious_cover))

# view sites only with scaled data
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_impervious_cover), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled impervious cover") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Developed Impervious Surface of Urban Areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 

# pretty high correlation here..
cor(grid_pop_dens$scaled_impervious_cover, grid_pop_dens$scaled_pop_den_km2)

## --------------------------------------------------
# Calculate land area of grid cells 

# save RDS for easy access later

#saveRDS(as.data.frame(grid_pop_dens), "./preprocessed_data/site_data_35km_200minpop_.RDS")
#saveRDS(df_id_dens, "preprocessed_data/occurrence_records_35km_200minpop_.RDS") 


# below in development
## --------------------------------------------------
# Calculate land TIN value of each grid cell
crs(tin)
plot(tin,
     #col=rev(terrain.colors(10)),
     legend=T,
     main="TIN")


# project the grid to the raster
crs_tin <- sf::st_crs(raster::crs(tin))
prj1 <- st_transform(grid, crs_tin)

tin2 <- crop(tin, CA)
tin3 <- mask(tin2, CA)
maxValue(tin3)

plot(tin3,
     col=rev(terrain.colors(10)),
     alpha=1,
     legend=T,
     main="NDVI")

my_window <- extent(-125, -112, 32, 42)
plot(my_window, col=NA, xlab="longitude", ylab = "latitude")
plot(tin3, add=T)


# this takes a long long time to plot
tin_df <- as.data.frame(tin3, xy = TRUE) # %>%
#mutate(logp1_pop_dens = log(pop_dens + 1))
plot <- ggplot(tin_df) +
  geom_sf(data=CA) 

plot <- plot + geom_tile(aes(x=x, y=y, fill=tin_2014w_ProjectRaster1)) 

