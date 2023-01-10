### Gather and prep spatial data for pollinator occurrence
# jcu; started oct 27, 2022

# Collect spatial data to describe where the observations occur i.e.
# We will construct a grid of cells across the spatial extent (5 Western States)
# Each cell will be a site, where observations of any species can
# be said to have occurred in any of the years of our timeline.
# Choosing the spatial size of the grid cells will be important,
# and we will reiterate the analysis to understand the sensitivity to this decision.
# and also importantly we will plan to limit the grid cells to those 
# that contain 'urban areas' with some minimum pop dense.

# THEN intersect the filtered pollinator occurrence data against the urban grid
# to retain only occurrence data from our sites. This is the data that we will 
# format and then feed to the model 

# This file walks through this process, with plots to visualize
# how the spatial data is gathered and how data 'collection' choices 
# like grid size and minimum urbanization intensity to include a site
# shape the data that we will analyze with our model.

# The partner file 'get_spatial_data.R' contains a function 'get_spatial_data' which
# can be called from a 'run_model_X_.R' file with specific data collection
# choices of grid size and minimum pop in the cell to be chosen or varied.
# In this way the data can be regathered over an interative loop to test sensitivity.

# load packages
library(tidyverse) # data carpentry
library(tigris) # get state shapefile
library(sf) # spatial data processing
library(raster) # raster data processing

## --------------------------------------------------
# Operation Functions
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
grid_size <- 30000 # e.g., 25000 = 25km x 25 km sites
# CRS for NAD83 / UTM Zone 10N
crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"
# minimum population size
# a grid cell must intersect a city with min_population_size or more
# people living there. This filters out grid cell sites that do not
# have a an urban center above some minimum threshold

# min_population_size of 38 (/km^2) is ~ 100/mile^2 which is a typical threshold for 
# considering an area to be 'urban'
# let's up the minimum a bit and go with 100 per sq km, which is about 260/sq mile
min_population_size <- 300 

## --------------------------------------------------
# Spatial extent and urban areas

# spatial data - California state shapefile
states <- tigris::states() %>%
  filter(NAME %in% c("California", "Oregon", "Washington", "Arizona", "Nevada"))
str(states)
st_crs(states)
crs(states)

## --------------------------------------------------
# ecological rasters 

# human pop density raster
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
# 2015 pop density at 1km resolution
pop_raster=raster("./data/spatial_data/population_density/gpw_v4_population_density_rev11_2015_30_sec.tif")
crs(pop_raster)

# MRLC data https://www.mrlc.gov/data

# DO NOT READ IF USING THE AGGREGATED FILE PRODUCED ONE TIME ONLY BELOW
# impervious surface raster
# 2016 urban imp surface cover https://www.mrlc.gov/data/nlcd-2016-percent-developed-imperviousness-conus
#imp=raster::raster("D:/urban_spatial_data/nlcd_2016_impervious_l48_20210604/nlcd_2016_impervious_l48_20210604.img")
#raster::crs(imp)

# TIN raster
# 2014 TIN surface cover (Time-Integrated NDVI)
# https://earthexplorer.usgs.gov/ 
# > Vegetation Monitoring > Phenology > eMODIS Phenology > 250m res - 2014 - Western NA V2
tin=raster::raster("D:/urban_spatial_data/PHEMUSW2014V02_TIF/TIN2014_wUSAeM250m_v2.tif")
#raster::crs(tin)

# MAX NDVI raster
# 2014MAX NDVI surface cover (Time-Integrated NDVI)
# https://earthexplorer.usgs.gov/ 
# > Vegetation Monitoring > Phenology > eMODIS Phenology > 250m res - 2014 - Western NA V2
maxn=raster::raster("D:/urban_spatial_data/PHEMUSW2014V02_TIF/MAXN2014_wUSAeM250m_v2.tif")
#raster::crs(maxn)

## --------------------------------------------------
# Data Aggregation IMPERVIOUS SURFACE

# reduce size of files by aggregating up from 30m x 30m to 300m x 300m pixels
# only need to do this once and then we have the file we are going to use for all variable extractions

# first the impervious surface cover

# try reducing extent first
# aggregate at scale x (to a coarser scale to speed the processing and save comp load)
# 1/30meters = x/500 meters -> x = 16.7
# 1/30meters = x/300 meters -> x = 10
#imp_cropped_agg <- aggregate(imp_cropped, 10, fun=mean)

# crop to smaller bounding box to make more manageable 
#extent(imp)
#new_extent <- extent(-2400000, -1100000, 950000, 3200000)
#imp_cropped <- crop(x = imp, y = new_extent)

#writeRaster(imp_agg_cropped, 
#            "D:/urban_spatial_data/nlcd_2016_impervious_l48_20210604/300m_aggregate_crop_nlcd_2016_impervious.tif",
#            overwrite=TRUE)

#rm(imp, imp_cropped, imp_agg_cropped)
#gc()

imp=raster::raster("D:/urban_spatial_data/nlcd_2016_impervious_l48_20210604/300m_aggregate_crop_nlcd_2016_impervious.tif")


## --------------------------------------------------
# Extract raster variables from sites

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

# let's just plot 100 points to get a picture of the data and shapefile
df_trans_100 <- df_trans %>%
  sample_n(100)

ggplot() +
  geom_sf(data=states) + 
  geom_sf(data = df_trans_100, 
          aes(fill = species), 
          size = 4, alpha = 0.5, shape = 23) +
  theme(legend.position="none")

## --------------------------------------------------
# Overlay spatial polygons with grid 

# create _km grid - here you can substitute by specifying grid_size above
grid <- st_make_grid(states_trans, cellsize = c(grid_size, grid_size)) %>% 
  st_sf(grid_id = 1:length(.))

# create labels for each grid_id
# grid_lab <- st_centroid(grid) %>% cbind(st_coordinates(.))

# view the grid on the polygons
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid, fill = 'transparent', lwd = 0.3) +
  # geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="none")


## --------------------------------------------------
# Extract pop density from each grid cell

pop_raster <- crop(pop_raster, states)
pop_raster <- mask(pop_raster, states)

# project the states to the raster and then crop so we cut off the ocean
crs_raster <- sf::st_crs(raster::crs(pop_raster))
prj_states <- st_transform(states_trans, crs_raster)

# project the grid to the raster
crs_raster <- sf::st_crs(raster::crs(pop_raster))
prj1 <- st_transform(grid, crs_raster)

plot(log(pop_raster+1),
     col=rev(terrain.colors(10)),
     alpha=1,
     legend=T,
     main=expression("(log) Population Density/km"^2),
     xlab="Longitude", ylab="Latitude")

plot(prj_states, colour = NA, add = TRUE)
# plot(prj1, colour = NA, add = TRUE)

# Extract raster values to list object
# this takes a while since there are many raster cells with their own values
# in each grid cell.
r.vals <- raster::extract(pop_raster, prj1)

# Use list apply to calculate mean raster value for each grid cell
r.mean <- lapply(r.vals, FUN=mean, na.rm=TRUE)

# this takes a long long time to plot
# much easier to just plot the raster directly rather than as df with ggplot
#pop_raster_df <- as.data.frame(pop_raster, xy = TRUE) %>%
#  rename("pop_dens" = "gpw_v4_population_density_rev11_2015_30_sec") %>%
#  mutate(logp1_pop_dens = log(pop_dens + 1))
#ggplot(pop_raster_df) +
#  geom_tile(aes(x=x, y=y, fill=gpw_v4_population_density_rev11_2015_30_sec))
#plot(pop_raster)

# Join mean values to the sites
grid_pop_dens <- cbind(grid, unlist(r.mean)) %>% 
  rename("pop_density_per_km2" = "unlist.r.mean.")

# now filter out the non-urban areas (areas below our pop density threshold)
# and make a scaled response variable
grid_pop_dens <- grid_pop_dens %>%
  filter(pop_density_per_km2 > min_population_size) %>%
  mutate(scaled_pop_den_km2 = center_scale(pop_density_per_km2))


## --------------------------------------------------
# Now tag records with site ID's based on spatial intersection

# which grid square is each point in?
df_id_dens <- df_trans %>% 
  st_join(grid_pop_dens, join = st_intersects) %>% as.data.frame %>%
  # filter out records from outside of the urban grid
  filter(!is.na(grid_id)) %>%
  left_join(., dplyr::select(
    df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") 


## --------------------------------------------------
# plot the spatial data

# Reproject the occurrence data if you want to plot random sample of NHC records
#df_w_dens_sf <- st_as_sf(df_id_dens,
#                         coords = c("decimalLongitude", "decimalLatitude"), 
#                         crs = 4326)

# and then transform it to the crs again
#df_w_dens_trans <- st_transform(df_w_dens_sf, crs = crs)

# create labels for each grid_id
#urban_grid_lab <- st_centroid(grid_pop_dens) %>% cbind(st_coordinates(.))

# view sampled points transposed to the grid on the polygon
#ggplot() +
#  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
#  geom_sf(data = grid_pop_dens, lwd = 0.3) +
#  geom_sf(data = sample_n(df_w_dens_trans, 100), 
#          aes(fill = species), 
#          size = 4, alpha = 0.5, shape = 23) +
  #geom_text(data = urban_grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
#  labs(x = "Longitude") +
#  labs(y = "Latitude") +
#  ggtitle("Random sample of 100 NHC records from urban site areas", 
#          subtitle = "(coloured by species)") +
# theme(legend.position = "none")

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
  mutate(area = st_area(.) %>% 
           as.numeric()) %>%
  # for each field, get area overlapping with admin area
  as_tibble() %>% 
  group_by(grid_id) %>% 
  summarize(area = sum(area)) %>%
  mutate(scaled_site_area = center_scale(area)) 


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

# project the raster to the first raster used (population density - "pop_raster")
crs_imp <- sf::st_crs(raster::crs(imp))
crs_raster <- sf::st_crs(raster::crs(pop_raster))

imp <- projectRaster(from=imp, to=pop_raster)

# make sure that the grid is still projected to the raster
crs_raster <- sf::st_crs(raster::crs(imp))
prj1 <- st_transform(grid_pop_dens, crs_raster)

# project the states to the raster and then crop so we cut off the ocean
crs_raster <- sf::st_crs(raster::crs(imp))
prj_states <- st_transform(states_trans, crs_raster)

# finally, mask the raster to the study area (prj_states)
imp <- raster::mask(imp, prj_states)

# make a plot of the sites on the state background with imp surface cover
plot(imp, 
     xlab = "Longitude", ylab = "Latitude")
plot(prj_states, colour = NA, add = TRUE)
plot(prj1, colour = NA, add = TRUE) +
  title("Impervious Surface Cover")

# then extract values and cbind with the grid_pop_dens
# extract raster values to list object
# this takes a while since there are many raster cells with their own values
# in each grid cell.
# Use list apply to calculate mean raster value for each grid cell

r.vals_imp <- raster::extract(imp, prj1)

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

## --------------------------------------------------
# Calculate tin surface cover value of each grid cell

# project the raster to the first raster used (population density - "pop_raster")
crs_tin <- sf::st_crs(raster::crs(tin))
crs_raster <- sf::st_crs(raster::crs(pop_raster))

tin <- projectRaster(from=tin, to=pop_raster)

# make sure that the grid is still projected to the raster
crs_raster <- sf::st_crs(raster::crs(tin))
prj1 <- st_transform(grid_pop_dens, crs_raster)

# project the states to the raster and then crop so we cut off the ocean
crs_raster <- sf::st_crs(raster::crs(tin))
prj_states <- st_transform(states_trans, crs_raster)

# finally, mask the raster to the study area (prj_states)
tin <- raster::crop(tin, prj_states)
tin <- raster::mask(tin, prj_states)

test <- tin

# make a plot of the sites on the state background with tin surface cover
maxValue(test)
test[test > 100] <- NA

plot(test, 
     xlab = "Longitude", ylab = "Latitude",
     legend.args = list(text = 'Index')) 
plot(prj_states, colour = NA, add = TRUE)
plot(prj1, colour = NA, add = TRUE) +
  title("\n Time-Integrated NDVI",
        adj = .5, line = 1)

rm(test)
tin[tin > 100] <- NA

# then extract values and cbind with the grid_pop_dens
# extract raster values to list object
# this takes a while since there are many raster cells with their own values
# in each grid cell.
# Use list apply to calculate mean raster value for each grid cell

r.vals_tin <- raster::extract(tin, prj1)

r.mean_tin <- lapply(r.vals_tin, FUN=mean, na.rm=TRUE)

grid_pop_dens <- cbind(grid_pop_dens, unlist(r.mean_tin)) %>%
  rename("tin_cover" = "unlist.r.mean_tin.") %>%
  mutate(scaled_tin_cover = center_scale(tin_cover))

# view sites only with scaled data
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_tin_cover), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled TIN") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Time-integrated NDVI in Urban Areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 

## --------------------------------------------------
# Calculate maxn surface cover value of each grid cell

# project the raster to the first raster used (population density - "pop_raster")
crs_maxn <- sf::st_crs(raster::crs(maxn))
crs_raster <- sf::st_crs(raster::crs(pop_raster))

maxn <- projectRaster(from=maxn, to=pop_raster)

# make sure that the grid is still projected to the raster
crs_raster <- sf::st_crs(raster::crs(maxn))
prj1 <- st_transform(grid_pop_dens, crs_raster)

# project the states to the raster and then crop so we cut off the ocean
crs_raster <- sf::st_crs(raster::crs(maxn))
prj_states <- st_transform(states_trans, crs_raster)

# finally, mask the raster to the study area (prj_states)
maxn <- raster::crop(maxn, prj_states)
maxn <- raster::mask(maxn, prj_states)

test <- maxn

# make a plot of the sites on the state background with tin surface cover
maxValue(test)
test[test > 254] <- NA

plot(test, 
     xlab = "Longitude", ylab = "Latitude",
     legend.args = list(text = 'Index')) 
plot(prj_states, colour = NA, add = TRUE)
plot(prj1, colour = NA, add = TRUE) +
  title("\n Maximum NDVI",
        adj = .5, line = 1)

rm(test)
maxn[maxn > 254] <- NA

# then extract values and cbind with the grid_pop_dens
# extract raster values to list object
# this takes a while since there are many raster cells with their own values
# in each grid cell.
# Use list apply to calculate mean raster value for each grid cell

r.vals_maxn <- raster::extract(maxn, prj1)

r.mean_maxn <- lapply(r.vals_maxn, FUN=mean, na.rm=TRUE)

grid_pop_dens <- cbind(grid_pop_dens, unlist(r.mean_maxn)) %>%
  rename("maxn" = "unlist.r.mean_maxn.") %>%
  mutate(scaled_maxn = center_scale(maxn))

# view sites only with scaled data
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_maxn), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled Max NDVI") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Maximum NDVI in Urban Areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 

## --------------------------------------------------
# Construct a correlation matrix for variables

df <- cbind(grid_pop_dens$scaled_pop_den_km2,
            grid_pop_dens$scaled_impervious_cover,
            grid_pop_dens$scaled_site_area,
            grid_pop_dens$scaled_maxn,
            grid_pop_dens$scaled_tin_cover)

(variable_correlations <- cor(df))

colnames(variable_correlations) <- c("scaled_pop_den_km2", 
                                     "scaled_impervious_cover", 
                                     "scaled_site_area",
                                     "scaled_maxn",
                                     "scaled_tin_cover")

rownames(variable_correlations) <- c("scaled_pop_den_km2", 
                                     "scaled_impervious_cover", 
                                     "scaled_site_area",
                                     "scaled_maxn",
                                     "scaled_tin_cover")
 

## --------------------------------------------------
# save RDS for easy access later

# saveRDS(as.data.frame(grid_pop_dens), "./preprocessed_data/site_data_30km_200minpop_.RDS")
# saveRDS(as.data.frame(df_id_dens), "preprocessed_data/occurrence_records_30km_200minpop_.RDS") 


# below in development
## --------------------------------------------------

test <- raster("D:/urban_spatial_data/Herbaceous_2009_2021/Herbaceous_2009_2021/rcmap_herbaceous_2009.tif")
test <- raster("D:/urban_spatial_data/Shrub_2009_2021/Shrub_2009_2021/rcmap_shrub_2009.tif")

plot(test)

# project the raster to the first raster used (population density - "pop_raster")
crs_test <- sf::st_crs(raster::crs(test))
crs_raster <- sf::st_crs(raster::crs(pop_raster))

test <- projectRaster(from=test, to=pop_raster)

# make sure that the grid is still projected to the raster
crs_raster <- sf::st_crs(raster::crs(test))
prj1 <- st_transform(grid_pop_dens, crs_raster)

# project the states to the raster and then crop so we cut off the ocean
crs_raster <- sf::st_crs(raster::crs(test))
prj_states <- st_transform(states_trans, crs_raster)

# finally, mask the raster to the study area (prj_states)
test <- raster::mask(test, prj_states)

# make a plot of the sites on the state background with test surface cover
plot(test, 
     xlab = "Longitude", ylab = "Latitude")
plot(prj_states, colour = NA, add = TRUE)
plot(prj1, colour = NA, add = TRUE) +
  title("Impervious Surface Cover")

# then extract values and cbind with the grid_pop_dens
# extract raster values to list object
# this takes a while since there are many raster cells with their own values
# in each grid cell.
# Use list apply to calculate mean raster value for each grid cell

r.vals_imp <- raster::extract(imp, prj1)

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

# do we need to pull vectors here or will they get pulled in the workflow?
## --------------------------------------------------


# and gather the site names and data that will be handled by the model
grid_id_names <- as.character(grid_pop_dens %>%
                                pull(grid_id))

scaled_pop_density <- grid_pop_dens %>%
  pull(scaled_pop_den_km2)

scaled_grid_area <- attArea %>% 
  pull(scaled_site_area)