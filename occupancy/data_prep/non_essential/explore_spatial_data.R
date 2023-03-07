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
library(exactextractr) # quick extraction of raster data

## --------------------------------------------------
# Operation Functions
## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}


## --------------------------------------------------
# Study design

## global options
# grid size
grid_size <- 15000 # e.g., 25000 = 25km x 25 km sites
# CRS for NAD83 / UTM Zone 10N
# crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"
# Albers equal area
crs <- 5070
# minimum population size
# a grid cell must intersect a city with min_population_size or more
# people living there. This filters out grid cell sites that do not
# have a an urban center above some minimum threshold

# min_population_size of 38 (/km^2) is ~ 100/mile^2 which is a typical threshold for 
# considering an area to be 'urban'
# let's up the minimum a bit and go with 100 per sq km, which is about 260/sq mile
min_population_size <- 1200 

# minimum site area 
# if sites are super tiny, the observation process could likely be very unstable
min_site_area = 0.10

#taxon = "syrphidae"
taxon = "bombus"


## --------------------------------------------------
# Spatial extent and urban areas

# spatial data - state shapefile
states <- tigris::states() %>%
  #filter(NAME %in% c("California", "Oregon", "Washington", "Arizona", "Nevada"))
  # lower 48 + DC
  filter(REGION != 9) %>%
  filter(!NAME %in% c("Alaska", "Hawaii"))
  
str(states)
st_crs(states)

## --------------------------------------------------
# Environmental rasters 

# human pop density raster
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
# 2015 pop density at 1km resolution
pop_raster=raster("./data/spatial_data/population_density/gpw_v4_population_density_rev11_2015_30_sec.tif")
crs(pop_raster)

## --------------------------------------------------
# Data Aggregation 

# impervious surface

# reduce size of files by aggregating up from 30m x 30m to 300m x 300m pixels
# only need to do this once and then we have the file we are going to use for all variable extractions

# first the impervious surface cover

# try reducing extent first
# aggregate at scale x (to a coarser scale to speed the processing and save comp load)
# 1/30meters = x/500 meters -> x = 16.7
# 1/30meters = x/300 meters -> x = 10
#imp_cropped_agg <- aggregate(imp_cropped, 10, fun=mean)

# land use
# project the states to the raster and then crop so we cut off the ocean
#crs_raster <- sf::st_crs(raster::crs(land))
#prj_states <- st_transform(states, crs_raster)

# crop to smaller bounding box to make more manageable 
#extent(land)
#new_extent <- extent(-2400000, -1800000, 1250000, 1800000)
#land_cropped <- crop(x = land, y = new_extent)

#extent(land)
#new_extent <- extent(-2075000, -1925000, 1400000, 1525000)
#land_cropped_zoom <- crop(x = land_cropped, y = new_extent)

# finally, mask the raster to the study area (prj_states)
#land <- raster::mask(land, prj_states)

#writeRaster(land_cropped_zoom, "./data/spatial_data/land_use/land_use_full_southCA_zoom.tif", overwrite=TRUE)

#rm(land, land_cropped)
#gc()

#land2 <- raster::raster("./data/spatial_data/land_use/land_use_full.gri")

#land=raster::raster("./data/spatial_data/land_use/land_use.tif")


## --------------------------------------------------
# Extract raster variables from sites

## --------------------------------------------------
# Overlay spatial polygons with grid 

## transform
# transform state shapefile to crs
states_trans <- st_transform(states, crs) # albers equal area

# create _km grid - here you can substitute by specifying grid_size above
grid <- st_make_grid(states_trans, cellsize = c(grid_size, grid_size)) %>% 
  st_sf(grid_id = 1:length(.))

# create labels for each grid_id
grid_lab <- st_centroid(grid) %>% cbind(st_coordinates(.))

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

#plot(log(pop_raster+1),
#     col=rev(terrain.colors(10)),
#     alpha=1,
#     legend=T,
#     main=expression("(log) Population Density/km"^2),
#     xlab="Longitude", ylab="Latitude",)

#plot(prj_states, colour = NA, add = TRUE)
# plot(prj1, colour = NA, add = TRUE)

# Extract raster values to list object and then summarize by the mean value
r.vals <- exactextractr::exact_extract(pop_raster, prj1, 'mean')

## --------------------------------------------------
# filter sites to urban areas using the r.vals for pop density

grid_pop_dens <- cbind(grid, r.vals) %>% 
  rename("pop_density_per_km2" = "r.vals")
# now filter out the non-urban areas (areas below our pop density threshold)
# and make a scaled response variable
grid_pop_dens <- grid_pop_dens %>%
  filter(pop_density_per_km2 > min_population_size) 

## --------------------------------------------------
# Try extracting median household income from each block

# census block income data
# 2020 income data 
# B19013. Median Household Income in the Past 12 Months (in 2020 Inflation-Adjusted Dollars)
#https://data2.nhgis.org/main
# Typically, Block Groups have a population of 600 to 3,000 people.
income=st_read("./data/spatial_data/socioeconomic_data/US_blck_grp_2020.shp")
gc()
income_data <- read.csv("./data/spatial_data/socioeconomic_data/nhgis0001_ds249_20205_blck_grp.csv")
income_data <- income_data %>%
  dplyr::select(GISJOIN, AMR8E001)

income <- left_join(income, income_data, by="GISJOIN")
rm(income_data)
gc()

# transform state shapefile to crs
income_trans <- st_transform(income, crs) # albers equal area

grid_pop_dens <- st_join(grid_pop_dens, income_trans)
grid_pop_dens <- grid_pop_dens %>%
  group_by(grid_id) %>%
  mutate(avg_income = mean(AMR8E001, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(grid_id, pop_density_per_km2, avg_income) %>%
  filter(!is.na(avg_income)) %>%
  mutate(scaled_avg_income = center_scale(avg_income))

rm(income, income_trans)
gc()

## --------------------------------------------------
# Extract environmental variables from each remaining site

# MRLC data https://www.mrlc.gov/data

# DO NOT READ IF USING THE AGGREGATED/CROPPED FILE PRODUCED ONE TIME ONLY BELOW
# land cover raster 30m x 30m
# 2016 land cover data https://www.mrlc.gov/data/nlcd-2016-land-cover-conus
land=raster::raster("./data/spatial_data/land_cover/land_cover/nlcd_2016_land_cover_l48_20210604.img")
raster::crs(land)

# make sure that the grid is still projected to the raster
crs_raster <- sf::st_crs(raster::crs(land))
prj1 <- st_transform(grid_pop_dens, crs_raster)

# project the states to the raster
crs_raster <- sf::st_crs(raster::crs(land))
prj_states <- st_transform(states_trans, crs_raster)

# project the grid labels
prj_grid_lab <- st_transform(grid_lab, crs_raster)

# make a plot of the sites on the state background with land cover
# 11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95

# natural habitat
#temp <- reclassify(land_cropped_zoom, cbind(c(41, 42, 43, 52, 71), 1))
#temp <- reclassify(temp, cbind(c(0, 11, 12, 21, 22, 23, 24, 31, 51, 72, 73, 74, 81, 82, 90, 95), NA))

# open developed habitat
#temp <- reclassify(land_cropped_zoom, cbind(c(21), 1))
#temp <- reclassify(temp, cbind(c(0, 11, 12, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95), NA))

# open dev
#temp <- land_cropped_zoom
#temp[temp != 21] <- NA

# natural
#temp <- land_cropped_zoom
#temp[temp != c( 41, 42, 43, 52, 71)] <- NA

#rasterVis::levelplot(temp, xlab = "Longitude", ylab = "Latitude", 
#                     colorkey=FALSE) 

plot(temp, 
     xlab = "Longitude", ylab = "Latitude",
     col = "darkblue")
plot(prj_states, colour = NA, add = TRUE)
plot(prj_grid_lab, add = TRUE)
plot(prj1, colour = NA, add = TRUE) +
  title("Land Use",
        adj = .5, line = 1)

# then extract values and cbind with the grid_pop_dens
# extract raster values to list object
# make sure that the grid is still projected to the raster
crs_raster <- sf::st_crs(raster::crs(land))
prj1 <- st_transform(grid_pop_dens, crs_raster)

# project the states to the raster
crs_raster <- sf::st_crs(raster::crs(land))
prj_states <- st_transform(states_trans, crs_raster)


r.vals_land <- exactextractr::exact_extract(land, prj1)
gc()

# want to reclassify open water as NA
r.vals_land_NA <- lapply(r.vals_land, function(x) na_if(x$value, 0))
r.vals_land_NA <- lapply(r.vals_land_NA, function(x) na_if(x,11))
gc()

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
gc()

# now drop NA values so that the below estimates are the proportion of cover
# of all land cover in the administrative area
r.vals_land_NA <- lapply(r.vals_land_NA, na.omit)
gc()

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
                                     (length(which(x %in% c(52,71,41,42,43))) / length(x))
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
   r.site_area, r.vals_land, r.vals_land_NA, states)
gc()

## --------------------------------------------------
# Visualize the data

# view extent only with population density
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  #geom_sf(data = grid_pop_dens, aes(fill = pop_density_per_km2), lwd = 0.3) +
  #scale_fill_gradient2(name = expression("Population/km"^2)) +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  #ggtitle("Population density in urban areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 


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

# scaled site area data
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

# view sites only with scaled income
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_avg_income), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled median household income") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Avg. household income in urban areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 

# scaled shrub and herb cover
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_herb_shrub_cover), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled herb/shrub cover") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Undeveloped Herb and Shrub Cover in Urban Areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 

# scaled open developed area
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_developed_open), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled open developed area") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Open Developed Land in Urban Areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 

# scaled medium/high developed
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_developed_med_high), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled med/high developed cover") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Med. and High Developed Land in Urban Areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 

# sclaed forest cover
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_forest), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled forest cover") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Forest Cover in Urban Areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 


# sclaed herb shrub forest cover
ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_herb_shrub_forest), lwd = 0.3) +
  scale_fill_gradient2(name = "Scaled cover") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  ggtitle("Undeveloped natural habitat in Urban Areas") +
  labs(x = "Longitude") +
  labs(y = "Latitude") 

## --------------------------------------------------
# Construct a correlation matrix for variables

df <- cbind(grid_pop_dens$scaled_pop_den_km2,
            grid_pop_dens$scaled_site_area,
            grid_pop_dens$scaled_developed_open,
            grid_pop_dens$scaled_developed_med_high,
            grid_pop_dens$scaled_herb_shrub_cover,
            grid_pop_dens$scaled_forest,
            grid_pop_dens$scaled_avg_income)

(variable_correlations <- cor(df))

colnames(variable_correlations) <- c("scaled_pop_den_km2", 
                                     "scaled_site_area", 
                                     "scaled_developed_open",
                                     "scaled_developed_med_high",
                                     "scaled_herb_shrub_cover",
                                     "scaled_forest",
                                     "scaled_avg_income")

rownames(variable_correlations) <- c("scaled_pop_den_km2", 
                                     "scaled_site_area", 
                                     "scaled_developed_open",
                                     "scaled_developed_med_high",
                                     "scaled_herb_shrub_cover",
                                     "scaled_forest",
                                     "scaled_avg_income")

View(as.data.frame(variable_correlations))

## --------------------------------------------------
# Occurrence data

# read either the syrphidae data or the bombus data
df <- read.csv(paste0("./data/occurrence_data/", taxon, "_data.csv"))

## --------------------------------------------------
# Prep the data

# make the df into a spatial file
(df_sf <- st_as_sf(df,
                   coords = c("decimalLongitude", "decimalLatitude"), 
                   crs = 4326))

# and then transform it to the crs
df_trans <- st_transform(df_sf, crs = crs)


## --------------------------------------------------
# plot the spatial data

# Reproject the occurrence data if you want to plot random sample of NHC records
#df_w_dens_sf <- st_as_sf(df_id_dens,
#                         coords = c("decimalLongitude", "decimalLatitude"), 
#                         crs = 4326)

# and then transform it to the crs again
#df_w_dens_trans <- st_transform(df_w_dens_sf, crs = crs)

# create labels for each grid_id
# urban_grid_lab <- st_centroid(grid_pop_dens) %>% cbind(st_coordinates(.))

# view sampled points transposed to the grid on the polygon
#ggplot() +
#  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
#  geom_sf(data = grid_pop_dens, lwd = 0.3) +
#  geom_sf(data = sample_n(df_w_dens_trans, 500), 
#          aes(fill = species), 
#          size = 4, alpha = 0.5, shape = 23) +
#geom_text(data = urban_grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
#  labs(x = "Longitude") +
#  labs(y = "Latitude") +
#  ggtitle("Random sample of 500 Syrphidae NHC records from urban site areas", 
#          subtitle = "(coloured by species)") +
# theme(legend.position = "none")

## --------------------------------------------------
# Now tag records with site ID's based on spatial intersection

# which grid square is each point in?
df_id_dens <- df_trans %>% 
  st_join(grid_pop_dens, join = st_intersects) %>% as.data.frame %>%
  # filter out records from outside of the urban grid
  filter(!is.na(grid_id)) %>%
  left_join(., dplyr::select(
    df, id, decimalLatitude, decimalLongitude), by="gbifID") 


## --------------------------------------------------
# save RDS for easy access later

#saveRDS(as.data.frame(grid_pop_dens), "./preprocessed_data/site_data_30km_200minpop_.RDS")

#saveRDS(as.data.frame(df_id_dens), paste0(
#  "preprocessed_data/",
#  taxon,
#  "_occurrence_records_30km_200minpop_.RDS")) 

# end file

## --------------------------------------------------
# Try quantifying landscape configuration
library(landscapemetrics)

check_landscape(land)

# make sure that the grid is still projected to the raster
crs_raster <- sf::st_crs(raster::crs(land))
prj1 <- st_transform(grid_pop_dens, crs_raster)

# natural habitat
temp <- reclassify(land, cbind(c(41, 42, 43, 52, 71), 1))
temp <- reclassify(temp, cbind(c(0, 11, 12, 21, 22, 23, 24, 31, 51, 72, 73, 74, 81, 82, 90, 95), NA))

plot(temp)

my_metrics = sample_lsm(temp, prj1, 
                        level = "landscape", metric = c("enn_mn"))


## --------------------------------------------------
# Connect sites for spatial autocorrelation smoothing
# https://github.com/ConnorDonegan/Stan-IAR

# cite: Donegan, Connor. Flexible Functions for ICAR, BYM, and BYM2 Models in Stan. 
# Code Repository. 2021. 
# Available online: https://github.com/ConnorDonegan/Stan-IAR (access date).

# should do this AFTER removing any small sites

C <- spdep::nb2mat(spdep::poly2nb(grid_pop_dens, queen = TRUE), style = "B", zero.policy = TRUE)

edges <- function (w) {
  lw <- apply(w, 1, function(r) {
    which(r != 0)
  })
  all.edges <- lapply(1:length(lw), function(i) {
    nbs <- lw[[i]]
    if (length(nbs)) 
      data.frame(node1 = i, node2 = nbs, weight = w[i, nbs])
  })
  all.edges <- do.call("rbind", all.edges)
  edges <- all.edges[which(all.edges$node1 < all.edges$node2), ]
  return(edges)
}

prep_icar_data <- function (C, inv_sqrt_scale_factor = NULL) {
  n <- nrow(C)
  E <- edges(C)
  G <- list(np = nrow(C), from = E$node1, to = E$node2, nedges = nrow(E))
  class(G) <- "Graph"
  nb2 <- spdep::n.comp.nb(spdep::graph2nb(G))
  k = nb2$nc
  if (inherits(inv_sqrt_scale_factor, "NULL")) inv_sqrt_scale_factor <- array(rep(1, k), dim = k)
  group_idx = NULL
  for (j in 1:k) group_idx <- c(group_idx, which(nb2$comp.id == j))
  group_size <- NULL
  for (j in 1:k) group_size <- c(group_size, sum(nb2$comp.id == j))
  # intercept per connected component of size > 1, if multiple.
  m <- sum(group_size > 1) - 1
  if (m) {
    GS <- group_size
    ID <- nb2$comp.id
    change.to.one <- which(GS == 1)
    ID[which(ID == change.to.one)] <- 1
    A = model.matrix(~ factor(ID))
    A <- as.matrix(A[,-1])
  } else {
    A <- model.matrix(~ 0, data.frame(C))
  }
  l <- list(k = k, 
            group_size = array(group_size, dim = k), 
            n_edges = nrow(E), 
            node1 = E$node1, 
            node2 = E$node2, 
            group_idx = array(group_idx, dim = n), 
            m = m,
            A = A,
            inv_sqrt_scale_factor = inv_sqrt_scale_factor, 
            comp_id = nb2$comp.id)
  return(l)
}

icar.data <- prep_icar_data(C)