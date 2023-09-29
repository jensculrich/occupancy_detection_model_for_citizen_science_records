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
grid_size <- 10000 # e.g., 25000 = 25km x 25 km sites
# Albers equal area projection
crs <- 5070
# minimum population size
min_population_size <- 1200 
# minimum site area 
min_site_area = 0.25

#taxon = "syrphidae"
taxon = "bombus"

## --------------------------------------------------
# Spatial extent and urban areas

# spatial data - state shapefile
states <- tigris::states() %>%
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
  theme(legend.position="none") + 
  ggtitle("")


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

# plot the raster (takes a long time to do this)
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

gc()

## --------------------------------------------------
# Join with Metro areas

# Metropolitan statistical areas
# https://catalog.data.gov/dataset/tiger-line-shapefile-2019-nation-u-s-current-metropolitan-statistical-area-micropolitan-statist
# updated 2018, based on 2010 census data
CBSA <- sf::read_sf("./data/spatial_data/tl_2019_us_cbsa/tl_2019_us_cbsa.shp")

## level 3 cluster (ecoregion3)
crs_CBSA <- sf::st_crs(raster::crs(CBSA))
prj1 <- st_transform(grid_pop_dens, crs_CBSA)

CBSA_names <- st_join(prj1, CBSA) %>%
  group_by(grid_id) %>%
  slice(which.max(ALAND)) %>% 
  pull(NAME) 

CBSA_lookup <- as.numeric(as.factor(CBSA_names))

n_CBSA <- unique(CBSA_lookup) %>%
  length()

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

# Sacramento county, california
filtered <- income_trans %>%
  filter(STATEFP == "06" & COUNTYFP == "067")


# Plot income by census block
income1 <- ggplot() + geom_sf(data = filtered, aes(fill = relative_AMR8E001)) +
  scale_fill_gradient2(high = ("goldenrod"), 
                       name="Relative median income",
                       limits=c(0,3)) +
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  ggtitle("Median household income in Sacramento County, California,\nrelative to median household income in California \n(by census block group) ")
extent(filtered)
income1 <- income1 +
  geom_sf(data = grid_pop_dens, fill = NA, color = "violetred4") +
  coord_sf(xlim = c(-2219577 , -2137015 ), ylim = c(1966995 , 2040501 ), expand = FALSE)

# Cook county, Illinois
filtered2 <- income_trans %>%
  filter(STATEFP == "17" & COUNTYFP == "031")

# Plot income by census block
income2 <- ggplot() + geom_sf(data = filtered2, aes(fill = relative_AMR8E001)) +
  scale_fill_gradient2(high = ("goldenrod"), 
                       name="Relative median income",
                       limits=c(0,3.5)) +
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  ggtitle("Median household income in Cook County, Illinois,\nrelative to median household income in Illinois \n(by census block group) ")
extent(filtered2)
income2 <- income2 +
  geom_sf(data = grid_pop_dens, fill = NA, color = "violetred4") +
  coord_sf(xlim = c(635130.2  ,701503.4  ), ylim = c(2080966  , 2157258 ), expand = FALSE)


# NA's: "It's common for data to be suppressed when there isn't a sufficient sample size. 
# The suppression takes place to protect information about the respondents in those areas
# and to limit unreliable statistics." Looking at those NA areas, they are mostly areas of
# natural habitat, farmland, unused open land, and airports

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

#land_cropped <- crop(land, prj1)
#land_masked <- mask(land_cropped, prj1)

# for plotting make some inset area
extent(land)
e <- c(-2060000, -1970000, 1380000, 1470000)
land_cropped_inset <- crop(land, e)

# can we mask it to the grid cells first?
# try cropping just grid cells and then plotting those to get the right area first?

# natural
temp <- land_cropped_inset
temp[temp != c( 41, 42, 43, 52, 71, 90, 95)] <- NA
temp <- as.data.frame(temp, xy=TRUE)
temp2 <- filter(temp, !is.na(NLCD.Land.Cover.Class_NLCD.Land.Cover.Class))

# open dev
temp <- land_cropped_inset
temp[temp != c( 21)] <- NA
temp <- as.data.frame(temp, xy=TRUE)
temp2 <- filter(temp, !is.na(NLCD.Land.Cover.Class_NLCD.Land.Cover.Class))

#ggplot(temp) +
 # geom_tile(aes(x=x, y=y, fill=factor(NLCD.Land.Cover.Class_NLCD.Land.Cover.Class),alpha=0.8)) 

# project the states to the raster
crs_raster <- sf::st_crs(raster::crs(land))
prj_states <- st_transform(states_trans, crs_raster)
prj1 <- st_transform(grid_pop_dens, crs_raster)

ggplot() +
  geom_sf(data = prj_states, fill = 'white', lwd = 0.05) +
  geom_sf(data = prj1, fill = NA, lwd = 0.3) +
  coord_sf(
    xlim = c(-2060000 , -1970000), 
    ylim = c(1380000, 1470000),
    expand = FALSE) +
  # is in the wrong proj?
  geom_tile(data = temp2, aes(x=x, y=y, fill=factor(NLCD.Land.Cover.Class_NLCD.Land.Cover.Class))) +
  guides(fill=guide_legend(title="Natural Habitat")) +
  xlab("Longitude") + 
  ylab("Latitude")
  
#rasterVis::levelplot(temp, xlab = "Longitude", ylab = "Latitude", 
 #                    colorkey=FALSE) 

#plot(temp, 
#     xlab = "Longitude", ylab = "Latitude",
#     col = "darkblue")
#plot(prj_states, colour = NA, add = TRUE)
#plot(prj_grid_lab, add = TRUE)
#plot(prj1, colour = NA, add = TRUE)
#  title("Land Use",
#        adj = .5, line = 1)

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
                                     (length(which(x %in% c(52,71,41,42,43, 90, 95))) / length(x))
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

grid_pop_dens <- grid_pop_dens %>%
  filter(scaled_pop_den_km2 < 5) %>%
  filter(scaled_pop_den_km2 > -5) %>%
  # rescale the variables (in case any sites removed)
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

library(cowplot)

# view extent only with population density
main_map <- ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = pop_density_per_km2), lwd = 0.3) +
  scale_fill_gradient2(name = expression("population/km"^2),
                       low="white", high="dodgerblue3") +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  #ggtitle("Population density in urban areas") +
  #labs(x = "Longitude") +
  #labs(y = "Latitude") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    # legend.justification defines the edge of the legend that the legend.position coordinates refer to
    legend.margin = margin(0, 12, 12, 12),
    legend.justification = c("right", "bottom"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.position = c(1, 0.01),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.background = element_rect(
      fill="white",
      size=1, linetype="solid", 
      colour ="black")
  ) 

# get lims
#layer_scales(main_population_map)$x$get_limits()
#layer_scales(main_population_map)$y$get_limits()

# southwestern inset
main_map +
  coord_sf(
    xlim = c(-2150000 , -1350000), 
    ylim = c(1000000, 1700000),
    expand = FALSE
  )


main_map <- main_map +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

ggdraw(main_map) +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(-2150000 , -1350000), 
          ylim = c(1000000, 1700000),
          expand = FALSE) +
        theme(legend.position = "none",
              plot.background=element_rect(fill="white", color="white"))
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = -0.075, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.01,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.55, 
    height = 0.55)

##
# view sites only with scaled population density
main_map <- ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_pop_den_km2), lwd = 0.3) +
  scale_fill_gradient2(name = expression("scaled pop. density"),
                       low="firebrick3", high="dodgerblue3") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    # legend.justification defines the edge of the legend that the legend.position coordinates refer to
    #legend.box.margin = margin(50,50,50,50),
    legend.margin = margin(0, 12, 12, 12),
    legend.justification = c("right", "bottom"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.position = c(1, 0.01),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.background = element_rect(
      fill="white",
      size=1, linetype="solid", 
      colour ="black")
  ) 

# get lims
#layer_scales(main_population_map)$x$get_limits()
#layer_scales(main_population_map)$y$get_limits()

main_map <- main_map +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

ggdraw(main_map) +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(-2150000 , -1350000), 
          ylim = c(1000000, 1700000),
          expand = FALSE) +
        theme(legend.position = "none",
              plot.background=element_rect(fill="white", color="white"))
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = -0.075, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.01,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.55, 
    height = 0.55)


##
# scaled site area data
main_map <- ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_site_area), lwd = 0.3) +
  scale_fill_gradient2(name = expression("scaled site area"),
                       low="firebrick3", high="dodgerblue3") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    # legend.justification defines the edge of the legend that the legend.position coordinates refer to
    #legend.box.margin = margin(50,50,50,50),
    legend.margin = margin(0, 12, 12, 12),
    legend.justification = c("right", "bottom"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.position = c(1, 0.01),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.background = element_rect(
      fill="white",
      size=1, linetype="solid", 
      colour ="black")
  ) 

# get lims
#layer_scales(main_population_map)$x$get_limits()
#layer_scales(main_population_map)$y$get_limits()

main_map <- main_map +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

ggdraw(main_map) +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(-2150000 , -1350000), 
          ylim = c(1000000, 1700000),
          expand = FALSE) +
        theme(legend.position = "none",
              plot.background=element_rect(fill="white", color="white"))
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = -0.075, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.01,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.55, 
    height = 0.55)

##
# scaled income
main_map <- ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_avg_income), lwd = 0.3) +
  scale_fill_gradient2(name = expression("scaled relative income"),
                       low="firebrick3", high="dodgerblue3") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    # legend.justification defines the edge of the legend that the legend.position coordinates refer to
    #legend.box.margin = margin(50,50,50,50),
    legend.margin = margin(0, 12, 12, 12),
    legend.justification = c("right", "bottom"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.position = c(1, 0.01),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.background = element_rect(
      fill="white",
      size=1, linetype="solid", 
      colour ="black")
  ) 

# get lims
#layer_scales(main_population_map)$x$get_limits()
#layer_scales(main_population_map)$y$get_limits()

main_map <- main_map +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

ggdraw(main_map) +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(-2150000 , -1350000), 
          ylim = c(1000000, 1700000),
          expand = FALSE) +
        theme(legend.position = "none",
              plot.background=element_rect(fill="white", color="white"))
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = -0.075, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.01,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.55, 
    height = 0.55)


##
# scaled natural habitat area
main_map <- ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, aes(fill = scaled_herb_shrub_forest), lwd = 0.3) +
  scale_fill_gradient2(name = expression("scaled natural habitat area"),
                       low="firebrick3", high="dodgerblue3") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    # legend.justification defines the edge of the legend that the legend.position coordinates refer to
    #legend.box.margin = margin(50,50,50,50),
    legend.margin = margin(0, 12, 12, 12),
    legend.justification = c("right", "bottom"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.position = c(1, 0.01),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.background = element_rect(
      fill="white",
      size=1, linetype="solid", 
      colour ="black")
  ) 

# get lims
#layer_scales(main_population_map)$x$get_limits()
#layer_scales(main_population_map)$y$get_limits()

main_map <- main_map +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

ggdraw(main_map) +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(-2150000 , -1350000), 
          ylim = c(1000000, 1700000),
          expand = FALSE) +
        theme(legend.position = "none",
              plot.background=element_rect(fill="white", color="white"))
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = -0.075, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.01,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.55, 
    height = 0.55)



##
# solid map with insets for southwest for both income and natural habitat
main_map <- ggplot() +
  geom_sf(data = states_trans, fill = 'grey95', lwd = 0.05) +
  geom_sf(data = grid_pop_dens, 
          #aes(fill = scaled_herb_shrub_forest), 
          lwd = 0.3,
          fill = "grey20",
          alpha = 0.8
  ) +
  #scale_fill_gradient2(name = expression("scaled natural habitat area"),
  #                     low="firebrick3", high="dodgerblue3") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    # legend.justification defines the edge of the legend that the legend.position coordinates refer to
    #legend.box.margin = margin(50,50,50,50),
    legend.margin = margin(0, 12, 12, 12),
    legend.justification = c("left", "bottom"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.position = c(0, 0),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.background = element_rect(
      fill="white",
      size=1, linetype="solid", 
      colour ="black")
  ) 

cities <- c("Los Angeles", "Phoenix", "Las Vegas", "Tuscon", "San Diego")
longitude <- c("-118.241736", "-112.072499", "-115.141447", "-110.978051", "-117.15")
latitude <- c("34.054091", "33.448200", "36.170334", "32.250236", "32.713097")

city_df <- as.data.frame(cbind(cities, longitude, latitude))
city_df <- st_as_sf(city_df, coords = c(2:3), crs = 4326)

crs <- sf::st_crs(raster::crs(states_trans))
prj_city_df <- st_transform(city_df, crs)

cities_labels <- c("Los Angeles", "Phoenix", "Las Vegas", "Tuscon", "San Diego")
longitude_labels <- c("-117.1", "-112.9", "-113.95", "-111.72", "-118.15")
latitude_labels <- c("34.575", "33.1", "36.0", "31.9", "32.25")

city_df_labels <- as.data.frame(cbind(cities_labels, longitude_labels, latitude_labels))
city_df_labels <- st_as_sf(city_df_labels, coords = c(2:3), crs = 4326)

crs <- sf::st_crs(raster::crs(states_trans))
prj_city_df_labels <- st_transform(city_df_labels, crs)

# southwestern inset
inset1 <- main_map +
  geom_sf(data = grid_pop_dens, 
          aes(fill = scaled_herb_shrub_forest), 
          lwd = 0.3,
          #fill = "grey20",
          alpha = 0.8
  ) +
  geom_sf(data = prj_city_df, size = 8, 
          shape = 1, fill = "black", alpha = 0.9) +
  geom_sf_label(data = prj_city_df_labels, aes(label = cities), size = 7) +
  scale_fill_gradient2(name = expression("scaled natural habitat area"),
                       low="firebrick3", high="dodgerblue3") +
  coord_sf(
    xlim = c(-2150000 , -1350000), 
    ylim = c(1000000, 1700000),
    expand = FALSE
  ) +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

inset2 <- main_map +
  geom_sf(data = grid_pop_dens, 
          aes(fill = scaled_avg_income), 
          lwd = 0.3,
          #fill = "grey20",
          alpha = 0.8
  ) +
  geom_sf(data = prj_city_df, size = 8, 
          shape = 1, fill = "black", alpha = 0.9) +
  geom_sf_label(data = prj_city_df_labels, aes(label = cities), size = 7) +
  scale_fill_gradient2(name = expression("scaled relative income"),
                       low="firebrick3", high="dodgerblue3") +
  coord_sf(
    xlim = c(-2150000 , -1350000), 
    ylim = c(1000000, 1700000),
    expand = FALSE
  ) +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

inset3 <- main_map +
  geom_sf(data = grid_pop_dens, 
          aes(fill = scaled_developed_open), 
          lwd = 0.3,
          #fill = "grey20",
          alpha = 0.8
  ) +
  geom_sf(data = prj_city_df, size = 8, 
          shape = 1, fill = "black", alpha = 0.9) +
  geom_sf_label(data = prj_city_df_labels, aes(label = cities), size = 7) +
  scale_fill_gradient2(name = expression("scaled developed greenspace"),
                       low="firebrick3", high="dodgerblue3") +
  coord_sf(
    xlim = c(-2150000 , -1350000), 
    ylim = c(1000000, 1700000),
    expand = FALSE
  ) +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )


main_map <- main_map +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

plot_grid(inset1, NULL, inset2, 
          #labels = c('b)', "", 'c)'),
          #label_size = 20,
          nrow=1,
          rel_widths = c(1,0.05,1))

plot_grid(inset1, NULL, inset3, 
          #labels = c('b)', "", 'c)'),
          #label_size = 20,
          nrow=1,
          rel_widths = c(1,0.05,1))

main_map

# get lims
#layer_scales(main_population_map)$x$get_limits()
#layer_scales(main_population_map)$y$get_limits()

main_map <- main_map +
  geom_rect(aes(
    xmin = -2150000,
    ymin = 1000000,
    xmax = -1350000,
    ymax = 1700000),
    fill = NA, 
    colour = "black",
    size = 1
  )

ggdraw(main_map) +
  draw_plot(
    {
      main_map + 
        coord_sf(
          xlim = c(-2150000 , -1350000), 
          ylim = c(1000000, 1700000),
          expand = FALSE) +
        theme(legend.position = "none",
              plot.background=element_rect(fill="white", color="white"))
    },
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = -0.075, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.01,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.55, 
    height = 0.55)