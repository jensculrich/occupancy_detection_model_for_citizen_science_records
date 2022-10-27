### Gather and prep spatial data for pollinator occurrence
# jcu; started oct 27, 2022

# Collect spatial data to describe where the observations occur i.e.
# We will construct a grid of cells across the spatial extent (California)
# Each cell will be a site, where observations of any species can
# be said to have occurred in any of the years of our timeline.
# Choosing the spatial size of the grid cells will be important,
# and also importantly we will plan to limit the grid cells to those 
# that contain 'urban areas'.

library(tidyverse)
library(tigris) # get state shapefile
library(sf) # spatial data processing

## --------------------------------------------------
# Get Data

CA <- tigris::states() %>%
  filter(NAME == "California")
str(CA)
st_crs(CA)

# occurrence data
df <- read.csv("./data/data_filtered.csv")

## --------------------------------------------------
# Overlay spatial polygon with grid 

crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"

# let's just plot 100 points to get a picture of the data and shapefile
df_100 <- df %>%
  sample_n(100)

ggplot() +
  geom_sf(data=CA) + 
  geom_point(data = df_100, 
             aes(x = decimalLongitude, y = decimalLatitude, fill = species), 
             size = 4, alpha = 0.5, shape = 23) +
  theme(legend.position="none")


## transform
CA_trans <- st_transform(CA, 26910) # NAD83 / UTM Zone 10N

(df_sf <- st_as_sf(df,
                     coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = 4326)
)

df_trans <- st_transform(df_sf, crs = crs)

# create 50km grid - here you can substitute 200 for 50000
grid_25 <- st_make_grid(CA_trans, cellsize = c(25000, 25000)) %>% 
  st_sf(grid_id = 1:length(.))

# create labels for each grid_id
grid_lab <- st_centroid(grid_25) %>% cbind(st_coordinates(.))

# view sampled points transposed to the grid on the polygon
ggplot() +
  geom_sf(data = CA_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = sample_n(df_trans, 500), 
             aes(fill = species), 
             size = 4, alpha = 0.5, shape = 23) +
  geom_sf(data = grid_25, fill = 'transparent', lwd = 0.3) +
  geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="none")

# which grid square is each point in?
df_id <- df_trans %>% st_join(grid_25, join = st_intersects) %>% as.data.frame
