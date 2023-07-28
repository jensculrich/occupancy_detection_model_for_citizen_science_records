library(tidyverse) # data carpentry
library(tigris) # get state shapefile
library(sf) # spatial data processing

# USA_Contiguous_Albers_Equal_Area_Conic
crs <- 5070

# spatial data - state shapefile
states <- tigris::states() %>%
  # lower 48 + DC
  filter(REGION != 9) %>%
  filter(!NAME %in% c("Alaska", "Hawaii"))

#str(states)
#st_crs(states)

states_trans <- states  %>% 
  st_transform(., crs) # USA_Contiguous_Albers_Equal_Area_Conic

level_three_cluster <- sf::read_sf("./data/spatial_data/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp")

level_three_cluster <- level_three_cluster %>% 
  st_transform(., crs)

level_three_crop <- st_crop(level_three_cluster, states_trans)
level_three_intersection <- st_intersection(level_three_crop, states_trans)

# view extent only with population density
main_map <- ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = level_three_intersection, aes(fill = NA_L3NAME), alpha=0.5, lwd = 0.3) +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  #ggtitle("Population density in urban areas") +
  #labs(x = "Longitude") +
  #labs(y = "Latitude") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none")

main_map

level_one_cluster <- sf::read_sf("./data/spatial_data/na_cec_eco_l1/NA_CEC_ECO_Level1.shp")

level_one_cluster <- level_one_cluster %>% 
  st_transform(., crs)

level_one_crop <- st_crop(level_one_cluster, states_trans)
level_one_intersection <- st_intersection(level_one_crop, states_trans)

# view extent only with population density
main_map2 <- ggplot() +
  geom_sf(data = states_trans, fill = 'white', lwd = 0.05) +
  geom_sf(data = level_one_intersection, aes(fill = NA_L1NAME), alpha=0.5, lwd = 0.3) +
  #geom_text(data = urban_grid_lab, 
  #          aes(x = X, y = Y, label = grid_id), size = 2) +
  #ggtitle("Population density in urban areas") +
  #labs(x = "Longitude") +
  #labs(y = "Latitude") +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none")

main_map2


