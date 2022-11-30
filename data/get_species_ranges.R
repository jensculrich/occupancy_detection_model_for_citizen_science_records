### Generate species ranges
# jcu; started nov 29, 2022

# Use data from historical to present day to determine the range of each species
# Then, identify which grid cells are intersecting the range for each species.
# We will pass back the range data to prep_data_x.R, and then sub in 0's
# in the binary indicator array for whether or not a species was sampled at a site
# in a given time period, so that max obs = 0 for a site outside the species range.

library(tidyverse)

get_species_ranges <- function(
  site_name_vector,
  species_vector,
  min_year_for_species_ranges
){
 
  # read occurrence data
  # this is occurrence data from all time records, not just from the study time span 
  df <- read.csv("./data/data_unfiltered.csv") 
  
  species_name <- species_vector[15]
  
  test <- df %>%
    filter(species == species_name,
           year > min_year_for_species_ranges) %>%
    dplyr::select(decimalLatitude, decimalLongitude)
    
    ch <- chull(test)
    
    coords <- test[c(ch, ch[1]), ]  # closed polygon
    plot(test, pch=19)
    title(main = "inferred range for Paragus haemorrhous")
    lines(coords, col="red")
    
    # now see which sites (that are being included) overlap with the range
    
}