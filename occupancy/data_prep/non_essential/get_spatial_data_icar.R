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

# Finally, extract spatial covariate data including:
# population density
# impervious surface cover
# site area

# This file contains a function 'get_spatial_data' which
# can be called from a 'run_model_X_.R' file with specific data collection
# choices of grid size and minimum pop in the cell to be chosen or varied

# The file 'explore_spatial_data.R' walks through this process, with plots to visualize
# how the spatial data is gathered and how data 'collection' choices 
# like grid size and minimum urbanization intensity to include a site
# shape the data that we will analyze with our model

library(tidyverse) # data carpentry
library(tigris) # get state shapefile
library(sf) # spatial data processing
library(raster) # read and format raster data
library(exactextractr) # quick extraction of raster data
library(spdep) # create site-wise adjacency matrix for spatial smoothing
library(INLA) # create site-wise adjacency matrix for spatial smoothing

get_spatial_data <- function(
  grid_size, # square edge dimensions of a site, in meters
  min_population_size, # minimum pop density of a site to be considered "urban"
  taxon, # prepare data for syrphidae or bombus
  min_site_area # min land area (in the state admin areas) of a site to be included
){
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  ## To save time and computing power, if study dimensions match preprocessing reqts.
  ## then use the preprocessed data to save time, 
  #if(grid_size == 30000 && min_population_size == 300 && min_site_area == 0.10){
    
    # read the occurrence data for the given taxon
    #df_id_dens <- readRDS(paste0(
      #"./preprocessed_data/",
      #taxon,
      #"_occurrence_records_30km_300minpop_.RDS"))
    
    #grid_pop_dens <- readRDS("./preprocessed_data/site_data_30km_300minpop_.RDS")
    
    ## --------------------------------------------------
    # Extract variables
    
    #scaled_pop_density <- grid_pop_dens %>% 
      #pull(scaled_pop_den_km2)
    
    #scaled_site_area <- grid_pop_dens %>% 
      #pull(scaled_site_area)
    
    #scaled_developed_open <- grid_pop_dens %>% 
      #pull(scaled_developed_open)
    
    #scaled_herb_shrub_cover <- grid_pop_dens %>% 
      #pull(scaled_herb_shrub_cover)
    
    #scaled_forest <- grid_pop_dens %>% 
      #pull(scaled_forest)
    
    #scaled_developed_med_high <- grid_pop_dens %>% 
      #pull(scaled_developed_med_high)
    
    ## --------------------------------------------------
    # Calculate correlations between variables 
    
    #my_variables <- as.data.frame(cbind(scaled_pop_density, 
                                        #scaled_site_area, 
                                        #scaled_developed_open,
                                        #scaled_developed_med_high,
                                        #scaled_herb_shrub_cover, 
                                        #scaled_forest))
    
    #correlation_matrix <- cor(my_variables)
    
  #} else{
    
    ## --------------------------------------------------
    # Envrionmental raster data
    
    # take this chunk out
    # CRS for NAD83 / UTM Zone 10N
    crs <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"
    
    # spatial data - California state shapefile
    states <- tigris::states() %>%
      filter(NAME %in% c("California", "Oregon", "Washington", "Arizona", "Nevada"))
    
    states_trans <- states  %>% 
      st_transform(., 26910) # NAD83 / UTM Zone 10N
    
    ## --------------------------------------------------
    # Environmental data rasters
    
    # pop density raster
    # https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
    # 2015 pop density at 1km resolution
    pop_raster=raster("./data/spatial_data/population_density/gpw_v4_population_density_rev11_2015_30_sec.tif")
    
    # land cover raster 30m x 30m
    # https://www.mrlc.gov/data/nlcd-2016-land-cover-conus
    # land cover data from 2016 
    land=raster::raster("./data/spatial_data/land_use/land_use.tif")
    
    ## --------------------------------------------------
    # Overlay the shapefile with a grid of sites of size == 'grid_size' 
    
    # create _km grid - here you can substitute by specifying grid_size above
    grid <- st_make_grid(states_trans, cellsize = c(grid_size, grid_size)) %>% 
      st_sf(grid_id = 1:length(.))
    
    ## --------------------------------------------------
    # Extract mean population density in each grid cell
    
    pop_raster <- crop(pop_raster, states)
    pop_raster <- mask(pop_raster, states)
    
    # project the grid to the raster
    crs_raster <- sf::st_crs(raster::crs(pop_raster))
    prj1 <- st_transform(grid, crs_raster)
    
    # Extract raster values to list object, and then summarize by the mean value
    r.vals <- exactextractr::exact_extract(pop_raster, prj1, 'mean')
    
    ## --------------------------------------------------
    # filter sites to urban areas using the r.vals for pop density
    
    grid_pop_dens <- cbind(grid, r.vals) %>% 
      rename("pop_density_per_km2" = "r.vals")
    # now filter out the non-urban areas (areas below our pop density threshold)
    # and make a scaled response variable
    grid_pop_dens <- grid_pop_dens %>%
      filter(pop_density_per_km2 > min_population_size) 
    
    # free unused space
    rm(pop_raster, prj1, r.vals, crs_raster)
    gc(verbose = FALSE)
    
    
    ## --------------------------------------------------
    # Extract environmental variables from each remaining site
    
    # make sure that the grid is still projected to the raster
    crs_raster <- sf::st_crs(raster::crs(land))
    prj1 <- st_transform(grid_pop_dens, crs_raster)
    
    # then extract values and cbind with the grid_pop_dens
    # extract raster values to list object
    
    r.vals_land <- exactextractr::exact_extract(land, prj1)
    
    # want to reclassify open water as NA
    r.vals_land_NA <- lapply(r.vals_land, function(x) na_if(x$value, 0))
    r.vals_land_NA <- lapply(r.vals_land_NA, function(x) na_if(x,11))
    
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
    
    
    # now drop NA values so that the below estimates are the proportion of cover
    # of all land cover in the administrative area
    r.vals_land_NA <- lapply(r.vals_land_NA, na.omit)
    
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
       r.site_area, r.vals_land, r.vals_land_NA, states, states_trans)
    gc()
    
    ## --------------------------------------------------
    # Connect sites for spatial autocorrelation smoothing
    # https://github.com/ConnorDonegan/Stan-IAR
    
    # cite: Donegan, Connor. Flexible Functions for ICAR, BYM, and BYM2 Models in Stan. 
    # Code Repository. 2021. 
    # Available online: https://github.com/ConnorDonegan/Stan-IAR (access date).
    
    C <- spdep::nb2mat(spdep::poly2nb(grid_pop_dens, queen = TRUE), style = "B", zero.policy = TRUE)
    
    #' convert connectivity matrix to unique pairs of connected nodes (graph structure)
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
    
    #' compute scaling factor for adjacency matrix
    #' accounts for differences in spatial connectivity 
    scale_c <- function(C) {
      
      #' compute geometric mean of a vector
      geometric_mean <- function(x) exp(mean(log(x))) 
      
      N = dim(C)[1]
      
      # Create ICAR precision matrix  (diag - C): this is singular
      # function Diagonal creates a square matrix with given diagonal
      Q =  Matrix::Diagonal(N, rowSums(C)) - C
      
      # Add a small jitter to the diagonal for numerical stability (optional but recommended)
      Q_pert = Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)
      
      # Function inla.qinv provides efficient way to calculate the elements of the
      # the inverse corresponding to the non-zero elements of Q
      Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N),e=0))
      
      # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
      scaling_factor <- geometric_mean(Matrix::diag(Q_inv)) 
      return(scaling_factor) 
    }
    
    #' prepare Stan data for ICAR model given a connectivity matrix
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
    
    ## calculate the scale factor for each of k connected group of nodes, using the scale_c function from M. Morris
    # k <- icar.data$k
    k <- 1
    scale_factor <- vector(mode = "numeric", length = k)
    for (j in 1:k) {
      g.idx <- which(icar.data$comp_id == j) 
      if (length(g.idx) == 1) {
        scale_factor[j] <- 1
        next
      }    
      Cg <- C[g.idx, g.idx] 
      scale_factor[j] <- scale_c(Cg) 
    }
    
    ## update the data list for Stan
    icar.data$inv_sqrt_scale_factor <- 1 / sqrt( scale_factor )
    
    rm(C, edges, scale_c, prep_icar_data)
    gc()
    
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
    # Now tag records with site ID's based on spatial intersection
    
    # which grid square is each point in?
    df_id_dens <- df_trans %>% 
      st_join(grid_pop_dens, join = st_intersects) %>% as.data.frame %>%
      # filter out records from outside of the urban grid
      filter(!is.na(grid_id)) %>%
      left_join(., dplyr::select(
        df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") 
    
    # free unused space
    rm(df, df_sf, df_trans)
    gc()
    
    
    ## --------------------------------------------------
    # Calculate correlations between variables
    
    correlation_matrix <- cor(as.data.frame(cbind(grid_pop_dens$scaled_pop_den_km2,
                                                  grid_pop_dens$scaled_site_area,
                                                  grid_pop_dens$scaled_developed_open,
                                                  grid_pop_dens$scaled_developed_med_high,
                                                  grid_pop_dens$scaled_herb_shrub_cover,
                                                  grid_pop_dens$scaled_forest,
                                                  grid_pop_dens$scaled_herb_shrub_forest)))
    
    colnames(correlation_matrix) <- c("scaled_pop_den_km2", 
                                      "scaled_site_area", 
                                      "scaled_developed_open",
                                      "scaled_developed_med_high",
                                      "scaled_herb_shrub_cover",
                                      "scaled_forest",
                                      "scaled_herb_shrub_forest")
    
    rownames(correlation_matrix) <- c("scaled_pop_den_km2", 
                                      "scaled_site_area", 
                                      "scaled_developed_open",
                                      "scaled_developed_med_high",
                                      "scaled_herb_shrub_cover",
                                      "scaled_forest",
                                      "scaled_herb_shrub_forest")
    
  #} # end else
  
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    df_id_urban_filtered = df_id_dens,
    urban_grid = grid_pop_dens,
    correlation_matrix = correlation_matrix,
    
    n_edges = icar.data$n_edges,
    node1 = icar.data$node1,
    node2 = icar.data$node2,
    inv_sqrt_scale_factor = icar.data$inv_sqrt_scale_factor,
    k = icar.data$k
    
  ))
}

