## ************************************************************
## simulate data for single season model 
## ************************************************************

library(stringr)
library(sf); library(rgeos);
library(concaveman); library(tidyverse)
library(reshape2); library(locfit)

# all code for simulation adapted from original developer: https://github.com/lmguzman/occ_historical/blob/main/analysis/simulation/src/simulate_ms.R
# following their outline for using opportunistic data: https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13896

# ----------------------------------------------------
# Function: GenPoly()
# Generates a specified number of species
# ranges with a given complexity which
# specifies the mean number of vertices
# for the set of polygons. Grid size is used to
# generate an artificial landscape on which to place
# randomly generated polygons. Grids must
# be square. Additionally, the number of
# polygons to generate are also specified.
# ----------------------------------------------------
GenPoly <- function(numpoly=100,
                    complexity=50,
                    gridsize=100){
  
  # make a placeholder polygon to create
  # a grid.
  sfc <- st_sfc(st_polygon(list(rbind(c(0,0),
                                      c(1,0),
                                      c(1,1),
                                      c(0,0)))))
  
  # create the grid using the specified grid size.
  grid <- st_make_grid(sfc, cellsize=1/gridsize, square=TRUE)
  
  # create a list to hold all of the generated polygons.
  # iterate while generating the number of specified polygons.
  poly_ls <- list()
  i <- 1
  while(i <= numpoly){
    
    # draw random number of vertices from normal distribution
    n_vert <- round(rnorm(1, mean=complexity, sd=4))
    if(n_vert >= 3){
      verts <- st_sample(grid, size=n_vert)
      verts <- st_as_sf(verts)
      poly <- concaveman(verts) # generate the concave polygon for these vertices
      
      # randomly scale and translate the polygons to have different size ranges.
      scale_flag <- 1-(0.1*rpois(1, 1))
      # print(paste("Scaling by:", scale_flag))
      
      shift_flag_x <- 0.1*runif(1, -2, 2) # translate by x factor
      shift_flag_y <- 0.1*runif(1, -2, 2) # translate by y factor
      # print(paste("Translating by:", shift_flag_x, shift_flag_y))
      
      poly_geom <- st_geometry(poly)
      poly_shift <- (poly_geom*scale_flag+matrix(data=c(shift_flag_x, shift_flag_y), ncol=2))
      
      # check to make sure that the polygon is intersecting something on the grid
      if(length(st_intersects(grid,poly_shift)) <= 1){
        poly_shift <- concaveman(verts)
      }
      
      poly <- poly_shift %>% st_crop(grid)
      poly_ls[[i]] <- poly
      i <- i+1
    }
    else{
      n_vert <- round(rnorm(1, mean=complexity, sd=4))
    }
  }
  
  return(list(poly_ls, grid))
}

test <- GenPoly()

# ----------------------------------------------------
# Function: PolyToMatrix()
# Takes a list of polygons and converts them to a
# binary matrix of grid intersections. 
# A list of matrices is returned. Takes only a 
# large list object from the function GenPoly.
# ----------------------------------------------------
PolyToMatrix <- function(polys){
  matrix_ls <- list()
  n_poly <- length(polys[[1]])
  grid_size <- sqrt(length(polys[[2]]))
  
  for(i in 1:n_poly){
    mtrx <- st_intersects(st_as_sf(polys[[2]]), 
                          st_as_sf(polys[[1]][[i]]))
    matrix_ls[[i]] <- matrix(as.matrix(mtrx), nrow=grid_size, ncol=grid_size)
  }
  
  return(matrix_ls)
}

## expit and logit functions
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

make.ranges <- function(nsp, nsite, type.range) {
  
  if(type.range == 'all'){
    nsite.by.sp <- rep(nsite, nsp)
  }
  if(type.range == 'equal') {
    nsite.by.sp <- sample.int(n=nsite, size=nsp, replace=TRUE)
  } else if(type.range == 'logn') {
    prop.sites <- rbeta(nsp, shape1=1, shape2=3)
    nsite.by.sp <- round(nsite * prop.sites)
  } else if(type.range == 'polys'){
    polys <- GenPoly(numpoly=nsp, gridsize=sqrt(nsite))
    res <- do.call(rbind, lapply(PolyToMatrix(polys), as.vector))
    return(res)
  }
  
  get.sites.within.range <- function(ii) {
    sites <- rep(0,nsite)
    sites[sample(x=1:nsite, size=ii, replace=FALSE)] <- 1
    sites
  }
  res <- t(sapply(nsite.by.sp, get.sites.within.range))
  names(dim(res)) <- c('nsp','nsite')
  res==1
}


visit.history <- function(nsp, nsite, nyr, nvisit, type.visit, mu.v.0, mu.v.yr, prop.visits.same){
  
  ## for later: add a spatial component that maybe visits cluster in space through time?
  
  if(type.visit == 'visit_all'){
    
    ## all of the visits for all of the species
    
    vis.arr <- array(1,
                     dim=c(nsp=nsp, 
                           nsite=nsite,
                           nyr=nyr,
                           nvisit=nvisit))
    visit_collector <- array(1,
                             dim=c(nsite=nsite,
                                   nyr=nyr,
                                   nvisit=nvisit))
    
  }else if(type.visit == 'visit_miss'){
    
    ## there are missing visits
    
    ## Here some of the visits are the same for all species or some of the visits are different for species
    ## Depends on the proportion of visits that are the same
    ## If prop.visits.same =0 then all of the species have separate visits and if prop.visits.same=1 all of the species have the same visits 
    
    ## visit_collector is like the "Collector ID" where we identify which visits were community visits = 1 and which visits were independent sampling visits = 0
    
    ## provide parameters for the change of visits through time
    
    sigma.v.0 <- 1
    sigma.v.yr <- 1
    
    mu.v.0.sp <- rnorm(nsp, mu.v.0, sd = sigma.v.0)
    mu.v.yr.sp <- rnorm(nsp, mu.v.yr, sd = sigma.v.yr)
    
    vis.arr <- array(NA,
                     dim=c(nsp=nsp, 
                           nsite=nsite,
                           nyr=nyr,
                           nvisit=nvisit))
    
    visit_collector <- array(NA,
                             dim=c(nsite=nsite,
                                   nyr=nyr,
                                   nvisit=nvisit))
    
    for(site in 1:nsite) {
      for(yr in 1:nyr) {
        for(visit in 1:nvisit) {
          
          # for each visit sample whether its a community or individual visit
          
          same_vis <- rbinom(1, 1, prob = prop.visits.same)
          
          visit_collector[site,yr,visit] <- same_vis
          
          if(same_vis == 1){
            
            # if it is a community visit, then all of the species get the same
            
            vis.arr[,site,yr,visit] <- rbinom(1,1,expit(mu.v.0 + mu.v.yr*(yr)))
            
          }else(
            for(sp in 1:nsp){
              
              # If it is an individual species visit then not all of the species get the same
              
              vis.arr[sp,site,yr,visit] <- rbinom(1,1, expit(mu.v.0.sp[sp] + mu.v.yr.sp[sp]*(yr)))
              
            }
          )
        }
      }
    }
  }
  
  return(list(vis.arr = vis.arr, visit_collector = visit_collector))
}


make.data <- function(## data structure set up 
  nsp=20,
  nsite=18, ## number of sites
  nyr=10, ## number of years
  nvisit=6, ## number of samples per year
  ## detection
  mu.p.0 = 0,
  p.yr = 0.5,
  sigma.p.site = 0.3,
  sigma.p.sp = 0.3,
  ## occupancy
  mu.psi.0 = 0,
  sigma.psi.sp = 0.5,
  mu.psi.yr = 0.5, 
  sigma.psi.yr = 0.2,
  ## visit
  mu.v.0 = 0, 
  mu.v.yr = -0.5,
  ## type sym
  type.range = "equal",
  type.visit = 'visit_miss',
  prop.visits.same = 1){
  
  
  ## ------------------------------------------------------------
  ## Simulate species ranges them for the
  ## specified number of species and sites.
  
  sp.range <- make.ranges(nsp, nsite, type.range)
  
  ## ------------------------------------------------------------
  ## Get visit history
  
  vis.info <- visit.history(nsp=nsp, nsite=nsite, 
                            nyr=nyr, nvisit=nvisit, type.visit=type.visit, mu.v.0, mu.v.yr,
                            prop.visits.same=prop.visits.same)
  vis.arr <- vis.info$vis.arr
  
  ## --------------------------------------------------
  ## specify species-specific colonization and persistence probabilities
  
  ## species-specific random intercepts
  psi.sp <- rnorm(n=nsp, mean=0, sd=sigma.psi.sp)
  
  ## effect of year on occupancy (species-specific random slopes)
  psi.yr <- rnorm(n=nsp, mean=mu.psi.yr, sd=sigma.psi.yr)
  
  ### detection
  
  ## effect of site on detection (year-specific)
  p.site <- matrix(rnorm(n=nsite*nyr, mean=0, sd=sigma.p.site),
                   nrow=nsite,
                   ncol=nyr)
  p.sp   <- rnorm(n=nsp, mean=0, sd=sigma.p.sp)
  
  
  ## Create matrices for psi and phi
  
  psi.mat <- array(NA, dim =c(nsp, nsite, nyr))
  
  p.mat <- array(NA, dim =c(nsp, nsite, nyr, nvisit))
  
  for(site in 1:nsite) {
    for(yr in 1:nyr) {
      for(sp in 1:nsp) {
        
        psi.mat[sp,site,yr] <- expit(mu.psi.0 +
                                       psi.sp[sp] +
                                       psi.yr[sp]*(yr))
        
        for(visit in 1:nvisit) {
          
          p.mat[sp,site,yr,visit] <- expit(mu.p.0 +
                                             p.sp[sp] +
                                             p.site[site,yr] +
                                             p.yr*(yr))
        }
      }
    }
  }
  
  
  ## --------------------------------------------------
  ## Generate presence and absence
  
  Z <- array(NA, dim=c(nsp=nsp,
                       nsite=nsite,
                       nyr=nyr))
  
  for(yr in 1:nyr){
    for(sp in 1:nsp){
      
      Z[sp, which(sp.range[sp,]) ,yr] <- rbinom(n = sum(sp.range[sp,]), size = 1, prob = psi.mat[sp,which(sp.range[sp,]),yr])
      
    }
  }
  
  
  ## --------------------------------------------------
  ## Generate detection non detection data 
  
  V <- array(NA, dim=c(nsp=nsp,
                       nsite=nsite,
                       nyr=nyr,
                       nvisit = nvisit))
  
  for(sp in 1:nsp){
    for(yr in 1:nyr){
      for(v in 1:nvisit){
        V[sp, which(sp.range[sp,]),yr,v] <- rbinom(n = sum(sp.range[sp,]), size = 1, prob = Z[sp, which(sp.range[sp,]),yr] * p.mat[sp,which(sp.range[sp,]),yr, v])
      }
    }
  }
  
  ## --------------------------------------------------
  ## Generate visit data 
  
  X <-V*vis.arr
  
  
  ## add dimension names to arrays
  
  arr.names <- list(sp=paste('sp',str_pad(1:nsp,4,pad='0'),sep='_'),
                    site=paste('site',str_pad(1:nsite,4,pad='0'),sep='_'),
                    yr=paste('y',1:nyr,sep='_'),
                    visit=paste('v',1:nvisit,sep='_'))
  dimnames(sp.range) <- arr.names[1:2]
  dimnames(Z)        <- arr.names[1:3]
  dimnames(X)        <- arr.names[1:4]
  
  
  return(list(sp.range=sp.range,
              vis.arr=vis.arr,
              vis.col=vis.info$visit_collector,
              Z=Z,
              X=X,
              nsp=nsp,
              nsite=nsite, ## number of sites
              nyr=nyr, ## number of years
              nvisit=nvisit, ## number of samples per year
              ## detection
              mu.p.0 = mu.p.0,
              p.yr = p.yr,
              sigma.p.site = sigma.p.site,
              sigma.p.sp = sigma.p.sp,
              ## occupancy
              mu.psi.0 = mu.psi.0,
              sigma.psi.sp = sigma.psi.sp,
              mu.psi.yr = mu.psi.yr, 
              sigma.psi.yr = sigma.psi.yr,
              ## visit
              mu.v.0 = mu.v.0, 
              mu.v.yr = mu.v.yr,
              ## species specific values
              psi.sp = psi.sp,
              psi.yr = psi.yr,
              ## type sym
              type.range = type.range,
              type.visit = type.visit,
              prop.visits.same = prop.visits.same))
  
}

my.data <- make.data()

#############################
# prep data #################
#############################
prep.data <- function(dd, limit.to.visits, limit.to.range, time.interval.yr, time.interval.visit) {
  
  ## keep only detected species:
  sp.keep <- which(apply(dd$X, 1, sum, na.rm=TRUE)>0)
  
  ## which sites to keep will depend on limit.to.visits
  ## keep all sites, even those without any detections
  if(limit.to.visits=='all') {
    site.keep <- 1:dd$nsite
  }
  ## keep only sites that were visited
  if(limit.to.visits=='visits') {
    site.keep <- which(apply(dd$vis.arr, 2, sum, na.rm = TRUE)>0)
  }
  ## keep only sites that yielded a detection of at least one species
  if(limit.to.visits=='detected') {
    site.keep <- which(apply(dd$X, 'site', sum, na.rm = TRUE)>0)
  }
  
  ## subset based on the above
  dd$nsite.old  <- dim(dd$X)[2]
  dd$Z        <- dd$Z[sp.keep,site.keep,,drop=FALSE]
  dd$X        <- dd$X[sp.keep,site.keep,,,drop=FALSE]
  dd$sp.range <- dd$sp.range[sp.keep,site.keep,drop=FALSE]
  dd$vis.arr  <- dd$vis.arr[sp.keep,site.keep,,,drop=FALSE]
  dd$nsp      <- length(sp.keep)
  dd$nsite    <- length(site.keep)
  dd$vis.col <- dd$vis.col[site.keep,,,drop=FALSE]
  
  
  ####### time interval ###
  
  ## convert data into long format
  
  presence.only <- which(dd$X==1, arr.ind=TRUE) %>% data.frame()
  
  visits.only <- which(dd$vis.arr==1, arr.ind=TRUE) %>% data.frame()
  
  names(presence.only) <- c("spn", "siten", "yr", "visit")
  
  names(visits.only) <- c("spn", "siten", "yr", "visit")
  
  new_year <- expand.grid(yr = 1:dd$nyr, visit = 1:dd$nvisit) %>% 
    arrange(yr) %>% 
    mutate(syr = 1:n())
  
  presence.new.year <- presence.only %>% 
    left_join(new_year)
  
  visits.new.year <- visits.only %>% 
    left_join(new_year)
  
  #### use given time interval to bin ###
  
  # time.interval.yr <- 2
  # time.interval.visit <- 3
  
  unique.syr <- unique(sort(presence.new.year$syr))
  
  time.cut.yr <- as.numeric(cut(unique.syr, time.interval.yr))
  
  time.cut.visit <- as.numeric(cut(unique.syr, time.interval.visit*time.interval.yr))
  
  time.cut.visit.2 <-  rep(rep(1:time.interval.visit, each = max(unique.syr)/(time.interval.yr*time.interval.visit)), 
                           time.interval.yr)
  
  new.time.interval <-data.frame(syr = unique.syr, time.cut.yr, time.cut.visit, time.cut.visit.2)
  
  sp_names <- data.frame(spn = 1:length(dimnames(dd$X)$sp), sp = dimnames(dd$X)$sp)
  site_names <-data.frame(siten = 1:length(dimnames(dd$X)$site), site = dimnames(dd$X)$site)
  
  presence.new.time.interval <- presence.new.year %>% 
    left_join(new.time.interval) %>% 
    group_by(spn, siten, time.cut.yr, time.cut.visit.2) %>% 
    summarise(presence = 1) %>% 
    left_join(sp_names) %>% left_join(site_names) %>% 
    mutate(yr = paste('yr',time.cut.yr,sep='_'),
           visit = paste('v',time.cut.visit.2,sep='_'))
  
  visits.new.time.interval <- visits.new.year %>% 
    left_join(new.time.interval) %>% 
    group_by(spn, siten, time.cut.yr, time.cut.visit.2) %>% 
    summarise(presence = 1) %>% 
    left_join(sp_names) %>% left_join(site_names) %>% 
    mutate(yr = paste('yr',time.cut.yr,sep='_'),
           visit = paste('v',time.cut.visit.2,sep='_'))
  
  
  # presence.new.time.interval$time.cut.yr %>% table()
  # presence.new.time.interval$time.cut.visit %>% table()
  
  occ.arr <- array(0, dim = c(dd$nsp, dd$nsite, time.interval.yr, time.interval.visit), 
                   dimnames = list(sp=dimnames(dd$X)$sp,
                                   site=dimnames(dd$X)$site,
                                   year= paste0("yr_", 1:time.interval.yr),
                                   visit=paste0("v_", 1:time.interval.visit)))
  
  occ.arr[cbind(match(presence.new.time.interval$sp, dimnames(dd$X)$sp), match(presence.new.time.interval$site, dimnames(dd$X)$site), 
                match(presence.new.time.interval$yr, paste0("yr_", 1:time.interval.yr)), 
                match(presence.new.time.interval$visit,paste0("v_", 1:time.interval.visit)))] <- 1 
  
  
  vis.arr <- array(0, dim = c(dd$nsp, dd$nsite, time.interval.yr, time.interval.visit), 
                   dimnames = list(sp=dimnames(dd$X)$sp,
                                   site=dimnames(dd$X)$site,
                                   year= paste0("yr_", 1:time.interval.yr),
                                   visit=paste0("v_", 1:time.interval.visit)))
  
  names(dim(occ.arr)) <- c("nsp", "nsite", 'nyr', 'nvisit')
  
  vis.arr[cbind(match(visits.new.time.interval$sp, dimnames(dd$X)$sp), match(visits.new.time.interval$site, dimnames(dd$X)$site), 
                match(visits.new.time.interval$yr, paste0("yr_", 1:time.interval.yr)), 
                match(visits.new.time.interval$visit,paste0("v_", 1:time.interval.visit)))] <- 1 
  
  dd$X2 <- occ.arr
  dd$vis.arr2 <- vis.arr
  
  
  ## generate master index (to improve model efficiency (this prevents
  ## unnecessary iterating through all irrelevant sites and visits)
  get.indices <- function(sp) {
    ## visited array
    vis.arr <- dd$vis.arr2[sp,,,]
    ## if modelling all visits, set visit array to 1 everywhere
    if(limit.to.visits=='all') 
      vis.arr[TRUE] <- 1
    ## if modelling sites with detections only, create new visit array
    if(limit.to.visits=='detected') {
      nsp.detected <- apply(dd$X2, 2:4, sum, na.rm = TRUE)
      vis.arr[TRUE] <- 1
      vis.arr[nsp.detected==0] <- 0
    }
    
    # either restrict or not restrict inference to inferred species ranges
    if(limit.to.range=='yes'){ # restrict to within species' ranges
      # infer species ranges from affirmative detections
      # make a placeholder polygon to create
      # a grid.
      sfc <- st_sfc(st_polygon(list(rbind(c(0,0),
                                          c(1,0),
                                          c(1,1),
                                          c(0,0)))))
      
      # create the grid using the specified number of simulated sites
      grid <- st_make_grid(sfc, cellsize=1/sqrt(dd$nsite.old), square=TRUE) %>%
        st_as_sf() %>%
        dplyr::mutate(GID=row_number())
      
      # iteratively add points for each species based on the centroids of
      # each grid cell in the matrix, then construct a convex polygon for each
      # species and add to a list
      
      
      sp_GID <- rowSums(dd$X2[sp,,,], c(2,3,4))
      sp_centroids <- c()
      for(j in 1:dd$nsite){
        if(sp_GID[j] >= 1){
          sp_centroids <- sp_centroids %>% base::append(st_centroid(grid[j,])$x)
          
        }
      }
      
      sp_poly <- st_convex_hull(st_union(sp_centroids))
      grid_size <- sqrt(dd$nsite.old)
      mtrx <- st_intersects(st_as_sf(grid),
                            st_as_sf(sp_poly))
      matrix_ls <- matrix(as.matrix(mtrx), nrow=grid_size, ncol=grid_size)
      
      vector_range <- as.vector(matrix_ls)
      vis.arr[!vector_range[site.keep],,] <- 0
      tmp <- which(vis.arr==1, arr.ind=TRUE)
      indices <- cbind(rep(sp,nrow(tmp)),tmp)
      
    }
    if(limit.to.range=='no'){ # unrestricted with respect to species' ranges
      vis.arr[!dd$sp.range[sp,],,] <- 1
      tmp <- which(vis.arr==1, arr.ind=TRUE)
      indices <- cbind(rep(sp,nrow(tmp)), tmp)
    }
    
    return(indices)
    
  }
  master.index <- do.call(rbind, lapply(1:dd$nsp, get.indices))
  colnames(master.index) <- c('sp','site','yr','visit')
  
  ## data structures to be returned
  
  my.data <- list(X=dd$X2[master.index])
  
  my.constants <- list(
    nsp=dim(dd$X2)['nsp'],
    nsite=dim(dd$X2)['nsite'],
    nyr=dim(dd$X2)['nyr'],
    nind=nrow(master.index),
    yrv=master.index[,'yr'],
    sitev=master.index[,'site'],
    spv=master.index[,'sp'])
  
  my.info <- list(
    time.interval.yr=time.interval.yr,
    time.interval.visit=time.interval.visit,
    limit.to.visits=limit.to.visits,
    limit.to.range=limit.to.range,
    sp.keep = names(sp.keep), 
    site.keep = names(site.keep)
  )
  
  return(list(my.constants = my.constants, my.data = my.data, my.info = my.info))
}
