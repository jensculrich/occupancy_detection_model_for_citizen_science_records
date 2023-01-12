## get K (max search values for binomial N-mixture model component)
# jcu; started jan 11, 2023

get_K <- function(V_citsci, K_addition, K_multiplier){

## --------------------------------------------------
# Max search range for each species*site*interval
K <- array(dim = c(n_species, n_sites, n_intervals))

for(i in 1:n_species){
  for(j in 1:n_sites){
    for(k in 1:n_intervals){
      
      # make K the same across all intervals in case detection changes too muh through time
      K[i,j,k] <- K_multiplier * (max(V_citsci[i,j,,]) + K_addition)
      
    }
  }
}

# preview V_museum
# K[1:n_species,1:n_sites,1]

  return(K)

}