# make sure we have everything needed to make a corr matrix w signifcance testing
library(Hmisc)
library(corrplot)

# read in the data
# manually
my_data <- readRDS(
  "./occupancy/analysis/prepped_data/bombus/10km_1200minpop_1minUniqueDetections_4ints_3visits.rds")
# or using analysis selections
my_data <- readRDS(paste0("./occupancy/analysis/prepped_data/",
                                     taxon, "/", grid_size / 1000, "km_",
                                     min_population_size, "minpop_",
                                     min_records_per_species, "minUniqueDetections", "_",
                                     n_intervals, "ints_",
                                     n_visits, "visits",
                                     #"_nonurban",
                                     ".rds"))


correlation_matrix <- my_data$correlation_matrix

rm(my_data)
gc()

correlation_matrix <- rcorr(correlation_matrix)

colnames(correlation_matrix$r) <- c("population density", 
                                    "site area", 
                                    "developed open space",
                                    "low development",
                                    "med/high development",
                                    "low/med/high development",
                                    "income",
                                    "natural habitat (herb/shrub)",
                                    "natural habitat (forest)",
                                    "natural habitat (any)"
                                    )

rownames(correlation_matrix$r) <- c("population density", 
                                    "site area", 
                                    "developed open space",
                                    "low development",
                                    "med/high development",
                                    "low/med/high development",
                                    "income",
                                    "natural habitat (herb/shrub)",
                                    "natural habitat (forest)",
                                    "natural habitat (any)"
                                    )

colnames(correlation_matrix$P) <- c("population density", 
                                    "site area", 
                                    "developed open space",
                                    "low development",
                                    "med/high development",
                                    "low/med/high development",
                                    "income",
                                    "natural habitat (herb/shrub)",
                                    "natural habitat (forest)",
                                    "natural habitat (any)"
                                    )

rownames(correlation_matrix$P) <- c("population density", 
                                    "site area", 
                                    "developed open space",
                                    "low development",
                                    "med/high development",
                                    "low/med/high development",
                                    "income",
                                    "natural habitat (herb/shrub)",
                                    "natural habitat (forest)",
                                    "natural habitat (any)"
                                    )

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

flattenCorrMatrix(correlation_matrix$r, correlation_matrix$P)

library(corrplot)
corrplot(correlation_matrix$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

# Insignificant correlations are leaved blank
corrplot(correlation_matrix$r, type="upper", 
         order="original", 
         diag = FALSE, 
         tl.col="black",
         tl.srt=45,
         p.mat = correlation_matrix$P, sig.level = 0.05, insig = "blank")

