# File to make sensitivity analysis figures and tables
# jcu, started mar 28, 2023.

library(tidyverse)

## --------------------------------------------------
## Sensitivity to temporal divisions

## --------------------------------------------------
## Read in model run results

ten_km <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits_.rds")
fifteen_km <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_15km_1000minpop_2minUniqueDetections_3ints_3visits_.rds")
twenty_km <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_20km_800minpop_2minUniqueDetections_3ints_3visits_.rds")
twentyfive_km <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_25km_600minpop_2minUniqueDetections_3ints_3visits_.rds")

fit_summary_10 <- rstan::summary(ten_km)
fit_summary_15 <- rstan::summary(fifteen_km)
fit_summary_20 <- rstan::summary(twenty_km)
fit_summary_25 <- rstan::summary(twentyfive_km)

View(cbind(1:nrow(fit_summary_10$summary), fit_summary_10$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_15$summary), fit_summary_15$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_20$summary), fit_summary_20$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_25$summary), fit_summary_25$summary)) # View to see which row corresponds to the parameter of interest

## --------------------------------------------------
## Make figure table 

# parameter means
params = 3 # mu_psi_0, mu_psi_natural_habitat[native], psi race
n_grains = 4 # 10, 15, 20, and 25km scales

x <- (rep(1:params, each=n_grains)) # parameter reference
y <- (rep(1:n_grains, times=params)) # spatial grain reference


estimate <-  c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[286,1],
  fit_summary_15$summary[288,1],
  fit_summary_20$summary[316,1],
  fit_summary_25$summary[330,1],

  # param 2 (psi race)
  fit_summary_10$summary[297,1],
  fit_summary_15$summary[299,1],
  fit_summary_20$summary[327,1],
  fit_summary_25$summary[341,1],
  
  # param 2 (psi natural (native species))
  fit_summary_10$summary[1335,1],
  fit_summary_15$summary[1140,1],
  fit_summary_20$summary[1168,1],
  fit_summary_25$summary[1193,1]
  
  
)

lower_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[286,4],
  fit_summary_15$summary[288,4],
  fit_summary_20$summary[316,4],
  fit_summary_25$summary[330,4],

  # param 2 (psi race)
  fit_summary_10$summary[297,4],
  fit_summary_15$summary[299,4],
  fit_summary_20$summary[327,4],
  fit_summary_25$summary[341,4],
  
  # param 2 (psi natural (native species))
  fit_summary_10$summary[1335,4],
  fit_summary_15$summary[1140,4],
  fit_summary_20$summary[1168,4],
  fit_summary_25$summary[1193,4]
) 

upper_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[286,8],
  fit_summary_15$summary[288,8],
  fit_summary_20$summary[316,8],
  fit_summary_25$summary[330,8],

  # param 2 (psi race)
  fit_summary_10$summary[297,8],
  fit_summary_15$summary[299,8],
  fit_summary_20$summary[327,8],
  fit_summary_25$summary[341,8],
  
  # param 2 (psi natural (native species))
  fit_summary_10$summary[1335,8],
  fit_summary_15$summary[1140,8],
  fit_summary_20$summary[1168,8],
  fit_summary_25$summary[1193,8]
) 

lower_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[286,5],
  fit_summary_15$summary[288,5],
  fit_summary_20$summary[316,5],
  fit_summary_25$summary[330,5],

  # param 2 (psi race)
  fit_summary_10$summary[297,5],
  fit_summary_15$summary[299,5],
  fit_summary_20$summary[327,5],
  fit_summary_25$summary[341,5],
  
  # param 2 (psi natural (native species))
  fit_summary_10$summary[1335,5],
  fit_summary_15$summary[1140,5],
  fit_summary_20$summary[1168,5],
  fit_summary_25$summary[1193,5]
) 

upper_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[286,7],
  fit_summary_15$summary[288,7],
  fit_summary_20$summary[316,7],
  fit_summary_25$summary[330,7],

  # param 2 (psi race)
  fit_summary_10$summary[297,7],
  fit_summary_15$summary[299,7],
  fit_summary_20$summary[327,7],
  fit_summary_25$summary[341,7],
  
  # param 2 (psi natural (native species))
  fit_summary_10$summary[1335,7],
  fit_summary_15$summary[1140,7],
  fit_summary_20$summary[1168,7],
  fit_summary_25$summary[1193,7]
) 

df <- as.data.frame(cbind(x, y, estimate, 
                          lower_95, upper_95,
                          lower_50, upper_50 
)) %>%
  mutate(x = as.factor(x),
         y = as.factor(y))



# converting the result to dataframe
#df <- as.data.frame(rev_data_frame)

(v <- ggplot(df, aes(x, y, width=1, height=1)) +
    geom_tile(aes(fill = estimate)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1, 2, 3),
                     labels=c(bquote(psi[0]),
                              bquote(psi[race]),
                              bquote(mu[psi["nat. green."]~"[native]"])

                     )) +
    scale_y_discrete(name="Spatial grain \n (landscape width (km))", breaks = rep(1:n_grains),
                     labels=c("10", "15", "20", "25")) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3"),
                         name="mean parameter \nestimate",
                         limits=c(-1,2), breaks = c(-1, 0, 1, 2)) +
    #geom_text(data = df, 
    #     aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
    
    geom_text(data = df, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "50% BCI: [", signif(lower_50,2), ", ", signif(upper_50,2), "]")),
              vjust = -1, size = 5) + 
    geom_text(data = df, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "95% BCI: [", signif(lower_95,2), ", ", signif(upper_95,2), "]")),
              vjust = 1, size = 5) +
    
    
    theme(legend.position = "right",
          legend.text=element_text(size=14),
          legend.title=element_text(size=16),
          axis.text.x = element_text(size = 24, angle = 45, hjust=1),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(size = 12),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank())
)

## --------------------------------------------------
## Sensitivity to temporal divisions 

## --------------------------------------------------
## Read in model run results

three <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits_.rds")
two <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1200minpop_2minUniqueDetections_5ints_2visits_.rds")
four <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1200minpop_2minUniqueDetections_2ints_4visits_.rds")

fit_summary_three <- rstan::summary(three)
fit_summary_four <- rstan::summary(four)
fit_summary_two <- rstan::summary(two)

View(cbind(1:nrow(fit_summary_three$summary), fit_summary_three$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_four$summary), fit_summary_four$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_two$summary), fit_summary_two$summary)) # View to see which row corresponds to the parameter of interest


## --------------------------------------------------
## Make figure table 

# parameter means
params = 3 # mu_psi_0, natural habitat native species, income 
n_grains = 3 # 3 4 or 5 intervals (divided by 4, 3, or 2 years)

x <- (rep(1:params, each=n_grains)) # parameter reference
y <- (rep(1:n_grains, times=params)) # spatial grain reference


estimate <-  c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[286,1],
  fit_summary_three$summary[286,1],
  fit_summary_four$summary[286,1],

  # param 3 (psi race)
  fit_summary_two$summary[297,1],
  fit_summary_three$summary[297,1],
  fit_summary_four$summary[297,1],
  
  # param 3 (psi nat greenspace)
  fit_summary_two$summary[1335,1],
  fit_summary_three$summary[1335,1],
  fit_summary_four$summary[1335,1]
)

lower_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[286,4],
  fit_summary_three$summary[286,4],
  fit_summary_four$summary[286,4],
  
  # param 3 (psi race)
  fit_summary_two$summary[297,4],
  fit_summary_three$summary[297,4],
  fit_summary_four$summary[297,4],
  
  # param 3 (psi nat greenspace)
  fit_summary_two$summary[1335,4],
  fit_summary_three$summary[1335,4],
  fit_summary_four$summary[1335,4]
) 

upper_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[286,8],
  fit_summary_three$summary[286,8],
  fit_summary_four$summary[286,8],
  
  # param 3 (psi race)
  fit_summary_two$summary[297,8],
  fit_summary_three$summary[297,8],
  fit_summary_four$summary[297,8],
  
  # param 3 (psi nat greenspace)
  fit_summary_two$summary[1335,8],
  fit_summary_three$summary[1335,8],
  fit_summary_four$summary[1335,8]
) 

lower_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[286,5],
  fit_summary_three$summary[286,5],
  fit_summary_four$summary[286,5],
  
  # param 3 (psi race)
  fit_summary_two$summary[297,5],
  fit_summary_three$summary[297,5],
  fit_summary_four$summary[297,5],
  
  # param 3 (psi nat greenspace)
  fit_summary_two$summary[1335,5],
  fit_summary_three$summary[1335,5],
  fit_summary_four$summary[1335,5]
) 

upper_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[286,7],
  fit_summary_three$summary[286,7],
  fit_summary_four$summary[286,7],
  
  # param 3 (psi race)
  fit_summary_two$summary[297,7],
  fit_summary_three$summary[297,7],
  fit_summary_four$summary[297,7],
  
  # param 3 (psi nat greenspace)
  fit_summary_two$summary[1335,7],
  fit_summary_three$summary[1335,7],
  fit_summary_four$summary[1335,7]
) 

df <- as.data.frame(cbind(x, y, estimate, 
                          lower_95, upper_95,
                          lower_50, upper_50 
)) %>%
  mutate(x = as.factor(x),
         y = as.factor(y))



# converting the result to dataframe
#df <- as.data.frame(rev_data_frame)

(w <- ggplot(df, aes(x, y, width=1, height=1)) +
    geom_tile(aes(fill = estimate)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1, 2, 3),
                     labels=c(bquote(psi[0]),
                              bquote(psi[race]),
                              bquote(mu[psi["nat. green."]~"[native]"])
                              
                     )) +
    scale_y_discrete(name="Years per occupancy interval", breaks = rep(1:n_grains),
                     labels=c("2", "3", "4")) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3"),
                         name="mean parameter \nestimate",
                         limits=c(-1,2), breaks = c(-1, 0, 1, 2)) +
    #geom_text(data = df, 
    #     aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
    
    geom_text(data = df, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "50% BCI: [", signif(lower_50,2), ", ", signif(upper_50,2), "]")),
              vjust = -1, size = 5) + 
    geom_text(data = df, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "95% BCI: [", signif(lower_95,2), ", ", signif(upper_95,2), "]")),
              vjust = 1, size = 5) +
    
    
    theme(legend.position = "right",
          legend.text=element_text(size=14),
          legend.title=element_text(size=16),
          axis.text.x = element_text(size = 24, angle = 45, hjust=1),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(size = 12),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank())
                              
)


## --------------------------------------------------
## Sensitivity to population density threshold

## --------------------------------------------------
## Read in model run results

low <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1000minpop_2minUniqueDetections_3ints_3visits_.rds")
medium <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits_.rds")
high <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1400minpop_2minUniqueDetections_3ints_3visits_.rds")

fit_summary_medium <- rstan::summary(medium)
fit_summary_high <- rstan::summary(high)
fit_summary_low <- rstan::summary(low)

View(cbind(1:nrow(fit_summary_medium$summary), fit_summary_medium$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_high$summary), fit_summary_high$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_low$summary), fit_summary_low$summary)) # View to see which row corresponds to the parameter of interest


## --------------------------------------------------
## Make figure table 

# parameter means
params = 3 # mu_psi_0, natural habitat native species 
n_grains = 3 # 3 4 or 6 intervals (divided by 4, 3, or 2 years)

x <- (rep(1:params, each=n_grains)) # parameter reference
y <- (rep(1:n_grains, times=params)) # spatial grain reference


estimate <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[308,1],
  fit_summary_medium$summary[286,1],
  fit_summary_high$summary[276,1],

  # param 3 (psi race)
  fit_summary_low$summary[319,1],
  fit_summary_medium$summary[297,1],
  fit_summary_high$summary[287,1],
  
  # param 3 (psi natural)
  fit_summary_low$summary[1563,1],
  fit_summary_medium$summary[1335,1],
  fit_summary_high$summary[1185,1]
 
)

lower_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[308,4],
  fit_summary_medium$summary[286,4],
  fit_summary_high$summary[276,4],
  
  # param 3 (psi race)
  fit_summary_low$summary[319,4],
  fit_summary_medium$summary[297,4],
  fit_summary_high$summary[287,4],
  
  # param 3 (psi natural)
  fit_summary_low$summary[1563,4],
  fit_summary_medium$summary[1335,4],
  fit_summary_high$summary[1185,4]
) 

upper_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[308,8],
  fit_summary_medium$summary[286,8],
  fit_summary_high$summary[276,8],
  
  # param 3 (psi race)
  fit_summary_low$summary[319,8],
  fit_summary_medium$summary[297,8],
  fit_summary_high$summary[287,8],
  
  # param 3 (psi natural)
  fit_summary_low$summary[1563,8],
  fit_summary_medium$summary[1335,8],
  fit_summary_high$summary[1185,8]
) 

lower_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[308,5],
  fit_summary_medium$summary[286,5],
  fit_summary_high$summary[276,5],
  
  # param 3 (psi race)
  fit_summary_low$summary[319,5],
  fit_summary_medium$summary[297,5],
  fit_summary_high$summary[287,5],
  
  # param 3 (psi natural)
  fit_summary_low$summary[1563,5],
  fit_summary_medium$summary[1335,5],
  fit_summary_high$summary[1185,5]
) 

upper_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[308,7],
  fit_summary_medium$summary[286,7],
  fit_summary_high$summary[276,7],
  
  # param 3 (psi race)
  fit_summary_low$summary[319,7],
  fit_summary_medium$summary[297,7],
  fit_summary_high$summary[287,7],
  
  # param 3 (psi natural)
  fit_summary_low$summary[1563,7],
  fit_summary_medium$summary[1335,7],
  fit_summary_high$summary[1185,7]
) 

df <- as.data.frame(cbind(x, y, estimate, 
                          lower_95, upper_95,
                          lower_50, upper_50 
)) %>%
  mutate(x = as.factor(x),
         y = as.factor(y))



# converting the result to dataframe
#df <- as.data.frame(rev_data_frame)

(w <- ggplot(df, aes(x, y, width=1, height=1)) +
    geom_tile(aes(fill = estimate)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1, 2, 3),
                     labels=c(bquote(psi[0]),
                              bquote(psi[race]),
                              bquote(mu[psi["nat. green."]~"[native]"])
                              
                     )) +
    scale_y_discrete(name="Minimum population density threshold \n(people per sq. km)", breaks = rep(1:n_grains),
                     labels=c("1000", "1200", "1400")) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3"),
                         name="mean parameter \nestimate",
                         limits=c(-1,2), breaks = c(-1, 0, 1, 2)) +
    #geom_text(data = df, 
    #     aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
    
    geom_text(data = df, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "50% BCI: [", signif(lower_50,2), ", ", signif(upper_50,2), "]")),
              vjust = -1, size = 5) + 
    geom_text(data = df, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "95% BCI: [", signif(lower_95,2), ", ", signif(upper_95,2), "]")),
              vjust = 1, size = 5) +
    
    
    theme(legend.position = "right",
          legend.text=element_text(size=14),
          legend.title=element_text(size=16),
          axis.text.x = element_text(size = 24, angle = 45, hjust=1),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(size = 12),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank())
  
)
