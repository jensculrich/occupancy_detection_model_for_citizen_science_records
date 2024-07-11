# File to make sensitivity analysis figures and tables
# jcu, started mar 28, 2023.

library(tidyverse)

## --------------------------------------------------
## Sensitivity to spatial grain

## --------------------------------------------------
## Read in model run results

ten_km <- readRDS("./occupancy/model_outputs/bombus/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits_.rds")
fifteen_km <- readRDS("./occupancy/model_outputs/bombus/bombus_15km_1000minpop_1minUniqueDetections_4ints_3visits_.rds")
twenty_km <- readRDS("./occupancy/model_outputs/bombus/bombus_20km_800minpop_1minUniqueDetections_4ints_3visits_.rds")
twentyfive_km <- readRDS("./occupancy/model_outputs/bombus/bombus_25km_600minpop_1minUniqueDetections_4ints_3visits_.rds")

fit_summary_10 <- rstan::summary(ten_km)
fit_summary_15 <- rstan::summary(fifteen_km)
fit_summary_20 <- rstan::summary(twenty_km)
fit_summary_25 <- rstan::summary(twentyfive_km)

fit_summary$summary[1,] # top row (first monitored paramter from the fit)
fit_summary$summary[1,1] # parameter mean
fit_summary$summary[1,3] # parameter sd
fit_summary$summary[1,4] # 2.5 % CI
fit_summary$summary[1,8] # 97.5 % CI

View(cbind(1:nrow(fit_summary_10$summary), fit_summary_10$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_15$summary), fit_summary_15$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_20$summary), fit_summary_20$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_25$summary), fit_summary_25$summary)) # View to see which row corresponds to the parameter of interest


## --------------------------------------------------
## Make figure table 

# parameter means
params = 3 # mu_psi_0, mu_psi_herb_shrub_forest, psi race
n_grains = 4 # 10, 15, 20, and 25km scales

x <- (rep(1:params, each=n_grains)) # parameter reference
y <- (rep(1:n_grains, times=params)) # spatial grain reference


estimate <-  c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[68,1],
  fit_summary_15$summary[70,1],
  fit_summary_20$summary[70,1],
  fit_summary_25$summary[70,1],
  # param 2 (psi race)
  fit_summary_10$summary[74,1],
  fit_summary_15$summary[76,1],
  fit_summary_20$summary[76,1],
  fit_summary_25$summary[76,1],
  # param 3 (psi natural)
  fit_summary_10$summary[75,1],
  fit_summary_15$summary[77,1],
  fit_summary_20$summary[77,1],
  fit_summary_25$summary[77,1]
)

lower_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[68,4],
  fit_summary_15$summary[70,4],
  fit_summary_20$summary[70,4],
  fit_summary_25$summary[70,4],
  # param 2 (psi race)
  fit_summary_10$summary[74,4],
  fit_summary_15$summary[76,4],
  fit_summary_20$summary[76,4],
  fit_summary_25$summary[76,4],
  # param 3 (psi natural)
  fit_summary_10$summary[75,4],
  fit_summary_15$summary[77,4],
  fit_summary_20$summary[77,4],
  fit_summary_25$summary[77,4]
) 

upper_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[68,8],
  fit_summary_15$summary[70,8],
  fit_summary_20$summary[70,8],
  fit_summary_25$summary[70,8],
  # param 2 (psi race)
  fit_summary_10$summary[74,8],
  fit_summary_15$summary[76,8],
  fit_summary_20$summary[76,8],
  fit_summary_25$summary[76,8],
  # param 3 (psi natural)
  fit_summary_10$summary[75,8],
  fit_summary_15$summary[77,8],
  fit_summary_20$summary[77,8],
  fit_summary_25$summary[77,8]
) 

lower_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[68,5],
  fit_summary_15$summary[70,5],
  fit_summary_20$summary[70,5],
  fit_summary_25$summary[70,5],
  # param 2 (psi race)
  fit_summary_10$summary[74,5],
  fit_summary_15$summary[76,5],
  fit_summary_20$summary[76,5],
  fit_summary_25$summary[76,5],
  # param 3 (psi natural)
  fit_summary_10$summary[75,5],
  fit_summary_15$summary[77,5],
  fit_summary_20$summary[77,5],
  fit_summary_25$summary[77,5]
) 

upper_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_10$summary[68,7],
  fit_summary_15$summary[70,7],
  fit_summary_20$summary[70,7],
  fit_summary_25$summary[70,7],
  # param 2 (psi income)
  fit_summary_10$summary[74,7],
  fit_summary_15$summary[76,7],
  fit_summary_20$summary[76,7],
  fit_summary_25$summary[76,7],
  # param 3 (psi natural)
  fit_summary_10$summary[75,7],
  fit_summary_15$summary[77,7],
  fit_summary_20$summary[77,7],
  fit_summary_25$summary[77,7]
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
                            bquote(psi["race"]),
                            bquote(mu[psi["nat. green."]])
                            #bquote(FTP[citsci]),
                            #bquote(FTP[museum])
                   )) +
    scale_y_discrete(name="Spatial grain \n (landscape width (km))", breaks = rep(1:n_grains),
                     labels=c("10", "15", "20", "25")) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3"),
                         name="mean parameter \nestimate",
                         limits=c(-2,1), breaks = c(-2, -1, 0, 1)) +
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

three <- readRDS("./occupancy/model_outputs/bombus/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits_.rds")
two <- readRDS("./occupancy/model_outputs/bombus/bombus_10km_1200minpop_1minUniqueDetections_6ints_2visits_.rds")
four <- readRDS("./occupancy/model_outputs/bombus/bombus_10km_1200minpop_1minUniqueDetections_3ints_4visits_.rds")

fit_summary_three <- rstan::summary(three)
fit_summary_two <- rstan::summary(two)
fit_summary_four <- rstan::summary(four)

View(cbind(1:nrow(fit_summary_two$summary), fit_summary_two$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_three$summary), fit_summary_three$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary_four$summary), fit_summary_four$summary)) # View to see which row corresponds to the parameter of interest

## --------------------------------------------------
## Make figure table 

# parameter means
params = 3 # mu_psi_0, mu_psi_herb_shrub_forest, psi race
n_grains = 3 # 3 4 or 6 intervals (divided by 4, 3, or 2 years)

x <- (rep(1:params, each=n_grains)) # parameter reference
y <- (rep(1:n_grains, times=params)) # spatial grain reference


estimate <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[68,1],
  fit_summary_three$summary[68,1],
  fit_summary_four$summary[68,1],
  # param 3 (psi race)
  fit_summary_two$summary[74,1],
  fit_summary_three$summary[74,1],
  fit_summary_four$summary[74,1],
  # param 2 (psi natural)
  fit_summary_two$summary[75,1],
  fit_summary_three$summary[75,1],
  fit_summary_four$summary[75,1]
)

lower_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[68,4],
  fit_summary_three$summary[68,4],
  fit_summary_four$summary[68,4],
  # param 3 (psi race)
  fit_summary_two$summary[74,4],
  fit_summary_three$summary[74,4],
  fit_summary_four$summary[74,4],
  # param 2 (psi natural)
  fit_summary_two$summary[75,4],
  fit_summary_three$summary[75,4],
  fit_summary_four$summary[75,4]
) 

upper_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[68,8],
  fit_summary_three$summary[68,8],
  fit_summary_four$summary[68,8],
  # param 3 (psi race)
  fit_summary_two$summary[74,8],
  fit_summary_three$summary[74,8],
  fit_summary_four$summary[74,8],
  # param 2 (psi natural)
  fit_summary_two$summary[75,8],
  fit_summary_three$summary[75,8],
  fit_summary_four$summary[75,8]
) 

lower_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[68,5],
  fit_summary_three$summary[68,5],
  fit_summary_four$summary[68,5],
  # param 3 (psi race)
  fit_summary_two$summary[74,5],
  fit_summary_three$summary[74,5],
  fit_summary_four$summary[74,5],
  # param 2 (psi natural)
  fit_summary_two$summary[75,5],
  fit_summary_three$summary[75,5],
  fit_summary_four$summary[75,5]
) 

upper_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_two$summary[68,7],
  fit_summary_three$summary[68,7],
  fit_summary_four$summary[68,7],
  # param 3 (psi race)
  fit_summary_two$summary[74,7],
  fit_summary_three$summary[74,7],
  fit_summary_four$summary[74,7],
  # param 2 (psi natural)
  fit_summary_two$summary[75,7],
  fit_summary_three$summary[75,7],
  fit_summary_four$summary[75,7]
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
                              bquote(psi["race"]),
                              bquote(mu[psi["nat. habitat"]])
                              #bquote(FTP[citsci]),
                              #bquote(FTP[museum])
                     )) +
    scale_y_discrete(name="Years per occupancy interval", breaks = rep(1:n_grains),
                     labels=c("4", "3", "2")) +
    scale_fill_gradient2(low = ("firebrick2"), 
                         name="mean parameter \nestimate",
                         limits=c(-1,1)) +
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
          axis.text.x = element_text(size = 18, angle = 45, hjust=1),
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

low <- readRDS("./occupancy/model_outputs/bombus/bombus_10km_1000minpop_1minUniqueDetections_4ints_3visits_.rds")
medium <- readRDS("./occupancy/model_outputs/bombus/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits_.rds")
high <- readRDS("./occupancy/model_outputs/bombus/bombus_10km_1400minpop_1minUniqueDetections_4ints_3visits_.rds")

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
  fit_summary_low$summary[70,1],
  fit_summary_medium$summary[68,1],
  fit_summary_high$summary[68,1],
  
  # param 3 (psi race)
  fit_summary_low$summary[76,1],
  fit_summary_medium$summary[74,1],
  fit_summary_high$summary[74,1],
  
  # param 3 (psi natural)
  fit_summary_low$summary[77,1],
  fit_summary_medium$summary[75,1],
  fit_summary_high$summary[75,1]
)

lower_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[70,4],
  fit_summary_medium$summary[68,4],
  fit_summary_high$summary[68,4],
  
  # param 3 (psi race)
  fit_summary_low$summary[76,4],
  fit_summary_medium$summary[74,4],
  fit_summary_high$summary[74,4],
  
  # param 3 (psi natural)
  fit_summary_low$summary[77,4],
  fit_summary_medium$summary[75,4],
  fit_summary_high$summary[75,4]
) 

upper_95 <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[70,8],
  fit_summary_medium$summary[68,8],
  fit_summary_high$summary[68,8],
  
  # param 3 (psi race)
  fit_summary_low$summary[76,8],
  fit_summary_medium$summary[74,8],
  fit_summary_high$summary[74,8],
  
  # param 3 (psi natural)
  fit_summary_low$summary[77,8],
  fit_summary_medium$summary[75,8],
  fit_summary_high$summary[75,8]
) 

lower_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[70,5],
  fit_summary_medium$summary[68,5],
  fit_summary_high$summary[68,5],
  
  # param 3 (psi race)
  fit_summary_low$summary[76,5],
  fit_summary_medium$summary[74,5],
  fit_summary_high$summary[74,5],
  
  # param 3 (psi natural)
  fit_summary_low$summary[77,5],
  fit_summary_medium$summary[75,5],
  fit_summary_high$summary[75,5]
) 

upper_50 <- c(
  # param 1 (mu_psi_0)
  fit_summary_low$summary[70,7],
  fit_summary_medium$summary[68,7],
  fit_summary_high$summary[68,7],
  
  # param 3 (psi race)
  fit_summary_low$summary[76,7],
  fit_summary_medium$summary[74,7],
  fit_summary_high$summary[74,7],
  
  # param 3 (psi natural)
  fit_summary_low$summary[77,7],
  fit_summary_medium$summary[75,7],
  fit_summary_high$summary[75,7]
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
                              bquote(mu[psi["nat. green."]])
                              
                     )) +
    scale_y_discrete(name="Minimum population density threshold \n(people per sq. km)", breaks = rep(1:n_grains),
                     labels=c("1000", "1200", "1400")) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3"),
                         name="mean parameter \nestimate",
                         limits=c(-1,1.25), breaks = c(-1, -0.5, 0, 0.5, 1)) +
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
