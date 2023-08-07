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

list_of_draws_10 <- as.data.frame(ten_km)
list_of_draws_15 <- as.data.frame(fifteen_km)
list_of_draws_20 <- as.data.frame(twenty_km)
list_of_draws_25 <- as.data.frame(twentyfive_km)

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
params = 3 # mu_psi_0, mu_psi_herb_shrub_forest, mu_psi_income
n_grains = 4 # 10, 15, 20, and 25km scales

x <- (rep(1:params, each=n_grains)) # parameter reference
y <- (rep(1:n_grains, times=params)) # spatial grain reference


estimate <-  c(
  # param 1 (psi_species)
  fit_summary_10$summary[68,1],
  fit_summary_15$summary[70,1],
  fit_summary_20$summary[70,1],
  fit_summary_25$summary[70,1],
  # param 2 (psi income)
  fit_summary_10$summary[73,1],
  fit_summary_15$summary[75,1],
  fit_summary_20$summary[75,1],
  fit_summary_25$summary[75,1],
  # param 3 (psi natural)
  fit_summary_10$summary[75,1],
  fit_summary_15$summary[77,1],
  fit_summary_20$summary[77,1],
  fit_summary_25$summary[77,1]
)

lower_95 <- c(
  # param 1 (psi_species)
  fit_summary_10$summary[68,4],
  fit_summary_15$summary[70,4],
  fit_summary_20$summary[70,4],
  fit_summary_25$summary[70,4],
  # param 2 (psi income)
  fit_summary_10$summary[73,4],
  fit_summary_15$summary[75,4],
  fit_summary_20$summary[75,4],
  fit_summary_25$summary[75,4],
  # param 3 (psi natural)
  fit_summary_10$summary[75,4],
  fit_summary_15$summary[77,4],
  fit_summary_20$summary[77,4],
  fit_summary_25$summary[77,4]
) 

upper_95 <- c(
  # param 1 (psi_species)
  fit_summary_10$summary[68,8],
  fit_summary_15$summary[70,8],
  fit_summary_20$summary[70,8],
  fit_summary_25$summary[70,8],
  # param 2 (psi income)
  fit_summary_10$summary[73,8],
  fit_summary_15$summary[75,8],
  fit_summary_20$summary[75,8],
  fit_summary_25$summary[75,8],
  # param 3 (psi natural)
  fit_summary_10$summary[75,8],
  fit_summary_15$summary[77,8],
  fit_summary_20$summary[77,8],
  fit_summary_25$summary[77,8]
) 

lower_50 <- c(
  # param 1 (psi_species)
  fit_summary_10$summary[68,5],
  fit_summary_15$summary[70,5],
  fit_summary_20$summary[70,5],
  fit_summary_25$summary[70,5],
  # param 2 (psi income)
  fit_summary_10$summary[73,5],
  fit_summary_15$summary[75,5],
  fit_summary_20$summary[75,5],
  fit_summary_25$summary[75,5],
  # param 3 (psi natural)
  fit_summary_10$summary[75,5],
  fit_summary_15$summary[77,5],
  fit_summary_20$summary[77,5],
  fit_summary_25$summary[77,5]
) 

upper_50 <- c(
  # param 1 (psi_species)
  fit_summary_10$summary[68,7],
  fit_summary_15$summary[70,7],
  fit_summary_20$summary[70,7],
  fit_summary_25$summary[70,7],
  # param 2 (psi income)
  fit_summary_10$summary[73,7],
  fit_summary_15$summary[75,7],
  fit_summary_20$summary[75,7],
  fit_summary_25$summary[75,7],
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
                            bquote(mu[psi["income"]]),
                            bquote(mu[psi["nat. habitat"]])
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
                "50% BCI: (", signif(lower_50,2), ", ", signif(upper_50,2), ")")),
              vjust = -1, size = 5) + 
    geom_text(data = df, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "95% BCI: (", signif(lower_95,2), ", ", signif(upper_95,2), ")")),
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

three <- readRDS("./occupancy/model_outputs/bombus/bombus_15km_1000minpop_10minpersp_3ints_4visits_.RDS")
four <- readRDS("./occupancy/model_outputs/bombus/bombus_15km_1000minpop_10minpersp_4ints_3visits_.RDS")
six <- readRDS("./occupancy/model_outputs/bombus/bombus_15km_1000minpop_10minpersp_6ints_2visits_.RDS")

fit_summary_three <- rstan::summary(three)
fit_summary_four <- rstan::summary(four)
fit_summary_six <- rstan::summary(six)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

## --------------------------------------------------
## Make figure table 

# parameter means
params = 3 # mu_psi_0, mu_psi_herb_shrub_forest, mu_psi_income
n_grains = 3 # 3 4 or 6 intervals (divided by 4, 3, or 2 years)

x <- (rep(1:params, each=n_grains)) # parameter reference
y <- (rep(1:n_grains, times=params)) # spatial grain reference


estimate <-  c(
  # param 1 (psi_species)
  fit_summary_three$summary[1,1],
  fit_summary_four$summary[1,1],
  fit_summary_six$summary[1,1],
  # param 2 (psi income)
  fit_summary_three$summary[6,1],
  fit_summary_four$summary[6,1],
  fit_summary_six$summary[6,1],
  # param 3 (psi natural)
  fit_summary_three$summary[8,1],
  fit_summary_four$summary[8,1],
  fit_summary_six$summary[8,1]
)

lower_95 <- c(
  # param 1 (psi_species)
  fit_summary_three$summary[1,4],
  fit_summary_four$summary[1,4],
  fit_summary_six$summary[1,4],
  # param 2 (psi income)
  fit_summary_three$summary[6,4],
  fit_summary_four$summary[6,4],
  fit_summary_six$summary[6,4],
  # param 3 (psi natural)
  fit_summary_three$summary[8,4],
  fit_summary_four$summary[8,4],
  fit_summary_six$summary[8,4]
) 

upper_95 <- c(
  # param 1 (psi_species)
  fit_summary_three$summary[1,8],
  fit_summary_four$summary[1,8],
  fit_summary_six$summary[1,8],
  # param 2 (psi income)
  fit_summary_three$summary[6,8],
  fit_summary_four$summary[6,8],
  fit_summary_six$summary[6,8],
  # param 3 (psi natural)
  fit_summary_three$summary[8,8],
  fit_summary_four$summary[8,8],
  fit_summary_six$summary[8,8]
) 

lower_50 <- c(
  # param 1 (psi_species)
  fit_summary_three$summary[1,5],
  fit_summary_four$summary[1,5],
  fit_summary_six$summary[1,5],
  # param 2 (psi income)
  fit_summary_three$summary[6,5],
  fit_summary_four$summary[6,5],
  fit_summary_six$summary[6,5],
  # param 3 (psi natural)
  fit_summary_three$summary[8,5],
  fit_summary_four$summary[8,5],
  fit_summary_six$summary[8,5]
) 

upper_50 <- c(
  # param 1 (psi_species)
  fit_summary_three$summary[1,7],
  fit_summary_four$summary[1,7],
  fit_summary_six$summary[1,7],
  # param 2 (psi income)
  fit_summary_three$summary[6,7],
  fit_summary_four$summary[6,7],
  fit_summary_six$summary[6,7],
  # param 3 (psi natural)
  fit_summary_three$summary[8,7],
  fit_summary_four$summary[8,7],
  fit_summary_six$summary[8,7]
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
                              bquote(psi[species[income]]),
                              bquote(psi[species[natural.hab]])
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
                "50% BCI: (", signif(lower_50,2), ", ", signif(upper_50,2), ")")),
              vjust = -1, size = 5) + 
    geom_text(data = df, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "95% BCI: (", signif(lower_95,2), ", ", signif(upper_95,2), ")")),
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
