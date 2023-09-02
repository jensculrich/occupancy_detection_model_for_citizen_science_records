# File to make results figures
# jcu, started dec 5, 2022.

library(tidyverse)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/large_files/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits.rds")
species_names <- readRDS("./figures/species_names/bombus_names_10km_urban.RDS")

#list_of_draws <- as.data.frame(stan_out)

fit_summary <- rstan::summary(stan_out)

fit_summary$summary[1,] # top row (first monitored paramter from the fit)
fit_summary$summary[1,1] # parameter mean
fit_summary$summary[1,3] # parameter sd
fit_summary$summary[1,4] # 2.5 % CI
fit_summary$summary[1,8] # 97.5 % CI

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

n_species <- length(species_names)

## --------------------------------------------------
## Plot ecological paramter means and variation

# parameter means
X_eco <- c(1, 2, 3) # 4 ecological params of interest
# mean of eco params
Y_eco <- c(fit_summary$summary[73,1], # mu psi income
           fit_summary$summary[75,1], # mu psi natural habitat
           fit_summary$summary[77,1]#, # psi site area
           #fit_summary$summary[1,1] # mu psi 0
          )

# confidence intervals
lower_95_eco <- c(fit_summary$summary[73,4], # mu psi income
                  fit_summary$summary[75,4], # mu psi natural habitat
                  fit_summary$summary[77,4]#, # psi site area
                  #fit_summary$summary[1,4] # mu psi 0
)

upper_95_eco <- c(fit_summary$summary[73,8], # mu psi income
                  fit_summary$summary[75,8], # mu psi natural habitat
                  fit_summary$summary[77,8]#, # psi site area
                  #fit_summary$summary[1,8] # mu psi 0
)

# confidence intervals
lower_50_eco <- c(fit_summary$summary[73,5], # mu psi income
                  fit_summary$summary[75,5], # mu psi natural habitat
                  fit_summary$summary[77,5]#, # psi site area
                  #fit_summary$summary[1,5] # mu psi 0
)

upper_50_eco <- c(fit_summary$summary[73,7], # mu psi income
                  fit_summary$summary[75,7], # mu psi natural habitat
                  fit_summary$summary[77,7]#, # psi site area
                  #fit_summary$summary[1,7] # mu psi 0
)


df_estimates_eco <- as.data.frame(cbind(X_eco, Y_eco, 
                                        lower_95_eco, upper_95_eco,
                                        lower_50_eco, upper_50_eco))

df_estimates_eco$X_eco <- as.factor(df_estimates_eco$X_eco)

## --------------------------------------------------
## Get species specific estimates

species_estimates <- data.frame()

for(i in 1:n_species){
  
  # row is one before the row of the first species estimate
  #species_estimates[1,i] <- NA # psi species
  #species_estimates[1,i] <- fit_summary$summary[23+i,1] # psi species
  species_estimates[1,i] <- NA # site area
  species_estimates[2,i] <- fit_summary$summary[152+i,1] # herb shrub forest
  species_estimates[3,i] <- fit_summary$summary[119+i,1] # income
  
}

## --------------------------------------------------
## Draw ecological parameter plot

(s <- ggplot(df_estimates_eco) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_discrete(name="", breaks = c(1, 2, 3),
                       labels=c(bquote(mu[psi["income"]]),
                                bquote(mu[psi["nat. habitat"]]),
                                bquote(psi["site area"])#, 
                                #bquote(mu[psi[0]])
                                )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-1.5, 1.75)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 24, angle=45, vjust=-0.5),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
)

df_estimates_eco_species <- cbind(df_estimates_eco, species_estimates)

for(i in 1:n_species){
  
  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,6+i])))
  #test[1,2] <- NA
  #test[2,2] <- NA
  #test[4,2] <- NA
  colnames(test) <- c("X_eco", "Y_eco")
  
  s <- s + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
                      col = "skyblue", size = 6, shape = "|", alpha = 0.75)
  
}

s <- s +
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
    geom_point(aes(x=X_eco, y=Y_eco),
               size = 5, alpha = 0.8) 


s

# do the syrphids too so we can plot them side by side
## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/large_files/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits.rds")
species_names <- readRDS("./figures/species_names/syrphidae_names_10km_urban.RDS")
nativity <- readRDS("./figures/species_names/syrphidae_nativity_10km_urban.RDS")


fit_summary <- rstan::summary(stan_out)

fit_summary$summary[1,] # top row (first monitored paramter from the fit)
fit_summary$summary[1,1] # parameter mean
fit_summary$summary[1,3] # parameter sd
fit_summary$summary[1,4] # 2.5 % CI
fit_summary$summary[1,8] # 97.5 % CI

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

n_species <- length(species_names)


## --------------------------------------------------
## Plot ecological paramter means and variation

# parameter means
X_eco <- c(1, 2, 3, 4, 5) # 3 ecological params of interest
# mean of eco params
Y_eco <- c(#fit_summary$summary[6,1], # mu psi income
  fit_summary$summary[10,1], # mu psi income
  fit_summary$summary[946,1], # mu psi natural habitat non native
  fit_summary$summary[945,1], # mu psi natural habitat native
  fit_summary$summary[947,1], # mu psi natural habitat all
  fit_summary$summary[9,1]#, # psi site area
  
)

# confidence intervals
lower_95_eco <- c(#fit_summary$summary[6,4], # mu psi income
  fit_summary$summary[10,4], # mu psi income
  fit_summary$summary[946,4], # mu psi natural habitat non native
  fit_summary$summary[945,4], # mu psi natural habitat native
  fit_summary$summary[947,4], # mu psi natural habitat all
  fit_summary$summary[9,4]#, # psi site area
)

upper_95_eco <- c(#fit_summary$summary[6,8], # mu psi income
  fit_summary$summary[10,8], # mu psi income
  fit_summary$summary[946,8], # mu psi natural habitat non native
  fit_summary$summary[945,8], # mu psi natural habitat native
  fit_summary$summary[947,8], # mu psi natural habitat all
  fit_summary$summary[9,8]#, # psi site area
)

# confidence intervals
lower_50_eco <- c(#fit_summary$summary[6,5], # mu psi income
  fit_summary$summary[10,5], # mu psi income
  fit_summary$summary[946,5], # mu psi natural habitat non native
  fit_summary$summary[945,5], # mu psi natural habitat native
  fit_summary$summary[947,5], # mu psi natural habitat all
  fit_summary$summary[9,5]#, # psi site area
)

upper_50_eco <- c(#fit_summary$summary[6,7], # mu psi income
  fit_summary$summary[10,7], # mu psi income
  fit_summary$summary[946,7], # mu psi natural habitat non native
  fit_summary$summary[945,7], # mu psi natural habitat native
  fit_summary$summary[947,7], # mu psi natural habitat all
  fit_summary$summary[9,7]#, # psi site area
)


df_estimates_eco <- as.data.frame(cbind(X_eco, Y_eco, 
                                        lower_95_eco, upper_95_eco,
                                        lower_50_eco, upper_50_eco))

df_estimates_eco$X_eco <- as.factor(df_estimates_eco$X_eco)

## --------------------------------------------------
## Get species specific estimates

subset_nonnative <- fit_summary$summary[158:298,] %>%
  cbind(., nativity) %>%
  as.data.frame(.) %>%
  filter(nativity == "0")

subset_native <- fit_summary$summary[158:298,] %>%
  cbind(., nativity) %>%
  as.data.frame(.) %>%
  filter(nativity == "1")

species_estimates <- data.frame()

for(i in 1:n_species){
  
  # row is one before the row of the first species estimate
  species_estimates[1,i] <- NA # site area
  species_estimates[2,i] <- fit_summary$summary[157+i,1] # all
  species_estimates[3,i] <- subset_native[i,1]
  species_estimates[4,i] <- subset_nonnative[i,1]
  species_estimates[5,i] <- NA # income
}

## --------------------------------------------------
## Draw ecological parameter plot

(s2 <- ggplot(df_estimates_eco) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2, 3, 4, 5),
                    labels=c(bquote(psi["income"]),
                      bquote(mu[psi["nat. habitat"]~"[nonnative]"]),
                      bquote(mu[psi["nat. habitat"]~"[native]"]),
                      bquote(mu[psi["nat. habitat"]~"[all]"]),
                      bquote(psi["site area"])#, 
                      #bquote(psi[0])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-1.5, 1.75)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 24, angle=45, vjust=-0.5),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip()
)

df_estimates_eco_species <- cbind(df_estimates_eco, species_estimates)

for(i in 1:n_species){
  
  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,6+i])))
  #test[1,2] <- NA
  #test[2,2] <- NA
  #test[4,2] <- NA
  colnames(test) <- c("X_eco", "Y_eco")
  
  s2 <- s2 + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
                      col = "skyblue", size = 6, shape = "|", alpha = 0.75)
  
}

s2 <- s2 +
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 

s2

library(cowplot)
plot_grid(s, s2, labels = c('a)', 'b)'),
          label_size = 18,
          rel_widths = c(1,1.2))

## --------------------------------------------------
## Plot detection parameter means and variation (community science)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/large_files/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits.rds")
species_names <- readRDS("./figures/species_names/bombus_names_10km_urban.RDS")

fit_summary <- rstan::summary(stan_out)

fit_summary$summary[1,] # top row (first monitored paramter from the fit)
fit_summary$summary[1,1] # parameter mean
fit_summary$summary[1,3] # parameter sd
fit_summary$summary[1,4] # 2.5 % CI
fit_summary$summary[1,8] # 97.5 % CI

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

n_species <- length(species_names)

# parameter means
X_detection <- c(1, 2, 3) # 4 detection params of interest

# mean of cit sci params and museum params
Y_detection <- c(
           fit_summary$summary[67,1],# rho
           # rc 
           #fit_summary$summary[18,1], # mu p 0
           # cs
           fit_summary$summary[83,1], # p pop dens
           fit_summary$summary[82,1]#, # p interval
           #fit_summary$summary[78,1] # mu p cs 0
          )

# confidence intervals
lower_95_detection <- c(
  fit_summary$summary[67,4], # rho
  # rc 
  #fit_summary$summary[18,4], # mu p 0
  # cs
  fit_summary$summary[83,4], # p pop dens
  fit_summary$summary[82,4] # p interval
  #fit_summary$summary[11,4] # mu p 0 
  ) 

upper_95_detection <- c(
  fit_summary$summary[67,8], # rho
  # rc 
  #fit_summary$summary[18,8], # mu p 0
  # cs
  fit_summary$summary[83,8], # p pop dens
  fit_summary$summary[82,8] # p interval
  #fit_summary$summary[11,8] # mu p 0 
) 

# confidence intervals
lower_50_detection <- c(
  fit_summary$summary[67,5], # rho
  # rc 
  #fit_summary$summary[18,5], # mu p 0
  # cs
  fit_summary$summary[83,5], # p pop dens
  fit_summary$summary[82,5] # p interval
  #fit_summary$summary[11,5] # mu p 0 
) 

upper_50_detection <- c(
  fit_summary$summary[67,7], # rho
  # rc 
  #fit_summary$summary[18,7], # mu p 0
  # cs
  fit_summary$summary[83,7], # p pop dens
  fit_summary$summary[82,7] # p interval
  #fit_summary$summary[11,7] # mu p 0 
) 

df_estimates_detection <- as.data.frame(cbind(X_detection, Y_detection, 
                                           lower_95_detection, upper_95_detection,
                                           lower_50_detection, upper_50_detection))

df_estimates_detection$X_detection <- as.factor(
  df_estimates_detection$X_detection)

## --------------------------------------------------
## Get species specific estimates

species_estimates_det <- data.frame()

for(i in 1:n_species){
  
  # row is one before the row of the first species estimate
  #species_estimates_det[1,i] <- fit_summary$summary[21+i,1] # mu p 0 citsci
  species_estimates_det[2,i] <- NA # p interval
  #species_estimates_det[3,i] <- fit_summary$summary[103+i,1] # p pop dens
  #species_estimates_det[4,i] <- fit_summary$summary[103+i,1] # mu p 0 mus
  species_estimates_det[5,i] <- NA # p total records
  
}


## --------------------------------------------------
## Draw parameter plot

(q <- ggplot(df_estimates_detection) +
   geom_errorbar(aes(x=X_detection, ymin=lower_95_detection, ymax=upper_95_detection),
                 color="black",width=0.1,size=1,alpha=0.5) +
   geom_errorbar(aes(x=X_detection, ymin=lower_50_detection, ymax=upper_50_detection),
                 color="black",width=0,size=3,alpha=0.8) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2, 3),
                    labels=c(
                      bquote(rho),
                      #bquote("p.rc"["total records"]),
                      #bquote(rc[0]),      
                      bquote(p.cs["pop. density"]),
                      bquote(p.cs["interval"^2]) 
                      #bquote(p.cs[0]))
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-0.5, 1.25), breaks = seq(-1, 1, by = 1)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 24, angle=45, vjust=-0.5),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         plot.title = element_text(size = 20),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip() + 
   ggtitle("")
)


df_estimates_eco_species <- cbind(df_estimates_citsci, species_estimates_det)

#for(i in 1:n_species){
  
#  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,4+i])))
  
#  colnames(test) <- c("X_citsci", "Y_citsci")
  
#  q <- q + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
#                      col = "skyblue", size = 6, shape = "|", alpha = 0.5)
  
#}

(q <- q +
    geom_point(aes(x=X_detection, y=Y_detection),
               size = 5, alpha = 0.8) 
)

## syrphidae
## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/large_files/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits.rds")
species_names <- readRDS("./figures/species_names/syrphidae_names_10km_urban.RDS")

fit_summary <- rstan::summary(stan_out)

fit_summary$summary[1,] # top row (first monitored paramter from the fit)
fit_summary$summary[1,1] # parameter mean
fit_summary$summary[1,3] # parameter sd
fit_summary$summary[1,4] # 2.5 % CI
fit_summary$summary[1,8] # 97.5 % CI

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

n_species <- length(species_names)

## --------------------------------------------------
## Plot detection paramter means and variation (citizen science)

# parameter means
X_detection <- c(1, 2) # 2 detection params of interest

# mean of cit sci params and museum params
Y_detection <- c(
  # cs
  fit_summary$summary[16,1], # p pop dens
  fit_summary$summary[15,1]#, # p interval
  #fit_summary$summary[78,1] # mu p cs 0
)

# confidence intervals
lower_95_detection <- c(
  # cs
  fit_summary$summary[16,4], # p pop dens
  fit_summary$summary[15,4] # p interval
  #fit_summary$summary[11,4] # mu p 0 
) 

upper_95_detection <- c(
  # cs
  fit_summary$summary[16,8], # p pop dens
  fit_summary$summary[15,8] # p interval
  #fit_summary$summary[11,8] # mu p 0 
) 

# confidence intervals
lower_50_detection <- c(
  # cs
  fit_summary$summary[16,5], # p pop dens
  fit_summary$summary[15,5] # p interval
  #fit_summary$summary[11,5] # mu p 0 
) 

upper_50_detection <- c(
  #fit_summary$summary[18,7], # mu p 0
  # cs
  fit_summary$summary[16,7], # p pop dens
  fit_summary$summary[15,7] # p interval
  #fit_summary$summary[11,7] # mu p 0 
) 

df_estimates_detection <- as.data.frame(cbind(X_detection, Y_detection, 
                                              lower_95_detection, upper_95_detection,
                                              lower_50_detection, upper_50_detection))

df_estimates_detection$X_detection <- as.factor(
  df_estimates_detection$X_detection)

## --------------------------------------------------
## Get species specific estimates

species_estimates_det <- data.frame()

for(i in 1:n_species){
  
  # row is one before the row of the first species estimate
  #species_estimates_det[1,i] <- fit_summary$summary[21+i,1] # mu p 0 citsci
  species_estimates_det[2,i] <- NA # p interval
  #species_estimates_det[3,i] <- fit_summary$summary[103+i,1] # p pop dens
  #species_estimates_det[4,i] <- fit_summary$summary[103+i,1] # mu p 0 mus
  species_estimates_det[5,i] <- NA # p total records
  
}


## --------------------------------------------------
## Draw parameter plot

(q2 <- ggplot(df_estimates_detection) +
   geom_errorbar(aes(x=X_detection, ymin=lower_95_detection, ymax=upper_95_detection),
                 color="black",width=0.1,size=1,alpha=0.5) +
   geom_errorbar(aes(x=X_detection, ymin=lower_50_detection, ymax=upper_50_detection),
                 color="black",width=0,size=3,alpha=0.8) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2),
                    labels=c(
                      bquote(p.cs["pop. density"]),
                      bquote(p.cs["interval"^2]) 
                      )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-0.5, 1.25), breaks = seq(-1, 1, by = 1)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 24, angle=45, vjust=-0.5),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         plot.title = element_text(size = 20),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip() + 
   ggtitle("")
)


#df_estimates_eco_species <- cbind(df_estimates_citsci, species_estimates_det)

#for(i in 1:n_species){
  
#  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,4+i])))
  
#  colnames(test) <- c("X_citsci", "Y_citsci")
  
#  q2 <- q2 + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
#                      col = "skyblue", size = 6, shape = "|", alpha = 0.5)
  
#}

(q2 <- q2 +
    geom_point(aes(x=X_detection, y=Y_detection),
               size = 5, alpha = 0.8) 
)

## --------------------------------------------------
## Panel p and q together

plot_grid(q, q2, labels = c('c)', 'd)'),
          label_size = 18,
          rel_widths = c(1,1.2))



#gridExtra::grid.arrange(q, p, ncol=2)

## --------------------------------------------------
## And now let's make a species-specific table
# "psi_species[8]",
# "psi_species[36]",
# "psi_species[41]",
# "psi_species[12]",
# "psi_species[23]"

stan_out <- readRDS("./occupancy/model_outputs/large_files/bombus_10km_1200minpop_1minUniqueDetections_4ints_3visits.rds")
species_names <- readRDS("./figures/species_names/bombus_names_10km_urban.RDS")
n_species <- length(species_names)

fit_summary <- rstan::summary(stan_out)

#stan_out2 <- readRDS("./occupancy/model_outputs/bombus/non_urban/bombus_40km_1000minpop_5minpersp_4ints_3visits_.RDS")
#fit_summary2 <- rstan::summary(stan_out2)
#species_names2 <- readRDS("./figures/species_names/bombus_names_40km_nonurban.RDS")

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
#View(cbind(1:nrow(fit_summary2$summary), fit_summary2$summary)) # View to see which row corresponds to the parameter of interest

species_names_df <- species_names %>%
  as.data.frame(.)
#species_names2_df <- species_names2 %>%
#  as.data.frame(.)

#mismatches <- anti_join(species_names2_df, species_names_df,  by = ".")
#mismatch_rows <- which(species_names2_df == "insularis")

# parameter means
params = 3

x <- (rep(1:params, each=n_species)) # parameter reference
y <- (rep(1:n_species, times=params)) # species reference

estimate <-  c(
  # param 1 (psi_species_rangewide)
  #fit_summary2$summary[87:103,1], fit_summary2$summary[105:119,1],
  # param 1 (psi_species_urban)
  rev(fit_summary$summary[88:119,1]),
  # param 2 (psi income)
  rev(fit_summary$summary[120:151,1]),
  # param 3 (psi natural)
  rev(fit_summary$summary[152:183,1])
)

lower <- c(
  # param 1 (psi_species_rangewide)
  #fit_summary2$summary[87:103,4], fit_summary2$summary[105:119,4],
  # param 1 (psi_species_urban)
  rev(fit_summary$summary[88:119,4]),
  # param 2 (psi income)
  rev(fit_summary$summary[120:151,4]),
  # param 3 (psi natural)
  rev(fit_summary$summary[152:183,4])
)

upper <- c(
  # param 1 (psi_species_rangewide)
  #fit_summary2$summary[87:103,8], fit_summary2$summary[105:119,8],
  # param 1 (psi_species_urban)
  rev(fit_summary$summary[88:119,8]),
  # param 2 (psi income)
  rev(fit_summary$summary[120:151,8]),
  # param 3 (psi natural)
  rev(fit_summary$summary[152:183,8])
) 

# species_names_label <- str_replace_all(species_names,'Bombus','B.')
species_names_label <- paste0("B. ", species_names) 

df <- as.data.frame(cbind(x, y, estimate, lower, upper 
                          )) %>%
  mutate(x = as.factor(x),
         y = as.factor(y))

# converting the result to dataframe
#df <- map_df(df, rev)
ordered_df <- filter(df, x == 3)
ordered_df <- ordered_df[order(ordered_df$estimate,decreasing=TRUE),]

df_intercepts <- filter(df, x == 1)
ordered_df <- rbind(ordered_df, df_intercepts)

species_names_label <- species_names_label %>%
  as.data.frame(.) %>%
  map_df(., rev) %>%
  pull(.)

p1 <- ggplot(df, aes(x, y, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(1, 2, 3),
                   labels=c(#bquote(psi[species - range]),
                            bquote(psi[species]),
                            bquote(psi[species[income]]),
                            bquote(psi[species["nat. habitat"]])
                            #bquote(FTP[citsci]),
                            #bquote(FTP[museum])
                            )) +
  scale_y_discrete(name="", breaks = rep(1:n_species),
                   labels=species_names_label) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df, 
   #     aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = df, 
            aes(x = x, y = y, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
                size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

# plot on separate colour scales
temp <- df %>% filter(x == 1)
p1.1 <- ggplot(temp, aes(x, y, width=1, height=1)) +
  geom_tile(aes(fill = estimate, width = 0.5)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(1),
                   labels=c(#bquote(psi[species - range]),
                     bquote(psi[species])
                   )) +
  scale_y_discrete(name="", breaks = rep(1:n_species),
                   labels=species_names_label) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df, 
  #     aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = temp, 
            aes(x = x, y = y, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

temp <- df %>% filter(x == 2)
p1.2 <- ggplot(temp, aes(x, y, width=1, height=1)) +
  geom_tile(aes(fill = estimate, width = 0.5)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(2),
                   labels=c(#bquote(psi[species - range]),
                     bquote(psi[species])
                     #bquote(psi[species[income]])
                   )) +
  scale_y_discrete(name="", breaks = rep(1:n_species),
                   labels=species_names_label) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df, 
  #     aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = temp, 
            aes(x = x, y = y, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

temp <- df %>% filter(x == 3)
p1.3 <- ggplot(temp, aes(x, y, width=1, height=1)) +
  geom_tile(aes(fill = estimate, width = 0.5)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(3),
                   labels=c(#bquote(psi[species - range]),
                     bquote(psi[species])
                     
                     #bquote(psi[species["nat. habitat"]])
                   )) +
  scale_y_discrete(name="", breaks = rep(1:n_species),
                   labels=species_names_label) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df, 
  #     aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = temp, 
            aes(x = x, y = y, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())



# parameter means
params = 2

x2 <- (rep(1:params, each=n_species)) # parameter reference
y2 = (rep(1:n_species, times=params)) # species reference


estimate2 <-  c(
  # param 4 (Freeman Tukey P cit sci)
  fit_summary$summary[760:791,1], 
  # param 5 (Freeman Tukey P rc)
  fit_summary$summary[856:887,1]
)

df2 = as.data.frame(cbind(x2,y2,estimate2)) %>%
  mutate(x2 = as.factor(x2),
         y2 = as.factor(y2))

p2 <- ggplot(df2, aes(x2, y2, width=1, height=1)) +
  geom_tile(aes(fill = estimate2)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(1, 2),
                   labels=c(
                            bquote(FTP[cs]),
                            bquote(FTP[rc])
                   )) +
  scale_y_discrete(name="", breaks = rep(1:n_species),
                   labels=species_names_label) +
  # or no species labels
  #scale_y_discrete(name="", breaks = "",
    #               labels="") +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  geom_text(data = df2, colour = "white",
       aes(x = x2, y = y2, label = signif(estimate2, 2)), size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

library(cowplot)
plot_grid(p1, p2, align = "h", axis = "bt", rel_widths = c(1, .5))


## --------------------------------------------------
# Make an ordered plot to match the hoverflies
# species_names_label <- str_replace_all(species_names,'Bombus','B.')
species_names_label <- paste0("B. ", species_names) 

df <- as.data.frame(cbind(x, y, estimate, lower, upper 
)) %>%
  mutate(x = as.factor(x),
         y = as.factor(y))

# converting the result to dataframe
#df <- map_df(df, rev)
ordered_df <- filter(df, x == 3)
ordered_df <- ordered_df[order(ordered_df$estimate,decreasing=FALSE),]

df_intercepts <- filter(df, x == 1)
ordered_df <- rbind(ordered_df, df_intercepts)

species_names_label <- species_names_label %>%
  as.data.frame(.) %>%
  #map_df(., rev) %>%
  pull(.)

species_names_df <- species_names_label %>%
  #rev(.) %>%
  as.data.frame(.) %>%
  mutate(row_id=row_number()) %>%
  mutate(row_id = as.factor(row_id)) %>%
  rename("y" = "row_id",
         "species_name" = ".")

ordered_df <- left_join(ordered_df, species_names_df, by = "y")

# make a column for the slopes
temp <- ordered_df %>% filter(x == 3) %>%
  mutate(row_id=row_number()) %>%
  mutate(row_id = as.factor(row_id))

# grab row id's for later to match rows
rows <- temp %>%
  select(species_name, row_id)

p1.2 <- ggplot(temp, aes(x, row_id, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(2),
                   labels=c(#bquote(psi[species - range]),
                     bquote(psi[species])
                     #bquote(psi[species["nat. habitat"]])
                     #bquote(FTP[citsci]),
                     #bquote(FTP[museum])
                   )) +
  scale_y_discrete(name="", breaks = rep(1:nrow(temp)),
                   labels=temp$species_name) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df_filtered, 
  #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = temp, 
            aes(x = x, y = row_id, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())


# make a column for the intercepts
temp2 <- ordered_df %>% filter(x == 1)
temp2 <- left_join(temp2, rows)

p1.1 <- ggplot(temp2, aes(x, row_id, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(2),
                   labels=c(#bquote(psi[species - range]),
                     bquote(psi[species])
                     #bquote(psi[species["nat. habitat"]])
                     #bquote(FTP[citsci]),
                     #bquote(FTP[museum])
                   )) +
  scale_y_discrete(name="", breaks = rep(1:nrow(temp)),
                   labels=temp$species_name) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df_filtered, 
  #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = temp2, 
            aes(x = x, y = row_id, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())


## --------------------------------------------------
## Prediction versus covariate

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

plot(NA, xlim=c(-3,3), ylim=c(0,1))
curve(ilogit(fit_summary$summary[77,1]*x), add=TRUE, col = "blue", lwd = 3)
curve(ilogit(fit_summary$summary[77,4]*x), add=TRUE)
curve(ilogit(fit_summary$summary[77,8]*x), add=TRUE)

## --------------------------------------------------
## Natural habitat

# effect and all others held at mean
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,75]
}

plot(NA, xlim=c(-3,3), ylim=c(0,1),
     xlab = "Natural habitat area (scaled)",
     ylab = "Pr(Occurrence)")

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[68,1] + # intercept
      # should add non-centered random effects
    fit_summary$summary[73,1] + # income  
    fit_summary$summary[77,1] + # site area 
      params[i,1]*x), 
        add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[68,1] + # intercept
    # should add non-centered random effects
    fit_summary$summary[73,1] + # income  
    fit_summary$summary[77,1] + # site area 
    fit_summary$summary[75,1]*x), 
      add=TRUE, col = "blue", lwd = 3)



## --------------------------------------------------
## Income

# effect and all others held at mean
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,73]
}

plot(NA, xlim=c(-3,3), ylim=c(0,1),
     xlab = "Household Income (scaled)",
     ylab = "Pr(Occurence)")

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[68,1] + # intercept
      # should add non-centered random effects
      fit_summary$summary[75,1] + # natural habitat  
      fit_summary$summary[77,1] + # site area 
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[68,1] + # intercept
    # should add non-centered random effects
    fit_summary$summary[75,1] + # nat hab  
    fit_summary$summary[77,1] + # site area 
    fit_summary$summary[73,1]*x), 
  add=TRUE, col = "blue", lwd = 3)


## --------------------------------------------------
## Site Area

# effect and all others held at mean
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,77]
}

plot(NA, xlim=c(-3,3), ylim=c(0,1),
     xlab = "Site area (scaled)",
     ylab = "Pr(Occupancy)")

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[68,1] + # intercept
      # should add non-centered random effects
      fit_summary$summary[75,1] + # natural habitat  
      fit_summary$summary[73,1] + # income
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[68,1] + # intercept
    # should add non-centered random effects
    fit_summary$summary[73,1] + # income  
    fit_summary$summary[75,1] + # nat hab
    fit_summary$summary[77,1]*x), 
  add=TRUE, col = "blue", lwd = 3)

# just effect
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,10]
}

plot(NA, xlim=c(-3,3), ylim=c(0,1),
     xlab = "Site area (scaled)",
     ylab = "Pr(Occupancy)")

for(i in 1:n_lines){
  curve(ilogit(params[i,1]*x), 
        add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(fit_summary$summary[10,1]*x), 
      add=TRUE, col = "blue", lwd = 3)


## --------------------------------------------------
## Detection - Time on detection


# just effect
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,82]
}

plot(NA, xlim=c(0,4), ylim=c(0,1),
     xlab = "Time interval",
     ylab = "Pr(Detection)")

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[78,1] + # intercept
      # should add non-centered random effects
      fit_summary$summary[83,1] + # pop dens
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[78,1] + # intercept
    # should add non-centered random effects
    fit_summary$summary[83,1] + # pop dens
    fit_summary$summary[82,1]*x),
  
      add=TRUE, col = "blue", lwd = 2)

## --------------------------------------------------
## Detection - Pop density on detection


# just effect
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,83]
}

plot(NA, xlim=c(-3,3), ylim=c(0,1),
     xlab = "Population density (scaled)",
     ylab = "Pr(Detection)")

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[78,1] + # intercept
      # should add non-centered random effects
      fit_summary$summary[82,1] + # time
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[78,1] + # intercept
    # should add non-centered random effects
    fit_summary$summary[82,1] + # time
    fit_summary$summary[83,1]*x),
  
  add=TRUE, col = "blue", lwd = 2)

