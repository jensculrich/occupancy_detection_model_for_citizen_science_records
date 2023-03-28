# File to make results figures
# jcu, started dec 5, 2022.

library(tidyverse)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_15km_1000minpop_5minpersp_4ints_3visits_.RDS")
species_names <- readRDS("./figures/species_names/syrphidae_names_15km_urban.RDS")

list_of_draws <- as.data.frame(stan_out)

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
X_eco <- c(1, 2, 3) # 3 ecological params of interest
# mean of eco params
Y_eco <- c(#fit_summary$summary[6,1], # mu psi income
           fit_summary$summary[7,1], # mu psi herb shrub forest
           fit_summary$summary[9,1], # psi site area
           fit_summary$summary[1,1] # mu psi 0
)

# confidence intervals
lower_95_eco <- c(#fit_summary$summary[6,4], # mu psi income
                  fit_summary$summary[7,4], # mu psi herb shrub forest
                  fit_summary$summary[9,4],
                  fit_summary$summary[1,4] # mu psi 0
)

upper_95_eco <- c(#fit_summary$summary[6,8], # mu psi income
                  fit_summary$summary[7,8], # mu psi herb shrub forest
                  fit_summary$summary[9,8],
                  fit_summary$summary[1,8] # mu psi 0
)

# confidence intervals
lower_50_eco <- c(#fit_summary$summary[6,5], # mu psi income
                  fit_summary$summary[7,5], # mu psi herb shrub forest
                  fit_summary$summary[9,5],
                  fit_summary$summary[1,5] # mu psi 0
)

upper_50_eco <- c(#fit_summary$summary[6,7], # mu psi income
                  fit_summary$summary[7,7], # mu psi herb shrub forest
                  fit_summary$summary[9,7],
                  fit_summary$summary[1,7] # mu psi 0
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
  species_estimates[1,i] <- fit_summary$summary[23+i,1] # psi species
  species_estimates[2,i] <- NA # site area
  species_estimates[3,i] <- fit_summary$summary[89+i,1] # herb shrub forest
  species_estimates[4,i] <- fit_summary$summary[56+i,1] # income
  
}

## --------------------------------------------------
## Draw ecological parameter plot

(s <- ggplot(df_estimates_eco) +
   geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                 color="black",width=0.1,size=1,alpha=0.5) +
   geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                 color="black",width=0,size=3,alpha=0.8) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2, 3),
                    labels=c(#bquote(psi[income]),
                             bquote(psi[natural.habitat]),
                             bquote(psi[site.area]), 
                             bquote(psi[0]))) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-0.5, 5)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip()
)

df_estimates_eco_species <- cbind(df_estimates_eco, species_estimates)

for(i in 1:n_species){
  
  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,4+i])))
  #test[1,2] <- NA
  #test[2,2] <- NA
  #test[4,2] <- NA
  colnames(test) <- c("X_eco", "Y_eco")
  
  s <- s + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
                      col = "skyblue", size = 6, shape = "|", alpha = 0.5)
  
}

s <- s +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 

s


## --------------------------------------------------
## Plot detection paramter means and variation (citizen science)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

# parameter means
X_citsci <- c(1, 2, 3, 4, 5, 6) # 6 ecological params of interest

# mean of cit sci params and museum params
Y_citsci <- c(
  # museum 
  fit_summary$summary[152,1], # p pop dens
  fit_summary$summary[151,1], # p interval
  fit_summary$summary[126,1], # mu p 0
  # cit sci
  fit_summary$summary[125,1], # p pop dens
  fit_summary$summary[124,1], # p interval
  fit_summary$summary[99,1]) # mu p 0

# confidence intervals
lower_95_citsci <- c(
  # museum 
  fit_summary$summary[152,4], # p pop dens
  fit_summary$summary[151,4], # p interval
  fit_summary$summary[126,4], # mu p 0
  # cit sci
  fit_summary$summary[125,4], # p pop dens
  fit_summary$summary[124,4], # p interval
  fit_summary$summary[99,4]) # mu p 0

upper_95_citsci <- c(
  # museum 
  fit_summary$summary[152,8], # p pop dens
  fit_summary$summary[151,8], # p interval
  fit_summary$summary[126,8], # mu p 0
  # cit sci
  fit_summary$summary[125,8], # p pop dens
  fit_summary$summary[124,8], # p interval
  fit_summary$summary[99,8]) # mu p 0

df_estimates_citsci <- as.data.frame(cbind(X_citsci, Y_citsci, 
                                           lower_95_citsci, upper_95_citsci))
df_estimates_citsci$X_citsci <- as.factor(df_estimates_citsci$X_citsci)


## --------------------------------------------------
## Draw parameter plot

(q <- ggplot(df_estimates_citsci) +
   geom_errorbar(aes(x=X_citsci, ymin=lower_95_citsci, ymax=upper_95_citsci),
                 color="black",width=0.1,size=1,alpha=0.5) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2, 3, 4, 5, 6),
                    labels=c(
                      bquote(p.museum[pop.density]),
                      bquote(p.museum[time]), 
                      bquote(p.museum[0]),      
                      bquote(p.citsci[pop.density]),
                      bquote(p.citsci[time]), 
                      bquote(p.citsci[0]))) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-7, 3)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         plot.title = element_text(size = 20),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip() + 
   ggtitle("Observation process")
)

(q <- q +
    geom_point(aes(x=X_citsci, y=Y_citsci),
               size = 3, alpha = 0.7) 
)

df_estimates_eco_species <- cbind(df_estimates_citsci, species_estimates)

for(i in 1:n_species){
  
  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,4+i])))
  
  colnames(test) <- c("X_citsci", "Y_citsci")
  
  q <- q + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
                      col = "skyblue", size = 6, shape = "|", alpha = 0.5)
  
}

q

## --------------------------------------------------
## Panel p and q together
#gridExtra::grid.arrange(q, p, ncol=2)

## --------------------------------------------------
## And now let's make a species-specific table
# "psi_species[8]",
# "psi_species[36]",
# "psi_species[41]",
# "psi_species[12]",
# "psi_species[23]"


fit_summary <- rstan::summary(stan_out)

# parameter means
params = 4

x <- (rep(1:params, each=n_species)) # parameter reference
y = (rep(1:n_species, times=params)) # species reference


estimate <-  c(
  # param 1 (psi_species)
  fit_summary$summary[21:202,1],
  # param 2 (psi natural)
  fit_summary$summary[203:384,1],
  # param 3 (Freeman Tukey P cit sci)
  fit_summary$summary[385:566,1], 
  # param 4 (Freeman Tukey P museum)
  fit_summary$summary[567:748,1]
)



lower <- c(
  # param 1 (psi_species)
  fit_summary$summary[21:202,4],
  # param 2 (psi natural)
  fit_summary$summary[203:384,4],
  # param 3 (Freeman Tukey P cit sci)
  fit_summary$summary[385:566,4], 
  # param 4 (Freeman Tukey P museum)
  fit_summary$summary[567:748,4]
) 

upper <- c(
  # param 1 (psi_species)
  fit_summary$summary[21:202,8],
  # param 2 (psi natural)
  fit_summary$summary[203:384,8],
  # param 3 (Freeman Tukey P cit sci)
  fit_summary$summary[385:566,8], 
  # param 4 (Freeman Tukey P museum)
  fit_summary$summary[567:748,8]
) 

df = as.data.frame(cbind(x,y,estimate, lower, upper)) %>%
  mutate(x = as.factor(x),
         y = as.factor(y))

df_filtered <- df %>%
  mutate(y_num = as.integer(y)) %>%
  filter(y_num > 19)%>%
  filter(y_num < 62)

# species_names <- str_replace_all(species_names,'Bombus','B.')

ggplot(df_filtered, aes(x, y, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(1, 2, 3, 4),
                   labels=c(bquote(psi[species]),
                            bquote(psi[species[natural]]),
                            bquote(FTP[citsci]),
                            bquote(FTP[museum])
                   )) +
  scale_y_discrete(name="", breaks = rep(1:n_species),
                   labels=species_names) +
  scale_fill_gradient2(low = ("firebrick2")) +
  geom_text(data = df_filtered, 
            aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  #geom_text(data = df, 
  #      aes(x = x, y = y, label = paste0(signif(estimate, 2),
  #          "\n(", signif(lower,2), ", ", signif(upper,2), ")")),
  #         size = 6) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())
