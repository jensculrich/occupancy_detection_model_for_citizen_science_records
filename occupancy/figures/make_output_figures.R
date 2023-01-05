# File to make results figures
# jcu, started dec 5, 2022.

library(tidyverse)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./model_outputs/stan_out_model_integrated_ranges_200_35km_25records.rds")

list_of_draws <- as.data.frame(stan_out)

fit_summary <- rstan::summary(stan_out)

fit_summary$summary[1,] # top row (first monitored paramter from the fit)
fit_summary$summary[1,1] # parameter mean
fit_summary$summary[1,3] # parameter sd
fit_summary$summary[1,4] # 2.5 % CI
fit_summary$summary[1,8] # 97.5 % CI

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
# mu_psi_0 [1,]
# mu_psi_pop_density [94,]
# mu_psi_interval [96,]
# psi_site_area [98,]
# mu_p_citsci_0 [99,]
# p_citsci_interval [147,]
# p_citsci_pop_density [148,]
# mu_p_museum_0 [149,]
# p_museum_interval [197,]
# p_museum_pop_density[198,]

n_species <- 45

## --------------------------------------------------
## Plot ecological paramter means and variation

# parameter means
X_eco <- c(1, 2, 3, 4) # 4 ecological params of interest
Y_eco <- c(fit_summary$summary[98,1], 
           fit_summary$summary[96,1],
           fit_summary$summary[94,1],
           fit_summary$summary[1,1]) # mean of eco params

# confidence intervals
lower_95_eco <- c(fit_summary$summary[98,4], 
              fit_summary$summary[96,4],
              fit_summary$summary[94,4],
              fit_summary$summary[1,4])

upper_95_eco <- c(fit_summary$summary[98,8], 
              fit_summary$summary[96,8],
              fit_summary$summary[94,8],
              fit_summary$summary[1,8]) 

df_estimates_eco <- as.data.frame(cbind(X_eco, Y_eco, lower_95_eco, upper_95_eco))
df_estimates_eco$X_eco <- as.factor(df_estimates_eco$X_eco)

## --------------------------------------------------
## Get species specific estimates

# Not going to plot these here but this is how you would grab them if wanted..
species_estimates <- data.frame()

for(i in 1:n_species){
  
  species_estimates[1,i] <- fit_summary$summary[1+i,1] # psi species
  species_estimates[2,i] <- fit_summary$summary[48+i,1] # psi pop density
  ## species_estimates[2,i] <- fit_summary$summary[_+i,1] # didn't track psi interval!
  
}

## --------------------------------------------------
## Draw ecological parameter plot

(s <- ggplot(df_estimates_eco) +
    geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                  color="black",width=0.1,size=1,alpha=0.5) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_discrete(name="", breaks = c(1, 2, 3, 4),
                       labels=c(bquote(psi[site.area]), 
                                bquote(mu[psi[time.interval]]),
                                bquote(mu[psi[population.density]]), 
                                bquote(psi[0]))) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-3, 3.5)) +
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

(s <- s +
    geom_point(aes(x=X_eco, y=Y_eco),
               size = 5, alpha = 0.8) 
)

df_estimates_eco_species <- cbind(df_estimates_eco, species_estimates)

for(i in 1:n_species){
  
  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[2,4+i])))
  test[1,2] <- NA
  test[2,2] <- NA
  test[4,2] <- NA
  colnames(test) <- c("X_eco", "Y_eco")
  
  s <- s + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
                      col = "skyblue", size = 6, shape = "|", alpha = 0.5)
  
}

s

## --------------------------------------------------
## Plot detection paramter means and variation (citizen science)

# mu_p_citsci_0 [99,]
# p_citsci_interval [147,]
# p_citsci_pop_density [148,]
# mu_p_museum_0 [149,]
# p_museum_interval [197,]
# p_museum_pop_density[198,]

# parameter means
X_citsci <- c(1, 2, 3) # 4 ecological params of interest
Y_citsci <- c(fit_summary$summary[148,1], 
           fit_summary$summary[147,1],
           fit_summary$summary[99,1]) # mean of cit sci params

# confidence intervals
lower_95_citsci <- c(fit_summary$summary[148,4], 
                     fit_summary$summary[147,4],
                     fit_summary$summary[99,4])

upper_95_citsci <- c(fit_summary$summary[148,8], 
                     fit_summary$summary[147,8],
                     fit_summary$summary[99,8]) 

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
   scale_x_discrete(name="", breaks = c(1, 2, 3),
                    labels=c(bquote(p[population.density]),
                             bquote(p[time.interval]), 
                             bquote(p[0]))) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-5.5, 3)) +
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
   ggtitle("Citizen science \n detection process")
)

(q <- q +
    geom_point(aes(x=X_citsci, y=Y_citsci),
               size = 3, alpha = 0.7) 
)

## --------------------------------------------------
## Plot detection paramter means and variation (museum records)

# mu_p_citsci_0 [99,]
# p_citsci_interval [147,]
# p_citsci_pop_density [148,]
# mu_p_museum_0 [149,]
# p_museum_interval [197,]
# p_museum_pop_density[198,]

# parameter means
X_museum <- c(1, 2, 3) # 4 ecological params of interest
Y_museum <- c(fit_summary$summary[198,1], 
              fit_summary$summary[197,1],
              fit_summary$summary[149,1]) # mean of cit sci params

# confidence intervals
lower_95_museum <- c(fit_summary$summary[198,4], 
                     fit_summary$summary[197,4],
                     fit_summary$summary[149,4])

upper_95_museum <- c(fit_summary$summary[198,8], 
                     fit_summary$summary[197,8],
                     fit_summary$summary[149,8]) 

df_estimates_museum <- as.data.frame(cbind(X_museum, Y_museum, 
                                           lower_95_museum, upper_95_museum))
df_estimates_museum$X_museum <- as.factor(df_estimates_museum$X_museum)

## --------------------------------------------------
## Draw parameter plot

(p <- ggplot(df_estimates_museum) +
   geom_errorbar(aes(x=X_museum, ymin=lower_95_museum, ymax=upper_95_museum),
                 color="black",width=0.1,size=1,alpha=0.5) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2, 3),
                    labels=c(bquote(p[population.density]),
                             bquote(p[time.interval]), 
                             bquote(p[0]))) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-5.5, 3)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(
     legend.text=element_text(size=10),
     axis.text.x = element_text(size = 18),
     axis.text.y = element_text(size = 18),
     axis.title.x = element_text(size = 18),
     axis.title.y = element_text(size = 18),
     plot.title = element_text(size = 20),
     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip() + 
   ggtitle("Museum record \n detection process")
)

(p <- p +
    geom_point(aes(x=X_museum, y=Y_museum),
               size = 3, alpha = 0.7) 
)

## --------------------------------------------------
## Panel p and q together
gridExtra::grid.arrange(q, p, ncol=2)

## --------------------------------------------------
## And now let's make a species-specific table
# "psi_species[8]",
# "psi_species[36]",
# "psi_species[41]",
# "psi_species[12]",
# "psi_species[23]"

# parameter means
species = 5
params = 2

x <- (rep(1:params, each=species)) # parameter reference
y = (rep(1:species, times=params)) # species reference

estimate <- c(# param 1 (psi_species)
       fit_summary$summary[9,1], 
       fit_summary$summary[37,1],
       fit_summary$summary[42,1],
       fit_summary$summary[13,1],
       fit_summary$summary[24,1],
       
       # param 2 (psi pop dens species)
       fit_summary$summary[56,1], 
       fit_summary$summary[84,1],
       fit_summary$summary[89,1],
       fit_summary$summary[60,1],
       fit_summary$summary[71,1]
       ) 
 
lower <- c(# param 1 (psi_species)
  fit_summary$summary[9,4], 
  fit_summary$summary[37,4],
  fit_summary$summary[42,4],
  fit_summary$summary[13,4],
  fit_summary$summary[24,4],
  
  # param 2 (psi pop dens species)
  fit_summary$summary[56,4], 
  fit_summary$summary[84,4],
  fit_summary$summary[89,4],
  fit_summary$summary[60,4],
  fit_summary$summary[71,4]
) 

upper <- c(# param 1 (psi_species)
  fit_summary$summary[9,8], 
  fit_summary$summary[37,8],
  fit_summary$summary[42,8],
  fit_summary$summary[13,8],
  fit_summary$summary[24,8],
  
  # param 2 (psi pop dens species)
  fit_summary$summary[56,8], 
  fit_summary$summary[84,8],
  fit_summary$summary[89,8],
  fit_summary$summary[60,8],
  fit_summary$summary[71,8]
) 

df = as.data.frame(cbind(x,y,estimate, lower, upper)) %>%
  mutate(x = as.factor(x),
         y = as.factor(y))

ggplot(df, aes(x, y, width=1, height=1)) +
  geom_tile(aes(fill = estimate), colour = "grey50") +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(1, 2),
                   labels=c(bquote(psi[species]),
                            bquote(psi[species[pop.density]]))) +
  scale_y_discrete(name="", breaks = rep(1:species),
                   labels=c("Copestylum mexicanum",
                            "Pseudodoros clavatus",
                            "Syritta flaviventris",
                            "Eoseristalis hirta",
                            "Helophilus fasciatus")) +
  scale_fill_gradient2(low = ("firebrick2")) +
  geom_text(data = df, 
            aes(x = x, y = y, label = paste0(signif(estimate, 2),
                "\n(", signif(lower,2), ", ", signif(upper,2), ")")),
                size = 6) +
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())
