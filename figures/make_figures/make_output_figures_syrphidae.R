# File to make results figures
# jcu, started dec 5, 2022.

library(tidyverse)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits_.rds")
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
X_eco <- c(1, 2) # 3 ecological params of interest
# mean of eco params
Y_eco <- c(#fit_summary$summary[6,1], # mu psi income
           fit_summary$summary[7,1], # mu psi herb shrub forest
           fit_summary$summary[9,1]#, # psi site area
           #fit_summary$summary[1,1] # mu psi 0
)

# confidence intervals
lower_95_eco <- c(#fit_summary$summary[6,4], # mu psi income
                  fit_summary$summary[7,4], # mu psi herb shrub forest
                  fit_summary$summary[9,4]#,
                  #fit_summary$summary[1,4] # mu psi 0
)

upper_95_eco <- c(#fit_summary$summary[6,8], # mu psi income
                  fit_summary$summary[7,8], # mu psi herb shrub forest
                  fit_summary$summary[9,8]#,
                  #fit_summary$summary[1,8] # mu psi 0
)

# confidence intervals
lower_50_eco <- c(#fit_summary$summary[6,5], # mu psi income
                  fit_summary$summary[7,5], # mu psi herb shrub forest
                  fit_summary$summary[9,5]#,
                  #fit_summary$summary[1,5] # mu psi 0
)

upper_50_eco <- c(#fit_summary$summary[6,7], # mu psi income
                  fit_summary$summary[7,7], # mu psi herb shrub forest
                  fit_summary$summary[9,7]#,
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
  species_estimates[2,i] <- fit_summary$summary[136+i,1] # herb shrub forest
  #species_estimates[4,i] <- fit_summary$summary[56+i,1] # income
  
}

## --------------------------------------------------
## Draw ecological parameter plot

(s <- ggplot(df_estimates_eco) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = c(1, 2),
                    labels=c(#bquote(psi[income]),
                             bquote(psi[natural.habitat]),
                             bquote(psi[site.area])#, 
                             #bquote(psi[0])
                             )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-1, 1)) +
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
  
  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,6+i])))
  #test[1,2] <- NA
  #test[2,2] <- NA
  #test[4,2] <- NA
  colnames(test) <- c("X_eco", "Y_eco")
  
  s <- s + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
                      col = "skyblue", size = 6, shape = "|", alpha = 0.5)
  
}

s <- s +
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 

s


## --------------------------------------------------
## Plot detection paramter means and variation (citizen science)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

# parameter means
X_detection <- c(1, 2, 3, 4, 5) # 5 detection params of interest

# mean of cit sci params and museum params
Y_detection <- c(
  # museum 
  fit_summary$summary[22,1], # p total records
  fit_summary$summary[17,1], # mu p 0
  # cit sci
  fit_summary$summary[16,1], # p pop dens
  fit_summary$summary[15,1], # p interval
  fit_summary$summary[10,1]) # mu p 0

# confidence intervals
lower_95_detection <- c(
  # museum 
  fit_summary$summary[22,4], # p total records
  fit_summary$summary[17,4], # mu p 0
  # cit sci
  fit_summary$summary[16,4], # p pop dens
  fit_summary$summary[15,4], # p interval
  fit_summary$summary[10,4]) # mu p 0

upper_95_detection <- c(
  # museum 
  fit_summary$summary[22,8], # p total records
  fit_summary$summary[17,8], # mu p 0
  # cit sci
  fit_summary$summary[16,8], # p pop dens
  fit_summary$summary[15,8], # p interval
  fit_summary$summary[10,8]) # mu p 0

# confidence intervals
lower_50_detection <- c(
  # museum 
  fit_summary$summary[22,5], # p total records
  fit_summary$summary[17,5], # mu p 0
  # cit sci
  fit_summary$summary[16,5], # p pop dens
  fit_summary$summary[15,5], # p interval
  fit_summary$summary[10,5]) # mu p 0

upper_50_detection <- c(
  # museum 
  fit_summary$summary[22,7], # p total records
  fit_summary$summary[17,7], # mu p 0
  # cit sci
  fit_summary$summary[16,7], # p pop dens
  fit_summary$summary[15,7], # p interval
  fit_summary$summary[10,7]) # mu p 0

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
   scale_x_discrete(name="", breaks = c(1, 2, 3, 4, 5),
                    labels=c(
                      bquote(p.museum[records.collected]),
                      bquote(p.museum[0]),      
                      bquote(p.citsci[pop.density]),
                      bquote(p.citsci[time.interval^2]), 
                      bquote(p.citsci[0]))) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-6, 1.5)) +
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
   ggtitle("")
)


df_estimates_eco_species <- cbind(df_estimates_citsci, species_estimates_det)

for(i in 1:n_species){
  
  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,4+i])))
  
  colnames(test) <- c("X_citsci", "Y_citsci")
  
  q <- q + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
                      col = "skyblue", size = 6, shape = "|", alpha = 0.5)
  
}

(q <- q +
    geom_point(aes(x=X_detection, y=Y_detection),
               size = 3, alpha = 0.7) 
)


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

stan_out <- readRDS("./occupancy/model_outputs/syrphidae/syrphidae_10km_1200minpop_5minpersp_3ints_3visits_.RDS")
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

#stan_out2 <- readRDS("./occupancy/model_outputs/syrphidae/non_urban/syrphidae_40km_1000minpop_5minpersp_3ints_3visits_.RDS")
fit_summary2 <- rstan::summary(stan_out2)
View(cbind(1:nrow(fit_summary2$summary), fit_summary2$summary)) # View to see which row corresponds to the parameter of interest

species_names <- readRDS("./figures/species_names/syrphidae_names_10km_urban.RDS")
#species_names2 <- readRDS("./figures/species_names/syrphidae_names_40km_nonurban.RDS")

n_species <- length(species_names)


# parameter means
params = 2

x <- (rep(1:params, each=n_species)) # parameter reference
y = (rep(1:n_species, times=params)) # species reference

species_names_df <- species_names %>%
  as.data.frame(.)
#species_names2_df <- species_names2 %>%
 # as.data.frame(.)

#mismatches1 <- anti_join(species_names2_df, species_names_df,  by = ".")
#mismatches2 <- anti_join(species_names_df, species_names2_df,  by = ".") %>%
#  cbind(mismatch = "Y") 

#test <- left_join(species_names_df, mismatches2) %>%
#  rename("species" = ".") 

# there are 48 species in the urban study that were not recovered in the random resampling
#mismatch_rows <- which(test$mismatch == "Y")

# add an empty row to fit_summary2 with NAs whenever species was not recovered (mismatch rows)
#newrow <- rep(NA, ncol(fit_summary2$summary))

#temp <- fit_summary2$summary

#for(i in 1:length(mismatch_rows)){
#  r <- mismatch_rows[i] + 12
#  temp <- rbind(temp[1:r,],newrow,temp[-(1:r),])
#}

#View(cbind(1:nrow(temp), temp)) # View to see which row corresponds to the parameter of interest

estimate <-  c(
  # nativity
  # param 1 (psi_species_rangewide)
  #temp[14:154,1], 
  # param 1 (psi_species)
  fit_summary$summary[18:158,1],
  # param 2 (psi natural)
  fit_summary$summary[159:299,1]#,
  # param 3 (Freeman Tukey P cit sci)
  #fit_summary$summary[253:367,1], 
  # param 4 (Freeman Tukey P museum)
  #fit_summary$summary[368:482,1]
)

lower <-  c(
  # param 1 (psi_species_rangewide)
  #temp[14:154,4], 
  # param 1 (psi_species)
  fit_summary$summary[18:158,4],
  # param 2 (psi natural)
  fit_summary$summary[159:299,4]#,
  # param 3 (Freeman Tukey P cit sci)
  #fit_summary$summary[253:367,1], 
  # param 4 (Freeman Tukey P museum)
  #fit_summary$summary[368:482,1]
)

upper <-  c(
  # param 1 (psi_species_rangewide)
  #temp[14:154,8], 
  # param 1 (psi_species)
  fit_summary$summary[18:158,8],
  # param 2 (psi natural)
  fit_summary$summary[159:299,8]#,
  # param 3 (Freeman Tukey P cit sci)
  #fit_summary$summary[253:367,1], 
  # param 4 (Freeman Tukey P museum)
  #fit_summary$summary[368:482,1]
)

 

df = as.data.frame(cbind(x,y,estimate, lower, upper)) %>%
  mutate(x = as.factor(x),
         y = as.factor(y))

#df_filtered <- df %>%
#  mutate(y_num = as.integer(y)) %>%
#  filter(y_num > 0) %>%
#  filter(y_num < 36)

# parameter means
params = 1

x2 <- (rep(1:params, each=n_species)) # parameter reference
y2 = (rep(1:n_species, times=params)) # species reference


estimate2 <-  c(
  # param 3 (Freeman Tukey P cit sci)
  fit_summary$summary[812:952,1]
)

df2 = as.data.frame(cbind(x2,y2,estimate2)) %>%
  mutate(x2 = as.factor(x2),
         y2 = as.factor(y2))

library(cowplot)

for(i in 1:5){
  
  df_filtered <- df %>%
    mutate(y_num = as.integer(y)) %>%
    filter(y_num > 30*i - 30) %>%
    filter(y_num < 30*i)
  
  p1 <- ggplot(df_filtered, aes(x, y, width=1, height=1)) +
    geom_tile(aes(fill = estimate)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1, 2),
                     labels=c(#bquote(psi[species - range]),
                              bquote(psi[species]),
                              bquote(psi[species["nat. habitat"]])
                              #bquote(FTP[citsci]),
                              #bquote(FTP[museum])
                     )) +
    scale_y_discrete(name="", breaks = rep(1:n_species),
                     labels=species_names) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
      #geom_text(data = df_filtered, 
    #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
    
    geom_text(data = df_filtered, 
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
  
  df_filtered2 <- df2 %>%
    mutate(y_num = as.integer(y2)) %>%
    filter(y_num > 30*i - 30) %>%
    filter(y_num < 30*i)
  
  p2 <- ggplot(df_filtered2, aes(x2, y2, width=.8, height=1)) +
    geom_tile(aes(fill = estimate2)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1),
                     labels=c(
                       bquote(FTP[cs])
                     )) +
    scale_y_discrete(name="", breaks = "",
                     labels="") +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
    geom_text(data = df_filtered2, colour = "black",
              aes(x = x2, y = y2, label = signif(estimate2, 2)), size = 3.5) +
    theme(legend.position = "none",
          #legend.text=element_text(size=14),
          #legend.title=element_text(size=16),
          axis.text.x = element_text(size = 16, angle = 45, hjust=1),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          plot.title = element_text(size = 12),
          plot.margin = unit(c(0,-1,0,-1), "cm"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank())
  
  nativity_temp <- nativity
  df3 <- as.data.frame(cbind(rep(1, length(nativity_temp)), seq(1:length(nativity_temp)), nativity_temp)) %>%
    mutate(new = ifelse(nativity_temp==1, "native", "non-native"))
  
  df_filtered3 <- df3 %>%
    mutate(y_num = as.integer(V2)) %>%
    filter(y_num > 30*i - 30) %>%
    filter(y_num < 30*i)
  
  df_filtered4 <- cbind(df_filtered2, df_filtered3$nativity_temp, df_filtered3$new) %>% 
    rename("nativity_temp" = "df_filtered3$nativity_temp",
           "new" = "df_filtered3$new") %>% 
    mutate(nativity_temp = as.factor(nativity_temp))
  
  str(df_filtered4)
  
  p3 <- ggplot(df_filtered4, aes(x2, y2, width=.8, height=1)) +
    geom_tile(aes(fill = as.factor(nativity_temp)), colour = "white") +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1),
                     labels=c(
                       bquote("nativity")
                     )) +
    scale_y_discrete(name="", breaks = "",
                     labels="") +
    scale_fill_manual(values = c("0"="firebrick3","1"="dodgerblue3")) +
      #values = c("0" = "white", "1" = "#EDE284")) +
    #geom_text(data = df_filtered3, 
    #          aes(x = V1, y = V2, label = new),
    #          size = 3.5) +
    theme(legend.position = "none",
          #legend.text=element_text(size=14),
          #legend.title=element_text(size=16),
          #legend.title=element_text(size=16),
          axis.text.x = element_text(size = 16, angle = 45, hjust=1),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          plot.title = element_text(size = 12),
          plot.margin = unit(c(0,-1.5,0,-1.5), "cm"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank())
  

  plot_grid(p1, p3, p2, nrow = 1, align = "h", axis = "bt", rel_widths = c(2, 0.5, 0.5))
  
}



## --------------------------------------------------
## Prediction versus covariate

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

list_of_draws <- as.data.frame(stan_out)

# 95% conf int
plot(NA, xlim=c(-3,3), ylim=c(0,1))
curve(ilogit(fit_summary$summary[7,1]*x), add=TRUE, col = "blue", lwd = 3)
curve(ilogit(fit_summary$summary[7,4]*x), add=TRUE)
curve(ilogit(fit_summary$summary[7,8]*x), add=TRUE)

## --------------------------------------------------
## Natural habitat

# effect and all others held at mean
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,953]
}

plot(NA, xlim=c(-3,3), ylim=c(0,1),
     xlab = "Natural habitat area (scaled)",
     ylab = "Pr(Occurrence (native hoverlfy species))")

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[1,1] + # intercept
      fit_summary$summary[10,1] + # site area 
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[1,1] + # intercept
    # should add non-centered random effects
    fit_summary$summary[10,1] + # site area 
    fit_summary$summary[953,1]*x), 
  add=TRUE, col = "blue", lwd = 3)
