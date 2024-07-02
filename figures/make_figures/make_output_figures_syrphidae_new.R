# File to make results figures
# jcu, started dec 5, 2022.

library(tidyverse)

## --------------------------------------------------
## Read in model run results

stan_out <- readRDS("./occupancy/model_outputs/large_files/syrphidae_10km_1200minpop_2minUniqueDetections_3ints_3visits_.rds")
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

#stan_out2 <- readRDS("./occupancy/model_outputs/syrphidae/non_urban/syrphidae_40km_1000minpop_5minpersp_3ints_3visits_.RDS")
#fit_summary2 <- rstan::summary(stan_out2)
#View(cbind(1:nrow(fit_summary2$summary), fit_summary2$summary)) # View to see which row corresponds to the parameter of interest

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

estimate <- c(
  # nativity
  # param 1 (psi_species_rangewide)
  #temp[14:154,1], 
  # param 1 (psi_species)
  rev(fit_summary$summary[309:449,1]),
  # param 2 (psi natural)
  rev(fit_summary$summary[450:590,1])#,
  # param 3 (Freeman Tukey P cit sci)
  #fit_summary$summary[253:367,1], 
  # param 4 (Freeman Tukey P museum)
  #fit_summary$summary[368:482,1]
)

lower <-  c(
  # param 1 (psi_species_rangewide)
  #temp[14:154,4], 
  # param 1 (psi_species)
  rev(fit_summary$summary[309:449,4]),
  # param 2 (psi natural)
  rev(fit_summary$summary[450:590,4])#,
  # param 3 (Freeman Tukey P cit sci)
  #fit_summary$summary[253:367,1], 
  # param 4 (Freeman Tukey P museum)
  #fit_summary$summary[368:482,1]
)

upper <-  c(
  # param 1 (psi_species_rangewide)
  #temp[14:154,8], 
  # param 1 (psi_species)
  rev(fit_summary$summary[309:449,8]),
  # param 2 (psi natural)
  rev(fit_summary$summary[450:590,8])#,
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
  rev(fit_summary$summary[805:945,1])
)

df2 = as.data.frame(cbind(x2,y2,estimate2)) %>%
  mutate(x2 = as.factor(x2),
         y2 = as.factor(y2))

library(cowplot)

# flip species names
species_names_label <- species_names %>%
  as.data.frame(.) %>%
  map_df(., rev) %>%
  pull(.)

num_per_page = 30 

for(i in 1:5){
  
  df_filtered <- df %>%
    mutate(y_num = as.integer(y)) %>%
    filter(y_num > num_per_page*i - num_per_page) %>%
    filter(y_num < num_per_page*i)
  
  p1 <- ggplot(df_filtered, aes(x, y, width=1, height=1)) +
    geom_tile(aes(fill = estimate)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1, 2),
                     labels=c(#bquote(psi[species - range]),
                              bquote(psi["species"]~"[species]"),
                              bquote(psi["nat. habitat"]~"[species]")
                              #bquote(FTP[citsci]),
                              #bquote(FTP[museum])
                     )) +
    scale_y_discrete(name="", breaks = rep(1:n_species),
                     labels=species_names_label) +
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
  
  temp <- df_filtered %>% filter(x == 1)
  p1.1 <- ggplot(temp, aes(x, y, width=1, height=1)) +
    geom_tile(aes(fill = estimate)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1),
                     labels=c(#bquote(psi[species - range]),
                       bquote(psi["species"]~"[species]")
                       #bquote(psi[species["nat. habitat"]])
                       #bquote(FTP[citsci]),
                       #bquote(FTP[museum])
                     )) +
    scale_y_discrete(name="", breaks = rep(1:n_species),
                     labels=species_names_label) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
    #geom_text(data = df_filtered, 
    #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
    
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
  
  temp <- df_filtered %>% filter(x == 2)
  p1.2 <- ggplot(temp, aes(x, y, width=1, height=1)) +
    geom_tile(aes(fill = estimate)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(2),
                     labels=c(#bquote(psi[species - range]),
                       bquote(psi["nat. habitat"]~"[species]")
                       #bquote(psi[species["nat. habitat"]])
                       #bquote(FTP[citsci]),
                       #bquote(FTP[museum])
                     )) +
    scale_y_discrete(name="", breaks = rep(1:n_species),
                     labels=species_names_label) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
    #geom_text(data = df_filtered, 
    #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
    
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
  
  temp <- df_filtered %>% filter(x == 2)
  p1.3 <- ggplot(temp, aes(x, y, width=1, height=1)) +
    geom_tile(aes(fill = estimate)) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(2),
                     labels=c(#bquote(psi[species - range]),
                       bquote(psi["nat. habitat"]~"[species]")
                       #bquote(psi[species["nat. habitat"]])
                       #bquote(FTP[citsci]),
                       #bquote(FTP[museum])
                     )) +
    #scale_y_discrete(name="", breaks = rep(1:n_species),
     #                labels=species_names_label) +
    scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
    #geom_text(data = df_filtered, 
    #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
    
    geom_text(data = temp, 
              aes(x = x, y = y, label = paste0(
                #signif(estimate, 2),"\n(", 
                "[", signif(lower,2), ", ", signif(upper,2), "]")),
              size = 3.5) +
    theme(legend.position = "none",
          #legend.text=element_text(size=14),
          #legend.title=element_text(size=16),
          axis.text.x = element_text(size = 16, angle = 45, hjust=1),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 12),
          axis.title.y =element_blank(),
          plot.title = element_text(size = 12),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank())
  
  df_filtered2 <- df2 %>%
    mutate(y_num = as.integer(y2)) %>%
    filter(y_num > num_per_page*i - num_per_page) %>%
    filter(y_num < num_per_page*i)
  
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
    geom_text(data = df_filtered2, colour = "white",
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
  
  nativity_temp <- rev(nativity)
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
  

  plot_grid(p1.1, p1.3, p3, p2, nrow = 1, align = "h", axis = "bt", rel_widths = c(1, 0.5, 0.3, 0.3))
  
}

## --------------------------------------------------
## make a table with top and bottom (winners and losers)

df_filtered_nathabitat <- filter(df, x == 2)
top <- top_n(df_filtered_nathabitat, 16, estimate)  
bottom <- top_n(df_filtered_nathabitat, -16, estimate) 

top <- top[order(top$estimate,decreasing=TRUE),]
bottom <- bottom[order(bottom$estimate,decreasing=TRUE),]

# join occurrence intercepts
df_intercepts <- filter(df, x == 1)
test <- filter(df_intercepts, y %in% top$y) 
top <- rbind(top, test)

test <- filter(df_intercepts, y %in% bottom$y) 
bottom <- rbind(bottom, test)

# join nativity
nativity_by_speciesID <- as.data.frame(cbind(nativity, rev(y)))
nativity_by_speciesID$y <- as.factor(nativity_by_speciesID$V2)

top <- left_join(top, nativity_by_speciesID, by = "y") 
top <- distinct(top) %>% select(-V2)

bottom <- left_join(bottom, nativity_by_speciesID, by = "y") 
bottom <- distinct(bottom) %>% select(-V2)

top_and_bottom <- rbind(top, bottom) %>% map_df(., rev)

# join species names
species_names_df <- species_names %>%
  rev(.) %>%
  as.data.frame(.) %>%
  mutate(row_id=row_number()) %>%
  mutate(row_id = as.factor(row_id)) %>%
  rename("y" = "row_id",
         "species_name" = ".")

# add nativity symbol
species_names_df <- left_join(species_names_df, nativity_by_speciesID) %>%
  distinct() %>%
  mutate(species_name = ifelse(nativity == 0, paste0(species_name, "*"), species_name))

top_and_bottom <- left_join(top_and_bottom, species_names_df, by = "y")

# make a single plot
p1 <- ggplot(top_and_bottom, aes(x, y, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(1, 2),
                   labels=c(#bquote(psi[species - range]),
                     bquote(psi["species"]~"[species]"),
                     bquote(psi["nat. habitat"]~"[species]")
                     #bquote(FTP[citsci]),
                     #bquote(FTP[museum])
                   )) +
  scale_y_discrete(name="", breaks = rep(1:nrow(species_names_df)),
                   labels=species_names_df$species_name) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df_filtered, 
  #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = top_and_bottom, 
            aes(x = x, y = y, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11, face="italic"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

# make a column for the slopes
temp <- top_and_bottom %>% filter(x == 2) %>%
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
                     #bquote(psi[species])
                     bquote(psi[italic("nat. green."~"[species]")])
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
              "[", sprintf("%.1f",lower), ", ",
              #"[", signif(lower,2), ", ", 
              sprintf("%.1f",upper), "]")),
            #signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

# row order needs to match row order of slope column
temp2 <- top_and_bottom %>% filter(x == 1)
temp2 <- left_join(temp2, rows)

p1.1 <- ggplot(temp2, aes(x, row_id, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(2),
                   labels=c(#bquote(psi[species - range]),
                     bquote(psi["species"]~"[species]")
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
              "[", sprintf("%.1f",lower), ", ",
              #"[", signif(lower,2), ", ", 
              sprintf("%.1f",upper), "]")),
              #signif(upper,2), "]")),
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

list_of_draws <- as.data.frame(stan_out)

# 95% conf int
plot(NA, xlim=c(-3,3), ylim=c(0,1))
curve(ilogit(fit_summary$summary[1335,1]*x), add=TRUE, col = "blue", lwd = 3)
curve(ilogit(fit_summary$summary[1335,4]*x), add=TRUE, col = "black", lwd = 3, lty="dashed")
curve(ilogit(fit_summary$summary[1335,8]*x), add=TRUE, col = "black", lwd = 3, lty="dashed")

## --------------------------------------------------
## Natural habitat

# effect and all others held at mean
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,1335]
}

library(scales)
percs <- runif(100)
yticks_val <- pretty_breaks(n=5)(percs)

plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "Natural greenspace area (scaled)",
     ylab = "Pr(Occurrence (native hoverlfy species))", yaxt="n")
axis(2, at=yticks_val, lab=percent(yticks_val))

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[286,1] + # intercept
      fit_summary$summary[295,1] + # site area 
      fit_summary$summary[298,1] + # mean dev open
      fit_summary$summary[296,1] + # mean income effect  
      fit_summary$summary[297,1] + # mean racial composition effect  
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[286,1] + # intercept
    fit_summary$summary[295,1] + # site area 
    fit_summary$summary[298,1] + # mean dev open
    fit_summary$summary[296,1] + # mean income effect  
    fit_summary$summary[297,1] + # mean racial composition effect  
    fit_summary$summary[1335,1]*x), 
  add=TRUE, col = "blue", lwd = 3)

low <- ilogit(
  fit_summary$summary[286,1] + # intercept
    fit_summary$summary[295,1] + # site area 
    fit_summary$summary[298,1] + # mean dev open
    fit_summary$summary[296,1] + # mean income effect  
    fit_summary$summary[297,1] + # mean racial composition effect  
    (fit_summary$summary[1335,1]*-2))

high <- ilogit(
  fit_summary$summary[286,1] + # intercept
    fit_summary$summary[295,1] + # site area 
    fit_summary$summary[298,1] + # mean dev open
    fit_summary$summary[296,1] + # mean income effect  
    fit_summary$summary[297,1] + # mean racial composition effect 
    (fit_summary$summary[1335,1]*2))

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[286,1] + # intercept
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}

curve(ilogit(
  fit_summary$summary[286,1] + # intercept
    fit_summary$summary[1335,1]*x), 
  add=TRUE, col = "blue", lwd = 3)

low <- ilogit(
  fit_summary$summary[286,1] + # intercept  
    (fit_summary$summary[1335,1]*-2))

high <- ilogit(
  fit_summary$summary[286,1] + # intercept
    (fit_summary$summary[1335,1]*2))

## --------------------------------------------------
## income

# effect and all others held at mean
n_lines <- 100
params <- matrix(nrow = n_lines, ncol = 1)

for(i in 1:n_lines){
  row <- sample(1:nrow(list_of_draws), 1)
  params[i,1] <- list_of_draws[row,10]
}

plot(NA, xlim=c(-2,2), ylim=c(0,1),
     xlab = "Household income (scaled)",
     ylab = "Pr(Occurrence(all hoverfly species))")

for(i in 1:n_lines){
  curve(ilogit(
    fit_summary$summary[1,1] + # intercept
      fit_summary$summary[9,1] + # site area 
      fit_summary$summary[947,1] + # mean natural habitat effect (all species)
      params[i,1]*x), 
    add=TRUE, col = "lightgrey", lwd = 1)
}


curve(ilogit(
  fit_summary$summary[1,1] + # intercept
    # should add non-centered random effects
    fit_summary$summary[9,1] + # site area 
    fit_summary$summary[947,1] + # mean income effect 
    fit_summary$summary[10,1]*x), 
  add=TRUE, col = "blue", lwd = 3)

