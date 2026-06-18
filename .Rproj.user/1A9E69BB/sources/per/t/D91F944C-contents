library(msm)
library(circular)
library(ggplot2)
library(zoo)
library(plyr)


#### 1. LOAD IN DATA ####

df <- read.csv("./Data_inputs/juv_depart_time_series.csv")
df$date <- as.Date(strptime(df$date, format = "%d/%m/%Y"))
df <- df[order(df$id, df$date), ]

#### 2. SETTING UP VARIABLES ####

# msm that requires states to take integer values (starting at 1)
df$state <- df$depart
df$state[df$state == 1] <- 2
df$state[df$state == 0] <- 1

# converting wind direction radians and getting cos and sin components
# from -180 to 180 to 0 to 360
hist(df$wd)
df$wd_360 <- df$wd
df$wd_360[df$wd_360 < 0] <- df$wd_360[df$wd_360 < 0] + 360
# converting to radians 
deg_to_rad <- function(deg) {
  return(deg * pi / 180)
}
df$wd_rad <- deg_to_rad(df$wd_360)
# getting cos and sin components
df$sin_rad <- sin(df$wd_rad)
df$cos_rad <- cos(df$wd_rad)

# transition matrix from data
statetable.msm(state,id,data=df) 
# manually giving initial values for the non-diagonal of the transition intensity matrix based on the values from the line above
Q_juv <- rbind(c(0,0.05),
               c(0,0)) 

#### 3. RUNNING BEST MARKOV MODELS FOR EACH SPECIES ####

##### 3A. BBA #####

bba <- subset(df, species == "BBA",)
bba  <- bba[order(bba$id, bba$age), ]
bba_best <-  msm(state ~ age, subject=id, data = bba,
               qmatrix = Q_juv, covariates = ~ age*mass_st)

##### 3B. GHA #####

gha <- subset(df, species == "GHA",)
gha  <- gha[order(gha$id, gha$age), ]
gha_best <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad+ ws + age*mass_st)

##### 3C. WA #####

wa <- subset(df, species == "WALB",)
wa  <- wa[order(wa$id, wa$age), ]
wa_best <-  msm(state ~ age, subject=id, data = wa,
              qmatrix = Q_juv, covariates = ~ age*mass_st)


##### 3D. WCP #####

wcp <- subset(df, species == "WCP",)
wcp  <- wcp[order(wcp$id, wcp$age), ]
wcp_best <-  msm(state ~ age, subject=id, data = wcp,
               qmatrix = Q_juv, covariates = ~ age)


#### 4. CREATING PREDICTION DATAFRAMES ####

##### 4A. BBA - AGE+MASS #####

hist(bba$mass)
mass_levels <- c(3.25, 3.5, 3.75) # predict at three levels of bird mass
mean_mass <- mean(bba$mass) 
sd_mass <- sd(bba$mass)
mass_st_levels <- (mass_levels - mean_mass) / sd_mass

# create prediction grid
ages <- seq(min(bba$age)-0.4, max(bba$age)+0.4, by = 0.1)
pred_grid <- expand.grid(
  age = ages,
  mass_st = mass_st_levels
)

# predict transition probabilities - takes a while!
pred_list <- lapply(1:nrow(pred_grid), function(i) {
  print(paste(i, nrow(pred_grid), sep = "..."))
  age_i <- pred_grid$age[i]
  mass_i <- pred_grid$mass_st[i]
  
  P <- pmatrix.msm(
    bba_best, t = 5,
    covariates = list(age = age_i, mass_st = mass_i),
    ci = "normal",
    cl = 0.95
  )
  
  c(prob  = P$est[1,2],
    lower = P$L[1,2],
    upper = P$U[1,2])
})

preds_bba <- do.call(rbind, pred_list)
preds_bba <- cbind(pred_grid, preds_bba)
# add original mass levels for plotting
preds_bba$mass <- rep(mass_levels, each = length(ages))

# average predictions
preds_bba$lower2 <- rollmean(preds_bba$lower, k = 15, fill = NA)
preds_bba$upper2 <- rollmean(preds_bba$upper, k = 15, fill = NA)
preds_bba$age2 <- round(preds_bba$age)
preds_bba_avg <- ddply(preds_bba,.(age2, mass), summarize, prob_mean = mean(prob), lower_mean = mean(na.omit(lower2)), upper_mean = mean(na.omit(upper2)))

# paste out predictions so we don't have to run again
out_name <- "./Data_outputs/bba_dep_prob_pred_mass.csv"
write.csv(preds_bba_avg, out_name)


##### 4B. GHA - AGE+MASS #####

hist(gha$mass)
mass_levels <- c(3.5, 4, 4.5) # predict at three levels of bird mass
mean_mass <- mean(gha$mass) 
sd_mass <- sd(gha$mass)
mass_st_levels <- (mass_levels - mean_mass) / sd_mass

# getting average values for wind speed and direction
hist(gha$ws)
mean(gha$ws)
hist(gha$wd)
mean_wd <- as.numeric(mean(circular(gha$wd, units = "degrees", template = "geographics"))) * pi / 180
sin_rad <- sin(mean_wd)
cos_rad <- cos(mean_wd)

# create prediction grid
ages <- seq(min(gha$age)-0.4, max(gha$age)+0.4, by = 0.1)

pred_grid <- expand.grid(
  age = ages,
  mass_st = mass_st_levels,
  ws = mean(gha$ws),
  sin_rad = sin(mean_wd),
  cos_rad = cos(mean_wd)
)

# predict transition probabilities - takes a while!
pred_list <- lapply(1:nrow(pred_grid), function(i) {
  print(paste(i, nrow(pred_grid), sep = "..."))
  
  P <- pmatrix.msm(
    gha_best, t = 5,
    covariates = list(age = pred_grid$age[i], mass_st = pred_grid$mass_st[i], 
                      sin_rad = pred_grid$sin_rad[i],
                      cos_rad = pred_grid$cos_rad[i], ws = pred_grid$ws[i]),
    ci = "normal",
    cl = 0.95
  )
  
  c(prob  = P$est[1,2],
    lower = P$L[1,2],
    upper = P$U[1,2])
})

preds_gha <- do.call(rbind, pred_list)
preds_gha <- cbind(pred_grid, preds_gha)
# add original mass levels for plotting
preds_gha$mass <- rep(mass_levels, each = length(ages))

# average predictions
preds_gha$lower2 <- rollmean(preds_gha$lower, k = 15, fill = NA)
preds_gha$upper2 <- rollmean(preds_gha$upper, k = 15, fill = NA)
preds_gha$age2 <- round(preds_gha$age)
preds_gha_avg <- ddply(preds_gha,.(age2, mass), summarize, prob_mean = mean(prob), lower_mean = mean(na.omit(lower2)), upper_mean = mean(na.omit(upper2)))

# paste out predictions so we don't have to run again
out_name <- "./Data_outputs/gha_dep_prob_pred_mass.csv"
write.csv(preds_gha_avg, out_name)


##### 4C. GHA - WIND #####

# wind direction levels
wd_levels <- c(0, 90, 180, 270) * pi / 180
sin_wd <- sin(wd_levels)
cos_wd <- cos(wd_levels)
# getting average mass
mass_avg_st <- 0
# get average age 
age_avg <- mean(gha$age[gha$state == 2])
# wind speed
ws_seq <- seq(min(gha$ws)-0.4, max(gha$ws)+0.4, 0.1)

pred_grid <- expand.grid(
  age = age_avg,
  mass_st = mass_avg_st,
  ws = ws_seq,
  sin_rad = sin_wd,
  cos_rad = cos_wd
)

# predict transition probabilities  - takes a while!
pred_list <- lapply(1:nrow(pred_grid), function(i) {
  print(paste(i, nrow(pred_grid), sep = "..."))
  
  P <- pmatrix.msm(
    gha_best, t = 5,
    covariates = list(age = pred_grid$age[i], mass_st = pred_grid$mass_st[i], 
                      sin_rad = pred_grid$sin_rad[i],
                      cos_rad = pred_grid$cos_rad[i], ws = pred_grid$ws[i]),
    ci = "normal",
    cl = 0.95
  )
  
  c(prob  = P$est[1,2],
    lower = P$L[1,2],
    upper = P$U[1,2])
})

preds_gha <- do.call(rbind, pred_list)
preds_gha <- cbind(pred_grid, preds_gha)
# add original mass levels for plotting
preds_gha$wd_levels <- rep(c(0, 90, 180, 270), each = length(ws_seq))

# average predictions
preds_gha$lower2 <- rollmean(preds_gha$lower, k = 15, fill = NA)
preds_gha$upper2 <- rollmean(preds_gha$upper, k = 15, fill = NA)
preds_gha$ws2 <- round(preds_gha$ws)
preds_gha_avg <- ddply(preds_gha,.(ws2, wd_levels), summarize, prob_mean = mean(prob), lower_mean = mean(na.omit(lower2)), upper_mean = mean(na.omit(upper2)))
preds_gha_avg$wd_levels <- factor(preds_gha_avg$wd_levels, levels = c("0", "90", "180", "270"))

# paste out predictions so we don't have to run again
out_name <- "./Data_outputs/gha_dep_prob_pred_wind.csv"
write.csv(preds_gha_avg, out_name)


##### 4D. WA - AGE+MASS #####

hist(wa$mass)
mass_levels <- c(10, 11, 12) # predict at three levels of bird mass
mean_mass <- mean(wa$mass)  
sd_mass <- sd(wa$mass)
mass_st_levels <- (mass_levels - mean_mass) / sd_mass
# create prediction grid
ages <- seq(min(wa$age)-0.4, max(wa$age)+0.4, by = 0.1)
pred_grid <- expand.grid(age = ages, mass_st = mass_st_levels)

# predict transition probabilities - takes a while!
pred_list <- lapply(1:nrow(pred_grid), function(i) {
  print(paste(i, nrow(pred_grid), sep = "..."))
  age_i <- pred_grid$age[i]
  mass_i <- pred_grid$mass_st[i]
  P <- pmatrix.msm(
    wa_best, t = 5,
    covariates = list(age = age_i, mass_st = mass_i),
    ci = "normal",
    cl = 0.95)
  c(prob  = P$est[1,2],
    lower = P$L[1,2],
    upper = P$U[1,2])
  })

preds_wa <- do.call(rbind, pred_list)
preds_wa <- cbind(pred_grid, preds_wa)
# add original mass levels for plotting
preds_wa$mass <- rep(mass_levels, each = length(ages))

# average predictions
preds_wa$lower2 <- rollmean(preds_wa$lower, k = 40, fill = NA)
preds_wa$upper2 <- rollmean(preds_wa$upper, k = 40, fill = NA)
preds_wa$age2 <- round(preds_wa$age)
preds_wa_avg <- ddply(preds_wa,.(age2, mass), summarize, prob_mean = mean(prob), lower_mean = mean(na.omit(lower2)), upper_mean = mean(na.omit(upper2)))

# paste out predictions so we don't have to run again
out_name <- "./Data_outputs/wa_dep_prob_pred_mass.csv"
write.csv(preds_wa_avg, out_name)


##### 4E. WCP - AGE #####

pred_grid <- data.frame(age = seq(min(wcp$age)-0.4, max(wcp$age)+0.4, 0.1))

# predict transition probabilities - takes a while!
pred_list <- lapply(1:nrow(pred_grid), function(i) {
  print(paste(i, nrow(pred_grid), sep = "..."))
  age_i <- pred_grid$age[i]
  P <- pmatrix.msm(
    wcp_best, t = 5,
    covariates = list(age = age_i),
    ci = "normal",
    cl = 0.95)
  c(prob  = P$est[1,2],
    lower = P$L[1,2],
    upper = P$U[1,2])
  })

preds_wcp <- do.call(rbind, pred_list)
preds_wcp <- cbind(pred_grid, preds_wcp)

# average predictions
preds_wcp$lower2 <- rollmean(preds_wcp$lower, k = 15, fill = NA)
preds_wcp$upper2 <- rollmean(preds_wcp$upper, k = 15, fill = NA)
preds_wcp$age2 <- round(preds_wcp$age)
preds_wcp_avg <- ddply(preds_wcp,.(age2), summarize, prob_mean = mean(prob), lower_mean = mean(na.omit(lower2)), upper_mean = mean(na.omit(upper2)))

# paste out predictions so we don't have to run again
out_name <- "./Data_outputs/wcp_dep_prob_pred_age.csv"
write.csv(preds_wcp_avg, out_name)



#### 5. PLOTTING AND PASTING OUT ####

plot_direx <- "./Plots/"

##### 5A. BBA - AGE+MASS- WEIRD RESULTS?!?! #####

# load in predictions
in_bba <- read.csv("./Data_outputs/bba_dep_prob_pred_mass.csv")

# change factor levels of mass
in_bba$mass <- as.character(in_bba$mass)
in_bba$mass[in_bba$mass == 3.5] <- "3.50"
in_bba$mass <- as.factor(in_bba$mass)

p_bba <- ggplot(in_bba, aes(x = age2, y = prob_mean, fill = mass, linetype = mass)) +
  geom_ribbon(aes(ymin = lower_mean, ymax = upper_mean), alpha = 0.2, fill = "#CC79A7") +
  geom_line(linewidth = 0.6, col = "#CC79A7") +
  labs(y = "Departure probability",
       x = "Chick age (days)", linetype = "Mass (kg)", fill = "Mass (kg)") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 11), legend.position = c(0.17, 0.775),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                     legend.text = element_text(size = 10),
                     panel.border = element_rect(colour = "black", linewidth = 0.8))+
  scale_x_continuous(breaks = seq(95, 125, 5), limits = c(93, 128))
print(p_bba)

# pasting out
out_name <- paste0(plot_direx, "figure1_bba_mass.png")
ggsave(out_name, width = 2000, height = 1800, units = "px", dpi = 600)
dev.off() 

##### 5B. GHA - AGE+MASS #####

# load in predictions
in_gha_mass <- read.csv("./Data_outputs/gha_dep_prob_pred_mass.csv")
in_gha_mass$mass <- as.factor(as.character(in_gha_mass$mass))

p_gha_mass <- ggplot(in_gha_mass, aes(x = age2, y = prob_mean, fill = mass, linetype = mass)) +
  geom_ribbon(aes(ymin = lower_mean, ymax = upper_mean), alpha = 0.2, fill = "#009E73") +
  geom_line(linewidth = 0.6, col = "#009E73") +
  labs(y = "Departure probability",
       x = "Chick age (days)", linetype = "Mass (kg)", fill = "Mass (kg)") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 11), legend.position = c(0.17, 0.773),
                     legend.text = element_text(size = 10),
                     panel.border = element_rect(colour = "black", linewidth = 0.8))+
  scale_x_continuous(breaks = seq(120, 160, 10), limits = c(120, 164))
print(p_gha_mass)

# pasting out
out_name <- paste0(plot_direx, "figure1_gha_mass.png")
ggsave(out_name, width = 2000, height = 1800, units = "px", dpi = 600)
dev.off() 


##### 5C. GHA - WIND #####

# load in predictions
in_gha_wind <- read.csv("./Data_outputs/gha_dep_prob_pred_wind.csv")
# remove 180 as its the same as 0
in_gha_wind$wd_levels <- as.character(in_gha_wind$wd_levels)
in_gha_wind <- subset(in_gha_wind, wd_levels != "180",)
in_gha_wind$wd_levels[in_gha_wind$wd_levels == "270"] <- "-90"
in_gha_wind$wd_levels <- as.factor(as.character(in_gha_wind$wd_levels))

p_gha_wind <-  ggplot(in_gha_wind, aes(x = ws2, y = prob_mean, fill = wd_levels, linetype = wd_levels)) +
  geom_ribbon(aes(ymin = lower_mean, ymax = upper_mean), alpha = 0.2, fill = "#009E73") +
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
  geom_line(size = 0.6, col = "#009E73") +
  labs(y = "Departure probability",
       x = "Wind speed (m/s)", linetype = "Wind direction", fill = "Wind direction") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 11), 
                     legend.text = element_text(size = 10),
                     panel.border = element_rect(colour = "black", linewidth = 0.8))+
  scale_x_continuous(breaks = seq(4, 14, 2), limits = c(4, 14))
print(p_gha_wind)

# pasting out
out_name <- paste0(plot_direx, "figure1_gha_wind.png")
ggsave(out_name, width = 2775, height = 1800, units = "px", dpi = 600)
dev.off() 


##### 5D. WA - AGE+MASS #####

# load in predictions
in_wa <- read.csv("./Data_outputs/wa_dep_prob_pred_mass.csv")
in_wa$mass <- as.factor(as.character(in_wa$mass))
in_wa <- na.omit(in_wa)

p_wa <- ggplot(in_wa, aes(x = age2, y = prob_mean, fill = factor(mass), linetype = factor(mass))) +
  geom_ribbon(aes(ymin = lower_mean, ymax = upper_mean), alpha = 0.2, fill = "#0072B2") +
  geom_line(size = 0.6, col = "#0072B2") +
  labs(y = "Departure probability",
       x = "Chick age (days)", linetype = "Mass (kg)", fill = "Mass (kg)") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 11), legend.position = c(0.17, 0.773),
                     legend.text = element_text(size = 10),
                     panel.border = element_rect(colour = "black", linewidth = 0.8))+
  scale_x_continuous(breaks = seq(250, 290, 10), limits = c(248, 290))
print(p_wa)

# pasting out
out_name <- paste0(plot_direx, "figure1_wa_mass.png")
ggsave(out_name, width = 2000, height = 1800, units = "px", dpi = 600)
dev.off() 


##### 5E. WCP - AGE #####

# load in predictions
in_wcp <- read.csv("./Data_outputs/wcp_dep_prob_pred_age.csv")

p_wcp <- ggplot(in_wcp, aes(x = age2, y = prob_mean)) +
  geom_line(linewidth = 0.6,  colour = "#E69F00") +
  geom_ribbon(aes(ymin = lower_mean, ymax = upper_mean), alpha = 0.2, fill = "#E69F00") +
  labs(y = "Departure probability",
       x = "Chick age (days)") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 11), 
                     panel.border = element_rect(colour = "black", linewidth = 0.8))+
  scale_x_continuous(breaks = seq(70, 110, 10), limits = c(70, 110))
print(p_wcp)
# not in figure 1