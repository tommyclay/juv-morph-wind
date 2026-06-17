library(msm)
library(bbmle) 
library(circular)
library(ggplot2)

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
#        to
# from   1   2
#     1 902  51
# running a null model
mod1 <- msm(state ~ age, subject=id, data = df,
            qmatrix = Q_juv)
print(mod1)
pmatrix.msm(mod1) # to get the transition matrix 
print(AIC(mod1))
# 399.8552


#### 3. RUN MARKOV MODELS FOR EACH SPECIES ####

##### 3A. BBA #####

bba <- subset(df, species == "BBA",)
bba  <- bba[order(bba$id, bba$age), ]

# 1. NULL MODEL WITH AGE
bba_m1 <-  msm(state ~ age, subject=id, data = bba,
               qmatrix = Q_juv)
# 2. AGE AND MASS
bba_m2 <-  msm(state ~ age, subject=id, data = bba,
               qmatrix = Q_juv, covariates = ~ mass_st)
# 3. AGE AND MASS INTERACTION
bba_m3 <-  msm(state ~ age, subject=id, data = bba,
               qmatrix = Q_juv, covariates = ~ age*mass_st)
# 4. WIND DIRECTION 
bba_m4 <-  msm(state ~ age , subject=id, data = bba,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad)
# 5. WIND SPEED
bba_m5 <-  msm(state ~ age, subject=id, data = bba,
               qmatrix = Q_juv, covariates = ~ ws)
# 6. WIND SPEED AND DIRECTION
bba_m6 <-  msm(state ~ age, subject=id, data = bba,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad + ws)
# 7. WIND SPEED AND DIRECTION, MASS AND AGE INTERACTION
bba_m7 <-  msm(state ~ age, subject=id, data = bba,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad+ ws + age*mass_st)
# 8. WIND SPEED, MASS AND AGE INTERACTION
bba_m8 <-  msm(state ~ age, subject=id, data = bba,
               qmatrix = Q_juv, covariates = ~ ws + age*mass_st)
# 9. WIND DIRECTION, MASS AND AGE INTERACTION
bba_m9 <-  msm(state ~ age, subject=id, data = bba,
                qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad + age*mass_st)

# COMPARE AIC
AICtab(bba_m1, bba_m2, bba_m3, bba_m4, bba_m5, bba_m6, bba_m7, bba_m8, bba_m9, 
       base = TRUE, delta = TRUE, weights = TRUE, sort = TRUE)
#        AIC  dAIC df weight
# bba_m3 75.4  0.0 4  0.3478
# bba_m9 75.9  0.6 6  0.2605
# bba_m7 76.4  1.0 7  0.2109
# bba_m8 76.8  1.4 5  0.1722
# bba_m4 84.1  8.7 3  0.0045
# bba_m6 85.2  9.9 4  0.0025
# bba_m1 87.2 11.8 1  <0.001
# bba_m5 89.0 13.6 2  <0.001
# bba_m2 89.2 13.8 2  <0.001

# BEST MODEL HAS AGE AND MASS INTERACTION
bba_best <- bba_m3


##### 3B. GHA #####

gha <- subset(df, species == "GHA",)
gha  <- gha[order(gha$id, gha$age), ]

# 1. NULL MODEL WITH AGE
gha_m1 <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv)
# 2. AGE AND MASS
gha_m2 <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ mass_st)
# 3. AGE AND MASS INTERACTION
gha_m3 <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ age*mass_st)
# 4. WIND DIRECTION 
gha_m4 <-  msm(state ~ age , subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad)
# 5. WIND SPEED
gha_m5 <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ ws)
# 6. WIND SPEED AND DIRECTION
gha_m6 <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad + ws)
# 7. WIND SPEED AND DIRECTION, MASS AND AGE INTERACTION
gha_m7 <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad+ ws + age*mass_st)
# 8. WIND SPEED, MASS AND AGE INTERACTION
gha_m8 <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ ws + age*mass_st)
# 9. WIND DIRECTION, MASS AND AGE INTERACTION
gha_m9 <-  msm(state ~ age, subject=id, data = gha,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad + age*mass_st)

# COMPARE AIC
AICtab(gha_m1, gha_m2, gha_m3, gha_m4, gha_m5, gha_m6, gha_m7, gha_m8, gha_m9, 
       base = TRUE, delta = TRUE, weights = TRUE, sort = TRUE)
#       AIC   dAIC  df weight
# gha_m7  97.6   0.0 7  0.7978
# gha_m8 101.2   3.6 5  0.1349
# gha_m9 102.7   5.1 6  0.0635
# gha_m3 108.5  10.9 4  0.0034
# gha_m6 114.5  16.8 4  <0.001
# gha_m4 114.9  17.3 3  <0.001
# gha_m1 127.4  29.7 1  <0.001
# gha_m5 128.6  30.9 2  <0.001
# gha_m2 129.4  31.7 2  <0.001

# BEST MODEL HAS AGE AND MASS INTERACTION AND WIND SPEED AND DIRECTION
gha_best <- gha_m7


##### 3C. WA #####

wa <- subset(df, species == "WALB",)
wa  <- wa[order(wa$id, wa$age), ]

# 1. NULL MODEL WITH AGE
wa_m1 <-  msm(state ~ age, subject=id, data = wa,
               qmatrix = Q_juv)
# 2. AGE AND MASS
wa_m2 <-  msm(state ~ age, subject=id, data = wa,
               qmatrix = Q_juv, covariates = ~ mass_st)
# 3. AGE AND MASS INTERACTION
wa_m3 <-  msm(state ~ age, subject=id, data = wa,
               qmatrix = Q_juv, covariates = ~ age*mass_st)
# 4. WIND DIRECTION 
wa_m4 <-  msm(state ~ age , subject=id, data = wa,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad)
# 5. WIND SPEED
wa_m5 <-  msm(state ~ age, subject=id, data = wa,
               qmatrix = Q_juv, covariates = ~ ws)
# 6. WIND SPEED AND DIRECTION
wa_m6 <-  msm(state ~ age, subject=id, data = wa,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad + ws)
# 7. WIND SPEED AND DIRECTION, MASS AND AGE INTERACTION
wa_m7 <-  msm(state ~ age, subject=id, data = wa,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad+ ws + age*mass_st)
# 8. WIND SPEED, MASS AND AGE INTERACTION
wa_m8 <-  msm(state ~ age, subject=id, data = wa,
               qmatrix = Q_juv, covariates = ~ ws + age*mass_st)
# 9. WIND DIRECTION, MASS AND AGE INTERACTION
wa_m9 <-  msm(state ~ age, subject=id, data = wa,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad + age*mass_st)

# COMPARE AIC
AICtab(wa_m1, wa_m2, wa_m3, wa_m4, wa_m5, wa_m6, wa_m7, wa_m8, wa_m9, 
       base = TRUE, delta = TRUE, weights = TRUE, sort = TRUE)
#       AIC   dAIC  df weight
# wa_m3  87.4   0.0 4  0.4978
# wa_m8  88.0   0.6 5  0.3602
# wa_m9  91.1   3.7 6  0.0775
# wa_m7  92.0   4.6 7  0.0500
# wa_m1  96.0   8.6 1  0.0067
# wa_m5  97.7  10.3 2  0.0028
# wa_m2  98.0  10.6 2  0.0025
# wa_m4  98.7  11.3 3  0.0017
# wa_m6 100.1  12.7 4  <0.001

# BEST MODEL HAS AGE AND MASS INTERACTION
wa_best <- wa_m3


##### 3D. WCP #####

wcp <- subset(df, species == "WCP",)
wcp  <- wcp[order(wcp$id, wcp$age), ]

# 1. NULL MODEL WITH AGE
wcp_m1 <-  msm(state ~ age, subject=id, data = wcp,
               qmatrix = Q_juv)
# 2. AGE AS COVARIATE 
wcp_m2 <-  msm(state ~ age, subject=id, data = wcp,
               qmatrix = Q_juv, covariates = ~ age)
# 3. WIND DIRECTION 
wcp_m3 <-  msm(state ~ age , subject=id, data = wcp,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad)
# 4. WIND SPEED
wcp_m4 <-  msm(state ~ age, subject=id, data = wcp,
               qmatrix = Q_juv, covariates = ~ ws)
# 5. WIND SPEED AND DIRECTION
wcp_m5 <-  msm(state ~ age, subject=id, data = wcp,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad + ws)
# 6. WIND SPEED AND DIRECTION, AND AGE
wcp_m6 <-  msm(state ~ age, subject=id, data = wcp,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad+ ws + age)
# 7. WIND SPEED AND AGE
wcp_m7 <-  msm(state ~ age, subject=id, data = wcp,
               qmatrix = Q_juv, covariates = ~ ws + age)
# 8. WIND DIRECTION AND AGE
wcp_m8 <-  msm(state ~ age, subject=id, data = wcp,
               qmatrix = Q_juv, covariates = ~ sin_rad + cos_rad + age)

# COMPARE
AICtab(wcp_m1, wcp_m2, wcp_m3, wcp_m4, wcp_m5, wcp_m6, wcp_m7, wcp_m8,
       base = TRUE, weights = TRUE, sort = TRUE)
#         AIC  dAIC df weight
# wcp_m2 85.2  0.0 2  0.4594
# wcp_m7 86.1  0.9 3  0.2976
# wcp_m8 87.3  2.1 4  0.1598
# wcp_m6 88.9  3.6 5  0.0748
# wcp_m4 95.0  9.8 2  0.0035
# wcp_m1 95.3 10.0 1  0.0031
# wcp_m3 97.5 12.3 3  <0.001
# wcp_m5 98.0 12.7 4  <0.001

# BEST MODEL HAS EFFECT OF AGE
wcp_best <- wcp_m2


##### 4. EXAMINING OBSERVED VERSUS PREDICTED TIME SERIES #####

##### 4A. BBA #####

prev <- prevalence.msm(bba_best)
prev_o_perc <- data.frame(prev$`Observed percentages`)
prev_e_perc <- data.frame(prev$`Expected percentages`)
obs_perc <- data.frame(Age = as.numeric(row.names(prev_o_perc)),
                       Perc = c(prev_o_perc$State.1, prev_o_perc$State.2),
                       State = rep(c("Colony", "Departed"), each = nrow(prev_o_perc)),
                       Type = "Observed")
exp_perc <- data.frame(Age = as.numeric(row.names(prev_e_perc)),
                       Perc = c(prev_e_perc$State.1, prev_e_perc$State.2),
                       State = rep(c("Colony", "Departed"), each = nrow(prev_e_perc)),
                       Type = "Predicted")
both_perc <- rbind(obs_perc, exp_perc)

# plotting
p_bba <- ggplot() + geom_line(data = both_perc, aes(Age, Perc, colour = State, linetype = Type))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 12), legend.position.inside = c(0.17, 0.55),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 10),
                     panel.border = element_rect(colour = "black", linewidth = 0.8))
p_bba


##### 4B. GHA #####

prev <- prevalence.msm(gha_best)
prev_o_perc <- data.frame(prev$`Observed percentages`)
prev_e_perc <- data.frame(prev$`Expected percentages`)
obs_perc <- data.frame(Age = as.numeric(row.names(prev_o_perc)),
                       Perc = c(prev_o_perc$State.1, prev_o_perc$State.2),
                       State = rep(c("Colony", "Departed"), each = nrow(prev_o_perc)),
                       Type = "Observed")
exp_perc <- data.frame(Age = as.numeric(row.names(prev_e_perc)),
                       Perc = c(prev_e_perc$State.1, prev_e_perc$State.2),
                       State = rep(c("Colony", "Departed"), each = nrow(prev_e_perc)),
                       Type = "Predicted")
both_perc <- rbind(obs_perc, exp_perc)

# plotting
p_gha <- ggplot() + geom_line(data = both_perc, aes(Age, Perc, colour = State, linetype = Type))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 12), legend.position.inside = c(0.17, 0.55),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 10),
                     panel.border = element_rect(colour = "black", linewidth = 0.8))
p_gha


##### 4C. WA #####

prev <- prevalence.msm(wa_best)
prev_o_perc <- data.frame(prev$`Observed percentages`)
prev_e_perc <- data.frame(prev$`Expected percentages`)
obs_perc <- data.frame(Age = as.numeric(row.names(prev_o_perc)),
                       Perc = c(prev_o_perc$State.1, prev_o_perc$State.2),
                       State = rep(c("Colony", "Departed"), each = nrow(prev_o_perc)),
                       Type = "Observed")
exp_perc <- data.frame(Age = as.numeric(row.names(prev_e_perc)),
                       Perc = c(prev_e_perc$State.1, prev_e_perc$State.2),
                       State = rep(c("Colony", "Departed"), each = nrow(prev_e_perc)),
                       Type = "Predicted")
both_perc <- rbind(obs_perc, exp_perc)

# plotting
p_wa <- ggplot() + geom_line(data = both_perc, aes(Age, Perc, colour = State, linetype = Type))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 12), legend.position.inside = c(0.17, 0.55),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 10),
                     panel.border = element_rect(colour = "black", linewidth = 0.8))
p_wa


##### 4D. WCP #####

prev <- prevalence.msm(wcp_best)
prev_o_perc <- data.frame(prev$`Observed percentages`)
prev_e_perc <- data.frame(prev$`Expected percentages`)
obs_perc <- data.frame(Age = as.numeric(row.names(prev_o_perc)),
                       Perc = c(prev_o_perc$State.1, prev_o_perc$State.2),
                       State = rep(c("Colony", "Departed"), each = nrow(prev_o_perc)),
                       Type = "Observed")
exp_perc <- data.frame(Age = as.numeric(row.names(prev_e_perc)),
                       Perc = c(prev_e_perc$State.1, prev_e_perc$State.2),
                       State = rep(c("Colony", "Departed"), each = nrow(prev_e_perc)),
                       Type = "Predicted")
both_perc <- rbind(obs_perc, exp_perc)

# plotting
p_wcp <- ggplot() + geom_line(data = both_perc, aes(Age, Perc, colour = State, linetype = Type))+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 12), legend.position.inside = c(0.17, 0.55),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 10),
                     panel.border = element_rect(colour = "black", linewidth = 0.8))
p_wcp

