
# library(geosphere)
# library(plyr)
library(ggplot2)
# library(fields)
# library(lme4)
library(MuMIn)
# library(visreg)
# library(psych)
# library(viridis)
library(glmmTMB)
# library(dplyr)
library(DHARMa)
library(psych)


#### 1. LOAD IN DATA ####

dat <- read.csv("C:/Users/44753/Documents/GitHub/juv-morph-wind/Data_inputs/juv_1st_month.csv")
# for displacement models
dat_dis <- read.csv("C:/Users/44753/Documents/GitHub/juv-morph-wind/Data_inputs/juv_displacement.csv")


#### 2. RUN MODELS ####

##### 2A. ORIENTATION WITH WINDS #####

wdrel_df <- na.omit(data.frame(species = dat$species, id = dat$id, days = dat$days, seq = dat$seq, wdrel = dat$wdrel))
wdrel_df$id <- as.factor(as.character(wdrel_df$id))
wdrel_df$species <- as.factor(as.character(wdrel_df$species))

# rescale to between 0 and 1
wdrel_df$wdrel_scale <- wdrel_df$wdrel/180
eps <- 1e-4  # small constant to avoid exact 0/1
wdrel_df <- wdrel_df %>%
  mutate(wdrel_scale2 = wdrel_scale * (1 - 2*eps) + eps)

# RUNNING FULL MODEL 

# sequence needs to be a factor
wdrel_df$seq_f <- factor(wdrel_df$seq)

all_wdrel <- glmmTMB(wdrel_scale2 ~ days*species + ar1(seq_f + 0 | id),  # AR1 for repeated measures
                     family = beta_family(link = "logit"),data = wdrel_df)
summary(all_wdrel)
# COMPARING CANDIDATE MODELS WITH AICc
wdrel_1 <-  glmmTMB(wdrel_scale2 ~ days + species + ar1(seq_f + 0 | id), 
                    family = beta_family(link = "logit"),data = wdrel_df)
wdrel_2 <-  glmmTMB(wdrel_scale2 ~ species + ar1(seq_f + 0 | id),  
                    family = beta_family(link = "logit"),data = wdrel_df)
wdrel_3 <-  glmmTMB(wdrel_scale2 ~ days + ar1(seq_f + 0 | id),  
                    family = beta_family(link = "logit"),data = wdrel_df)
wdrel_0 <-  glmmTMB(wdrel_scale2 ~ ar1(seq_f + 0 | id),  
                    family = beta_family(link = "logit"),data = wdrel_df)
AIC_df <- data.frame(AICc(all_wdrel, wdrel_1, wdrel_2, wdrel_3, wdrel_0))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
#           df      AICc     dAICc
# wdrel_1    8 -970.7268  0.000000
# all_wdrel 11 -970.5030  0.223854
# wdrel_2    7 -967.8226  2.904260
# wdrel_3    5 -939.6196 31.107234
# wdrel_0    4 -937.3598 33.367058

# BEST MODEL HAS DAYS AND SPECIES

all_wdrel_fin <- glmmTMB(wdrel_scale2 ~ days+species + ar1(seq_f + 0 | id),  
                         dispformula = ~ 1,    
                         family = beta_family(),data = wdrel_df)
summary(all_wdrel_fin)
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.04706    0.13442  -7.789 6.74e-15 ***
# days        -0.01227    0.00552  -2.223   0.0262 *  
# speciesGHA   0.17222    0.13865   1.242   0.2142    
# speciesWALB  0.77648    0.14287   5.435 5.48e-08 ***
# speciesWCP   0.37830    0.14973   2.526   0.0115 *  

# RESIDUAL CHECKS
sim <- simulateResiduals(all_wdrel_fin)
plot(sim)
# residuals not uniform
testDispersion(sim)
testTemporalAutocorrelation(sim, time = dat$time)
# fine


##### 2B. WINDS EXPERIENCED #####

# SETTING UP DATAFRAME 
wind_df <- na.omit(data.frame(species = dat$species, id = dat$id, ws = dat$ws, days = dat$days, seq = dat$seq))
wind_df$id <- as.factor(as.character(wind_df$id))
wind_df$species <- as.factor(as.character(wind_df$species))

# RUNNING FULL MODEL 

# sequence needs to be a factor
wind_df$seq_f <- factor(wind_df$seq)

all_wind <- glmmTMB(ws ~ days*species + ar1(seq_f + 0 | id), family = gaussian(), data = wind_df)
summary(all_wind)
# Model selection table 
#    cnd((Int)) dsp((Int))  cnd(dys) cnd(spc) cnd(dys:spc) df    logLik    AICc delta weight
# 8      9.796          + -0.134300        +            + 11 -5306.668 10635.5  0.00      1
# 3      8.075          +                  +               7 -5326.485 10667.0 31.56      0
# 4      8.181          + -0.008300        +               8 -5326.335 10668.7 33.28      0
# 1      8.220          +                                  4 -5346.308 10700.6 65.17      0
# 2      8.344          + -0.009421                        5 -5346.137 10702.3 66.84      0

# COMPARING CANDIDATE MODELS WITH AICc
wind_1 <-  glmmTMB(ws ~ days + species + ar1(seq_f + 0 | id), 
                   family = gaussian(), data = wind_df)
wind_2 <-  glmmTMB(ws ~ species + ar1(seq_f + 0 | id),  
                   family = gaussian(), data = wind_df)
wind_3 <-  glmmTMB(ws ~ days + ar1(seq_f + 0 | id),  
                   family = gaussian(), data = wind_df)
wind_0 <-  glmmTMB(ws ~ ar1(seq_f + 0 | id),  
                   family = gaussian(), data = wind_df)
AIC_df <- data.frame(AICc(all_wind, wind_1, wind_2, wind_3, wind_0))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
#         df     AICc    dAICc
# all_wind 11 10635.46  0.00000
# wind_2    7 10667.02 31.56201
# wind_1    8 10668.74 33.27738
# wind_0    4 10700.63 65.17403
# wind_3    5 10702.30 66.84199


# BEST MODEL HAD INTERACTION BETWEEN DAYS AND SPECIES

all_wind_fin <- all_wind
summary(all_wind_fin)

# RESIDUAL CHECKS
sim <- simulateResiduals(all_wind_fin)
plot(sim)
testDispersion(sim)
testTemporalAutocorrelation(sim, time = dat$time)
# fine


##### 2C. TRAVEL SPEEDS AS A FUNTION OF WIND SPEEDS #####


# SETTING UP DATAFRAME 
speed_df <- na.omit(data.frame(species = dat$species, id = dat$id, ws = dat$ws, days = dat$days,
                               speed = dat$speed, seq = dat$seq))
speed_df$id <- as.factor(as.character(speed_df$id))
speed_df$species <- as.factor(as.character(speed_df$species))

# RUNNING FULL MODEL WITH THREE-WAY INTERACTION

# sequence needs to be a factor
speed_df$seq_f <- factor(speed_df$seq)
# scaling days to help convergence
speed_df$days_s <- scale(speed_df$days)

all_speed <- glmmTMB(speed ~ days_s*ws*species + ar1(seq_f + 0 | id), family = tweedie(), data = speed_df)
summary(all_speed)
# COMPARING CANDIDATE MODELS TO DETERMINE SIGNIFICANCE OF SPECIES INTERACTIONS
speed_1 <-  glmmTMB(speed ~ days*ws + species + ar1(seq_f + 0 | id), family = tweedie(), data = speed_df)
speed_2 <-  glmmTMB(speed ~ days*species + ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_df)
speed_3 <-  glmmTMB(speed ~ days + ws*species + ar1(seq_f + 0 | id), family = tweedie(), data = speed_df)
speed_4 <-  glmmTMB(speed ~ days + ws + species + ar1(seq_f + 0 | id), family = tweedie(), data = speed_df)
AIC_df <- data.frame(AICc(all_speed, speed_1, speed_2, speed_3, speed_4))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
# #    df     AICc      dAICc
# speed_2   13 6826.229  0.0000000
# all_speed 20 6827.033  0.8043302
# speed_3   13 6909.834 83.6045313
# speed_4   10 6917.728 91.4990399
# speed_1   11 6918.801 92.5717432

# SPECIES x WIND SPEED INTERACTIONS ARE RETAINED - SPLIT BY SPECIES

# BBA #

speed_bba_df <- subset(speed_df, species == "BBA",)
speed_bba_df$id <- as.factor(as.character(speed_bba_df$id))
bba_speed <- glmmTMB(speed ~ days_s*ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_bba_df)
# COMPARING CANDIDATE MODELS WITH AICc
bba_1 <-  glmmTMB(speed ~ days_s + ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_bba_df)
bba_2 <-  glmmTMB(speed ~ days_s + ar1(seq_f + 0 | id), family = tweedie(), data = speed_bba_df)
bba_3 <-  glmmTMB(speed ~ ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_bba_df)
bba_0 <-  glmmTMB(speed ~ ar1(seq_f + 0 | id), family = tweedie(), data = speed_bba_df)
AIC_df <- data.frame(AICc(bba_speed, bba_1, bba_2, bba_3, bba_0))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
# df     AICc     dAICc
# bba_speed  8 1443.851  0.000000
# bba_1      7 1447.112  3.261900
# bba_3      6 1453.151  9.300146
# bba_2      6 1456.842 12.991408
# bba_0      5 1462.189 18.338087

# BEST MODEL HAD WIND SPEED, DAYS AND THEIR INTERACTION

bba_speed_fin <- glmmTMB(speed ~ days_s*ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_bba_df)
summary(bba_speed_fin)
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.97116    0.10967   8.856  < 2e-16 ***
# days_s      -0.45402    0.11339  -4.004 6.23e-05 ***
# ws           0.04249    0.01101   3.860 0.000113 ***
# days_s:ws    0.02713    0.01153   2.353 0.018624 *  

# RESIDUAL CHECKS
sim <- simulateResiduals(bba_speed_fin)
plot(sim)
# quantile deviations
testDispersion(sim)
testTemporalAutocorrelation(sim, time = speed_bba_df$time)
# fine


# GHA ##### RESULTS CHANGED ######

speed_gha_df <- subset(speed_df, species == "GHA",)
speed_gha_df$id <- as.factor(as.character(speed_gha_df$id))
gha_speed <- glmmTMB(speed ~ days_s*ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_gha_df)
# COMPARING CANDIDATE MODELS WITH AICc
gha_1 <-  glmmTMB(speed ~ days_s + ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_gha_df)
gha_2 <-  glmmTMB(speed ~ days_s + ar1(seq_f + 0 | id), family = tweedie(), data = speed_gha_df)
gha_3 <-  glmmTMB(speed ~ ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_gha_df)
gha_0 <-  glmmTMB(speed ~ ar1(seq_f + 0 | id), family = tweedie(), data = speed_gha_df)
AIC_df <- data.frame(AICc(gha_speed, gha_1, gha_2, gha_3, gha_0))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
#          df     AICc      dAICc
# gha_1      7 2510.766  0.0000000
# gha_speed  8 2511.574  0.8083083
# gha_3      6 2512.369  1.6026062
# gha_2      6 2521.676 10.9097152
# gha_0      5 2522.384 11.6181801

# BEST MODEL HAS WIND SPEED AND DAYS  ##### RESULTS CHANGED ######

gha_speed_fin <- glmmTMB(speed ~ days+ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_gha_df)
summary(gha_speed_fin)
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.944466   0.173184   5.454 4.94e-08 ***
# days        -0.020684   0.010337  -2.001 0.045404 *  
# ws           0.023897   0.006624   3.608 0.000309 ***

# RESIDUAL CHECKS
sim <- simulateResiduals(gha_speed_fin)
plot(sim)
# quantile deviations
testDispersion(sim)
# dispersion test significant
testTemporalAutocorrelation(sim, time = speed_gha_df$time)
# fine


# WCP #

speed_wcp_df <- subset(speed_df, species == "WCP",)
speed_wcp_df$id <- as.factor(as.character(speed_wcp_df$id))
wcp_speed <- glmmTMB(speed ~ days_s*ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wcp_df)
# COMPARING CANDIDATE MODELS WITH AICc
wcp_1 <-  glmmTMB(speed ~ days_s + ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wcp_df)
wcp_2 <-  glmmTMB(speed ~ days_s + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wcp_df)
wcp_3 <-  glmmTMB(speed ~ ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wcp_df)
wcp_0 <-  glmmTMB(speed ~ ar1(seq_f + 0 | id), family = tweedie(), data = speed_wcp_df)
AIC_df <- data.frame(AICc(wcp_speed, wcp_1, wcp_2, wcp_3, wcp_0))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
# df     AICc    dAICc
# wcp_3      6 1509.772 0.000000
# wcp_1      7 1510.021 0.249130
# wcp_speed  8 1512.074 2.301967
# wcp_2      6 1517.339 7.566992
# wcp_0      5 1518.097 8.324782

# BEST MODEL HAS JUST WIND SPEED

wcp_speed_fin <- glmmTMB(speed ~ ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wcp_df)
summary(wcp_speed_fin)
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.45137    0.12274   3.678 0.000236 ***
# ws           0.03550    0.01099   3.231 0.001234 ** 

# RESIDUAL CHECKS
sim <- simulateResiduals(wcp_speed_fin)
plot(sim)
# quantile deviations
testDispersion(sim)
testTemporalAutocorrelation(sim, time = speed_wcp_df$time)
# fine


# WA #

speed_wa_df <- subset(speed_df, species == "WALB",)
speed_wa_df$id <- as.factor(as.character(speed_wa_df$id))
wa_speed <- glmmTMB(speed ~ days_s*ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wa_df)
# COMPARING CANDIDATE MODELS WITH AICc
wa_1 <-  glmmTMB(speed ~ days_s + ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wa_df)
wa_2 <-  glmmTMB(speed ~ days_s + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wa_df)
wa_3 <-  glmmTMB(speed ~ ws + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wa_df)
wa_0 <-  glmmTMB(speed ~ ar1(seq_f + 0 | id), family = tweedie(), data = speed_wa_df)
AIC_df <- data.frame(AICc(wa_speed, wa_1, wa_2, wa_3, wa_0))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
#         df     AICc      dAICc
# wa_1      7 1231.690  0.0000000
# wa_speed  8 1232.580  0.8902969
# wa_2      6 1255.158 23.4685452
# wa_3      6 1296.440 64.7499777
# wa_0      5 1317.764 86.0741969

# BEST MODEL HAS DAYS AND WIND SPEED

wa_speed_fin <- glmmTMB(speed ~ ws + days_s + ar1(seq_f + 0 | id), family = tweedie(), data = speed_wa_df)
summary(wa_speed_fin)
# Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.65148    0.11582  -5.625 1.85e-08 ***
#   ws           0.05618    0.01101   5.101 3.38e-07 ***
#   days_s       0.81122    0.07866  10.313  < 2e-16 ***

# RESIDUAL CHECKS
sim <- simulateResiduals(wa_speed_fin)
plot(sim)
testDispersion(sim)
# dispersion significant
testTemporalAutocorrelation(sim, time = speed_wa_df$time)
# autocorrelation significant



##### 2D. DISPLACEMENT DURING FIRST WEEK #####


# CHECKING CORRELATION BETWEEN VARIABLES
MyVar <- c("mass_scale", "tw_wk")
psych::corr.test(dat_dis[,MyVar], method = "spearman")
#             mass_scale tw_wk
# mass_scale       1.00 -0.03
# tw_wk           -0.03  1.00

hist(dat_dis$tw_wk)
hist(dat_dis$mass_scale)
hist(dat_dis$dist_wk_scale)

wk_mod <- lm(dist_wk_scale ~ mass_scale*species + tw_wk, data = dat_dis)
summary(wk_mod)
m1 <- lm(dist_wk_scale ~ tw_wk, data = dat_dis)
m2 <- lm(dist_wk_scale ~ mass_scale, data = dat_dis)
m3 <- lm(dist_wk_scale ~ species, data = dat_dis)
m4 <- lm(dist_wk_scale ~ mass_scale + tw_wk, data = dat_dis)
m5 <- lm(dist_wk_scale ~ species + tw_wk, data = dat_dis)
m6 <- lm(dist_wk_scale ~ mass_scale + species, data = dat_dis)
m7 <- lm(dist_wk_scale ~ mass_scale + species + tw_wk,data = dat_dis)
m8 <- lm(dist_wk_scale ~ mass_scale * species,data = dat_dis)
AIC_df <- data.frame(AICc(wk_mod, m1, m2, m3, m4, m5, m6, m7, m8))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
#         df      AICc     dAICc
# m4      4  90.87196 0.0000000
# m7      6  91.19347 0.3215094
# m2      3  91.72960 0.8576402
# m1      3  94.09544 3.2234857
# m5      5  94.43011 3.5581465
# wk_mod  8  96.43020 5.5582370
# m6      5  97.12424 6.2522762
# m3      4  99.93328 9.0613194
# m8      7 100.53184 9.6598847

# BEST MODEL JUST HAS MASS SCALED

wk_fin <- lm(dist_wk_scale ~ mass_scale, data = dat_dis)
summary(wk_fin)
# Coefficients:
#   Estimate              Std. Error t value Pr(>|t|)  
# (Intercept)  0.00000000000000002291  0.15730982752032832139   0.000   1.0000  
# mass_scale  -0.39515614564927048491  0.16498793900744068708  -2.395   0.0228 *


# RESIDUAL CHECKS
sim <- simulateResiduals(wk_fin)
plot(sim)
testDispersion(sim)
# all good



##### 2E. DISPLACEMENT DURING FIRST MONTH #####

# CHECKING CORRELATION BETWEEN VARIABLES
MyVar <- c("mass_scale", "tw_mo")
psych::corr.test(dat_dis[,MyVar], method = "spearman")
#             mass_scale tw_mo
# mass_scale       1.00 -0.04
# tw_mo           -0.04  1.00

hist(dat_dis$tw_mo)
hist(dat_dis$dist_mo_scale)

mo_mod <- lm(dist_mo_scale ~ mass_scale*species + tw_mo, data = dat_dis)
summary(mo_mod)
m1 <- lm(dist_mo_scale ~ tw_mo, data = dat_dis)
m2 <- lm(dist_mo_scale ~ mass_scale, data = dat_dis)
m3 <- lm(dist_mo_scale ~ species, data = dat_dis)
m4 <- lm(dist_mo_scale ~ mass_scale + tw_mo, data = dat_dis)
m5 <- lm(dist_mo_scale ~ species + tw_mo, data = dat_dis)
m6 <- lm(dist_mo_scale ~ mass_scale + species, data = dat_dis)
m7 <- lm(dist_mo_scale ~ mass_scale + species + tw_mo,data = dat_dis)
m8 <- lm(dist_mo_scale ~ mass_scale * species,data = dat_dis)
AIC_df <- data.frame(AICc(mo_mod, m1, m2, m3, m4, m5, m6, m7, m8))
AIC_df <- AIC_df[order(AIC_df$AICc),]
AIC_df$dAICc <- AIC_df$AICc-AIC_df$AICc[1]
AIC_df
#        df      AICc     dAICc
# m7      6  88.79835  0.000000
# m5      5  92.56062  3.762267
# m4      4  92.65420  3.855851
# m2      3  94.13295  5.334598
# mo_mod  8  94.13707  5.338718
# m1      3  94.16794  5.369583
# m6      5  99.52759 10.729234
# m3      4  99.93328 11.134925
# m8      7 105.00797 16.209612

# BEST MODEL JUST HAS SCALED MASS, TW AND SPECIES


mo_fin <- m7
summary(mo_fin)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -2.3531     0.6889  -3.416 0.001962 ** 
# mass_scale   -0.3720     0.1473  -2.525 0.017505 *  
# speciesGHA   -0.1790     0.3737  -0.479 0.635617    
# speciesWALB   1.4043     0.5305   2.647 0.013182 *  
# tw_mo         0.4642     0.1221   3.802 0.000712 ***

# RESIDUAL CHECKS
sim <- simulateResiduals(mo_fin)
plot(sim)
testDispersion(sim)
# all good