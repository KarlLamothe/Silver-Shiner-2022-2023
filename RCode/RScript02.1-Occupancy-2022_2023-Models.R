################################################################################
################################################################################
# Occupancy model analysis of Silver Shiner data using spOccupancy package
# RScript02.1-Occupancy-2022_2023-Models.R 
# This is Script 6 of occupancy modelling
# This script is used to build single species occupancy models for Silver Shiner
# using data collected from 16 Mile Creek in 2022 and 2023. 
################################################################################
################################################################################
# load packages
library(pacman)     # downloads and loads packages simultaneously
p_load(ggplot2)     # nice plots
p_load(spOccupancy) # occupancy models
p_load(xlsx)        # read in xlsx files
p_load(sf)          # coordinates to projections
p_load(oce)         # decimal degress to utm
p_load(beepr)       # alarm when models finish

# Personal ggplot theme
theme_set(theme_bw() + 
            theme(axis.title   = element_text(size=11, family="sans", colour="black"),
                  axis.text.x  = element_text(size=10, family="sans", colour="black"),
                  axis.text.y  = element_text(size=10, family="sans", colour="black"),
                  legend.text  = element_text(size=10, family="sans", colour="black"),
                  plot.title   = element_text(size=11, family="sans", colour="black"),
                  strip.text   = element_text(size=11, family="sans", colour="black"),
                  panel.border = element_rect(colour = "black", fill=NA),
                  axis.ticks   = element_line(colour = "black"),
                  legend.background = element_blank()))

################################################################################
# Prepare data for modelling
################################################################################
# read in detection and covariate data collected in 2022, 2023
SS_data_2022 <- read.xlsx("Data/Silver-Shiner-16Mile-2011-2022.xlsx", header=T, sheetName = "SS_SMC_2022")
SS_data_2022 <- SS_data_2022[order(SS_data_2022$SiteNo),]
SS_data_2023 <- read.xlsx("Data/Silver-Shiner-16Mile-2011-2022.xlsx", header=T, sheetName = "SS_SMC_2023")
SS_data_2023 <- SS_data_2023[order(SS_data_2023$SiteNo),]
SS_data_New  <- rbind(SS_data_2022, SS_data_2023)
colnames(SS_data_New)

# confirm structure of data
str(SS_data_New)

# convert year to factor
SS_data_New$Year <- as.character(SS_data_New$Year)
head(SS_data_New)

# Occupancy model covariates - year and depth
Covariates.std <- cbind.data.frame(Depth     = scale(SS_data_New$Depth, center=T),
                                   Discharge = scale(SS_data_New$Discharge_m3_s, center=T),
                                   Year      = SS_data_New$Year)
ggplot(Covariates.std, aes(x=Discharge, group=Year, fill = Year))+geom_density()

# Standardized site coordinates
Coordinates <- cbind.data.frame(Latitude = SS_data_New$Latitude,
                                Longitude = SS_data_New$Longitude)

# UTM
Coordinates.UTM<-data.frame(utm_x = lonlat2utm(SS_data_New$Longitude, SS_data_New$Latitude, 17, km = FALSE)$easting,
                            utm_y = lonlat2utm(SS_data_New$Longitude, SS_data_New$Latitude, 17, km = FALSE)$northing,
                            year  = SS_data_New$Year)

# extract capture data
SS.frame.abundance<-cbind.data.frame(Haul.1 = SS_data_New$Haul1, 
                                     Haul.2 = SS_data_New$Haul2, 
                                     Haul.3 = SS_data_New$Haul3)

# copy data to make presence/absence data frame
SS.frame.PA <- SS.frame.abundance
SS.frame.PA[SS.frame.PA>0]<-1 # convert to presence absence
head(SS.frame.PA)

# format data to consider depletion design.
for(z in 2:3){ # stop sampling after a detection has occurred
  SS.frame.PA[,z][SS.frame.PA[,(z-1)]>0]<-NA
  SS.frame.PA[,z][is.na(SS.frame.PA[,(z-1)])]<-NA
}
head(SS.frame.PA)

# naive occupancy per year
SS.PA <- SS.frame.abundance
SS.PA <- rowSums(SS.PA)
SS.PA[SS.PA > 0] <- 1        # convert to presence absence

# make data frame
Covariates <- cbind.data.frame(Covariates.std, 
                               Occur = SS.PA,
                               Count = rowSums(SS.frame.abundance))
aggregate(Covariates$Occur, list(Covariates$Year), sum)

# sites per year with detections
sum(Covariates$Occur[Covariates$Year=='2023'])/length(Covariates$Occur[Covariates$Year=='2023']) #0.4747475
sum(Covariates$Occur[Covariates$Year=='2022'])/length(Covariates$Occur[Covariates$Year=='2022']) #0.7425743

################################################################################
# Prepare and run occupancy models
################################################################################
# data for models
data.list <- c(y        = list(SS.frame.PA), 
               occ.covs = list(Covariates.std), 
               det.covs = list(Covariates.std),
               coords   = list(Coordinates.UTM[1:2]))

#####################################################################
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
###   non-spatial single species single season occupancy models   ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#####################################################################
set.seed(0528)
# Models
det.formula1 <- ~ 1
occ.formula1 <- ~ 1

det.formula2 <- ~ Discharge
occ.formula2 <- ~ 1

det.formula3 <- ~ Year
occ.formula3 <- ~ 1

det.formula4 <- ~ Depth
occ.formula4 <- ~ 1

# ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 1: occ ~ 1, det ~ 1 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~ # 
## model initial parameters, priors, and run parameters
out1.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))

## uninformative normal priors for the intercepts
out1.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = 0, var = 2.72))

## model run parameters
n.samples <- 50000; n.burn <- 25000; n.report <- 10000; n.chains <- 3; n.thin <- 25

# run model
(ptm <- Sys.time()) # time the run
out1 <- PGOcc(occ.formula = occ.formula1, det.formula = det.formula1, data = data.list, 
              inits = out1.inits, priors = out1.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains, n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 57.37981 secs
save(out1, file = paste("Results/2022-2023/Models/out1-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out1)
Model.coefs1 <- cbind.data.frame(Intercept.psi = c(0.5397, 0.1589, 0.2378, 0.5371, 0.8469, 1.0008, 2796),
                                 Intercept.p   = c(0.7628, 0.2149, 0.3326, 0.7629, 1.1707, 1.0024, 3000),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
# reshape the data into more usable format
Model.coefs1.m <- melt(Model.coefs1, id.vars=c("Summary"))
Model.coefs1.c <- dcast(Model.coefs1.m, variable~Summary, value.var="value")
Model.coefs1.c <- Model.coefs1.c[2:8]
Model.coefs1   <- cbind.data.frame(ModelNo = "Model 1", Model.coefs1.c,
                                   Variable = "Int",    Submodel = c("Occ","Det"),
                                   Type     = "Fixed")
Model.coefs1

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## MODEL 2: occ ~ 1, det ~ Discharge #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out2.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out2.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = 0, var = 2.72))

(ptm <- Sys.time()) # time the run
out2 <- PGOcc(occ.formula = occ.formula2, det.formula = det.formula2, data = data.list, 
              inits = out2.inits, priors = out2.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains, n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 58.02938  secs
save(out2, file = paste("Results/2022-2023/Models/out2-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out2)
Model.coefs2 <- cbind.data.frame(Intercept.psi = c( 0.7069, 0.1957,  0.3569,  0.6983,  1.1283, 1.0015, 3000),
                                 Intercept.p   = c( 0.5094, 0.2518,  0.0191,  0.5128,  0.9950, 1.0098, 3000),
                                 Discharge.p   = c(-0.9285, 0.2887, -1.5656, -0.9074, -0.4498, 1.0020, 2843),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs2.m <- melt(Model.coefs2, id.vars=c("Summary"))
Model.coefs2.c <- dcast(Model.coefs2.m, variable~Summary, value.var="value")
Model.coefs2.c <- Model.coefs2.c[2:8]
Model.coefs2   <- cbind.data.frame(ModelNo = "Model 2", Model.coefs2.c,
                                   Variable = c("Int","Int","Discharge"),
                                   Submodel = c("Occ","Det","Det"),
                                   Type     = "Fixed")
Model.coefs2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 3: occ ~ 1, det ~ Year #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out3.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out3.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = 0, var = 2.72))

(ptm <- Sys.time()) # time the run
out3 <- PGOcc(occ.formula = occ.formula3, det.formula = det.formula3, data = data.list, 
              inits = out3.inits, priors = out3.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains, n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 1.105818 mins
save(out3, file = paste("Results/2022-2023/Models/out3-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out3)
Model.coefs3 <- cbind.data.frame(Intercept.psi = c( 0.9704, 0.2390,  0.5241,  0.9588,  1.4564, 1.0019, 3000),
                                 Intercept.p   = c( 1.2907, 0.2833,  0.7198,  1.2892,  1.8512, 1.0036, 3000),
                                 Year2023.p    = c(-1.9295, 0.3927, -2.6652, -1.9327, -1.1221, 0.9997, 2777),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs3.m <- melt(Model.coefs3,   id.vars=c("Summary"))
Model.coefs3.c <- dcast(Model.coefs3.m, variable~Summary, value.var="value")
Model.coefs3.c <- Model.coefs3.c[2:8]
Model.coefs3   <- cbind.data.frame(ModelNo = "Model 3", Model.coefs3.c,
                                   Variable = c("Int","Int", "Year 2023"),
                                   Submodel = c("Occ","Det","Det"),
                                   Type     = "Fixed")
Model.coefs3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 3: occ ~ 1, det ~ Year #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out4.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out4.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = 0, var = 2.72))

(ptm <- Sys.time()) # time the run
out4 <- PGOcc(occ.formula = occ.formula4, det.formula = det.formula4, data = data.list, 
              inits = out4.inits, priors = out4.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains, n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 55.96973 secs
save(out4, file = paste("Results/2022-2023/Models/out4-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out4)

Model.coefs4 <- cbind.data.frame(Intercept.psi = c( 0.5391, 0.1612,  0.2389,  0.5374, 0.8559, 0.9997, 3000),
                                 Intercept.p   = c( 0.8124, 0.2267,  0.3505,  0.8132, 1.2533, 0.9999, 2757),
                                 Depth.p       = c(-0.2003, 0.1952, -0.5701, -0.2012, 0.1873, 1.0004, 3000),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs4.m <- melt(Model.coefs4,   id.vars=c("Summary"))
Model.coefs4.c <- dcast(Model.coefs4.m, variable~Summary, value.var="value")
Model.coefs4.c <- Model.coefs4.c[2:8]
Model.coefs4   <- cbind.data.frame(ModelNo = "Model 4", Model.coefs4.c,
                                   Variable = c("Int","Int", "Depth"),
                                   Submodel = c("Occ","Det","Det"),
                                   Type     = "Fixed")
Model.coefs4

################################################################################
# Detection model comparison
################################################################################
# compare models
n.mods <- 4
model.comp.variables <- function() {
  df <- data.frame(); res <- data.frame()
  
  for (i in 1:n.mods) {
    # standard model
    model_name <- paste("out", i, sep = "")
    model <- get(model_name)
    WAIC  <- waicOcc(model)[3]
    pD    <- waicOcc(model)[2]
    elpd  <- waicOcc(model)[1]
    dev   <- model$k.fold.deviance
    
    res   <- cbind.data.frame(WAIC = WAIC, pD = pD, 
                              elpd = elpd, deviance = dev)
    df    <- rbind(df, res)
  }
  
  df$ModelNo   <- c(paste("Model", seq(1:n.mods)))
  colnames(df) <- c("WAIC","pD","elpd","deviance","model")
  rownames(df) <- NULL
  return(df)
}

mod.comp<-model.comp.variables()
(mod.comp<-mod.comp[order(mod.comp$WAIC),])

# create model coefficient plot
ModelCoefs <- rbind(Model.coefs1, Model.coefs2, Model.coefs3, Model.coefs4)
ModelCoefs$spatial <- "Non-spatial" # add column to identify non-spatial mods

#########################################################
#########################################################
# Occupancy models
#########################################################
#########################################################
# Models
det.formula5 <- ~ Year
occ.formula5 <- ~ Depth

det.formula6 <- ~ Year
occ.formula6 <- ~ Year

det.formula7 <- ~ Year
occ.formula7 <- ~ Depth + Year

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## MODEL 5: occ ~ Depth, det ~ Year #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out5.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out5.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = c(0, 1), var = c(2.72, 1.50)))

(ptm <- Sys.time()) # time the run
out5 <- PGOcc(occ.formula = occ.formula5, det.formula = det.formula5, data = data.list, 
              inits = out5.inits, priors = out5.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains, n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 55.96973 secs
save(out5, file = paste("Results/2022-2023/Models/out5-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out5)
Model.coefs5 <- cbind.data.frame(Intercept.psi = c( 1.1502, 0.2764,  0.6323,  1.1398,  1.7254, 1.0011, 3000),
                                 Depth.psi     = c( 0.5480, 0.2777,  0.0364,  0.5364,  1.1242, 1.0014, 3000),
                                 Intercept.p   = c( 1.2961, 0.2829,  0.7504,  1.2977,  1.8539, 1.0025, 2890),
                                 Year2023.p    = c(-2.0547, 0.3767, -2.8036, -2.0565, -1.2949, 1.0005, 3000),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs5.m <- melt(Model.coefs5,   id.vars=c("Summary"))
Model.coefs5.c <- dcast(Model.coefs5.m, variable~Summary, value.var="value")
Model.coefs5.c <- Model.coefs5.c[2:8]
Model.coefs5   <- cbind.data.frame(ModelNo = "Model 5", Model.coefs5.c,
                                  Variable = c("Int","Depth","Int", "Year 2023"),
                                  Submodel = c("Occ","Occ","Det","Det"),
                                  Type     = "Fixed")
Model.coefs5

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## MODEL 6: occ ~ Year, det ~ Year #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out6.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out6.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = 0, var = 2.72))

(ptm <- Sys.time()) # time the run
out6 <- PGOcc(occ.formula = occ.formula6, det.formula = det.formula6, data = data.list, 
              inits = out6.inits, priors = out6.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains, n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 59.64559  secs
save(out6, file = paste("Results/2022-2023/Models/out6-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out6)
Model.coefs6 <- cbind.data.frame(Intercept.psi = c( 1.0922, 0.2410,  0.6359,  1.0863,  1.5721, 1.0090, 3000),
                                 Year2023.psi  = c(-0.6481, 0.6597, -1.5152, -0.7668,  1.1229, 1.0101, 1327),
                                 Intercept.p   = c( 1.3014, 0.2895,  0.7297,  1.3048,  1.8676, 1.0044, 2928),
                                 Year2023.p    = c(-1.4523, 0.5339, -2.5940, -1.4116, -0.5199, 1.0102, 2162),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs6.m <- melt(Model.coefs6,   id.vars=c("Summary"))
Model.coefs6.c <- dcast(Model.coefs6.m, variable~Summary, value.var="value")
Model.coefs6.c <- Model.coefs6.c[2:8]
Model.coefs6   <- cbind.data.frame(ModelNo = "Model 6",
                                   Model.coefs6.c,
                                   Variable = c("Int","Year 2023","Int","Year 2023"),
                                   Submodel = c("Occ","Occ","Det","Det"),
                                   Type     = "Fixed")
Model.coefs6

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 7: occ ~ Depth + Year, det ~ Year #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out7.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out7.priors <- list(alpha.normal  = list(mean = 0, var = 2.72), 
                    beta.normal   = list(mean = c(0, 1, 0), var = c(2.72, 1.50, 2.72)))

(ptm <- Sys.time()) # time the run
out7 <- PGOcc(occ.formula = occ.formula7, det.formula = det.formula7, data = data.list, 
              inits = out7.inits, priors = out7.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains, n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 1.176672 mins
save(out7, file = paste("Results/2022-2023/Models/out7-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out7)

Model.coefs7 <- cbind.data.frame(Intercept.psi = c( 1.2106, 0.2744,  0.7124,  1.1988,  1.7704, 1.0120, 2213),
                                 Depth.psi     = c( 0.7196, 0.3836,  0.0347,  0.7159,  1.4920, 1.0062, 1101),
                                 Year2023.psi  = c( 1.0618, 1.3499, -1.1800,  1.0537,  3.7617, 1.0108,  788),
                                 Intercept.p   = c( 1.2817, 0.2807,  0.7372,  1.2765,  1.8364, 1.0002, 3137),
                                 Year2023.p    = c(-2.2203, 0.5354, -3.0825, -2.3029, -0.9214, 1.0023, 1007),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs7.m <- melt(Model.coefs7, id.vars=c("Summary"))
Model.coefs7.c <- dcast(Model.coefs7.m, variable~Summary, value.var="value")
Model.coefs7.c <- Model.coefs7.c[2:8]
Model.coefs7   <- cbind.data.frame(ModelNo = "Model 7",
                                   Model.coefs7.c,
                                   Variable = c("Int","Depth","Year 2023","Int","Year 2023"),
                                   Submodel = c("Occ","Occ","Occ","Det","Det"),
                                   Type     = c("Fixed","Fixed","Fixed","Fixed","Fixed"))
Model.coefs7

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# create data frame with all model summaries #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
ModelCoefs.1 <- rbind(Model.coefs1, Model.coefs2, Model.coefs3, Model.coefs4, 
                      Model.coefs5, Model.coefs6, Model.coefs7)
ModelCoefs.1$spatial <- "Non-spatial" # add column to identify non-spatial mods

n.mods <- 7
model.comp.variables <- function() {
  df <- data.frame(); res <- data.frame()
  
  for (i in 1:n.mods) {
    # standard model
    model_name <- paste("out", i, sep = "")
    model <- get(model_name)
    WAIC  <- waicOcc(model)[3]
    pD    <- waicOcc(model)[2]
    elpd  <- waicOcc(model)[1]
    dev   <- model$k.fold.deviance
    res   <- cbind.data.frame(WAIC = WAIC, pD = pD, 
                              elpd = elpd, deviance = dev)
    df    <- rbind(df, res)
  }
  
  df$ModelNo   <- c(paste("Model", seq(1:n.mods)))
  colnames(df) <- c("WAIC","pD","elpd","deviance","model")
  rownames(df) <- NULL
  return(df)
}

mod.comp.fin<-model.comp.variables()
(mod.comp.fin<-mod.comp.fin[order(mod.comp.fin$WAIC),])

# export model selection table
write.csv(mod.comp.fin, "Results/2022-2023/Models/Model-selection_nonspatial-2022-2023.csv")

          #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
          ####     spatial single species single season model    ####
          #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

batch.length <- 25; n.burn <- 5000; n.report <- 100; n.chains <- 3; n.thin <- 20; n.batch <- 1000

# set up priors and inits
dist.SS   <- dist(Coordinates.UTM) # Pair-wise distances between all sites
cov.model <- "exponential"     # Exponential covariance model

# initial values
SS.inits  <- list(alpha = 0, beta  = 0, 
                  z = apply(data.list$y, 1, max, na.rm = TRUE), 
                  sigma.sq = 2, 
                  phi = 3 / mean(dist.SS), 
                  w = rep(0, nrow(data.list$y)))

# priors
min.dist  <- min(dist.SS)
max.dist  <- max(dist.SS)
SS.priors1 <- list(beta.normal  = list(mean = 0, var = 2.72), # uninformative
                   alpha.normal = list(mean = 0, var = 2.72), # uninformative
                   sigma.sq.ig  = c(2, 1), 
                   phi.unif     = c(3/max.dist, 3/min.dist))

################################################################################
# Occupancy models
################################################################################
det.formula.sp1 <- ~ Year
occ.formula.sp1 <- ~ 1

det.formula.sp2 <- ~ Year
occ.formula.sp2 <- ~ Depth

det.formula.sp3 <- ~ Year
occ.formula.sp3 <- ~ Year

det.formula.sp4 <- ~ Year
occ.formula.sp4 <- ~ Depth + Year

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## spatial model 1: occ ~ 1, det ~ Year ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
SS.priors1 <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal  = list(mean = 0, var = 2.72),
                   sigma.sq.ig  = c(2, 1), 
                   phi.unif     = c(3/max.dist, 3/min.dist))

(start <- Sys.time())
out.sp1 <- spPGOcc(occ.formula = occ.formula.sp1, det.formula = det.formula.sp1, cov.model=cov.model, NNGP = TRUE,
                   data = data.list, n.batch = n.batch, batch.length = batch.length, n.omp.threads = 3,
                   n.chains = n.chains, n.burn = n.burn, n.report = n.report, n.neighbors = 5,
                   n.thin = n.thin, priors = SS.priors1, inits = SS.inits,
                   k.fold = 10, k.fold.threads = 10, verbose = TRUE)
(Sys.time()-start) #1.647771 
save(out.sp1, file = paste("Results/2022-2023/Models/out.sp1-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out.sp1)

Model.coefs.sp1 <- cbind.data.frame(Intercept.psi = c( 1.1369, 0.3193,  0.5983,  1.1087  ,1.8657, 1.0078, 671),
                                    Intercept.p   = c( 1.2962, 0.2821,  0.7376,  1.3002  ,1.8411, 1.0009, 3000),
                                    Year2023.p    = c(-1.9207, 0.3957, -2.6657, -1.9331, -1.1135, 1.0017, 2819),
                                    sigma.sq      = c( 1.1895, 1.9340,  0.1923,  0.6757  ,5.0343, 1.1612, 252),
                                    phi           = c( 1.0707, 0.7071,  0.0471,  1.0531  ,2.2760, 1.0453, 768),
                                    Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs.sp1m <- melt(Model.coefs.sp1,   id.vars=c("Summary"))
Model.coefs.sp1c <- dcast(Model.coefs.sp1m, variable~Summary, value.var="value")
Model.coefs.sp1c <- Model.coefs.sp1c[2:8]
Model.coefs.sp1  <- cbind.data.frame(ModelNo  = "Model 1 - spatial", Model.coefs.sp1c,
                                     Variable = c("Int","Int","Year 2023","sigma^2","phi"),
                                     Submodel = c("Occ","Det","Det","Spatial", "Spatial"),
                                     Type     = "Fixed")
Model.coefs.sp1

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## spatial model 2: occ ~ Depth, det ~ Year ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
SS.priors2 <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal  = list(mean = c(0, 1), var = c(2.72, 1.50)),
                   sigma.sq.ig  = c(2, 1), 
                   phi.unif     = c(3/max.dist, 3/min.dist))

(start <- Sys.time())
out.sp2 <- spPGOcc(occ.formula = occ.formula.sp2, det.formula = det.formula.sp2, cov.model=cov.model, NNGP = TRUE,
                   data = data.list, n.batch = n.batch, batch.length = batch.length, n.omp.threads = 3,
                   n.chains = n.chains, n.burn = n.burn, n.report = n.report, n.neighbors = 5,
                   n.thin = n.thin, priors = SS.priors2, inits = SS.inits,
                   k.fold = 10, k.fold.threads = 10, verbose = TRUE)
(Sys.time()-start); beep(sound = "coin", expr = NULL) # 1.878988  mins
save(out.sp2, file = paste("Results/2022-2023/Models/out.sp2-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out.sp2)

Model.coefs.sp2 <- cbind.data.frame(Intercept.psi = c( 1.3390, 0.3556,  0.7299,  1.3134,  2.1274, 1.0121, 1219),
                                   Depth.psi      = c( 0.6302, 0.3238,  0.0613,  0.6060,  1.3180, 1.0040, 2135),
                                   Intercept.p    = c( 1.3048, 0.2890,  0.7340,  1.2975,  1.8807, 1.0006, 3000),
                                   Year2023.p     = c(-2.0537, 0.3671, -2.7764, -2.0608, -1.3139, 1.0024, 3151),
                                   sigma.sq       = c( 1.1324, 1.5500,  0.1983,  0.6771,  5.0338, 1.1929, 401),
                                   phi            = c( 1.0237, 0.7193,  0.0184,  0.9527,  2.2783, 1.0048, 431),
                                   Summary        = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs.sp2m <- melt(Model.coefs.sp2,   id.vars=c("Summary"))
Model.coefs.sp2c <- dcast(Model.coefs.sp2m, variable~Summary, value.var="value")
Model.coefs.sp2c <- Model.coefs.sp2c[2:8]
Model.coefs.sp2  <- cbind.data.frame(ModelNo  = "Model 2 - spatial", Model.coefs.sp2c,
                                     Variable = c("Int","Depth","Int","Year 2023","sigma^2","phi"),
                                     Submodel = c("Occ","Occ","Det","Det","Spatial", "Spatial"),
                                     Type     = "Fixed")
Model.coefs.sp2

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## spatial model 3: occ ~ Year,  det ~ Year ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
SS.priors3 <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal  = list(mean = 0, var = 2.72),
                   sigma.sq.ig  = c(2, 1), 
                   sigma.sq.p.ig = list(shape = 0.1, scale = 0.1),
                   phi.unif     = c(3/max.dist, 3/min.dist))

(start <- Sys.time())
out.sp3 <- spPGOcc(occ.formula = occ.formula.sp3, det.formula = det.formula.sp3, cov.model=cov.model, NNGP = TRUE,
                   data = data.list, n.batch = n.batch, batch.length = batch.length, n.omp.threads = 3,
                   n.chains = n.chains, n.burn = n.burn, n.report = n.report, n.neighbors = 5,
                   n.thin = n.thin, priors = SS.priors3, inits = SS.inits,
                   k.fold = 10, k.fold.threads = 10, verbose = TRUE)
(Sys.time()-start); beep(sound = "coin", expr = NULL)# 1.652926 mins
save(out.sp3, file = paste("Results/2022-2023/Models/out.sp3-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out.sp3)
Model.coefs.sp3 <- cbind.data.frame(Intercept.psi = c( 1.2855, 0.3481,  0.7145,  1.2569,  2.0972, 1.0020,  557),
                                    Year2023.psi  = c(-0.7830, 0.7139, -1.8554, -0.8732,  1.1004, 1.0024, 1335),
                                    Intercept.p   = c( 1.3178, 0.2820,  0.7740,  1.3178,  1.8922, 1.0018, 3000),
                                    Year2023.p    = c(-1.4599, 0.5230, -2.5334, -1.4206, -0.5097, 1.0001, 1970),
                                    sigma.sq      = c( 1.2780, 2.0736,  0.1899,  0.7062,  6.0240, 1.1030, 245),
                                    phi           = c( 1.0150, 0.7118,  0.0218,  0.9563,  2.2593, 1.0101, 376),
                                    Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs.sp3m <- melt(Model.coefs.sp3,   id.vars=c("Summary"))
Model.coefs.sp3c <- dcast(Model.coefs.sp3m, variable~Summary, value.var="value")
Model.coefs.sp3c <- Model.coefs.sp3c[2:8]
Model.coefs.sp3  <- cbind.data.frame(ModelNo = "Model 3 - spatial", Model.coefs.sp3c,
                                     Variable = c("Int","Year 2023","Int","Year 2023","sigma^2","phi"),
                                     Submodel = c("Occ","Occ","Det","Det","Spatial", "Spatial"),
                                     Type = "Fixed")
Model.coefs.sp3

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## spatial model 4: occ ~ Depth + Year, det ~ Year ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
SS.priors4 <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal  = list(mean = c(0, 1, 0), var = c(2.72, 1.50, 2.72)),
                   sigma.sq.ig  = c(2, 1), 
                   sigma.sq.p.ig = list(shape = 0.1, scale = 0.1),
                   phi.unif     = c(3/max.dist, 3/min.dist))

(start <- Sys.time())
out.sp4 <- spPGOcc(occ.formula = occ.formula.sp4, det.formula = det.formula.sp4, cov.model=cov.model, NNGP = TRUE,
                   data = data.list, n.batch = n.batch, batch.length = batch.length, n.omp.threads = 3,
                   n.chains = n.chains, n.burn = n.burn, n.report = n.report, n.neighbors = 5,
                   n.thin = n.thin, priors = SS.priors4, inits = SS.inits,
                   k.fold = 10, k.fold.threads = 10, verbose = TRUE)
(Sys.time()-start); beep(sound = "coin", expr = NULL)# 1.870804  mins
save(out.sp4, file = paste("Results/2022-2023/Models/out.sp4-2022_2023-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out.sp4)
Model.coefs.sp4 <- cbind.data.frame(Intercept.psi = c( 1.3935, 0.3727,  0.7745,  1.3550,  2.2078, 1.0035, 718),
                                    Depth         = c( 0.7238, 0.4358, -0.0044,  0.7005,  1.6251, 1.0030, 955),
                                    Year2023.psi  = c( 0.6251, 1.3880, -1.5420,  0.4728,  3.5119, 1.0059, 627),
                                    Intercept.p   = c( 1.2868, 0.2869,  0.7446,  1.2774,  1.8467, 1.0015, 3000),
                                    Year2023.p    = c(-2.0721, 0.5719, -3.0243, -2.1607 ,-0.8377, 1.0011,  856),
                                    sigma.sq      = c( 1.2242, 1.8040,  0.1906,  0.7214,  5.4297, 1.0236,  357),
                                    phi           = c( 1.1051, 0.6836,  0.0523,  1.0782,  2.2764, 1.0061, 1073),
                                    Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs.sp4m <- melt(Model.coefs.sp4,   id.vars=c("Summary"))
Model.coefs.sp4c <- dcast(Model.coefs.sp4m, variable~Summary, value.var="value")
Model.coefs.sp4c <- Model.coefs.sp4c[2:8]
Model.coefs.sp4  <- cbind.data.frame(ModelNo = "Model 4 - spatial", Model.coefs.sp4c,
                                     Variable = c("Int","Depth","Year 2023","Int","Year 2023","sigma^2","phi"),
                                     Submodel = c("Occ","Occ","Occ","Det","Det","Spatial","Spatial"),
                                     Type = "Fixed")
Model.coefs.sp4

################################################################################
################################################################################
n.mods <- 4
model.comp.variables.sp <- function() {
  df <- data.frame(); res <- data.frame()
  
  for (i in 1:n.mods) {
    
    # spatial model
    model_name <- paste("out.sp", i, sep = "")
    model <- get(model_name)
    WAIC  <- waicOcc(model)[3]
    pD    <- waicOcc(model)[2]
    elpd  <- waicOcc(model)[1]
    dev   <- model$k.fold.deviance
    res   <- cbind.data.frame(WAIC = WAIC, pD = pD, 
                              elpd = elpd, deviance = dev)
    df    <- rbind(df, res)
  }
  
  df$ModelNo   <- c(paste("Spatial model", seq(1:n.mods)))
  colnames(df) <- c("WAIC","pD","elpd","deviance","model")
  rownames(df) <- NULL
  return(df)
}

mod.comp.sp.fin<-model.comp.variables.sp()
(mod.comp.sp.fin<-mod.comp.sp.fin[order(mod.comp.sp.fin$WAIC),])

write.csv(mod.comp.sp.fin, "Results/2022-2023/Models/Model_selection_spatial_2022_2023.csv")

Model.comparison<-rbind.data.frame(mod.comp.sp.fin, mod.comp.fin)
(Model.comparison<-Model.comparison[order(Model.comparison$WAIC),])

write.csv(Model.comparison, "Results/2022-2023/Models/Model_selection_total_2022_2023.csv")

################################################################################
################################################################################
ModelCoefs.tot2 <- rbind(Model.coefs.sp1, Model.coefs.sp2, Model.coefs.sp3, Model.coefs.sp4)
ModelCoefs.tot2$spatial <- "Spatial"

# export model selection table
ModelCoefs <- rbind(ModelCoefs.1, ModelCoefs.tot2)
write.csv(ModelCoefs, "Results/2022-2023/Models/Model_coefficients_2022_2023.csv")
