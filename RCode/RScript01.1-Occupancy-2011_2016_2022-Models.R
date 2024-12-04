################################################################################
################################################################################
# Occupancy model analysis of Silver Shiner data using spOccupancy package
# RScript01.1-Occupancy-2011_2016_2022-Models.R 
# This is Script 1 of occupancy modelling
# This script is used to build single species occupancy models for Silver Shiner
# using data collected from 16 Mile Creek in 2011, 2016, and 2022. 
################################################################################
################################################################################
# load packages
library(pacman)     # downloads and loads packages simultaneously
p_load(spOccupancy) # occupancy models
p_load(oce)         # convert coords to UTM
p_load(xlsx)        # import xlsx files
p_load(beepr)       # alarm when models finish

# Personal ggplot theme
theme_set(theme_bw() + 
            theme(axis.title   = element_text(size=9, family="sans", colour="black"),
                  axis.text.x  = element_text(size=9, family="sans", colour="black"),
                  axis.text.y  = element_text(size=9, family="sans", colour="black"),
                  legend.text  = element_text(size=9, family="sans", colour="black"),
                  plot.title   = element_text(size=10, family="sans", colour="black"),
                  strip.text   = element_text(size=9, family="sans", colour="black"),
                  panel.border = element_rect(colour = "black", fill=NA),
                  axis.ticks   = element_line(colour = "black"),
                  legend.title = element_blank()))

################################################################################
# Prepare data for modelling
################################################################################
# read in detection and covariate data collected in 2011, 2016, 2022
SS_data <- read.xlsx("Data/Silver-Shiner-16Mile-2011-2022.xlsx", header=T, sheetName = "2011_2016_2022")
str(SS_data)

# convert year to factor
SS_data$Year <- as.character(SS_data$Year)

# extract depth and year - primary variables of interest
colnames(SS_data)

# Occupancy model covariates - year and depth
Covariates <- cbind.data.frame(Depth = SS_data$Mean_Depth_m,
                               Year  = SS_data$Year)
mean(Covariates$Depth)

# Standardized site coordinates
Coordinates <- cbind.data.frame(Latitude = SS_data$Latitude,
                                Longitude = SS_data$Longitude)

# UTM
Coordinates.UTM<-data.frame(utm_x = lonlat2utm(SS_data$Longitude, SS_data$Latitude, 17, km = FALSE)$easting,
                            utm_y = lonlat2utm(SS_data$Longitude, SS_data$Latitude, 17, km = FALSE)$northing,
                            year  = SS_data$Year)

# Standardized depth
Covariates.std <- cbind.data.frame(Depth = scale(Covariates$Depth, center=T),
                                   Year  = Covariates$Year)

# extract capture data
SS.frame.abundance<-cbind.data.frame(Haul.1 = SS_data$Haul1, 
                                     Haul.2 = SS_data$Haul2, 
                                     Haul.3 = SS_data$Haul3)

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
SS.PA      <- cbind.data.frame(SS.PA, Year = Covariates$Year)
Covariates <- cbind.data.frame(Covariates, SS.PA)

# sites per year with detections
aggregate(SS.PA$SS.PA, list(SS.PA$Year), sum)
aggregate(SS.PA$SS.PA, list(SS.PA$Year), length)

#proportion
sum(SS.PA$SS.PA[SS.PA$Year=='2022'])/length(SS.PA$SS.PA[SS.PA$Year=='2022']) #0.7425743
sum(SS.PA$SS.PA[SS.PA$Year=='2016'])/length(SS.PA$SS.PA[SS.PA$Year=='2016']) #0.3636364
sum(SS.PA$SS.PA[SS.PA$Year=='2011'])/length(SS.PA$SS.PA[SS.PA$Year=='2011']) #0.5

################################################################################
################################################################################
#                   Prepare and run occupancy models                           #
################################################################################
################################################################################
# data for models
data.list <- c(y        = list(SS.frame.PA), 
               occ.covs = list(Covariates.std), 
               det.covs = list(Covariates.std),
               coords   = list(Coordinates.UTM[1:2]))

##############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# regular single species single season model #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
##############################################
# define models
occ.formula1 <- ~ 1
det.formula1 <- ~ 1

occ.formula2 <- ~ 1
det.formula2 <- ~ Depth

occ.formula3 <- ~ 1
det.formula3 <- ~ Year

## ~~~~~~~~~~~~~~~~~~~~~~~~~ #
## MODEL 1: det ~ 1; occ ~ 1 #
## ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# model initial parameters, priors, and run parameters
# Initial values 0 for detection (alpha) and occupancy (beta) covariates, and 
# detection history for occurrence
out1.inits <- list(alpha = 0, 
                   beta  = 0, 
                   z     = apply(data.list$y, 1, max, na.rm = TRUE))

# model priors - uninformed priors for detection and occupancy other than for
# depth, which has a reduced variance and effect size of 1. Normal priors
out1.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = 0, var = 2.72))
n.samples <- 50000; n.burn <- 25000; n.report <- 10000; n.chains <- 3; n.thin <- 25

(ptm <- Sys.time()) # time the run
out1 <- PGOcc(occ.formula = occ.formula1, det.formula = det.formula1, data = data.list, 
              inits = out1.inits, priors = out1.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains,n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 40.90425 secs
save(out1, file = paste("Results/2011-2016-2022/Models/out1-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out1)

# create data frame with model estimates. I just copied and pasted the values in
# manually from the console, and then added commas and cleaned it up. There has to
# be a better way your values may differ slightly from the values below
Model.coefs1 <- cbind.data.frame(Intercept.psi = c(0.6269, 0.1757, 0.2696, 0.6277, 0.9782, 1.0024, 3000),
                                 Intercept.p   = c(1.3831, 0.2685, 0.8822, 1.3764, 1.9138, 1.0006, 3000),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs1m <- melt(Model.coefs1,   id.vars=c("Summary"))
Model.coefs1c <- dcast(Model.coefs1m, variable~Summary, value.var="value")
Model.coefs1c <- Model.coefs1c[2:8]
Model.coefs1  <- cbind.data.frame(ModelNo  = "Model 1", Model.coefs1c,
                                  Variable = "Int", Submodel = c("Occ","Det"), 
                                  Type = "Fixed")
Model.coefs1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 2: det ~ Depth; occ ~ 1  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out2.inits <- list(alpha = 0, beta = 0, z = apply(data.list$y, 1, max, na.rm = TRUE))
out2.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = 0, var = 2.72))

(ptm <- Sys.time()) # time the run
out2 <- PGOcc(occ.formula = occ.formula2, det.formula = det.formula2, data = data.list, 
              inits = out2.inits, priors = out1.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains,n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 39.0303 secs
save(out2, file = paste("Results/2011-2016-2022/Models/out2-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out2)
Model.coefs2 <- cbind.data.frame(Intercept.psi = c( 0.6149, 0.1735,  0.2832,  0.6138 ,0.9527 ,1.0052 ,3000),
                                 Intercept.p   = c( 1.6147, 0.3023,  1.0400,  1.6054 ,2.2242 ,1.0039 ,2842),
                                 Depth.p       = c(-0.3806, 0.2563, -0.8727, -0.3849, 0.1218, 1.0026, 2531),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs2m <- melt(Model.coefs2,   id.vars=c("Summary"))
Model.coefs2c <- dcast(Model.coefs2m, variable~Summary, value.var="value")
Model.coefs2c <- Model.coefs2c[2:8]
Model.coefs2  <- cbind.data.frame(ModelNo  = "Model 2", Model.coefs2c,
                                  Variable = c("Int","Int", "Depth"),
                                  Submodel = c("Occ","Det","Det"),
                                  Type     = c("Fixed", "Fixed", "Fixed"))
Model.coefs2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 3: det ~ Year; occ ~ 1  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out3.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out3.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = 0, var = 2.72))

(ptm <- Sys.time()) # time the run
out3 <- PGOcc(occ.formula = occ.formula3, det.formula = det.formula3, data = data.list, 
              inits = out3.inits, priors = out3.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains,n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 41.55321 secs
save(out3, file = paste("Results/2011-2016-2022/Models/out3-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out3)
Model.coefs3 <- cbind.data.frame(Intercept.psi = c( 0.7914, 0.2235,  0.3941,  0.7822, 1.2553, 1.0053, 3000),
                                 Intercept.p   = c( 0.8981, 0.7960, -0.6346,  0.9070, 2.3941, 1.0065, 3000),
                                 Year2016.p    = c(-1.1501, 1.0951, -3.1001, -1.2459, 1.2579, 1.0071, 3000),
                                 Year2022.p    = c( 0.5186, 0.8132, -1.0156,  0.5088, 2.0860, 1.0089, 3000),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs3m <- melt(Model.coefs3,   id.vars=c("Summary"))
Model.coefs3c <- dcast(Model.coefs3m, variable~Summary, value.var="value")
Model.coefs3c <- Model.coefs3c[2:8]
Model.coefs3  <- cbind.data.frame(ModelNo  = "Model 3", Model.coefs3c,
                                  Variable = c("Int",
                                               "Int", "Year 2016","Year 2022"),
                                  Submodel = c("Occ","Det","Det","Det"),
                                  Type     = c("Fixed", "Fixed", "Fixed","Fixed"))
Model.coefs3

################################################################################
################################################################################
# Compare model output
################################################################################
################################################################################
n.mods <- 3
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
    res   <- cbind.data.frame(WAIC = WAIC, pD = pD, elpd = elpd, deviance = dev)
    df    <- rbind(df, res)
  }
  
  df$ModelNo   <- c(paste("Model", seq(1:n.mods)))
  colnames(df) <- c("WAIC","pD","elpd","deviance","model")
  rownames(df) <- NULL
  return(df)
}

mod.comp<-model.comp.variables()
mod.comp[order(mod.comp$WAIC),] 
# Inclusion of covariates does not improve model performance

################################################################################
################################################################################

# Specify occupancy models
occ.formula4 <- ~ Depth
det.formula4 <- ~ 1

occ.formula5 <- ~ Year
det.formula5 <- ~ 1

occ.formula6 <- ~ Depth + Year 
det.formula6 <- ~ 1

############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 4: det ~ 1; occ ~ Depth #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out4.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out4.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = c(0, 1), var  = c(2.72, 1.50)))

# run model
(ptm <- Sys.time()) # time the run
out4 <- PGOcc(occ.formula = occ.formula4, det.formula = det.formula4, data = data.list, 
              inits = out4.inits, priors = out4.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains,n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 43.0126 secs
save(out4, file = paste("Results/2011-2016-2022/Models/out4-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out4)
Model.coefs4 <- cbind.data.frame(Intercept.psi = c(0.7982, 0.2074, 0.4127, 0.7907, 1.2234, 1.0003, 3140),
                                 Depth         = c(1.0556, 0.2593, 0.5576, 1.0440, 1.5784, 1.0013, 3000),
                                 Intercept.p   = c(1.3633, 0.2684, 0.8473, 1.3562, 1.8982, 1.0017, 3000),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs4m <- melt(Model.coefs4,   id.vars=c("Summary"))
Model.coefs4c <- dcast(Model.coefs4m, variable~Summary, value.var="value")
Model.coefs4c <- Model.coefs4c[2:8]
Model.coefs4  <- cbind.data.frame(ModelNo = "Model 4", Model.coefs4c,
                                  Variable = c("Int","Depth","Int"),
                                  Submodel = c("Occ","Occ","Det"),
                                  Type     = "Fixed")
Model.coefs4

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 5: det ~ 1; occ ~ Year #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out5.inits <- list(alpha = 0, beta = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out5.priors <- list(alpha.normal  = list(mean = 0, var = 2.72), 
                    beta.normal   = list(mean = 0, var = 2.72))

(ptm <- Sys.time()) # time the run
out5 <- PGOcc(occ.formula = occ.formula5, det.formula = det.formula5, data = data.list, 
              inits = out5.inits, priors = out5.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains,n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 38.81981 secs
save(out5, file = paste("Results/2011-2016-2022/Models/out5-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out5)
Model.coefs5 <- cbind.data.frame(Intercept.psi = c( 0.0506, 0.3943, -0.7035,  0.0449, 0.8530, 1.0076, 2904),
                                 Year2016.psi  = c(-0.5683, 0.5731, -1.7375, -0.5507, 0.4834, 1.0053, 2442),
                                 Year2022.psi  = c( 1.0430, 0.4491,  0.1410,  1.0443, 1.9229, 1.0077, 3234),
                                 Intercept.p   = c( 1.3753, 0.2701,  0.8533,  1.3636, 1.9158, 1.0081, 2225),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs5m <- melt(Model.coefs5,   id.vars=c("Summary"))
Model.coefs5c <- dcast(Model.coefs5m, variable~Summary, value.var="value")
Model.coefs5c <- Model.coefs5c[2:8]
Model.coefs5  <- cbind.data.frame(ModelNo  = "Model 5", Model.coefs5c,
                                  Variable = c("Int", "Year 2016", "Year 2022", "Int"),
                                  Submodel = c("Occ","Occ","Occ","Det"),
                                  Type     = "Fixed")
Model.coefs5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 6: det ~ 1; occ ~ Depth + Year 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
out6.inits <- list(alpha = 0, beta  = 0, 
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
out6.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal  = list(mean = c(0, 1, 0, 0), 
                                        var  = c(2.72, 1.50, 2.72, 2.72)))

(ptm <- Sys.time()) # time the run
out6 <- PGOcc(occ.formula = occ.formula6, det.formula = det.formula6, data = data.list, 
              inits = out6.inits, priors = out6.priors, n.samples = n.samples,
              n.omp.threads = 3, n.chains = n.chains, n.burn = n.burn, n.report = n.report, 
              n.thin = n.thin, verbose = TRUE, k.fold = 10, k.fold.threads = 10)
Sys.time() - ptm # 39.95565 secs
save(out6, file = paste("Results/2011-2016-2022/Models/out6-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out6)
Model.coefs6 <- cbind.data.frame(Intercept.psi = c(-0.0220, 0.4434, -0.9259, -0.0128, 0.8302, 1.0009, 3000),
                                 Depth         = c( 1.0162, 0.2557,  0.5446,  1.0057, 1.5585, 1.0001, 3000),
                                 Year2016.psi  = c(-0.1896, 0.6332, -1.4193, -0.1988, 1.0585, 1.0008, 2710),
                                 Year2022.psi  = c( 1.2112, 0.4868,  0.2734,  1.2100, 2.1768, 1.0035, 3000),
                                 Intercept.p   = c( 1.3694, 0.2643,  0.8505,  1.3663, 1.9024, 1.0047, 3000),
                                 Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs6m <- melt(Model.coefs6,   id.vars=c("Summary"))
Model.coefs6c <- dcast(Model.coefs6m, variable~Summary, value.var="value")
Model.coefs6c <- Model.coefs6c[2:8]
Model.coefs6  <- cbind.data.frame(ModelNo  = "Model 6", Model.coefs6c,
                                  Variable = c("Int", "Depth", "Year 2016", "Year 2022",
                                               "Int"),
                                  Submodel = c("Occ","Occ","Occ","Occ","Det"),
                                  Type     = "Fixed")
Model.coefs6

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# create data frame with all model summaries #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
ModelCoefs.1 <- rbind(Model.coefs1, Model.coefs2, Model.coefs3, 
                      Model.coefs4, Model.coefs5, Model.coefs6)
ModelCoefs.1$spatial <- "Non-spatial"

#################################################################################
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
##                         Model Comparisons                                    #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#################################################################################

n.mods <- 6
mod.comp<-model.comp.variables()
(mod.comp<-mod.comp[order(mod.comp$WAIC),])

# Model 6 with lowest WAIC
write.csv(mod.comp,"Results/2011-2016-2022/Models/Model.comparisons.nonspatial.csv")

(ggplot(mod.comp, aes(y=WAIC, x = as.character(model)))+
    geom_point() + 
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()))/
  (ggplot(mod.comp, aes(y=deviance, x = as.character(model)))+
     geom_point() + 
     labs(y="Deviance", x="Model Number")+
     theme(legend.position = 'none',
           axis.text.x = element_text(angle=45, vjust=1, hjust=1),
           axis.title.x=element_blank()))

################################################################################
################################################################################
          ###########################################################
          ###########################################################
          #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
          ####     spatial single species single season model    ####
          #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
          ###########################################################
          ###########################################################
################################################################################
################################################################################
# Occupancy models; no need to redo the detection models from before
################################################################################
################################################################################
det.formula.sp1 <- ~ 1
occ.formula.sp1 <- ~ 1

det.formula.sp2 <- ~ 1
occ.formula.sp2 <- ~ Depth

det.formula.sp3 <- ~ 1
occ.formula.sp3 <- ~ Year

det.formula.sp4 <- ~ 1
occ.formula.sp4 <- ~ Depth + Year

# run parameters
batch.length <- 25; n.burn <- 5000; n.report <- 100; 
n.chains <- 3; n.thin <- 20; n.batch <- 1000

# set up priors and inits
dist.SS   <- dist(Coordinates.UTM) # Pair-wise distances between all sites
cov.model <- "exponential"         # Exponential covariance model

# initial values
SS.inits  <- list(alpha = 0, beta  = 0, 
                  z = apply(data.list$y, 1, max, na.rm = TRUE), 
                  sigma.sq = 2, 
                  phi = 3 / mean(dist.SS), 
                  w = rep(0, nrow(data.list$y)))

# priors
min.dist  <- min(dist.SS)
max.dist  <- max(dist.SS)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## spatial model 1: occ ~ 1, det ~ 1 ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
SS.priors1 <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal  = list(mean = 0, var = 2.72),
                   sigma.sq.ig  = c(2, 1), 
                   phi.unif     = c(3/max.dist, 3/min.dist))

(start <- Sys.time())
out.sp1 <- spPGOcc(occ.formula = occ.formula.sp1, det.formula = det.formula.sp1, 
                   cov.model=cov.model,NNGP = TRUE, data = data.list, 
                   n.batch = n.batch, batch.length = batch.length, 
                   n.omp.threads = 3, n.neighbors = 5, n.chains = n.chains, 
                   n.burn = n.burn, n.report = n.report, n.thin = n.thin, 
                   priors = SS.priors1, inits = SS.inits,
                   k.fold = 10, k.fold.threads = 10, verbose = TRUE)
(Sys.time()-start); beep(sound = "coin", expr = NULL) #1.235045 mins
save(out.sp1, file = paste("Results/2011-2016-2022/Models/outsp1-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out.sp1)
Model.coefs.sp1 <- cbind.data.frame(Intercept.psi = c(0.7353, 0.2312, 0.3196, 0.7175, 1.2309, 1.0016, 1652),
                                    Intercept.p   = c(1.3861, 0.2626, 0.8821, 1.3830, 1.8999, 1.0049, 3000),
                                    sigma.sq      = c(0.9257, 1.1327, 0.1802, 0.5862, 3.9191, 1.0076,  594),
                                    phi           = c(0.2674, 0.1452, 0.0161, 0.2773, 0.4967, 1.0112, 1134),
                                    Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs.sp1m <- melt(Model.coefs.sp1,   id.vars=c("Summary"))
Model.coefs.sp1c <- dcast(Model.coefs.sp1m, variable~Summary, value.var="value")
Model.coefs.sp1c <- Model.coefs.sp1c[2:8]
Model.coefs.sp1  <- cbind.data.frame(ModelNo  = "Model 1 - spatial", Model.coefs.sp1c,
                                     Variable = c("Int","Int","sigma^2","phi"),
                                     Submodel = c("Occ","Det","Spatial", "Spatial"),
                                     Type     = "Fixed")
Model.coefs.sp1

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## spatial model 2: occ ~ Depth, det ~ 1 ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
SS.priors2 <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal  = list(mean = c(0, 1), var = c(2.72, 1.50)),
                   sigma.sq.ig  = c(2, 1), 
                   phi.unif     = c(3/max.dist, 3/min.dist))

(start <- Sys.time())
out.sp2 <- spPGOcc(occ.formula = occ.formula.sp2, det.formula = det.formula.sp2, cov.model=cov.model,NNGP = TRUE,
                   data = data.list, n.batch = n.batch, batch.length = batch.length,n.omp.threads = 3,n.neighbors = 5,
                   n.chains = n.chains, n.burn = n.burn, n.report = n.report,
                   n.thin = n.thin, priors = SS.priors2, inits = SS.inits,
                   k.fold = 10, k.fold.threads = 10, verbose = TRUE)
(Sys.time()-start); beep(sound = "coin", expr = NULL) # 1.143345 mins
save(out.sp2, file = paste("Results/2011-2016-2022/Models/outsp2-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out.sp2)
Model.coefs.sp2 <- cbind.data.frame(Intercept.psi = c(0.9269, 0.2821, 0.4466, 0.9034, 1.5447, 1.0085, 1206),
                                    Depth.psi     = c(1.2259, 0.3319, 0.6665, 1.1984, 1.9902, 1.0340, 1157),
                                    Intercept.p   = c(1.3612, 0.2641, 0.8471, 1.3588, 1.8927, 1.0011, 2624),
                                    sigma.sq      = c(1.1441, 1.6387, 0.1883, 0.6489, 5.4124, 1.0491,  526),
                                    phi           = c(0.2506, 0.1509, 0.0086, 0.2513, 0.4975, 1.0080, 1084),
                                    Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs.sp2m <- melt(Model.coefs.sp2,   id.vars=c("Summary"))
Model.coefs.sp2c <- dcast(Model.coefs.sp2m, variable~Summary, value.var="value")
Model.coefs.sp2c <- Model.coefs.sp2c[2:8]
Model.coefs.sp2  <- cbind.data.frame(ModelNo  = "Model 2 - spatial", Model.coefs.sp2c,
                                     Variable = c("Int","Depth","Int","sigma^2","phi"),
                                     Submodel = c("Occ","Occ","Det","Spatial", "Spatial"),
                                     Type     = "Fixed")
Model.coefs.sp2

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## spatial model 3: occ ~ Year,  det ~ 1 ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
SS.priors3 <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal  = list(mean = 0, var = 2.72),
                   sigma.sq.ig  = c(2, 1), 
                   sigma.sq.p.ig = list(shape = 0.1, scale = 0.1),
                   phi.unif     = c(3/max.dist, 3/min.dist))

(start <- Sys.time())
out.sp3 <- spPGOcc(occ.formula = occ.formula.sp3, det.formula = det.formula.sp3, cov.model=cov.model,NNGP = TRUE,
                   data = data.list, n.batch = n.batch, batch.length = batch.length,n.omp.threads = 3,n.neighbors = 5,
                   n.chains = n.chains, n.burn = n.burn, n.report = n.report,
                   n.thin = n.thin, priors = SS.priors3, inits = SS.inits,
                   k.fold = 10, k.fold.threads = 10, verbose = TRUE)
(Sys.time()-start); beep(sound = "coin", expr = NULL)# 1.306716 mins
save(out.sp3, file = paste("Results/2011-2016-2022/Models/outsp3-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out.sp3)
Model.coefs.sp3 <- cbind.data.frame(Intercept.psi = c( 0.0783, 0.4485, -0.8137,  0.0757, 0.9538, 1.0013, 3000),
                                    Year2016.psi  = c(-0.6778, 0.6585, -2.0266, -0.6815, 0.6192, 1.0020, 3000),
                                    Year2022.psi  = c( 1.2355, 0.5447,  0.2029,  1.2146, 2.3832, 1.0057, 1797),
                                    Intercept.p   = c( 1.3829, 0.2641,  0.8810,  1.3771, 1.9047, 1.0000, 3000),
                                    sigma.sq      = c( 1.2515, 1.8523,  0.1960,  0.6974, 5.8522, 1.0095, 453),
                                    phi           = c( 0.2543, 0.1511,  0.0074,  0.2549, 0.4974, 1.0198, 845),
                                    Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs.sp3m <- melt(Model.coefs.sp3,   id.vars=c("Summary"))
Model.coefs.sp3c <- dcast(Model.coefs.sp3m, variable~Summary, value.var="value")
Model.coefs.sp3c <- Model.coefs.sp3c[2:8]
Model.coefs.sp3  <- cbind.data.frame(ModelNo = "Model 3 - spatial", Model.coefs.sp3c,
                                     Variable = c("Int","Year 2016","Year 2022","Int","sigma^2","phi"),
                                     Submodel = c("Occ","Occ","Occ","Det","Spatial", "Spatial"),
                                     Type = "Fixed")
Model.coefs.sp3

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## spatial model 4: occ ~ Depth + Year, det ~ 1 ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
SS.priors4 <- list(alpha.normal = list(mean = 0, var = 2.72), 
                   beta.normal  = list(mean = c(0, 1, 0, 0), var = c(2.72, 1.50, 2.72 , 2.72)),
                   sigma.sq.ig  = c(2, 1), 
                   sigma.sq.p.ig = list(shape = 0.1, scale = 0.1),
                   phi.unif     = c(3/max.dist, 3/min.dist))

(start <- Sys.time())
out.sp4 <- spPGOcc(occ.formula = occ.formula.sp4, det.formula = det.formula.sp4, cov.model=cov.model,NNGP = TRUE,
                   data = data.list, n.batch = n.batch, batch.length = batch.length,n.omp.threads = 3,n.neighbors = 5,
                   n.chains = n.chains, n.burn = n.burn, n.report = n.report,
                   n.thin = n.thin, priors = SS.priors4, inits = SS.inits,
                   k.fold = 10, k.fold.threads = 10, verbose = TRUE)
(Sys.time()-start); beep(sound = "coin", expr = NULL)# 1.200091 mins
save(out.sp4, file = paste("Results/2011-2016-2022/Models/outsp4-2011_2016_2022-", n.chains, "-", Sys.Date(), ".R", sep = ''))

summary(out.sp4)
Model.coefs.sp4 <- cbind.data.frame(Intercept.psi = c( 0.0082, 0.5351, -1.0047, -0.0027, 1.0579, 1.0000, 2735),
                                    Depth         = c( 1.2147, 0.3594,  0.6329,  1.1799, 2.0417, 1.0062,  734),
                                    Year2016.psi  = c(-0.2617, 0.7390, -1.6812, -0.2577, 1.1508, 1.0028, 3000),
                                    Year2022.psi  = c( 1.4230, 0.6247,  0.2834,  1.3905, 2.7388, 1.0008, 1354),
                                    Intercept.p   = c( 1.3790, 0.2620,  0.8622,  1.3781, 1.8931, 1.0019, 3000),
                                    sigma.sq      = c( 1.4485, 2.2339,  0.1938,  0.7411, 7.0508, 1.0265, 333),
                                    phi           = c( 0.2141, 0.1609,  0.0011,  0.2018, 0.4921, 1.0425, 201),
                                    Summary       = c("Mean", "SD", "CI2.5", "CI50", "CI97.5", "Rhat", "ESS"))
Model.coefs.sp4m <- melt(Model.coefs.sp4,   id.vars=c("Summary"))
Model.coefs.sp4c <- dcast(Model.coefs.sp4m, variable~Summary, value.var="value")
Model.coefs.sp4c <- Model.coefs.sp4c[2:8]
Model.coefs.sp4  <- cbind.data.frame(ModelNo = "Model 4 - spatial", Model.coefs.sp4c,
                                     Variable = c("Int","Depth","Year 2016","Year 2022","Int","sigma^2","phi"),
                                     Submodel = c("Occ","Occ","Occ","Occ","Det","Spatial","Spatial"),
                                     Type = "Fixed")
Model.coefs.sp4

################################################################################
# Compare models
################################################################################
n.mods<-4
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

write.csv(mod.comp.sp.fin,"Results/2011-2016-2022/Models/Model.comparisons.spatial.csv")

Mod.comparisons <- rbind.data.frame(mod.comp, mod.comp.sp.fin)
(Mod.comparisons <- Mod.comparisons[order(Mod.comparisons$WAIC),])

write.csv(Mod.comparisons,"Results/2011-2016-2022/Models/Model.comparisons.total.csv")

################################################################################
################################################################################
ModelCoefs.tot2 <- rbind(Model.coefs.sp1, Model.coefs.sp2, Model.coefs.sp3, Model.coefs.sp4)
ModelCoefs.tot2$spatial <- "Spatial"

# export model selection table
ModelCoefs <- rbind(ModelCoefs.1, ModelCoefs.tot2)
write.csv(ModelCoefs, "Results/2011-2016-2022/Models/Model.coefficients.2011-2016-2022.csv")
