################################################################################
################################################################################
# Investigating model fit for 2011, 2016, and 2022 models
################################################################################
################################################################################
# load packages
library(pacman)     # downloads and loads packages simultaneously
p_load(ggplot2)     # nice plots
p_load(spOccupancy) # occupancy models
p_load(patchwork)   # combine ggplots
p_load(coda)        # traceplots
p_load(reshape2)    # reshape data frames
p_load(MCMCvis)     # quick plot of posteriors
p_load(ggh4x)       # easily resize ggplot facets
p_load(sf)          # coordinates to projections
p_load(oce)         # decimal degress to utm

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

# include year as a detection covariate in models
detcovs <- cbind.data.frame(Year=Covariates.std$Year, 
                            Depth=Covariates.std$Depth)

#####################################################################
# Load results
#####################################################################
load("Results/2011-2016-2022/Models/out1-2011_2016_2022-3-2024-12-03.R")
load("Results/2011-2016-2022/Models/out2-2011_2016_2022-3-2024-12-03.R")
load("Results/2011-2016-2022/Models/out3-2011_2016_2022-3-2024-12-03.R")
load("Results/2011-2016-2022/Models/out4-2011_2016_2022-3-2024-12-03.R")
load("Results/2011-2016-2022/Models/out5-2011_2016_2022-3-2024-12-03.R")
load("Results/2011-2016-2022/Models/out6-2011_2016_2022-3-2024-12-03.R")

load("Results/2011-2016-2022/Models/outsp1-2011_2016_2022-3-2024-12-03.R")
load("Results/2011-2016-2022/Models/outsp2-2011_2016_2022-3-2024-12-03.R")
load("Results/2011-2016-2022/Models/outsp3-2011_2016_2022-3-2024-12-03.R")
load("Results/2011-2016-2022/Models/outsp4-2011_2016_2022-3-2024-12-03.R")

###################### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ######################
###################### Evaluate convergence and model fit ######################
###################### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ######################

# tukey freeman and chi-squared
# posterior predictive check when grouping data both across sites (group = 1) 
# as well as across replicates (group = 2), as they may reveal (or fail to 
# reveal).Binning the data across sites (group = 1) may help reveal 
# whether the model fails to adequately represent variation in occurrence and 
# detection probability across space, while binning the data across replicates 
# (group = 2) may help reveal whether the model fails to adequately represent 
# variation in detection probability across the different replicate surveys. 

###################### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ######################
n.mods <- 6

# Nonspatial models first
for (i in 1:n.mods) {
  model_name <- paste("out", i, sep = "")

  model <- get(model_name)
  ppc.out1 <- ppcOcc(model, fit.stat = 'freeman-tukey', group = 1)
  ppc.out2 <- ppcOcc(model, fit.stat = 'freeman-tukey', group = 2)
  ppc.out3 <- ppcOcc(model, fit.stat = 'chi-squared',   group = 1)
  ppc.out4 <- ppcOcc(model, fit.stat = 'chi-squared',   group = 2)
  
  print(summary(ppc.out1))
  print(summary(ppc.out2))
  print(summary(ppc.out3))
  print(summary(ppc.out4))
  
  # plotting dataframe base r
  ppc.df2 <- data.frame(fit     = c(ppc.out1$fit.y,      ppc.out2$fit.y, 
                                    ppc.out3$fit.y,      ppc.out4$fit.y), 
                        fit.rep = c(ppc.out1$fit.y.rep,  ppc.out2$fit.y.rep, 
                                    ppc.out3$fit.y.rep,  ppc.out4$fit.y.rep), 
                        group   = c(rep("Freeman-Tukey 1", length(ppc.out1$fit.y)),
                                    rep("Freeman-Tukey 2", length(ppc.out2$fit.y)),
                                    rep("Chi-square 1",    length(ppc.out3$fit.y)),
                                    rep("Chi-square 2",    length(ppc.out4$fit.y))),
                        color   = "black")
  ppc.df2$color[ppc.df2$fit.rep > ppc.df2$fit] <- 'grey'
  
  #plot 
  temp.plot<-ggplot(ppc.df2, aes(y=fit, x=fit.rep, group=group))+
    geom_point(aes(color=color))+
    geom_line(aes(y=fit, x=fit), color="black")+
    scale_color_manual(values=c("grey","black"))+
    facet_wrap(~group, scales="free")+
    labs(y="Fit statistic of the model generated data",
         x="Fit statistic of the true data")+
    theme(legend.position='none') 
  ggsave(temp.plot, file=paste0("Results/2011-2016-2022/Model-Fit/Model-Fit-M", i, ".tiff"),
         height=6, width=7.5, units='in')
  
  # check out what sites may be limiting the fit
  diff.fit2 <- cbind.data.frame(fit   = c(c(ppc.out1$fit.y.rep.group.quants[3,] - ppc.out1$fit.y.group.quants[3,]),
                                          c(ppc.out3$fit.y.rep.group.quants[3,] - ppc.out3$fit.y.group.quants[3,])),
                                Site  = rep(SS_data$Field.Number, 2),
                                group = c(rep("Freeman-Tukey 1", length(ppc.out1$fit.y.rep.group.quants[3,])),
                                          rep("Chi-square 1",    length(ppc.out3$fit.y.rep.group.quants[3,]))),
                                Year  = rep(c(rep("2011", 24), 
                                              rep("2016", 22),
                                              rep("2022", 101)), 2))
  
  temp.plot2<-ggplot(diff.fit2, aes(y=fit, x=Site, group=group, color=Year))+
    geom_point()+
    facet_wrap(~group)+
    scale_color_manual(values = c("pink", "darkblue","darkgrey")) +
    labs(y="Replicate - true discrepancy", x = "Site")+
    theme(legend.position = 'none',
          axis.text.x=element_blank())+ 
    facet_grid(group~Year, scales="free_x", space = "free_x")
  ggsave(temp.plot2, file=paste0("Results/2011-2016-2022/Model-Fit/Outliers_M", i, ".tiff"),
         height=6, width=7.5, units='in')
}

# I manually entered the significance values into this data frame
Significance <- cbind.data.frame(Model = c(rep("Non-spatial model 1",4), rep("Non-spatial model 2",4), 
                                           rep("Non-spatial model 3",4), rep("Non-spatial model 4",4), 
                                           rep("Non-spatial model 5",4), rep("Non-spatial model 6",4)),
                                 Appr  = rep(c("FT-1", "FT-2", "CS-1", "CS-2"), 12),
                                 Sign. = c(0.9930, 0.5340, 0.9510, 0.5727,
                                           0.9940, 0.5687, 0.9667, 0.5683,
                                           0.9230, 0.3397, 0.5160, 0.3560,
                                           0.9903, 0.5317, 0.9547, 0.5140,
                                           0.9927, 0.5487, 0.9530, 0.5503,
                                           0.9900, 0.5503, 0.9593, 0.5493))

################################################################################
# spatial
################################################################################
n.mods<-4
for (i in 1:n.mods) {
  model_name <- paste("out.sp", i, sep = "")
  model <- get(model_name)
  ppc.out1 <- ppcOcc(model, fit.stat = 'freeman-tukey', group = 1)
  ppc.out2 <- ppcOcc(model, fit.stat = 'freeman-tukey', group = 2)
  ppc.out3 <- ppcOcc(model, fit.stat = 'chi-squared',   group = 1)
  ppc.out4 <- ppcOcc(model, fit.stat = 'chi-squared',   group = 2)
  
  print(summary(ppc.out1))
  print(summary(ppc.out2))
  print(summary(ppc.out3))
  print(summary(ppc.out4))
  
  # plotting dataframe base r
  ppc.df2 <- data.frame(fit     = c(ppc.out1$fit.y,      ppc.out2$fit.y, 
                                    ppc.out3$fit.y,      ppc.out4$fit.y), 
                        fit.rep = c(ppc.out1$fit.y.rep,  ppc.out2$fit.y.rep, 
                                    ppc.out3$fit.y.rep,  ppc.out4$fit.y.rep), 
                        group   = c(rep("Freeman-Tukey 1", length(ppc.out1$fit.y)),
                                    rep("Freeman-Tukey 2", length(ppc.out2$fit.y)),
                                    rep("Chi-square 1",    length(ppc.out3$fit.y)),
                                    rep("Chi-square 2",    length(ppc.out4$fit.y))),
                        color   = "black")
  ppc.df2$color[ppc.df2$fit.rep > ppc.df2$fit] <- 'grey'
  
  #plot 
  temp.plot<-ggplot(ppc.df2, aes(y=fit, x=fit.rep, group=group))+
    geom_point(aes(color=color))+
    geom_line(aes(y=fit, x=fit), color="black")+
    scale_color_manual(values=c("grey","black"))+
    facet_wrap(~group, scales="free")+
    labs(y="Fit statistic of the model generated data",
         x="Fit statistic of the true data")+
    theme(legend.position='none') 
  ggsave(temp.plot, file=paste0("Results/2011-2016-2022/Model-Fit/Model-Fit-spM", i, ".tiff"),
         height=6, width=7.5, units='in')
  
  # check out what sites may be limiting the fit
  diff.fit2 <- cbind.data.frame(fit   = c(c(ppc.out1$fit.y.rep.group.quants[3,] - ppc.out1$fit.y.group.quants[3,]),
                                          c(ppc.out3$fit.y.rep.group.quants[3,] - ppc.out3$fit.y.group.quants[3,])),
                                Site  = rep(SS_data$Field.Number, 2),
                                group = c(rep("Freeman-Tukey 1", length(ppc.out1$fit.y.rep.group.quants[3,])),
                                          rep("Chi-square 1",    length(ppc.out3$fit.y.rep.group.quants[3,]))),
                                Year  = rep(c(rep("2011", 24), 
                                              rep("2016", 22),
                                              rep("2022", 101)), 2))
  
  temp.plot2<-ggplot(diff.fit2, aes(y=fit, x=Site, group=group, color=Year))+
    geom_point()+
    facet_wrap(~group)+
    scale_color_manual(values = c("pink", "darkblue","darkgrey")) +
    labs(y="Replicate - true discrepancy", x = "Site")+
    theme(legend.position = 'none',
          axis.text.x=element_blank())+ 
    facet_grid(group~Year, scales="free_x", space = "free_x")
  ggsave(temp.plot2, file=paste0("Results/2011-2016-2022/Model-Fit/Outliers_spM", i, ".tiff"),
         height=6, width=7.5, units='in')
}

# I manually entered the significance values into this data frame
Significance2 <- cbind.data.frame(Model = c(rep("Spatial model 1",4), rep("Spatial model 2",4), 
                                           rep("Spatial model 3",4), rep("Spatial model 4",4)),
                                 Appr  = rep(c("FT-1", "FT-2", "CS-1", "CS-2"), 12),
                                 Sign. = c(0.9950, 0.5683, 0.9560, 0.5607,
                                           0.9937, 0.5243, 0.9550, 0.5503,
                                           0.9923, 0.5467, 0.9507, 0.5610,
                                           0.9940, 0.5577, 0.9580, 0.5820))
Significance.Total <- rbind.data.frame(Significance,Significance2)
write.csv(Significance.Total, "Results/2011-2016-2022/Model-Fit/Fit.statistics.csv")

##########################################################
########### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###########
########### Evaluate convergence and model fit ###########
########### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###########
##########################################################
# Rhat
rhat.df<-cbind.data.frame(Rhat = c(out1$rhat$beta, out2$rhat$beta, out3$rhat$beta,
                                   out4$rhat$beta, out5$rhat$beta, out6$rhat$beta, 
                                   out.sp1$rhat$beta, out.sp2$rhat$beta, out.sp3$rhat$beta, 
                                   out.sp4$rhat$beta),
                          Model = c(rep("Model 1",  length(out1$rhat$beta)),  
                                    rep("Model 2",  length(out2$rhat$beta)),
                                    rep("Model 3",  length(out3$rhat$beta)),  
                                    rep("Model 4",  length(out4$rhat$beta)),
                                    rep("Model 5",  length(out5$rhat$beta)),
                                    rep("Model 6",  length(out6$rhat$beta)),
                                    rep("Spatial 1",  length(out.sp1$rhat$beta)),
                                    rep("Spatial 2",  length(out.sp2$rhat$beta)),
                                    rep("Spatial 3",  length(out.sp3$rhat$beta)),
                                    rep("Spatial 4",  length(out.sp4$rhat$beta))))
                                    
tiff("Results/2011-2016-2022/Model-Fit/Rhat-2011-2016-2022.tiff", height=3, width=4, units='in', res=800)
ggplot(rhat.df, aes(x=Model, y=Rhat))+
  geom_point(size=1)+
  labs(y="Rhat")+
  coord_flip()+
  theme(axis.title.y=element_blank())
dev.off()

write.csv(rhat.df, "Results/2011-2016-2022/Model-Fit/Rhat-2011-2016-2022.csv")

min(rhat.df$Rhat)
max(rhat.df$Rhat)
