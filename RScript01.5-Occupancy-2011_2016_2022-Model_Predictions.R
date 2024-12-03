################################################################################
################################################################################
# Check out model predictions
################################################################################
################################################################################
# load packages
library(pacman)     # downloads and loads packages simultaneously
p_load(ggplot2)     # nice plots
p_load(patchwork)   # combine ggplots

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

# load models
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

#################################################################################
## Posterior estimates
#################################################################################
n.mods <- 6
df <- data.frame(); Occ.Det.results <- data.frame()
  
for (i in 1:n.mods) {
  model_name <- paste("out", i, sep = "")
  model <- get(model_name)

  ##########################################################################
  #### Occupancy probability
  ##########################################################################
  Psi.results <- cbind.data.frame(Psi  = c(model$psi.samples),
                                  Year = c(rep("2011",72000), rep("2016",66000), rep("2022",303000)))

  Psi.res.df <- cbind.data.frame(Model    = model_name, 
                                 Submodel = "Occupancy",
                                 Year     = c("2011", "2016" , "2022"),
                                 Mean     = c(aggregate(Psi.results$Psi, list(Psi.results$Year), mean)$x),
                                 Med      = c(aggregate(Psi.results$Psi, list(Psi.results$Year), median)$x),
                                 SD       = c(aggregate(Psi.results$Psi, list(Psi.results$Year), sd)$x),
                                 Low.025  = c(aggregate(Psi.results$Psi, list(Psi.results$Year), quantile, probs = c(0.025, 0.975))$x[1:3]),
                                 High.975 = c(aggregate(Psi.results$Psi, list(Psi.results$Year), quantile, probs = c(0.025, 0.975))$x[4:6]))

  # create plot of annual occupancy probability
  Psi.results.gg<-ggplot(Psi.results, aes(x=Psi, group=Year, fill=Year))+
    geom_density(alpha=0.5, adjust=1.5) + 
    xlim(0,1) + labs(y = "Density", x = "Occupancy probability") +
    scale_fill_manual(values = c("pink", "darkblue","darkgrey")) +
    theme(legend.position = 'none')
  
  #### Site-level occupancy estimates
  psi.mean.samples <- cbind.data.frame(mean     = apply(model$psi.samples, 2, mean),
                                       low2.5   = apply(model$psi.samples, 2, quantile, 0.025),
                                       high97.5 = apply(model$psi.samples, 2, quantile, 0.975),
                                       Site     = seq(1:147), PA = SS.PA$SS.PA)
  psi.mean.1 <- psi.mean.samples[psi.mean.samples$PA == 1,]
  
  ### Export site level predictions
  write.csv(psi.mean.1, paste0("Results/2011-2016-2022/Model-Predictions/Site_Pred_M", i, "-2011-2016-2022.csv"))

  # Create site level predictions plot
  site.predsgg<-ggplot(psi.mean.samples, aes(x=Site, y = mean, ymin=low2.5, ymax=high97.5)) +
    geom_linerange(linewidth = 0.5) +
    geom_pointrange(size=0.3) +
    geom_point(aes(y= PA, x=Site), inherit.aes = F, color="red")

  # save predictions plot
  ggsave(site.predsgg, file=paste0("Results/2011-2016-2022/Model-Predictions/Site_Pred_M", i, ".tiff"),
         height=3, width=6, units='in')

  ##########################################################################
  # Detection probability
  ##########################################################################
  # extract fitted detection probability estimates 
  p.estimates<-cbind.data.frame(p = c(fitted(model)$p.samples[,,1]))
  
  # include year and plot
  p.estimates$Year <- c(rep("2011",72000),rep("2016",66000),rep("2022",303000))

   P.res.df <- cbind.data.frame(Model    = model_name, 
                               Submodel = "Detection",
                               Year     = c("2011", "2016" , "2022"),
                               Mean     = c(aggregate(p.estimates$p, list(p.estimates$Year), mean)$x),
                               Med      = c(aggregate(p.estimates$p, list(p.estimates$Year), median)$x),
                               SD       = c(aggregate(p.estimates$p, list(p.estimates$Year), sd)$x),
                               Low.025  = c(aggregate(p.estimates$p, list(p.estimates$Year), quantile, probs = c(0.025, 0.975))$x[1:3]),
                               High.975 = c(aggregate(p.estimates$p, list(p.estimates$Year), quantile, probs = c(0.025, 0.975))$x[4:6]))
  
  # detection probability plot
  det.prob.gg <- ggplot(p.estimates, aes(x = p, group = Year, fill = Year)) +
    geom_density(color = "black", adjust = 2, alpha = 0.7, lwd=0.5) +
    scale_fill_manual(values = c("pink", "darkblue","darkgrey")) +
    xlim(0, 1) +
    guides(fill=guide_legend(position='inside'))+
    labs(x = "Detection probability", y = "Density")+
    theme(legend.position.inside = c(0.25, 0.75),
          legend.background = element_blank())
  det.prob.gg

  # save the detection and occupancy probability plots together
  ggsave((det.prob.gg + Psi.results.gg), file=paste0("Results/2011-2016-2022/Model-Predictions/Occ_det_M", i, ".tiff"),
         height=2.5, width=6, units='in')
  
  ## Extract detection and occupancy results
  df <- rbind.data.frame(P.res.df,Psi.res.df)
  write.csv(df, paste0("Results/2011-2016-2022/Model-Predictions/Occ_det_M", i, ".csv"))

  Occ.Det.results <- rbind(Occ.Det.results, df)
}

#################################################################################
#################################################################################
# Spatial
#################################################################################
#################################################################################
n.mods <- 4
df <- data.frame(); Occ.Det.results <- data.frame()
  
for (i in 1:n.mods) {
  model_name <- paste("out.sp", i, sep = "")
  model <- get(model_name)

  ##########################################################################
  #### Occupancy probability
  ##########################################################################
  Psi.results <- cbind.data.frame(Psi  = c(model$psi.samples),
                                  Year = c(rep("2011",72000), rep("2016",66000), rep("2022",303000)))

  Psi.res.df <- cbind.data.frame(Model    = model_name, 
                                 Submodel = "Occupancy",
                                 Year     = c("2011", "2016" , "2022"),
                                 Mean     = c(aggregate(Psi.results$Psi, list(Psi.results$Year), mean)$x),
                                 Med      = c(aggregate(Psi.results$Psi, list(Psi.results$Year), median)$x),
                                 SD       = c(aggregate(Psi.results$Psi, list(Psi.results$Year), sd)$x),
                                 Low.025  = c(aggregate(Psi.results$Psi, list(Psi.results$Year), quantile, probs = c(0.025, 0.975))$x[1:3]),
                                 High.975 = c(aggregate(Psi.results$Psi, list(Psi.results$Year), quantile, probs = c(0.025, 0.975))$x[4:6]))

  # create plot of annual occupancy probability
  Psi.results.gg<-ggplot(Psi.results, aes(x=Psi, group=Year, fill=Year))+
    geom_density(alpha=0.5, adjust=1.5) + 
    xlim(0,1) + labs(y = "Density", x = "Occupancy probability") +
    scale_fill_manual(values = c("pink", "darkblue","darkgrey")) +
    theme(legend.position = 'none')
  
  #### Site-level occupancy estimates
  psi.mean.samples <- cbind.data.frame(mean     = apply(model$psi.samples, 2, mean),
                                       low2.5   = apply(model$psi.samples, 2, quantile, 0.025),
                                       high97.5 = apply(model$psi.samples, 2, quantile, 0.975),
                                       Site     = seq(1:147), PA = SS.PA$SS.PA)
  psi.mean.1 <- psi.mean.samples[psi.mean.samples$PA == 1,]
  
  ### Export site level predictions
  write.csv(psi.mean.1, paste0("Results/2011-2016-2022/Model-Predictions/Site_Pred_M", i, "sp-2011-2016-2022.csv"))
  
  # Create site level predictions plot
  site.predsgg<-ggplot(psi.mean.samples, aes(x=Site, y = mean, ymin=low2.5, ymax=high97.5)) +
    geom_linerange(linewidth = 0.5) +
    geom_pointrange(size=0.3) +
    geom_point(aes(y= PA, x=Site), inherit.aes = F, color="red")

  # save predictions plot
  ggsave(site.predsgg, file=paste0("Results/2011-2016-2022/Model-Predictions/Site_Pred_M", i, "sp.tiff"),
         height=3, width=6, units='in')
  
  ##########################################################################
  # Detection probability
  ##########################################################################
  # extract fitted detection probability estimates 
  p.estimates<-cbind.data.frame(p = c(fitted(model)$p.samples[,,1]))
  
  # include year and plot
  p.estimates$Year <- c(rep("2011",72000),rep("2016",66000),rep("2022",303000))

   P.res.df <- cbind.data.frame(Model    = model_name, 
                               Submodel = "Detection",
                               Year     = c("2011", "2016" , "2022"),
                               Mean     = c(aggregate(p.estimates$p, list(p.estimates$Year), mean)$x),
                               Med      = c(aggregate(p.estimates$p, list(p.estimates$Year), median)$x),
                               SD       = c(aggregate(p.estimates$p, list(p.estimates$Year), sd)$x),
                               Low.025  = c(aggregate(p.estimates$p, list(p.estimates$Year), quantile, probs = c(0.025, 0.975))$x[1:3]),
                               High.975 = c(aggregate(p.estimates$p, list(p.estimates$Year), quantile, probs = c(0.025, 0.975))$x[4:6]))
  
  # detection probability plot
  det.prob.gg <- ggplot(p.estimates, aes(x = p, group = Year, fill = Year)) +
    geom_density(color = "black", adjust = 2, alpha = 0.7, lwd=0.5) +
    scale_fill_manual(values = c("pink", "darkblue","darkgrey")) +
    xlim(0, 1) +
    guides(fill=guide_legend(position='inside'))+
    labs(x = "Detection probability", y = "Density")+
    theme(legend.position.inside = c(0.25, 0.75),
          legend.background = element_blank())
  det.prob.gg

  # save the detection and occupancy probability plots together
  ggsave((det.prob.gg + Psi.results.gg), file=paste0("Results/2011-2016-2022/Model-Predictions/Occ_det_M", i, "sp.tiff"),
         height=2.5, width=6, units='in')
  
  ## Extract detection and occupancy results
  df <- rbind.data.frame(P.res.df,Psi.res.df)
  write.csv(df, paste0("Results/2011-2016-2022/Model-Predictions/Occ_det_M", i, "sp.csv"))
  
  Occ.Det.results <- rbind(Occ.Det.results, df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# create data frame with annual estimates of detection and occupancy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

Occ.det <- rbind.data.frame(
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M1.csv", header=T),
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M2.csv", header=T),
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M3.csv", header=T),
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M4.csv", header=T),
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M5.csv", header=T),
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M6.csv", header=T),
  
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M1sp.csv", header=T),
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M2sp.csv", header=T),
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M3sp.csv", header=T),
  read.csv("Results/2011-2016-2022/Model-Predictions/Occ_det_M4sp.csv", header=T))
Occ.det

Occ.det$Year <- as.character(Occ.det$Year)
Occ.det$spatial <- c(rep("Non-spatial", 36), rep("Spatial",24))

# revise variable text
Occ.det$Submodel <- factor(Occ.det$Submodel,
                           levels = c("Detection", "Occupancy"),
                           labels = c('Detection Pr.', "Occupancy Pr."))

Occ.det$Year <- factor(Occ.det$Year,
                           levels = c("2011", "2016", "2022"),
                           labels = c('Year 2011', "Year 2016", "Year 2022"))
Occ.det1 <- Occ.det[Occ.det$Model=='out6' |
                      Occ.det$Model=='out.sp4',]

Occ.vars <- Occ.det1[Occ.det1$Submodel=="Occupancy Pr.",]

# plot
dat_text <- data.frame(label = "b)",Submodel = 'Occupancy Pr.', Year = 0.97, Mean = 0.6)
Occ.gg<-ggplot(Occ.vars, aes(y=Mean, x=Year))+
  geom_pointrange(aes(ymin = Low.025, ymax = High.975, color=spatial), size=0.25,
                  position = position_dodge2(width=0.5, preserve="single"),
                  show.legend = TRUE, linewidth=0.5)+
  geom_text(data = dat_text, mapping = aes(y = Year, x = Mean, label = label), size=3)+
  scale_color_manual(values= c("#000000", "darkgrey"),
                     labels= c("Non-spatial", "Spatial"))+
  labs(y = "Occupancy probability", y="Year") +
  theme(legend.position   = 'none',
        axis.title.x = element_blank()) 
Occ.gg

aggregate(Occ.det$Mean, list(Occ.det$Submodel), min)
aggregate(Occ.det$Mean, list(Occ.det$Submodel), max)
aggregate(Occ.det$Mean, list(Occ.det$Submodel), mean)

aggregate(Occ.det$Mean, list(Occ.det$Submodel, Occ.det$Year), min)
aggregate(Occ.det$Mean, list(Occ.det$Submodel, Occ.det$Year), max)

tiff("Results/2011-2016-2022/Model-Predictions/Posterior_Occupancy.tiff",
     height=2.25, width=7, units='in',res=800)
(Occ.var.gg+ Occ.gg)
dev.off()
