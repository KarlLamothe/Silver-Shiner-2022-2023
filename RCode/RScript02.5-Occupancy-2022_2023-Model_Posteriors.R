################################################################################
################################################################################
# Plot posterior probabilities
# RScript02.5-Occupancy-2022_2023-Model-Posteriors.R 
# This is Script 10 of occupancy modelling
# This script is used to plot posterior probabilities for single-species occupancy 
# models built for Silver Shiner using data collected from 16 Mile Creek in 2022
# and 2023.
################################################################################
################################################################################
# load packages
library(pacman)     # downloads and loads packages simultaneously
p_load(ggplot2)     # nice plots
p_load(patchwork)   # combine ggplots
p_load(spOccupancy) # occupancy models
p_load(xlsx)        # read in xlsx files

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
                  legend.background = element_blank()))

# read in detection and covariate data collected in 2022, 2023
SS_data_2022 <- read.xlsx("Data/Silver-Shiner-16Mile-2011-2022.xlsx", header=T, sheetName = "SS_SMC_2022")
SS_data_2022 <- SS_data_2022[order(SS_data_2022$SiteNo),]
SS_data_2023 <- read.xlsx("Data/Silver-Shiner-16Mile-2011-2022.xlsx", header=T, sheetName = "SS_SMC_2023")
SS_data_2023 <- SS_data_2023[order(SS_data_2023$SiteNo),]
SS_data_New  <- rbind(SS_data_2022, SS_data_2023)
SS_data_New$Year <- as.character(SS_data_New$Year)
head(SS_data_New)

# load the models
load("Results/2022-2023/Models/out1-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out2-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out3-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out4-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out5-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out6-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out7-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out.sp1-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out.sp2-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out.sp3-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out.sp4-2022_2023-3-2024-12-03.R")

################################################################################
# Posterior estimates
################################################################################
n.mods <- 7
df <- data.frame(); Occ.Det.results <- data.frame()
  
for (i in 1:n.mods) {
  model_name <- paste("out", i, sep = "")
  model <- get(model_name)
  
  ##########################################################################
  # Occupancy probability
  ##########################################################################
  Psi.results <- cbind.data.frame(Psi  = c(model$psi.samples),
                                  Year = c(rep("2022", 303000),rep("2023", 297000)))

  Psi.res.df <- cbind.data.frame(Model    = model_name, 
                                 Submodel = "Occupancy",
                                 Year     = c("2022","2023"),
                                 Mean     = c(aggregate(Psi.results$Psi, list(Psi.results$Year), mean)$x),
                                 Med      = c(aggregate(Psi.results$Psi, list(Psi.results$Year), median)$x),
                                 SD       = c(aggregate(Psi.results$Psi, list(Psi.results$Year), sd)$x),
                                 Low.025  = c(aggregate(Psi.results$Psi, list(Psi.results$Year), quantile, probs = c(0.025, 0.975))$x[1:2]),
                                 High.975 = c(aggregate(Psi.results$Psi, list(Psi.results$Year), quantile, probs = c(0.025, 0.975))$x[3:4]))

  # create plot
  Psi.results.gg<-ggplot(Psi.results, aes(x=Psi, group=Year, fill=Year))+
    geom_density(alpha=0.5, adjust=1.5) + 
    xlim(0,1) + labs(y = "Density", x = "Occupancy probability") +
  scale_fill_manual(values = c("pink", "darkblue","darkgrey")) +
  theme(legend.position = 'none')

  # Site-level occupancy estimates
  psi.mean.samples <- cbind.data.frame(mean     = apply(model$psi.samples, 2, mean),
                                       low2.5   = apply(model$psi.samples, 2, quantile, 0.025),
                                       high97.5 = apply(model$psi.samples, 2, quantile, 0.975),
                                       Site     = SS_data_New$SiteNo, 
                                       PA = SS.PA)
  psi.mean.1 <- psi.mean.samples[psi.mean.samples$PA == 1,]
  write.csv(psi.mean.1, paste0("Results/2022-2023/Predictions/Site_Pred_M", i, ".csv"))

  # Create site level estimates plot
  site.predsgg<-ggplot(psi.mean.samples, aes(x=Site, y = mean, ymin=low2.5, ymax=high97.5)) +
    geom_linerange(linewidth = 0.5) +
    geom_pointrange(size=0.3) +
    geom_point(aes(y= PA, x=Site), inherit.aes = F, color="red")
  ggsave(site.predsgg, file=paste0("Results/2022-2023/Predictions/Site_Pred_M", i, ".tiff"),
         height=3, width=6, units='in')
  
  ##########################################################################
  # Detection probability
  ##########################################################################
  p.estimates<-cbind.data.frame(p = c(fitted(model)$p.samples[,,1]))
  p.estimates$Year <- c(rep("2022", 303000),rep("2023", 297000))

  P.res.df <- cbind.data.frame(Model    = model_name, 
                               Submodel = "Detection",
                               Year     = c("2022","2023"),
                               Mean     = c(aggregate(p.estimates$p, list(p.estimates$Year), mean)$x),
                               Med      = c(aggregate(p.estimates$p, list(p.estimates$Year), median)$x),
                               SD       = c(aggregate(p.estimates$p, list(p.estimates$Year), sd)$x),
                               Low.025  = c(aggregate(p.estimates$p, list(p.estimates$Year), quantile, probs = c(0.025, 0.975))$x[1:2]),
                               High.975 = c(aggregate(p.estimates$p, list(p.estimates$Year), quantile, probs = c(0.025, 0.975))$x[3:4]))
  
  # detection probability plot
  det.prob.gg <- ggplot(p.estimates, aes(x = p, group = Year, fill = Year)) +
    geom_density(color = "black", adjust = 2, alpha = 0.7, lwd=0.5) +
    scale_fill_manual(values = c("pink", "darkblue","darkgrey")) +
    xlim(0, 1) +
    guides(fill=guide_legend(position='inside'))+
    labs(x = "Detection probability", y = "Density")+
    theme(legend.position.inside = c(0.25, 0.75),
          legend.background = element_blank())
  ggsave((det.prob.gg + Psi.results.gg ), file=paste0("Results/2022-2023/Predictions/Occ_det_M", i, ".tiff"),
         height=2.5, width=6, units='in')
  
  ## Extract detection and occupancy results
  df <- rbind.data.frame(P.res.df, Psi.res.df)
  write.csv(df, paste0("Results/2022-2023/Predictions/Occ_det_M", i, ".csv"))

  Occ.Det.results <- rbind(Occ.Det.results, df)
}

##################################
# Spatial
##################################
n.mods <- 4
df <- data.frame(); Occ.Det.results <- data.frame()
  
for (i in 1:n.mods) {
  model_name <- paste("out.sp", i, sep = "")
  model <- get(model_name)
  Psi.results <- cbind.data.frame(Psi  = c(model$psi.samples),
                                  Year = c(rep("2022", 303000),rep("2023", 297000)))
  Psi.res.df <- cbind.data.frame(Model    = model_name, 
                                 Submodel = "Occupancy",
                                 Year     = c("2022","2023"),
                                 Mean     = c(aggregate(Psi.results$Psi, list(Psi.results$Year), mean)$x),
                                 Med      = c(aggregate(Psi.results$Psi, list(Psi.results$Year), median)$x),
                                 SD       = c(aggregate(Psi.results$Psi, list(Psi.results$Year), sd)$x),
                                 Low.025  = c(aggregate(Psi.results$Psi, list(Psi.results$Year), quantile, probs = c(0.025, 0.975))$x[1:2]),
                                 High.975 = c(aggregate(Psi.results$Psi, list(Psi.results$Year), quantile, probs = c(0.025, 0.975))$x[3:4]))

  # create plot
  Psi.results.gg<-ggplot(Psi.results, aes(x=Psi, group=Year, fill=Year))+
    geom_density(alpha=0.5, adjust=1.5) + 
    xlim(0,1) + labs(y = "Density", x = "Occupancy probability") +
    scale_fill_manual(values = c("pink", "darkblue","darkgrey")) +
    theme(legend.position = 'none')
  
  # Site-level occupancy estimates
  psi.mean.samples <- cbind.data.frame(mean     = apply(model$psi.samples, 2, mean),
                                       low2.5   = apply(model$psi.samples, 2, quantile, 0.025),
                                       high97.5 = apply(model$psi.samples, 2, quantile, 0.975),
                                       Site     = SS_data_New$SiteNo, 
                                       PA = SS.PA)
  psi.mean.1 <- psi.mean.samples[psi.mean.samples$PA == 1,]
  write.csv(psi.mean.1, paste0("Results/2022-2023/Predictions/Site_Pred_M", i, "sp.csv"))

  # Create site level plot
  site.predsgg<-ggplot(psi.mean.samples, aes(x=Site, y = mean, ymin=low2.5, ymax=high97.5)) +
    geom_linerange(linewidth = 0.5) +
    geom_pointrange(size=0.3) +
    geom_point(aes(y= PA, x=Site), inherit.aes = F, color="red")
  ggsave(site.predsgg, file=paste0("Results/2022-2023/Predictions/Site_Pred_M", i, "sp.tiff"),
         height=3, width=6, units='in')

  # Detection probability
  p.estimates<-cbind.data.frame(p = c(fitted(model)$p.samples[,,1]))
  p.estimates$Year <- c(rep("2022", 303000),rep("2023", 297000))
  
  P.res.df <- cbind.data.frame(Model    = model_name, 
                               Submodel = "Detection",
                               Year     = c("2022","2023"),
                               Mean     = c(aggregate(p.estimates$p, list(p.estimates$Year), mean)$x),
                               Med      = c(aggregate(p.estimates$p, list(p.estimates$Year), median)$x),
                               SD       = c(aggregate(p.estimates$p, list(p.estimates$Year), sd)$x),
                               Low.025  = c(aggregate(p.estimates$p, list(p.estimates$Year), quantile, probs = c(0.025, 0.975))$x[1:2]),
                               High.975 = c(aggregate(p.estimates$p, list(p.estimates$Year), quantile, probs = c(0.025, 0.975))$x[3:4]))
  
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
  ggsave((det.prob.gg + Psi.results.gg ), file=paste0("Results/2022-2023/Predictions/Occ_det_M", i, "sp.tiff"),
         height=2.5, width=6, units='in')
  
  ## Extract detection and occupancy results
  df <- rbind.data.frame(P.res.df, Psi.res.df)
  write.csv(df, paste0("Results/2022-2023/Predictions/Occ_det_M", i, "sp.csv"))
  Occ.Det.results <- rbind(Occ.Det.results, df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# create data frame with annual estimates of detection and occupancy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
Occ.det <- rbind.data.frame(
  read.csv("Results/2022-2023/Predictions/Occ_det_M1.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M2.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M3.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M4.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M5.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M6.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M7.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M1sp.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M2sp.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M3sp.csv", header=T),
  read.csv("Results/2022-2023/Predictions/Occ_det_M4sp.csv", header=T)
)
Occ.det$Year <- as.character(Occ.det$Year)
Occ.det$spatial <- c(rep("Non-spatial", 28),rep("Spatial",16))

# remove intercept models
Occ.det2 <- Occ.det[!c(Occ.det$spatial=="Non-spatial"),]
Occ.det2$Submodel <- factor(Occ.det2$Submodel,
                            levels = c("Detection","Occupancy"),
                            labels = c('Detection probability', "Occupancy probability"))
Occ.det2$Year <- factor(Occ.det2$Year,
                            levels = c("2022","2023"),
                            labels = c('Year 2022', "Year 2023"))
Occ.det2 <- Occ.det2[Occ.det2$Model == "out.sp4" |
                       Occ.det2$Model == "out.sp3" |
                       Occ.det2$Model == "out.sp2",]

# plot
Occ.Pred<-ggplot(Occ.det2, aes(y=Mean, x=Year))+
  geom_pointrange(aes(ymin = Low.025, ymax = High.975, group=Model, color=Model), size=0.25,
                  position = position_dodge2(width=0.5, preserve="single"),
                  show.legend = TRUE, linewidth=0.5)+
  geom_point(aes(color=Model), position = position_dodge2(width=0.5, preserve="single"), size=0.25)+
  scale_color_manual(values= c("#000000", "pink", "#009E73"),
                     labels= c(expression(italic('p')~'~Year;'~italic("ψ")~'~Depth'),
                               expression(italic('p')~'~Year;'~italic("ψ")~'~Year'), 
                               expression(italic('p')~'~Year;'~italic("ψ")~'~Depth + Year')))+  
  facet_wrap(~Submodel)+
  labs(x = "Posterior probability", y="Year") +
  theme(legend.title      = element_blank(),
        legend.position   = 'top',
        axis.title.y = element_blank(),
        legend.key        = element_blank(),
        legend.background = element_blank(),
        legend.margin     = margin(0, unit='cm'),
        legend.text = element_text(hjust = 0)) 

tiff("Results/2022-2023/Models/Occ.det.gg.tiff",height=2.5, width=7, units='in',res=800)
Occ.Pred
dev.off()

################################################################################
################################################################################
################################################################################
################################################################################
load("Results/2022-2023/Models/out.sp2-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out.sp3-2022_2023-3-2024-12-03.R")
load("Results/2022-2023/Models/out.sp4-2022_2023-3-2024-12-03.R")

Psi.posts <- cbind.data.frame(Psi  = c(out.sp2$psi.samples,
                                       out.sp3$psi.samples,
                                       out.sp4$psi.samples),
                              Year = rep(c(rep("2022", 303000),
                                           rep("2023", 297000)), 3),
                              Model = rep(c("p ~ Year; ψ ~ Depth", 
                                            "p ~ Year; ψ ~ Year", 
                                            "p ~ Year; ψ ~ Depth + Year"), each=600000))
Psi.posts$Model <- factor(Psi.posts$Model,
                          levels = c("p ~ Year; ψ ~ Depth", 
                                     "p ~ Year; ψ ~ Year", 
                                     "p ~ Year; ψ ~ Depth + Year"))

post.distribution.gg<-ggplot(Psi.posts, aes(x=Psi, group=Year, fill=Year))+
  geom_density(alpha = 0.5, adjust=1)+
  guides(fill=guide_legend(position='inside'))+
  facet_wrap(~Model)+
  labs(x = "Occupancy probability")+
  scale_fill_manual(values=c("grey","blue"))+
  theme(legend.position.inside = c(0.08, 0.75),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

tiff("Results/2022-2023/Predictions/Post.distribution.occ.tiff",
     height=2.5, width=7.5, units='in',res=800)
post.distribution.gg
dev.off()

aggregate(Psi.posts$Psi, list(Psi.posts$Year, Psi.posts$Model), mean)
aggregate(Psi.posts$Psi, list(Psi.posts$Year, Psi.posts$Model), quantile, probs=c(0.025, 0.975))

# Model 5 
# 2022: 0.75 (0.28-0.99)
# 2023: 0.73 (0.26-0.99)
# 
# Model 6
# 2022: 0.74 (0.31-0.98)
# 2023: 0.59 (0.13-0.97)
# 
# Model 7 
# 2022: 0.75 (0.26-0.99)
# 2023: 0.77 (0.21-1.00)
