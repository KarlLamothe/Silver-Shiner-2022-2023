################################################################################
################################################################################
# Investigating model fit for 2022 and 2023 models
# RScript02.2-Occupancy-2022_2023-Model-Fit.R 
# This is Script 7 of occupancy modelling
# This script is used to evaluate the fit and convergence of single-species 
# occupancy models for Silver Shiner using data collected from 16 Mile Creek in 
# 2022 and 2023.
################################################################################
################################################################################
library(pacman)     # downloads and loads packages simultaneously
p_load(ggplot2)     # nice plots
p_load(spOccupancy) # occupancy models

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

#####################################################################
# Load results
#####################################################################
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

###################### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ######################
###################### Evaluate convergence and model fit ######################
###################### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ######################
# Rhat
rhat.df<-cbind.data.frame(Rhat = c(out1$rhat$beta, out2$rhat$beta, out3$rhat$beta, out4$rhat$beta,
                                   out5$rhat$beta, out6$rhat$beta, out7$rhat$beta,
                                   out.sp1$rhat$beta, out.sp2$rhat$beta, out.sp3$rhat$beta, out.sp4$rhat$beta),
                          Model = c(rep("Non-spatial 1",  length(out1$rhat$beta)),  
                                    rep("Non-spatial 2",  length(out2$rhat$beta)),
                                    rep("Non-spatial 3",  length(out3$rhat$beta)),  
                                    rep("Non-spatial 4",  length(out4$rhat$beta)),
                                    rep("Non-spatial 5",  length(out5$rhat$beta)),
                                    rep("Non-spatial 6",  length(out6$rhat$beta)),  
                                    rep("Non-spatial 7",  length(out7$rhat$beta)),
                                    rep("Spatial 1",  length(out.sp1$rhat$beta)),  
                                    rep("Spatial 2",  length(out.sp2$rhat$beta)),
                                    rep("Spatial 3",  length(out.sp3$rhat$beta)),  
                                    rep("Spatial 4",  length(out.sp4$rhat$beta))))
min(rhat.df$Rhat)
max(rhat.df$Rhat)

tiff("Results/2022-2023/Model-Fit/Rhat_2022_2023.tiff", 
     height=3, width=5, units='in', res=800)
ggplot(rhat.df, aes(x=Model, y=Rhat))+
  geom_point(size=1)+
  labs(y="Rhat")+
  coord_flip()+
  theme(axis.title.y=element_blank())
dev.off()

write.csv(rhat.df, "Results//2022-2023/Model-Fit/Rhat_2022_2023.csv")
str(rhat.df)

################################################################################
# Evaluate model fit
# Standard, non-spatial models first
################################################################################
n.mods<-7
for (i in 1:n.mods) {
  # standard model
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
  ggsave(temp.plot, file=paste0("Results/2022-2023/Model-Fit/Model_fit_M", i, ".tiff"),
         height=6, width=7.5, units='in')
  
  # check out what sites may be limiting the fit
  diff.fit2 <- cbind.data.frame(fit   = c(c(ppc.out1$fit.y.rep.group.quants[3,] - ppc.out1$fit.y.group.quants[3,]),
                                          c(ppc.out3$fit.y.rep.group.quants[3,] - ppc.out3$fit.y.group.quants[3,])),
                                Site  = rep(SS_data_New$Field_Number, 2),
                                group = c(rep("Freeman-Tukey 1", length(ppc.out1$fit.y.rep.group.quants[3,])),
                                          rep("Chi-square 1",    length(ppc.out3$fit.y.rep.group.quants[3,]))),
                                Year  = rep(c(rep("2022", 101),
                                              rep("2023", 99)), 2))
  
  temp.plot2<-ggplot(diff.fit2, aes(y=fit, x=Site, group=group, color=Year))+
    geom_point()+
    facet_wrap(~group)+
    scale_color_manual(values = c("pink", "darkblue","darkgrey")) +
    labs(y="Replicate - true discrepancy", x = "Site")+
    theme(legend.position = 'none',
          axis.text.x=element_blank())+ 
    facet_grid(group~Year, scales="free_x", space = "free_x")
  ggsave(temp.plot2, file=paste0("Results/2022-2023/Model-Fit/Outliers_M", i, ".tiff"),
         height=6, width=7.5, units='in')
}

# Fill in the significance values manually after the models have run
Significance <- cbind.data.frame(Model = c(rep("Model 1",4), rep("Model 2",4),
                                           rep("Model 3",4), rep("Model 4",4),
                                           rep("Model 5",4), rep("Model 6",4),
                                           rep("Model 7",4)),
                                 Appr  = c("FT-1", "FT-2", "CS-1", "CS-2"),
                                 Sign. = c(0.9893 ,0.538 ,0.976 ,0.5303 ,0.9377 ,
                                           0.5763 ,0.3603 ,0.574 ,0.7183 ,
                                           0.4367 ,0.2157 ,0.4523 ,0.9877 ,
                                           0.5557 ,0.9673 ,0.5417 ,0.6517 ,
                                           0.377 ,0.118 ,0.398 ,0.8817 ,
                                           0.5693 ,0.682 ,0.544 ,0.5387 ,
                                           0.2897 ,0.1477 ,0.3047))

################################################################################
# Spatial Models
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
  ggsave(temp.plot, file=paste0("Results/2022-2023/Model-Fit/Model_fit_M", i, "sp.tiff"),
         height=6, width=7.5, units='in')
  
  # check out what sites may be limiting the fit
  diff.fit2 <- cbind.data.frame(fit   = c(c(ppc.out1$fit.y.rep.group.quants[3,] - ppc.out1$fit.y.group.quants[3,]),
                                          c(ppc.out3$fit.y.rep.group.quants[3,] - ppc.out3$fit.y.group.quants[3,])),
                                Site  = rep(SS_data_New$Field_Number, 2),
                                group = c(rep("Freeman-Tukey 1", length(ppc.out1$fit.y.rep.group.quants[3,])),
                                          rep("Chi-square 1",    length(ppc.out3$fit.y.rep.group.quants[3,]))),
                                Year  = rep(c(rep("2022", 101),
                                              rep("2023", 99)), 2))
  
  temp.plot2<-ggplot(diff.fit2, aes(y=fit, x=Site, group=group, color=Year))+
    geom_point()+
    facet_wrap(~group)+
    scale_color_manual(values = c("pink", "darkblue","darkgrey")) +
    labs(y="Replicate - true discrepancy", x = "Site")+
    theme(legend.position = 'none',
          axis.text.x=element_blank())+ 
    facet_grid(group~Year, scales="free_x", space = "free_x")
  ggsave(temp.plot2, file=paste0("Results/2022-2023/Model-Fit/Outliers_M", i, "sp.tiff"),
         height=6, width=7.5, units='in')
}

# Fill in the significance values manually after the models have run
Significance2 <- cbind.data.frame(Model = c(rep("Spatial model 1",4), rep("Spatial model 2",4),
                                           rep("Spatial model 3",4), rep("Spatial model 4",4)),
                                 Appr  = c("FT-1", "FT-2", "CS-1", "CS-2"),
                                 Sign. = c(0.7087 ,0.4543 ,0.2207 ,0.4457 ,
                                           0.6547 ,0.3837 ,0.118 ,0.384 ,
                                           0.8923 ,0.5707 ,0.6943 ,0.5713 ,
                                           0.6023 ,0.3473 ,0.2297 ,0.3443 ))
Fit.statistics <- rbind.data.frame(Significance, Significance2)
write.csv(Fit.statistics, "Results/2022-2023/Model-Fit/Fit-statistics-2022-2023.csv")
