################################################################################
################################################################################
# Plot model coefficients
# RScript02.4-Occupancy-2022_2023-Model_Coefficients.R 
# This is Script 9 of occupancy modelling
# This script is used to plot model coefficients for single-species occupancy 
# models built for Silver Shiner using data collected from 16 Mile Creek in 2022 
# and 2022.
################################################################################
################################################################################
library(ggplot2)

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

# Import data
ModelCoefs.tot<- read.csv("Results/2022-2023/Models/Model_coefficients_2022_2023.csv",header=T)

# ~~~~~~~~~~~~~~~~~~~~~~ #
# PLOT #
# ~~~~~~~~~~~~~~~~~~~~~~ #
ModelCoefs.tot1 <- ModelCoefs.tot[ModelCoefs.tot$ModelNo == "Model 1" |
                                     ModelCoefs.tot$ModelNo == "Model 2" |
                                     ModelCoefs.tot$ModelNo == "Model 3" |
                                     ModelCoefs.tot$ModelNo == "Model 4" |
                                     ModelCoefs.tot$ModelNo == "Model 1 - spatial" |
                                     ModelCoefs.tot$ModelNo == "Model 2 - spatial" |
                                     ModelCoefs.tot$ModelNo == "Model 3 - spatial" |
                                     ModelCoefs.tot$ModelNo == "Model 4 - spatial",]
ModelCoefs.tot1 <- ModelCoefs.tot1[!ModelCoefs.tot1$Submodel == 'Spatial',]

ModelCoefs.tot1$Submodel <- factor(ModelCoefs.tot1$Submodel,
                                   levels=c("Det","Occ"), label=c("Detection probability", "Occupancy probability"))
ModelCoefs.tot1$Variable <- factor(ModelCoefs.tot1$Variable,
                                   levels=c("Discharge","Depth","Year 2023","Int"),
                                   labels=c("Discharge","Depth","Year 2023","Intercept"))

occ.det.facet<-ggplot(ModelCoefs.tot1, aes(x=Mean, y=Variable))+
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_pointrange(aes(xmin = CI2.5, xmax = CI97.5, color=ModelNo), size=0.15,
                  position = position_dodge2(width=0.7, preserve="single"),
                  show.legend = TRUE, linewidth=0.5)+
  facet_wrap(~Submodel, scales = "free_x")+
 scale_color_manual(values= c("black","grey","red","orange","gold","darkgreen", "blue","purple"),
                    labels= c(expression(italic('p')~'~1;'~italic("ψ")~'~1'),
                              expression(italic('p')~'~Discharge;'~italic("ψ")~'~1'),
                              expression(italic('p')~'~Year;' ~italic("ψ")~'~1'),
                              expression(italic('p')~'~Depth;'~italic("ψ")~'~1'),
                              expression(italic('p')~'~Year;' ~italic("ψ")~'~sp'),
                              expression(italic('p')~'~Year;' ~italic("ψ")~'~Depth + sp'),
                              expression(italic('p')~'~Year;' ~italic("ψ")~'~Year + sp'),
                              expression(italic('p')~'~Year;' ~italic("ψ")~'~Depth + Year + sp')))+  
  labs(x = "Posterior estimate", y = "Occupancy variable")+
  theme(legend.title    = element_blank(),
        axis.title.y = element_blank(),
        legend.key      = element_blank(),
        legend.text = element_text(hjust=0)) 

tiff("Results/2022-2023/Models/Occ_det_2022_2023_120224.tiff",
     height=2.25, width=7, units='in', res=800)
occ.det.facet
dev.off()
