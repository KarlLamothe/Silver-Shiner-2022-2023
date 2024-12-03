################################################################################
################################################################################
# Check out model coefficients
################################################################################
################################################################################
# load packages
library(pacman)     # downloads and loads packages simultaneously
p_load(ggplot2)     # nice plots
p_load(patchwork)   # combine ggplots

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

################################################################################
# Create plot
ModelCoefs.tot<- read.csv("Results/2011-2016-2022/Models/Model.coefficients.2011-2016-2022.csv", header=T)
ModelCoefs.tot

# Remove spatial variables
ModelCoefs.tot1 <- ModelCoefs.tot[ModelCoefs.tot$ModelNo=="Model 6" |
                                   ModelCoefs.tot$ModelNo=="Model 4 - spatial",]
ModelCoefs.tot1 <- ModelCoefs.tot1[!ModelCoefs.tot1$Variable=="sigma^2",]
ModelCoefs.tot1 <- ModelCoefs.tot1[!ModelCoefs.tot1$Variable=="phi",]

# ordering variables for plotting
ModelCoefs.tot1$Variable <- factor(ModelCoefs.tot1$Variable,
                                   levels = c("Depth", "Year 2022", "Year 2016","Int"),
                                   labels = c("Depth", "Year 2022", "Year 2016","Intercept"))

# creating new label for facet headers
ModelCoefs.tot1$Submodel <- factor(ModelCoefs.tot1$Submodel,
                                   levels = c("Det", "Occ"),
                                   labels = c('Detection Pr.', "Occupancy Pr."))

# only occupancy
Occ.Coefs.tot1 <- ModelCoefs.tot1[ModelCoefs.tot1$Submodel == "Occupancy Pr.",]

dat_text2 <- data.frame(label =  "a)", Submodel = "Occupancy Pr.", Mean = -1.75, Variable = 4.3)
Occ.var.gg<-ggplot(Occ.Coefs.tot1, aes(x=Mean, y=Variable))+
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual(values= c("#000000", "darkgrey"),
                     labels= c("Non-spatial", "Spatial"))+  
  guides(color=guide_legend(position='inside'))+
  geom_pointrange(aes(xmin = CI2.5, xmax = CI97.5, color=ModelNo), size=0.25,
                  position = position_dodge2(width=0.5, preserve="single"),
                  show.legend = TRUE, linewidth=0.5)+
  geom_text(data = dat_text2, mapping = aes(x = Mean, y = Variable, label = label), size=3)+
  labs(x = "Posterior estimate", 
       y = "Occupancy variable")+
  theme(legend.title    = element_blank(),
        legend.position.inside = c(0.2,0.2),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.key      = element_blank(),
        legend.margin   = margin(0, unit='cm'),
        legend.background = element_blank(),
        legend.text = element_text(hjust=0),
        legend.justification.top = "left",
        legend.justification.left = "top",
        legend.justification.bottom = "right") 
Occ.var.gg
# we'll export this figure with a figure from the next script
