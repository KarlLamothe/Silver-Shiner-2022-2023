################################################################################
################################################################################
# Inspecting traceplots
# RScript02.3-Occupancy-2022_2023-Traceplots.R 
# This is Script 8 of occupancy modelling
# This script is used to print traceplots for single-species occupancy models 
# built for Silver Shiner using data collected from 16 Mile Creek in 2022 and 
# 2023.
################################################################################
################################################################################
library(coda)

# load the results
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

# ~~~~~~~~~~ #
# Traceplots #
# ~~~~~~~~~~ #
create_trace_lists <- function(out, n.chains) {
  n.post.samples     <- out$n.post
  beta.samples.list  <- list()
  alpha.samples.list <- list()
  
  for (i in 1:n.chains) {
    beta.samples.list[[i]]  <- mcmc(out$beta.samples[((i - 1)  * n.post.samples + 1):(i * n.post.samples), ])
    alpha.samples.list[[i]] <- mcmc(out$alpha.samples[((i - 1) * n.post.samples + 1):(i * n.post.samples), ])
  }
  
  beta.samples.list  <- as.mcmc.list(beta.samples.list)
  alpha.samples.list <- as.mcmc.list(alpha.samples.list)
  
  return(list(beta.samples = beta.samples.list, alpha.samples = alpha.samples.list))
}

################################################
# Run function for each model
################################################
n.chains <- 3

out1.trace <- create_trace_lists(out1, n.chains)
out2.trace <- create_trace_lists(out2, n.chains)
out3.trace <- create_trace_lists(out3, n.chains)
out4.trace <- create_trace_lists(out4, n.chains)
out5.trace <- create_trace_lists(out5, n.chains)
out6.trace <- create_trace_lists(out6, n.chains)
out7.trace <- create_trace_lists(out7, n.chains)

#######################################################################
# Spatial occupancy models                                            #
# Function to create alpha and beta traceplots individually by chains #
#######################################################################
create_trace_lists.sp <- function(out, n.chains) {
  n.post.samples     <- out$n.post
  beta.samples.list  <- list()
  alpha.samples.list <- list()
  theta.samples.list <- list()
  
  for (i in 1:n.chains) {
    beta.samples.list[[i]]  <- mcmc(out$beta.samples[((i - 1)  * n.post.samples + 1):(i * n.post.samples), ])
    alpha.samples.list[[i]] <- mcmc(out$alpha.samples[((i - 1) * n.post.samples + 1):(i * n.post.samples), ])
    theta.samples.list[[i]] <- mcmc(out$theta.samples[((i - 1) * n.post.samples + 1):(i * n.post.samples), ])
  }
  beta.samples.list  <- as.mcmc.list(beta.samples.list)
  alpha.samples.list <- as.mcmc.list(alpha.samples.list)
  theta.samples.list <- as.mcmc.list(theta.samples.list)
  return(list(beta.samples = beta.samples.list, alpha.samples = alpha.samples.list,
              theta.samples = theta.samples.list))
}

################################################
# Run function for each model
################################################
outsp1.trace <- create_trace_lists.sp(out.sp1, n.chains)
outsp2.trace <- create_trace_lists.sp(out.sp2, n.chains)
outsp3.trace <- create_trace_lists.sp(out.sp3, n.chains)
outsp4.trace <- create_trace_lists.sp(out.sp4, n.chains)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 1 # occ ~ 1, det ~ 1 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
par(mar = c(4, 3, 1.5, 1.5)) 
tiff("Results/2022-2023/Traceplots/M1_1_det_2022_2023.tiff", res=800, height=3, width=7.5, units='in')
plot(out1.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/M1_1_occ_2022_2023.tiff", res=800, height=3, width=7.5, units='in')
plot(out1.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 2 # occ ~ 1, det ~ Discharge       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/M2_D_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(out2.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/M2_1_occ_2022_2023.tiff", res=800, height=3, width=7.5, units='in')
plot(out2.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 3 # occ ~ 1, det ~ Year                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/M3_Y_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(out3.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/M3_1_occ_2022_2023.tiff", res=800, height=3, width=7.5, units='in')
plot(out3.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 4 # occ ~ 1, det ~ Depth            #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/M4_Y_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(out4.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/M4_DY_occ_2022_2023.tiff", res=800, height=3, width=7.5, units='in')
plot(out4.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 5 # occ ~ Depth, det ~ Year     
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/M5_Y_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(out5.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/M5_DY_occ_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(out5.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 6 # occ ~ Year, det ~ Year      
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/M6_Y_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(out6.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/M6_DY_occ_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(out6.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 7 # occ ~ Depth + Year, det ~ Year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/M7_Y_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(out7.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/M7_DY_occ_2022_2023.tiff", res=800, height=9, width=7.5, units='in')
plot(out7.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#             # Spatial Models #           #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# MODEL 1 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/MSP1_1_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp1.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/MSP1_1_occ_2022_2023.tiff", res=800, height=3, width=7.5, units='in')
plot(outsp1.trace$beta.samples,  density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/MSP1_1_theta_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp1.trace$theta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 2       
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/MSP2_D_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp2.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/MSP2_1_occ_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp2.trace$beta.samples,  density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/MSP2_1_theta_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp2.trace$theta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 3 #                        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/MSP3_Y_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp3.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/MSP3_1_occ_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp3.trace$beta.samples,  density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/MSP3_1_theta_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp3.trace$theta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 4 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2022-2023/Traceplots/MSP4_Y_det_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp4.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/MSP4_dY_occ_2022_2023.tiff", res=800, height=9, width=7.5, units='in')
plot(outsp4.trace$beta.samples,  density = TRUE, las=1)
dev.off()

tiff("Results/2022-2023/Traceplots/MSP4_1_theta_2022_2023.tiff", res=800, height=6, width=7.5, units='in')
plot(outsp4.trace$theta.samples,  density = TRUE, las=1)
dev.off()
