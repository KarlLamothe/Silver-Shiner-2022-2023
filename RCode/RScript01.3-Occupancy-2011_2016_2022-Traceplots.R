################################################################################
################################################################################
# Inspecting traceplots
# RScript01.3-Occupancy-2011_2016_2022-Traceplots.R 
# This is Script 3 of occupancy modelling
# This script is used to print traceplots for single-species occupancy models 
# built for Silver Shiner using data collected from 16 Mile Creek in 2011, 2016, 
# and 2022.
################################################################################
################################################################################
library(coda)

# Load results
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

#### Function to create alpha and beta traceplots individually by chains ####
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

# Run function for each nonspatial model
n.chains <- 3
out1.trace <- create_trace_lists(out1, n.chains)
out2.trace <- create_trace_lists(out2, n.chains)
out3.trace <- create_trace_lists(out3, n.chains)
out4.trace <- create_trace_lists(out4, n.chains)
out5.trace <- create_trace_lists(out5, n.chains)
out6.trace <- create_trace_lists(out6, n.chains)

#### spatial models ####
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

out1.sp.trace <- create_trace_lists.sp(out.sp1, n.chains)
out2.sp.trace <- create_trace_lists.sp(out.sp2, n.chains)
out3.sp.trace <- create_trace_lists.sp(out.sp3, n.chains)
out4.sp.trace <- create_trace_lists.sp(out.sp4, n.chains)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## MODEL 1 # occ ~ 1, det ~ 1               #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/M1_det_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out1.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/M1_occ_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out1.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 2 # occ ~ 1, det ~ Depth                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/M2_det_2011_2016_2022.tiff", res=800, height=6, width=7.5, units='in')
plot(out2.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/M2_occ_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out2.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 3 # occ ~ 1, det ~ Year #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/M3_det_2011_2016_2022.tiff", res=800, height=9, width=7.5, units='in')
plot(out3.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/M3_occ_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out3.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 4 # occ ~ Depth, det ~ 1 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/M4_det_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out4.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/M4_occ_2011_2016_2022.tiff", res=800, height=6, width=7.5, units='in')
plot(out4.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 5 # occ ~ Year, det ~ 1     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/M5_det_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out5.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/M5_occ_2011_2016_2022.tiff", res=800, height=9, width=7.5, units='in')
plot(out5.trace$beta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MODEL 6 # occ ~ Depth + Year, det ~ 1       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/M6_det_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out6.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/M6_occ_2011_2016_2022.tiff", res=800, height=12, width=7.5, units='in')
plot(out6.trace$beta.samples,  density = TRUE, las=1)
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## SPATIAL MODEL 1
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/SPM1_det_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out1.sp.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/SPM1_occ_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out1.sp.trace$beta.samples,  density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/SPM1_theta_2011_2016_2022.tiff", res=800, height=6, width=7.5, units='in')
plot(out1.sp.trace$theta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# SPATIAL MODEL 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/SPM2_det_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out2.sp.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/SPM2_occ_2011_2016_2022.tiff", res=800, height=6, width=7.5, units='in')
plot(out2.sp.trace$beta.samples,  density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/SPM2_theta_2011_2016_2022.tiff", res=800, height=6, width=7.5, units='in')
plot(out2.sp.trace$theta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# SPATIAL MODEL 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/SPM3_det_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out3.sp.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/SPM3_occ_2011_2016_2022.tiff", res=800, height=9, width=7.5, units='in')
plot(out3.sp.trace$beta.samples,  density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/SPM3_theta_2011_2016_2022.tiff", res=800, height=6, width=7.5, units='in')
plot(out3.sp.trace$theta.samples,  density = TRUE, las=1)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# SPATIAL MODEL 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
tiff("Results/2011-2016-2022/Traceplots/SPM4_det_2011_2016_2022.tiff", res=800, height=3, width=7.5, units='in')
plot(out4.sp.trace$alpha.samples, density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/SPM4_occ_2011_2016_2022.tiff", res=800, height=12, width=7.5, units='in')
plot(out4.sp.trace$beta.samples,  density = TRUE, las=1)
dev.off()

tiff("Results/2011-2016-2022/Traceplots/SPM4_theta_2011_2016_2022.tiff", res=800, height=6, width=7.5, units='in')
plot(out4.sp.trace$theta.samples,  density = TRUE, las=1)
dev.off()
