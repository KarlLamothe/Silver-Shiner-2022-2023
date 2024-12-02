################################################################################
################################################################################
# Prospective power analysis for occupancy models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Presence-absence data are simulated to reflect proportional declines in 
# occupancy probability between two periods. For each simulation, an occupancy
# model is developed with sampling period as an occupancy covariate. Power to 
# detect a change is estimated as the proportion of simulations with a 
# significant sampling period effect at alpha = 0.05 and 0.20. Depletion 
# sampling can be performed and differences in detection probability can be
# implemented. 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
################################################################################
################################################################################
# load packages
library(spOccupancy)
library(doParallel)  
library(foreach)     
library(beepr)       
library(data.table)
library(ggplot2)
library(xlsx)

################################################################################
# Power analysis function
################################################################################
power_analysis_decline <- function(n.years   = 2,                # number of years in survey (Can only be 2 years right now)
                                   n.sites.per.year = 100,       # number of sites surveyed per year
                                   n.surveys = 3,                # number of surveys performed at each site
                                   declines  = c(0.3, 0.5, 0.7), # proportional decline
                                   psi.init  = 0.8,              # initial occupancy probability (period 1)
                                   Year.p    = FALSE,            # Do you want to include year as a covariate for detection?
                                   p.init    = 0.8,              # initial detection probability
                                   p.post    = 0.8,              # detection probability in time 2
                                   Depletion = TRUE,             # depletion sampling (TRUE) or not (FALSE)
                                   n.sims    = 500,              # number of simulations to run
                                   n.chains  = 3,                # number of chains in model
                                   n.burn    = 1000,             # number of model burn in samples
                                   n.thin    = 25,               # model thinning interval
                                   n.samples = 3000,             # number of posterior samples
                                   n.cores   = 4,                # number of computer cores to run across
                                   verbose   = FALSE) {          # if TRUE, spits out model output from spOccupancy
  
  # Calculate total number of sites across years
  n.sites <- n.years * n.sites.per.year
  
  # Set up parallel backend to use n_cores
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  # Parallelize the loop over different decline values using foreach
  power_results <- foreach(decline = declines, .combine = rbind, .packages = "spOccupancy") %:%
    foreach(i = 1:n.sims, .combine = rbind) %dopar% {
      
      # Define occupancy probabilities based on the specified decline for each year
      psi_post <- psi.init * (1 - decline)  # Year 2: reduced occupancy
      
      # Assign occupancy states based on the calculated probabilities
      # This is a challenge here. How do you reduce the population. Here we are
      # assuming that the sample is random each year
      z <- numeric(n.sites)
      z[1:n.sites.per.year] <- rbinom(n.sites.per.year, 1, psi.init)
      z[(n.sites.per.year + 1):n.sites] <- rbinom(n.sites.per.year, 1, psi_post)
      
      # Simulate detections based on year-specific detection probabilities
      # Year 1 detections based on initial detection probability
      y <- matrix(0, nrow = n.sites, ncol = n.surveys)  # detection matrix
      y[1:n.sites.per.year,] <- rbinom(n.sites.per.year * n.surveys, 1, z[1:n.sites.per.year] * p.init)
      
      # Period 2 detections based on Year 2 (post) detection probability
      y[(n.sites.per.year + 1):n.sites,] <- rbinom(n.sites.per.year * n.surveys, 1, z[(n.sites.per.year + 1):n.sites] * p.post)
      
      # convert data frame for multinomial estimation if depletion sampling performed
      if(Depletion == TRUE){
        y = t(sapply(1:nrow(y), function(i) {
          y <- y[i,]
          if(sum(y == 0) < n.surveys) { # no detection - do nothing
            x = min(which(y == 1))      # On which samples where species detected
            if(x < n.surveys) {         # If detection only on last survey - do nothing
              y[(x+1):n.surveys] <- NA  # if detection before last survey - assign surveys after to NA
            }
          }
          y # output vector of survey detections
        }))
        }
      
      # Create year assignment for each site
      year <- rep(1:n.years, each = n.sites.per.year)  # Assign each year its own sites
      
      # Prepare data for spOccupancy
      occ_data <- list(y = y, 
                       occ.covs = data.frame(year = year),  # Year as a covariate
                       det.covs = data.frame(year = year))
      
      # fit models
      if(Year.p == FALSE){
        fit <- PGOcc(occ.formula = ~ factor(year),  # Year as a factor
                     det.formula = ~ 1,
                     n.samples = n.samples, n.chains = n.chains, data = occ_data, 
                     n.burn = n.burn, n.thin = n.thin, verbose = verbose,
                     n.omp.threads = 3)
      } else {
        fit <- PGOcc(occ.formula = ~ factor(year),
                     det.formula = ~ factor(year),
                     n.samples = n.samples, n.chains = n.chains, data = occ_data, 
                     n.burn = n.burn, n.thin = n.thin, verbose = verbose,
                     n.omp.threads = 3)   
      }

    # Extracts posterior samples for year variable
    post_samples <- fit$beta.samples[, "factor(year)2"]  # Adjust this based on output

    # Calculate confidence intervals
    ci_95 <- quantile(post_samples, probs = c(0.025, 0.975))
    ci_80 <- quantile(post_samples, probs = c(0.1, 0.9))
    
    # Determine significance
    significant_95 <- ifelse(ci_95[1] > 0 | ci_95[2] < 0, 1, 0)
    significant_80 <- ifelse(ci_80[1] > 0 | ci_80[2] < 0, 1, 0)
    
    # Return results for each simulation
    return(cbind(Decline = decline, 
                 significant_95 = significant_95, 
                 significant_80 = significant_80))
    }
    
  # Stop parallel backend
  stopCluster(cl)
  
  # Summarize power results by decline level
  power_summary <- data.table(power_results)[, .(power95 = mean(significant_95), 
                                                 power80 = mean(significant_80)), by = Decline]
  
  return(power_summary)
}

################################################################################
################ ~~~~~~~~~~~~~~~~~ Scenarios ~~~~~~~~~~~~~~~~~ #################
################################################################################
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 75 sites, constant p = 0.81
## 1000 simulations, 5000 samples
## 3 chains, 1000 burn ins, 10 thinning interval
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 75,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50, 0.70),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 38.06436   mins
#print(power_result)
##30% ### power95 = 0.383   ; power80 = 0.677
##50% ### power95 = 0.810   ; power80 = 0.936
##70% ### power95 = 0.986   ; power80 = 0.999
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 100 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result1 <- power_analysis_decline(
#  n.sites.per.year = 100,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50, 0.70),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 49.60283 mins
#print(power_result1)
##30% ### power95 = 0.473; power80 = 0.744
##50% ### power95 = 0.923; power80 = 0.974
##70% ### power95 = 0.998; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 125 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result2 <- power_analysis_decline(
#  n.sites.per.year = 125,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 41.60741  mins
#print(power_result2)
##30% ### power95 = 0.554; power80 = 0.799
##50% ### power95 = 0.962; power80 = 0.992
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 150 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result3 <- power_analysis_decline(
#  n.sites.per.year = 150,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 53.47947 mins
#print(power_result3)
##30% ### power95 = 0.656; power80 = 0.860
##50% ### power95 = 0.982; power80 = 0.999
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 175 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result4 <- power_analysis_decline(
#  n.sites.per.year = 175,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 1.18646 hours
#print(power_result4)
##30% ### power95 = 0.733; power80 = 0.902
##50% ### power95 = 0.994; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 200 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time())
#power_result5 <- power_analysis_decline(
#  n.sites.per.year = 200,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 33.73389  mins
#print(power_result5)
##30% ### power95 = 0.766; power80 = 0.914
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 75 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result6 <- power_analysis_decline(
#  n.sites.per.year = 75,
#  psi.init   = 0.70, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50, 0.70),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 45.12707   mins
#print(power_result6)
##30% ### power95 = 0.707; power80 = 0.893
##50% ### power95 = 0.992; power80 = 1.000
##70% ### power95 = 1.000; power80 = 1.000 
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 100 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result7 <- power_analysis_decline(
#  n.sites.per.year = 100,
#  psi.init   = 0.70, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 44.36945  mins
#print(power_result7)
##30% ### power95 = 0.859; power80 = 0.967
##50% ### power95 = 1.000; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 125 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result8 <- power_analysis_decline(
#  n.sites.per.year = 125,
#  psi.init   = 0.70, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 20.3695  mins
#print(power_result8)
##30% ### power95 = 0.921; power80 = 0.984
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 150 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result9 <- power_analysis_decline(
#  n.sites.per.year = 150,
#  psi.init   = 0.70, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 25.20975  mins
#print(power_result9)
##30% ### power95 = 0.962; power80 = 0.993
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 175 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result10 <- power_analysis_decline(
#  n.sites.per.year = 175,
#  psi.init   = 0.70, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 28.9305   mins
#print(power_result10)
##30% ### power95 = 0.976; power80 = 0.998
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 200 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time())
#power_result11 <- power_analysis_decline(
#  n.sites.per.year = 200,
#  psi.init   = 0.70, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  mins
#print(power_result11)
##30% ### power95 = ; power80 = 
#
################################################################################
# Changing the number of surveys per site, n = 4
################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# psi = 0.43, Depletion, 75 sites, constant p = 0.81
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 75,
#  psi.init   = 0.43, 
#  n.surveys  = 4,
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50, 0.70),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  33.66132  mins
#print(power_result)
##30% ### power95 = 0.374; power80 = 0.662
##50% ### power95 = 0.832; power80 = 0.936   
##70% ### power95 = 0.993; power80 = 0.997
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 100 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 100,
#  psi.init   = 0.43, 
#  n.surveys  = 4,
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 37.3201 mins
#print(power_result)
##30% ### power95 = 0.475; power80 = 0.736
##50% ### power95 = 0.903; power80 = 0.971
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 125 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 125,
#  psi.init   = 0.43, 
#  n.surveys  = 4,
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 43.25369   mins
#print(power_result)
##30% ### power95 = 0.546   ; power80 = 0.798
##50% ### power95 = 0.970   ; power80 = 0.993
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 150 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 150,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.surveys  = 4,
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 52.01973   mins
#print(power_result)
##30% ### power95 = 0.674; power80 = 0.867
##50% ### power95 = 0.980; power80 = 0.996
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 175 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 175,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.surveys  = 4,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 54.56436    mins
#print(power_result)
##30% ### power95 = 0.715; power80 = 0.890
##50% ### power95 = 0.990; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 200 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time())
#power_result <- power_analysis_decline(
#  n.sites.per.year = 200,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.surveys  = 4,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  58.0598  mins
#print(power_result)
##30% ### power95 = 0.764    ; power80 = 0.916
##50% ### power95 = 0.995   ; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 75 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 75,
#  psi.init   = 0.70, 
#  n.surveys  = 4,
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50, 0.70),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #35.28265  mins
#print(power_result)
##30% ### power95 = 0.728; power80 = 0.908
##50% ### power95 = 0.989; power80 = 0.998
##70% ### power95 = 1.000; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 100 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 100,
#  psi.init   = 0.70, 
#  n.surveys  = 4,
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  31.33799  mins
#print(power_result)
##30% ### power95 = 0.856; power80 = 0.957
##50% ### power95 = 0.997; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 125 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 125,
#  psi.init   = 0.70, 
#  n.surveys  = 4,
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 20.95398   mins
#print(power_result)
##30% ### power95 = 0.929; power80 = 0.985
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 150 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 150,
#  psi.init   = 0.70, 
#  n.surveys  = 4,
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 24.7405   mins
#print(power_result)
##30% ### power95 = 0.962; power80 = 0.989
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 175 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 175,
#  psi.init   = 0.70, 
#  n.sims     = 1000, 
#  n.surveys  = 4,
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 28.47759  mins
#print(power_result)
##30% ### power95 = 0.978; power80 = 0.997
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 200 sites, constant p = 0.81
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time())
#power_result <- power_analysis_decline(
#  n.sites.per.year = 200,
#  psi.init   = 0.70, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.surveys  = 4,
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = 0.30,   
#  p.init     = 0.81, 
#  p.post     = 0.81,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 32.66046  mins
#print(power_result)
##30% ### power95 = 0.987; power80 = 0.998
#
#################################################################################
## Changing p
#################################################################################
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 75 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 75,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.surveys  = 3, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50, 0.70),   
#  p.init     = 0.40, 
#  p.post     = 0.40,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 48.91112 mins
#print(power_result)
##30% ### power95 = 0.216; power80 = 0.464 
##50% ### power95 = 0.572; power80 = 0.825
##70% ### power95 = 0.932; power80 = 0.993
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 100 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 100,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.surveys  = 3, 
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50, 0.70),   
#  p.init     = 0.40, 
#  p.post     = 0.40,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 1.100657  hours
#print(power_result)
##30% ### power95 = 0.307; power80 = 0.555
##50% ### power95 = 0.754; power80 = 0.916
##70% ### power95 = 0.978; power80 = 0.999    
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 125 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 125,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.surveys  = 3, 
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50, 0.70),   
#  p.init     = 0.40, 
#  p.post     = 0.40,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 1.387695 hours
#print(power_result)
##30% ### power95 = 0.357; power80 = 0.623
##50% ### power95 = 0.837; power80 = 0.954
##70% ### power95 = 0.995; power80 = 1.000    
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 150 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 150,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.surveys  = 3, 
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.40, 
#  p.post     = 0.40,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  1.120055 hours
#print(power_result)
##30% ### power95 = 0.449; power80 = 0.708
##50% ### power95 = 0.928; power80 = 0.997 
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 175 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 175,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.surveys  = 3, 
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.40, 
#  p.post     = 0.40,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  1.346156 hours
#print(power_result)
##30% ### power95 = 0.550   ; power80 = 0.787
##50% ### power95 = 0.946   ; power80 = 0.992 
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 200 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time())
#power_result <- power_analysis_decline(
#  n.sites.per.year = 200,
#  psi.init   = 0.43, 
#  n.sims     = 1000, 
#  n.samples  = 5000, 
#  n.chains   = 3,
#  n.surveys  = 3, 
#  n.burn     = 1000,
#  n.thin     = 10,
#  n.years    = 2,     
#  declines   = c(0.30, 0.50),   
#  p.init     = 0.40, 
#  p.post     = 0.40,
#  Year.p     = FALSE,
#  Depletion  = TRUE,
#  verbose    = FALSE,
#  n.cores    = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   1.628137  hours
#print(power_result)
##30% ### power95 = 0.588; power80 = 0.820
##50% ### power95 = 0.974; power80 = 0.996
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 75 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 75, psi.init = 0.70, n.sims = 1000, n.samples = 5000, 
#  n.chains = 3, n.surveys = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50, 0.70), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  45.84577  mins
#print(power_result)
##30% ### power95 = 0.401; power80 = 0.649 
##50% ### power95 = 0.882; power80 = 0.967
##70% ### power95 = 1.000; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 100 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 100, psi.init = 0.70,  n.sims = 1000, n.samples = 5000, 
#  n.chains = 3, n.surveys = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   mins
#print(power_result)
##30% ### power95 = 0.535; power80 = 0.785
##50% ### power95 = 0.956; power80 = 0.992
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 125 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 125, psi.init = 0.70, n.sims = 1000, n.samples = 5000, 
#  n.chains = 3, n.surveys = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   mins
#print(power_result)
##30% ### power95 = 0.670; power80 = 0.870
##50% ### power95 = 0.989; power80 = 0.997 
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 150 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 150, psi.init = 0.70, n.sims = 1000, n.samples = 5000, 
#  n.chains = 3, n.surveys = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   mins
#print(power_result)
##30% ### power95 = 0.739; power80 = 0.901
##50% ### power95 = 0.999; power80 = 0.999
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 175 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 175, psi.init = 0.70, n.sims = 1000, n.samples = 5000, 
#  n.chains = 3, n.surveys = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   mins
#print(power_result)
##30% ### power95 = 0.840; power80 = 0.946
##50% ### power95 = 1.000; power80 = 1.000 
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 200 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time())
#power_result <- power_analysis_decline(
#  n.sites.per.year = 200, psi.init = 0.70, n.sims = 1000, n.samples = 5000, 
#  n.chains = 3, n.surveys = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  mins
#print(power_result)
##30% ### power95 = 0.901; power80 = 0.974 
##50% ### power95 = 1.000; power80 = 1.000
#
#################################################################################
## Changing the number of surveys per site, n = 4
#################################################################################
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 75 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 75, psi.init = 0.43, n.surveys = 4, n.sims = 1000, 
#  n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50, 0.70), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 40.22985 mins
#print(power_result)
##30% ### power95 = 0.306; power80 = 0.576
##50% ### power95 = 0.721; power80 = 0.900
##70% ### power95 = 0.970; power80 = 0.994 
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 100 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 100, psi.init = 0.43, n.surveys = 4, n.sims = 1000, 
#  n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50, 0.70), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   1.114829 hours
#print(power_result)
##30% ### power95 = 0.379; power80 = 0.669
##50% ### power95 = 0.836; power80 = 0.950
##70% ### power95 = 0.994; power80 = 0.999
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 125 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 125, psi.init = 0.43, n.surveys = 4, n.sims = 1000, 
#  n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  57.61772  mins
#print(power_result)
##30% ### power95 = 0.474; power80 = 0.741
##50% ### power95 = 0.900; power80 = 0.980
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 150 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 150, psi.init = 0.43, n.sims = 1000, n.surveys = 4,
#  n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  51.34351 mins
#print(power_result)
##30% ### power95 = 0.567; power80 = 0.800
##50% ### power95 = 0.965; power80 = 0.995
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 175 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 175, psi.init = 0.43, n.sims = 1000, n.samples = 5000, 
#  n.chains = 3, n.surveys = 4, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #    1.002899  hours
#print(power_result)
##30% ### power95 = 0.613; power80 = 0.839
##50% ### power95 = 0.975; power80 = 0.997
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.43, Depletion, 200 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time())
#power_result <- power_analysis_decline(
#  n.sites.per.year = 200, psi.init = 0.43, n.sims = 1000, n.samples = 5000, 
#  n.chains = 3, n.surveys = 4, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   1.233032 hours
#print(power_result)
##30% ### power95 = 0.699; power80 = 0.887
##50% ### power95 = 0.983; power80 = 0.998
#
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 75 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 75, psi.init = 0.70, n.surveys = 4, n.sims = 1000, 
#  n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10,
#  n.years = 2, declines = c(0.30, 0.50, 0.70), p.init = 0.40, p.post = 0.40,
#  Year.p = FALSE, Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 46.916 mins
#print(power_result)
##30% ### power95 = 0.544; power80 = 0.780
##50% ### power95 = 0.958; power80 = 0.994
##70% ### power95 = 1.000; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 100 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
# n.sites.per.year = 100, psi.init = 0.70, n.surveys = 4, n.sims = 1000, 
# n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
# declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40, Year.p = FALSE,
# Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   42.28654 mins
#print(power_result)
##30% ### power95 = 0.711; power80 = 0.883
##50% ### power95 = 0.986; power80 = 1.000
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 125 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#(ptm <- Sys.time()) # time the run
#power_result <- power_analysis_decline(
#  n.sites.per.year = 125, psi.init = 0.70, n.surveys = 4, n.sims = 1000, 
#  n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30, 0.50), p.init = 0.40, p.post = 0.40, Year.p = FALSE,
#  Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  50.97366  mins
#print(power_result)
##30% ### power95 = 0.797; power80 = 0.940
##50% ### power95 = 0.998; power80 = 0.999
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 150 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 150, psi.init = 0.70, n.surveys = 4, n.sims = 1000, 
#  n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2, 
#  declines = c(0.30), p.init = 0.40, p.post = 0.40, Year.p = FALSE, 
#  Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #  31.63667  mins
#print(power_result)
##30% ### power95 = 0.881; power80 = 0.961
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 175 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time()) 
#power_result <- power_analysis_decline(
#  n.sites.per.year = 175, psi.init = 0.70, n.sims = 1000, n.surveys = 4,
#  n.samples = 5000, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30), p.init = 0.40, p.post = 0.40, Year.p = FALSE,
#  Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm # 36.9117  mins
#print(power_result)
##30% ### power95 = 0.916    ; power80 = 0.98
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## psi = 0.70, Depletion, 200 sites, constant p = 0.4
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#(ptm <- Sys.time())
#power_result <- power_analysis_decline(
#  n.sites.per.year = 200, psi.init = 0.70, n.sims = 1000, n.samples = 5000, 
#  n.surveys = 4, n.chains = 3, n.burn = 1000, n.thin = 10, n.years = 2,     
#  declines = c(0.30), p.init = 0.40, p.post = 0.40, Year.p = FALSE,
#  Depletion = TRUE, verbose = FALSE, n.cores = 4)
#beep(sound = 'mario')
#Sys.time() - ptm #   mins
#print(power_result)
##30% ### power95 = ; power80 = 
#
#################################################################################
################################################################################
# plot results
#################################################################################
################################################################################

# Personal ggplot theme
theme_set(theme_bw() + 
            theme(axis.title   = element_text(size=9, family="sans", colour="black"),
                  axis.text.x  = element_text(size=9, family="sans", colour="black"),
                  axis.text.y  = element_text(size=9, family="sans", colour="black"),
                  legend.text  = element_text(size=9, family="sans", colour="black"),
                  legend.title = element_text(size=9,  family='sans', colour='black'),
                  plot.title   = element_text(size=11, family="sans", colour="black"),
                  strip.text   = element_text(size=9,  family="sans", colour="black"),
                  panel.border = element_rect(colour = "black", fill=NA),
                  axis.ticks   = element_line(colour = "black"),
                  legend.background = element_blank(),
                  legend.margin = margin(0, 0, 0, 0),
                  legend.spacing.x = unit(0, "mm"),
                  legend.spacing.y = unit(0, "mm")))

# read in power data
Power.df <- read.xlsx("Power-simulations-results.xlsx", sheetName = "Sheet1", header=T)

# create data frame for plotting
head(Power.df)
Powergg.df <- cbind.data.frame(Psi.0 = as.factor(rep(Power.df$Initial.Psi,2)),
                               Sites = rep(Power.df$Sites, 2),
                               Surveys = rep(Power.df$Surveys, 2),
                               Reduction = rep(Power.df$Reduction, 2),
                               alpha = rep(c("α = 0.05","α = 0.20"), each = 144),
                               Power = c(Power.df$alpha...0.05, Power.df$alpha...0.20),
                               p     = rep(Power.df$p, 2))

# Revise the levels and labels for psi, nsurveys, and % reduction
Powergg.df$Psi.0     <- factor(Powergg.df$Psi.0,
                               levels=c("0.43","0.7"),
                               labels=c("ψ = 0.43","ψ = 0.70"))
Powergg.df$Surveys   <- factor(Powergg.df$Surveys,
                               levels=c(3,4),
                               labels=c("k = 3","k = 4"))
Powergg.df$Reduction <- factor(Powergg.df$Reduction,
                               levels=c(0.3,0.5,0.7),
                               labels=c("30%","50%","70%"))
Powergg.df$p     <- factor(Powergg.df$p,
                               levels=c("0.81","0.4"),
                               labels=c("p = 0.81","p = 0.40"))

# plot
ggplot(Powergg.df, aes(x=as.factor(Reduction), y=Power, group=as.factor(Sites)))+
  geom_hline(yintercept = 0.81, lty = 'dashed')+
  geom_hline(yintercept = 0.95, lty = 'dashed')+
  geom_point(aes(color=as.factor(Sites)))+
  scale_color_manual(values=c("grey","darkred","orange","darkgreen","blue","purple"))+
  geom_line(aes(color=as.factor(Sites)))+
  facet_grid(Psi.0+p~Surveys+alpha)+
  labs(x="Percent reduction", color="Sites")

alpha05df <- Powergg.df[Powergg.df$alpha=="α = 0.05",]

ggplot(alpha05df, aes(x=as.factor(Reduction), y=Power, group=as.factor(Sites)))+
  geom_hline(yintercept = 0.81, lty = 'dashed')+
  geom_hline(yintercept = 0.95, lty = 'dashed')+
  geom_point(aes(color=as.factor(Sites)))+
  scale_color_manual(values=c("grey","darkred","orange","darkgreen","blue","purple"))+
  geom_line(aes(color=as.factor(Sites)))+
  facet_grid(p~Psi.0+Surveys)+
  labs(x="Percent reduction", color="Sites")

alpha20df <- Powergg.df[Powergg.df$alpha=="α = 0.20",]

ggplot(alpha20df, aes(x=as.factor(Reduction), y=Power, group=as.factor(Sites)))+
  geom_hline(yintercept = 0.81, lty = 'dashed')+
  geom_hline(yintercept = 0.95, lty = 'dashed')+
  geom_point(aes(color=as.factor(Sites)))+
  scale_color_manual(values=c("grey","darkred","orange","darkgreen","blue","purple"))+
  geom_line(aes(color=as.factor(Sites)))+
  facet_grid(p~Psi.0+Surveys)+
  labs(x="Percent reduction", color="Sites")
