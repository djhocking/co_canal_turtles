######### Load Libraries #########
##

library(parallel)
library(dplyr)
library(jagsUI)

######### Load Data from Previous script #########

testing <- TRUE
Species <- "PRUB"
run_date <- Sys.Date()

if(testing) {
  ni = 300
  nt = 1
  nc = 2
  nb = 200
} else {
  nb = 2000
  ni = 12000
  nc = 8
  nt = 4
}

######### Load Data from Previous script #########

if(testing) {
  load(file = paste0("Data/Derived/", Species, "_all_site_testing_", run_date, ".RData"))
} else {
  load(file = paste0("Data/Derived/", Species, "_all_site_reg_", run_date, ".RData"))
}

# make M variable by site to speed code
# EDF <- read.csv(file = "Data/EDF.csv", stringsAsFactors = FALSE)
# 
# #Take out sites H and I
# EDF_CPIC <- EDF %>%
#   filter(site != "H" & site != "I" & species == "CPIC")
# str(EDF_CPIC)


# restrict M to being the maximum size of the data array so the jags loops doesn't go out of bounds
M <- if_else(df_M$M > M, M, trunc(df_M$M))
# M <- if_else(df_M$M < max_ind_sp, max_ind_sp, M) # real captures can be up to row 750 at any site, all after are definitely augmentation - not anymore

# zeros <- matrix(NA, max(M), n_sites)
# for(z in 1:n_sites) {
# zoo = matrix(0, n_ind_site$n[z], 1)
# zeros[ , z] = rbind(zoo, matrix(NA, max(M) - n_ind_site$n[z], 1))
# }

zeros <- matrix(0, max(M), n_sites)

######### JAGS Model with Zero Trick ##########


scr_zeros <- cat("
model{
  
  alpha2 ~ dnorm(0, pow(1.5, -2)) # Trap behavior universal distribution across sites
  mu_0 ~ dnorm(0, pow(1.5, -2))
  sd_0 ~ dunif(0, 5)
  # alpha_1_sex ~ dnorm(0, pow(1.5, -2))
  beta_0 ~ dnorm(0, pow(10, -2))
  beta_1 ~ dnorm(0, pow(10, -2))
  beta_2 ~ dnorm(0, pow(10, -2))
  # beta_3 ~ dnorm(0, pow(10, -2))
  # alpha_1_int ~ dnorm(0, pow(3, -2))
  psi_sex ~ dunif(0, 1)
  
  sigma ~ dnorm(3, pow(2, -2))T(0, ) # half normal with mean = 1 -> roughly 1 ha home range. Same for both sexes b/c lack of data to fit separately
  alpha1 <- 1 / (2 * sigma * sigma)
  
  # for(t in 1:2) {
  # sigma[t] <- pow(1 / (2 * alpha1[t]), 0.5) # sd of half normal
  # # log(alpha1[t]) <- alpha_1_int + alpha_1_sex * (t-1) # affect of being female on home range
  # alpha1[t] ~ dnorm(1.5, pow(2, -2))T(0.02, ) # prior home range mean = 250m, median = 200m, and sd = 180
  #   }
    
  for(g in 1:n_sites) {
    psi[g] ~ dunif(0, 1) # prob of individual being in the population
    
    for(k in 1:K) {
      alpha0[g, k] ~ dnorm(mu_0, sd_0)
    }
                 
     for(i in 1:M[g]) { 
     Sex[g, i] ~ dbern(psi_sex)
     Sex2[g, i] <- Sex[g, i] + 1
     z[g, i] ~ dbern(psi[g])
     s[g, i] ~ dunif(xlim[g, 1], xlim[g, 2])
      
      for(j in 1:max_trap[g]) { 
        d[g,i,j] <- abs(s[g, i] - trap_locs[g, j])
                 
        for(k in 1:K) {
          p[g, i, j, k] <- p0[g, i, j, k] * exp(-1 * alpha1 * d[g, i, j] * d[g, i, j]) # alpha1[Sex2[g, i]] 
            } # k
          } # j
        } # i
    
    for(i in 1:n0[g]) {
      for(j in 1:max_trap[g]) { 
        for(k in 1:K) {
          y[i, j, k, g] ~ dbern(p[g, i, j, k] * z[g, i])
          logit(p0[g, i, j, k]) <- alpha0[g, k] + (alpha2 * C[i, k, g])
        }
      }
    }
                 
    for(i in (n0[g] + 1):M[g]) {
      for(k in 1:K) {
        for(j in 1:max_trap[g]) {
          logit(p0[g, i, j, k]) <- alpha0[g, k]
        } # j
      } # k
      zeros[i, g] ~ dbern(1 - prod(1 - p[g, i, 1:max_trap[g], 1:K] * z[g, i]))
    } # i
    
    ##### Derived parameters #####
    N[g] <- sum(z[g , 1:M[g]])
    density[g] <- sum(z[g , 1:M[g]]) / (xlim[g, 2] - xlim[g, 1]) # divided distances by 100 so calculates turtles per 100 m of canal
    
    density_ha[g] <- sum(z[g , 1:M[g]]) / ((xlim[g, 2] - xlim[g, 1]) * 100 * width_m[g]) * 10000 # density per hectare
    
    # p_cap_site[g] <- n0[g] / (n0[g] + sum(z[g , (n0[g]+1):M[g]]))
    for(i in 1:M[g]) {
      p_site_ind[g, i] <- 1 - prod(1 - p[g, i, 1:max_trap[g], 1:K] * z[g, i])
    }
    p_cap_site[g] <- mean(p_site_ind[g, 1:M[g]])
    
#     for(i in 1:M[g]) {
# p_site_ind[g, i] <- 1 - prod(1 - p[g, i, 1:max_trap[g], 1:K] * z[g, i])
#     } # i

    # site_zeros[g] ~ dnorm(density[g] - (beta_0 + beta_1 * forest[g] + beta_2 * depth[g] + beta_3 * width[g]), pow(3, -2))
    site_zeros[g] ~ dnorm(density_ha[g] - (beta_0 + beta_1 * forest[g] + beta_2 * depth[g]), pow(3, -2))
    
  } # g
  
  for(k in 1:K) {
    p_cap_day[k] <- caps_day[k] / sum(N[1:n_sites])
  }
  
  # sigma_mean <- mean(sigma[1:2])
  
  # linear home range
    home_50 <- 2 * 0.675 * sigma * 100
    home_95 <- 2 * 1.96 * sigma * 100
  
}
", file = "Code/JAGS/zero_test.txt")

######### Set data and MCMC Conditions ########
forest <- read.csv(file = "Data/LandUse/Forest_Cover_Spatial.csv", header = FALSE)
depth <- read.csv(file = "Data/LandUse/Avg_Depth_Spatial.csv", header = TRUE)
width <- read.csv(file = "Data/LandUse/Width_Spatial.csv", header = TRUE)

forest_std <- as.numeric(scale(forest))
depth_std <- as.numeric(scale(depth))
width_std <- as.numeric(scale(width$width))

jags_data_site <- list(y = EM_array, 
                       Sex = sex, 
                       trap_locs = trap_locs, 
                       K = n_days, 
                       M = M, 
                       xlim = xlim, 
                       max_trap = n_traps_site$max_trap, 
                       forest = forest_std,
                       depth = depth_std,
                       # width = width_std,
                       width_m = as.numeric(width$width),
                       site_zeros = rep(0, n_sites),
                       C = recaptured, 
                       zeros = zeros,
                       n_sites = n_sites,
                       n0 = n_ind_site$n,
                       caps_day = caps_day$caps)

initsf <- function() {
  psi_st <- runif(1, 0.05, 0.2)
  Z_st <- matrix(NA_integer_, n_sites, max(M))
  for(n in 1:n_sites) {
    Z_st[n, 1:M[n]] <- rbinom(M[n], 1, psi_st)
    Z_st[n, 1:n_ind_site$n[n]] <- 1
  }
  s_st[is.na(Z_st)] <- NA
  return(
    list( 
         z = Z_st, 
         psi = (n_ind_site$n + 1) / rowSums(Z_st, na.rm = TRUE),
         psi_sex = (sum(sex*Z_st) + 1) / sum(Z_st, na.rm = TRUE),
         s = s_st)
  )
}

# rowSums(Z_st, na.rm = TRUE)

## not working: needs testing - all z should be 1 or 0 and psi should be a vector of values between 0 and 1 that are based on z and n0.
# initsf <- function() {
#   psi_st <- runif(1, 0.05, 0.2)
#   Z_st <- matrix(rbinom(n_sites * max(M, 1, psi_st)), n_sites, max(M))
#   for(l in 1:n_sites) {
#     zoo = matrix(1, 1, n_ind_site$n[l])
#     Z_st[l, ] = cbind(zoo, matrix(NA, 1, M[l] - n_ind_site$n[l]))
#   }
#   
#   return(
#   list(s = s_st, 
#        z = Z_st, 
#        psi = n_ind_site$n / rowSums(Z_st), 
#        psi_sex = runif(1, 0.3, 0.8))
#   )
# }

parameters <- c("density", "N", "alpha2", "alpha0", "alpha1", "mu_0", "sd_0", "mu_1", "sd_1", "beta_0", "beta_1", "beta_2", "psi_sex", "p_cap_day", "mu_psi", "sd_psi", "sigma_mean", "sigma", "p_cap_site", "home_50", "home_95", "density_ha") ## "beta_3", "p_site_ind", "sigma", # "C", maybe C or a summary stat, might blow up if saving each activity center "s".

# no p-cap-site for species without caps at all sites
# parameters <- c("density", "N", "alpha2", "alpha0", "alpha1", "mu_0", "sd_0", "mu_1", "sd_1", "beta_0", "beta_1", "beta_2", "psi_sex", "p_cap_day", "mu_psi", "sd_psi", "sigma_mean", "sigma", "p_cap_site", "home_50", "home_95", "density_ha") 

start_zeros <- Sys.time()

out <- jagsUI(data = jags_data_site, 
              inits = initsf, 
              model.file = "Code/JAGS/zero_test.txt", 
              n.chains = nc,
              n.iter = ni, 
              n.burnin = nb, 
              n.thin = nt, 
              parallel = TRUE, 
              codaOnly = parameters, 
              parameters.to.save = parameters)

end_zeros <- Sys.time()

end_zeros - start_zeros

if(testing) {
  saveRDS(out, file = paste0("Results/JAGS/", Species, "_all_sites_reg_testing_", run_date, ".rds"))
} else {
  saveRDS(out, file = paste0("Results/JAGS/", Species, "_all_sites_reg_", run_date, ".rds"))
}



