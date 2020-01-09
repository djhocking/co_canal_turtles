######### Load Libraries #########
##

library(dplyr)
library(jagsUI)

######### Load Data from Previous script #########

testing <- TRUE
run_date <- Sys.Date()

if(testing) {
  ni = 101
  nt = 1
  nc = 2
  nb = 1
} else {
  na = 5000
  nb = 5000
  ni = 25000
  nc = min(6, detectCores())
  nt = 3
}

######### Load Data from Previous script #########

if(testing) {
  load(file = paste0("Data/Derived/all_site_testing_", run_date, ".RData"))
} else {
  load(file = "Data/Derived/all_site.RData")
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


######### JAGS Model with Zero Trick ##########


scr_zeros <- cat("
model{
  
  alpha2 ~ dnorm(0, pow(1.5, -2)) # Trap behavior universal distribution across sites
  mu_0 ~ dnorm(0, pow(1.5, -2))
  sd_0 ~ dunif(0, 3)
  # alpha_1_sex ~ dnorm(0, pow(1.5, -2))
  beta_0 ~ dnorm(0, pow(10, -2))
  beta_1 ~ dnorm(0, pow(10, -2))
  beta_2 ~ dnorm(0, pow(10, -2))
  beta_3 ~ dnorm(0, pow(10, -2))
  # alpha_1_int ~ dnorm(0, pow(3, -2))
  psi_sex ~ dbeta(1, 1)
  
  for(t in 1:2) {
  sigma[t] <- pow(1 / (2 * alpha1[t]), 0.5) # sd of half normal
  # log(alpha1[t]) <- alpha_1_int + alpha_1_sex * (t-1) # affect of being female on home range
  alpha1[t] ~ dunif(0, 20)
    }
    
  for(g in 1:n_sites) {
    psi[g] ~ dbeta(1, 1) # prob of individual being in the population
    
    for(k in 1:K) {
      alpha0[g, k] ~ dnorm(mu_0, sd_0)
    }
                 
     for(i in 1:M[g]) {  ## here
     Sex[g, i] ~ dbern(psi_sex)
     Sex2[g, i] <- Sex[g, i] + 1
     z[g, i] ~ dbern(psi[g])
     s[g, i] ~ dunif(xlim[g, 1], xlim[g, 2])
      
      for(j in 1:max_trap[g]) { 
        d[g,i,j] <- abs(s[g, i] - trap_locs[g, j])
                 
        for(k in 1:K) {
          p[g, i, j, k] <- p0[g, i, j, k] * exp(-1 * alpha1[Sex2[g, i]] * d[g, i, j] * d[g, i, j])
            } # k
          } # j
        } # i
    
    for(i in 1:n0[g]) {
      for(j in 1:max_trap[g]) { 
        for(k in 1:K) {
          y[i, j, k, g] ~ dbern(p[g, i, j, k])
          logit(p0[g, i, j, k]) <- alpha0[g, k] + (alpha2 * C[i, k, g])
        }
      }
    }
    
    for(i in 1:n0[g]) {
      for(j in 1:max_trap[g]) { 
        for(k in 1:K) {
          y[i, j, k, g] ~ dbern(p[g, i, j, k])
          logit(p0[g, i, j, k]) <- alpha0[g, k] + (alpha2 * C[i, k, g])
        }
      }
    }
                 
    for(i in (n0[g] + 1):M[g]) {
      for(k in 1:K) {
        for(j in 1:max_trap[g]) {
          logit(p0[g, i, j, k]) <- alpha0[g, k]
        } # j
    #  } # k
      zeros[i, g] ~ dbern(1 - prod(p[g, 1:M[g], 1:max_trap[g], k] * z[g, i]))
    } # i
}
    
    # Derived parameters
    N[g] <- sum(z[g , 1:M[g]])
    density[g] <- sum(z[g , 1:M[g]]) / (xlim[g, 2] - xlim[g, 1]) # divided distances by 100 so calculates turtles per 100 m of canal
    
   for(k in 1:K) {
     for(g in 1:n_sites) { 
        p_cap_site_day[g, k] <- mean(p[g, 1:M[g], 1:max_trap[g], k])
        }
      p_cap_day[k] <- mean(p_cap_site_day[ , k])
   }

for(t in 1:2) {
      for(g in 1:n_sites) {
        for(i in 1:M[g]) {
          p_cap_ind_site[i, g] <- mean(p0[g, i, 1:max_trap[g], 1:K])
        }
        p_cap_site[g] <- mean(p_cap_ind_site[1:M[g], g])
      }
      p_cap_sex[t] <- mean(p_cap_sex[t])
      sigma_sex[t] <- mean(sigma[t])
}

sigma_mean <- mean(sigma[ ])
    
    # for(g in 1:n_sites) {
    #   for(i in 1:M[g]) {
    # #     p_cap_site_ind[g, i] <- sum(p[g, i, 1:max_trap[g], 1:K])
    # #   }
    # #   p_cap_site[g] <- mean(p_cap_site_ind[g, 1:M[g]])
    # # }

    site_zeros[g] ~ dnorm(density[g] - (beta_0 + beta_1 * forest[g] + beta_2 * depth[g] + beta_3 * width[g]), pow(3, -2))
  } # g
  
}
", file = "Code/JAGS/zero_test.txt")

######### Set data and MCMC Conditions ########
forest <- read.csv(file = "Data/LandUse/Forest_Cover_Spatial.csv", header = FALSE)
depth <- read.csv(file = "Data/LandUse/Avg_Depth_Spatial.csv", header = TRUE)
width <- read.csv(file = "Data/LandUse/Width_Spatial.csv", header = TRUE)

forest_std <- as.numeric(scale(forest))
depth_std <- as.numeric(scale(depth))
width_std <- as.numeric(scale(width))

jags_data_site <- list(y = EM_array, 
                       Sex = sex, 
                       trap_locs = trap_locs, 
                       K = n_days, 
                       M = M, 
                       xlim = xlim, 
                       max_trap = n_traps_site$max_trap, 
                       forest = forest_std,
                       depth = depth_std,
                       width = width_std,
                       site_zeros = rep(0, n_sites),
                       C = recaptured, 
                       zeros = matrix(0, max(M), n_sites),
                       n_sites = n_sites,
                       n0 = n_ind_site$n)

initsf <- function() {
  list(s = s_st, 
       z = Z_st, 
       psi = runif(n_sites, psi_st * 0.5, psi_st*2), 
       psi_sex = runif(1, 0.3, 0.8))
}

parameters <- c("density", "N", "alpha2", "alpha0", "alpha1", "mu_0", "sd_0", "mu_1", "sd_1", "alpha_1_sex", "beta_0", "beta_1", "beta_2", "beta_3", "sigma_mean", "psi_sex", "p_cap_day", "p_cap_sex", "mu_psi", "sd_psi", "sigma_mean", "sigma_sex", "p_cap_site") ## "sigma", # "C", maybe C or a summary stat, might blow up if saving each activity center "s".

start_zeros <- Sys.time()
# cl <- makeCluster(nc)                       # Request # cores
# clusterExport(cl, c("scr_zeros", "jags_data_site", "initsf", "parameters", "EM_array", "sex", "trap_locs", "n_days", "M", "xlim", "n_traps_site", "forest_std", "depth_std", "n_sites", "recaptured", "n_ind_site", "Z_st", "s_st", "psi_st", "ni", "nb", "nt", "nc")) # Make these available
# 
# clusterSetRNGStream(cl = cl, 54354354)


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

# stopCluster(cl)
