######### Load Libraries #########
##

library(dplyr)
library(nimble)

######### Load Data from Previous script #########

testing <- TRUE
run_date <- "2019-12-12"

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
EDF <- read.csv(file = "Data/EDF.csv", stringsAsFactors = FALSE)

#Take out sites H and I
EDF_CPIC <- EDF %>%
  filter(site != "H" & site != "I" & species == "CPIC")
str(EDF_CPIC)

# get number of unique individuals caught per site
inds <- EDF_CPIC %>%
  group_by(site) %>%
  select(site, ind) %>%
  distinct() %>%
  summarise(individuals = n())

# assume capture rate of min_cap_rate (~0.03) to get max individuals to augment (M[g])
min_cap_rate <- 0.02
df_M <- inds %>%
  mutate(M = individuals / min_cap_rate)
df_M

# restrict M to being the maximum size of the data array so the jags loops doesn't go out of bounds
M <- if_else(df_M$M > M, M, trunc(df_M$M))
# M <- if_else(df_M$M < max_ind_sp, max_ind_sp, M) # real captures can be up to row 750 at any site, all after are definitely augmentation - not anymore


######### JAGS Model with Zero Trick ##########


scr_zeros <- cat("
model{
  
  alpha2 ~ dnorm(0, pow(1.5, -3)) # Trap behavior universal distribution across sites
  mu_0 ~ dnorm(0, pow(1.5, -2))
  sd_0 ~ dunif(0, 3)
  mu_1 ~ dnorm(0, pow(1.5, -2))
  sd_1 ~ dunif(0, 3)
  alpha_1_sex ~ dnorm(0, pow(1.5, -2))
  beta_0 ~ dnorm(0, pow(3, -2))
  beta_1 ~ dnorm(0, pow(3, -2))
  beta_2 ~ dnorm(0, pow(3, -2))
  
  for(g in 1:n_sites) {
    psi[g] ~ dbeta(1, 1) # prob of individual being in the population
    psi.sex[g] ~ dbeta(1, 1) 
    alpha_1_int[g] ~ dnorm(mu_1, pow(sd_1, -2))
    
    for(t in 1:2) {
      sigma[g, t] <- pow(1 / (2*alpha1[g, t]), 0.5) # sd of half normal
      log(alpha1[g, t]) <- alpha_1_int[g] + alpha_1_sex * (t-1) # affect of being female on home range
    }
    
    for(k in 1:K) {
      alpha0[g, k] ~ dnorm(mu_0, sd_0)
    }
    
    for(i in 1:M[g]) {
      Sex[g, i] ~ dbern(psi.sex[g])
      Sex2[g, i] <- Sex[g, i] + 1
      z[g, i] ~ dbern(psi[g])
      s[g, i] ~ dunif(xlim[g, 1], xlim[g, 2]) ##??
      
      for(j in 1:max_trap[g]) { 
        d[g,i,j] <- abs(s[g, i] - trap_locs[g, j])
        
        for(k in 1:K) {
          p[g, i, j, k] <- p0[g, i, j, k] * exp(-1 * alpha1[g, Sex2[g, i]] * d[g, i, j] * d[g, i, j])
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
    
    for(i in (n0[g] + 1):M[g]) {
      zeros[i, g] ~ dbern((1 - prod(1 - p[g, i, 1:max_trap[g], 1:K])) * z[g, i])
      for(k in 1:K) {
        logit(p0[g, i, 1:max_trap[g], k]) <- alpha0[g, k]
      } # k
    } # i
    
    # Derived parameters
    N[g] <- sum(z[g , 1:M[g]])
    density[g] <- sum(z[g , 1:M[g]]) / (xlim[g, 2] - xlim[g, 1]) # divided distances by 100 so calculates turtles per 100 m of canal
    
    site_zeros[g] ~ dnorm(density[g] - (beta_0 + beta_1 * forest[g] + beta_2 * depth[g]), pow(3, -2))
  } # g
  
}
", file = "Code/JAGS/zero_test.txt")

######### Set data and MCMC Conditions ########
# Real covariate data
forest <- read.csv(file = "Data/LandUse/Forest_Cover_SingleColumn.csv", header = FALSE)
depth <- read.csv(file = "Data/LandUse/Avg_Depth_m.csv", header = TRUE)

forest_std <- as.numeric(scale(forest))
depth_std <- as.numeric(scale(depth))

jags_data_site <- list(y = EM_array, 
                       Sex = sex, 
                       trap_locs = trap_locs, 
                       K=n_days, 
                       M=M, 
                       xlim=xlim, 
                       max_trap = n_traps_site$max_trap, 
                       forest = forest_std,
                       depth = depth_std,
                       site_zeros = rep(0, n_sites),
                       C = recaptured, 
                       zeros = matrix(0, max(M), n_sites),
                       n_sites = n_sites,
                       n0 = n_ind_site$n) #, n_ind = n_ind)
# "initial values for the observed data have to be specified as NA"
initsf <- function() {
  list(s = s_st, 
       z = Z_st, 
       psi = runif(n_sites, psi_st * 0.5, psi_st*2), 
       psi.sex = runif(n_sites, 0.3, 0.8))
}

parameters <- c("sigma", "density", "N", "alpha2", "alpha0", "alpha1", "mu_0", "sd_0", "mu_1", "sd_1", "alpha_1_sex", "beta_0", "beta_1", "beta_2") #



out <- jagsUI(data = jags_data_site, 
              inits = initsf, 
              model.file = "Code/JAGS/zero_test.txt", 
              n.chains = 2,
              n.adapt = 100,
              n.iter = 101, 
              n.burnin = 1, 
              n.thin = 1, 
              parallel = TRUE, 
              codaOnly = parameters, 
              parameters.to.save = parameters)