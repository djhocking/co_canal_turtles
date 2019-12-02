######### Load Libraries #########
##

library(dplyr)
library(nimble)

######### Load Data from Previous script #########

testing <- TRUE
run_date <- "2019-12-01"

######### Load Data from Previous script #########

if(testing) {
  load(file = paste0("Data/Derived/all_site_testing_", run_date, ".RData"))
} else {
  load(file = "Data/Derived/all_site.RData")
}

######### NIMBLE Model ##########

scr_reg <- nimbleCode( {
    
  # # sigma0 ~ T(dnorm(0, 1 / 25), 0, 10000)
  # sigma0 ~ dunif(0, 5)
  mu_a0 ~ dnorm(0, 0.1)
  # sd_a0 ~ dt(0, pow(1.5, -2), 1)T(0, )
  sd_a0 ~ dunif(0, 100)
  # sd_a1 ~ dt(0, pow(1, -2), 1)T(0, 1000)
  sd_a1 ~ dunif(0, 100)
  mu_psi ~ dnorm(0, 0.5)
  # sd_psi ~ dt(0, pow(1.5, -2), 1)T(0, )
  sd_psi ~ dunif(0, 100)
  alpha2 ~ dnorm(0, 0.5)
  alpha3 ~ dnorm(0, 0.5)
  # beta1 ~ dnorm(0, 0.5)
  # beta2 ~ dnorm(0, 0.5)
  
  for(t in 1:2){
    mu_a1[t] ~ dnorm(0, 0.5) # fixed intercept differing by sex
  }
  
  for(g in 1:n_sites) {
    # mu_psi_site[g] ~ dnorm(mu_psi, 1 / sd_psi / sd_psi) # random intercept by site inseparable from detection intercept
    mu_psi_site[g] ~ dnorm(0, 0.5) # independent prior on augmentation intercept - fixed by site
    
    for(k in 1:K) {
      alpha0[g, k] ~ dnorm(mu_a0, 1 / sd_a0 / sd_a0)
    }
    
    for(i in 1:M[g]) {
      Sex[g, i] ~ dbern(psi.sex[g])
      Sex2[g, i] <- Sex[g, i] + 1
    }
    
    for(t in 1:2){
      sigma[g, t] <- pow(1 / (2*alpha1[g, t]), 0.5) # sd of half normal - derived parameter
      # alpha1[g, t] ~ dt(0, pow(1.5, -2), 1)T(0, )
      alpha1[g, t] ~ dunif(0, 100)
      # lalpha1[g, t] ~ dnorm(mu_a1[t], 1 / sd_a1 / sd_a1) # random intercept
      #  alpha1[g, t] <- exp(lalpha1[g, t])
    } # t
    
    # prob of individual being in the population (for augmentation since N unknown)
    logit(psi[g]) <- mu_psi_site[g] # + beta1 * depth[g] + beta2 * forest[g] 
    psi.sex[g] ~ dunif(0, 1)
    
    # for(i in 1:(max(M[1:n_sites]))) {
    for(i in 1:M[g]) {
      z[g, i] ~ dbern(psi[g])
      s[g, i] ~ dunif(xlim[g, 1], xlim[g, 2]) ##??
      
      for(j in 1:max_trap[g]) { 
        d[g,i,j] <- abs(s[g, i] - trap_locs[g, j])
        
        for(k in 1:K) {
          # for(t in 1:2) {
            logit(p0[g, i, j, k]) <- alpha0[g, k] + alpha3 * Sex[g, i] + (alpha2 * C[i, k, g])  
          # } # t
        } # k
      } # j
    } # i    
    
    for(i in 1:M[g]) {
      for (j in 1:max_trap[g]) {
        for (k in 1:K) {
          y[i, j, k, g] ~ dbern(p[g, i, j, k])
          p[g, i, j, k] <- z[g, i] * p0[g, i, j, k] * exp(- alpha1[g, Sex2[g, i]] * d[g, i, j] * d[g, i, j])
        } # i
      } # j
    } # k
    
    # Derived parameters
    N[g] <- sum(z[g , 1:M[g]])
  #   density[g] <- sum(z[g , 1:M[g]]) / (xlim[g, 2] - xlim[g, 1]) # divided distances by 100 so calculates turtles per 100 m of canal
  #   sigma_site[g] <- mean(sigma[g, 1:2])    
  } # g
  # 
  # for(k in 1:K) {
  #   for(g in 1:n_sites) { 
  #     p_cap_site_day[g, k] <- mean(p[g, 1:M[g], 1:max_trap[g], k])
  #   }
  #   p_cap_day[k] <- mean(p_cap_site_day[ , k])
  # }
  # 
  # for(t in 1:2) {
  #   for(g in 1:n_sites) { 
  #     for(i in 1:M[g]) {
  #       p_cap_ind_site_sex[i, g, t] <- mean(p0[g, i, 1:max_trap[g], 1:K, t])
  #     }
  #     p_cap_site_sex[g, t] <- mean(p_cap_ind_site_sex[1:M[g], g, t])
  #   }
  #   p_cap_sex[t] <- mean(p_cap_site_sex[1:n_sites, t])
  #   sigma_sex[t] <- mean(sigma[ , t])
  # }
  # 
  # sigma_mean <- mean(sigma[ , ])
  # 
  # for(g in 1:n_sites) {
  #   for(i in 1:M[g]) {
  #     p_cap_site_ind[g, i] <- sum(p[g, i, 1:max_trap[g], 1:K])
  #   }
  #   p_cap_site[g] <- mean(p_cap_site_ind[g, 1:M[g]])
  # }
    
  })
  
######### Set data and MCMC Conditions ########

############ TEMP #################
# Fake covariate data for testing
forest_std <- rnorm(12, 0, 2)
depth_std <- rnorm(12, 0, 2)
#############################


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
M <- if_else(df_M$M < max_ind_sp, max_ind_sp, M) # real captures can be up to row 750 at any site, all after are definitely augmentation

# in previous code add max(EM_array[ , , , ], na.rm = TRUE) for each species and save as max_ind_sp (683 for CPIC)


############ TEMP #################
# Real covariate data
forest <- read.csv(file = "Data/LandUse/Forest_Cover_SingleColumn.csv", header = FALSE)
depth <- read.csv(file = "Data/LandUse/Avg_Depth_m.csv", header = TRUE)

forest_std <- as.numeric(scale(forest))
depth_std <- as.numeric(scale(depth))

jags_data_site <- list(y = EM_array, 
                       Sex = Sex, 
                       trap_locs = trap_locs, 
                       K=K, 
                       M=M, 
                       xlim=xlim, 
                       max_trap = max_trap, 
                       # forest = forest_std,
                       # depth = depth_std,
                       C = C, 
                       n_sites = G) #, n_ind = n_ind)
# "initial values for the observed data have to be specified as NA"
initsf <- function() {
  list(# alpha0 = rnorm(n_sites, -2, 0.5), 
    # alpha1 = matrix(abs(rnorm(n_sites * 2, 1, 2)), n_sites, 2),
    # alpha2 = matrix(rnorm(n_sites * M, 1, 2), n_sites, M),
    # alpha2 = rnorm(1, 1, 1)),
    # s = t(sst), 
    # z = z[ , 1:max(M)]) #, 
    mu_psi_site = rep(exp(3)/(3 + exp(3)), n_sites),
    mu_a0 = rnorm(1, -5, 1),
    sd_a0 = runif(1, 1, 5)) # , 
  # psi.sex = runif(n_sites)) #, Sex = c(rep(NA, n_ind))) ## Error = "Invalid parameters for chain 1: non-numeric intial values supplied for variable(s) Sex"   #### ALPHA2????
}


parameters <- c("sigma", "N", "alpha2", "alpha0", "alpha1", "sigma") #, "beta1", "beta2") # "C", maybe C or a summary stat, might blow up if saving each activity center "s". "N", 

samples <- nimbleMCMC(
  code = scr_reg,
  constants = jags_data_site, ## provide the combined data & constants as constants
  inits = initsf,
  monitors = parameters,
  niter = 300,
  nburnin = 200,
  thin = 1)

inits1 <- initsf()
Rmodel <- nimbleModel(code = scr_reg, name = 'scr_reg', inits = inits1) # , constants = pumpConsts, data = pumpData,
Rmcmc <- buildMCMC(Rmodel)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(100)
samples <- as.matrix(Cmcmc$mvSamples)

samples <- runMCMC(CpumpMCMC, niter = niter)



if(run_par) {
  cl <- makeCluster(nc)                        # Request # cores
  clusterExport(cl, c("jags_data_site", "inits", "parameters", "z", "sst", "Sex", "ni", "na", "nt", "K", "C", "M", "G", "n_sites")) # Make these available
  clusterSetRNGStream(cl = cl, 54354354)
  
  system.time({ # no status bar (% complete) when run in parallel
    out <- clusterEvalQ(cl, {
      library(rjags)
      jm <- jags.model("Code/JAGS/scr_test_varM.txt", jags_data_site, inits = inits, n.adapt = na, n.chains = 1) # Compile model and run burnin - consider settingburnin vs adapt
      out <- coda.samples(jm, parameters, n.iter = ni, thin = nt) # Sample from posterior distribution
      return(as.mcmc(out))
    })
  }) #
  
  stopCluster(cl)
  
  samples <- mcmc.list(out)
}
}




