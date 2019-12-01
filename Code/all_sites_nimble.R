######### Load Libraries #########
##

library(dplyr)
library(nimble)

######### Load Data from Previous script #########

testing <- TRUE
run_date <- "2019-11-30"

######### Load Data from Previous script #########

if(testing) {
  load(file = paste0("Data/Derived/all_site_testing_", run_date, ".RData"))
} else {
  load(file = "Data/Derived/all_site.RData")
}

######### NIMBLE Model ##########

scr_reg <- nimbleCode( {
    
    mu0 ~ dnorm(0, 0.01)
   # sigma0 ~ dt(0, pow(5, -2), 1)T(0, )
  # sigma0 ~ T(dnorm(0, 1 / 25), 0, 10000)
  sigma0 ~ dunif(0, 5)
    
    # mu_a1 ~ dnorm(0, 1 / (5^2))I(0, ) ## half normal
    # sd_a1 ~ dunif(0, 5)
    mu_a2 ~ dnorm(0, 0.01)
    # sd_a2 ~ dt(0, pow(5, -2), 1)T(0, ) # half cauchy prior with scale = 5 (25?)
    sd_a2 ~ dunif(0, 5)
    
    beta1 ~ dnorm(0, 0.1)
    # beta2 ~ dnrom(0, 0.1)
    
    for(t in 1:2){
      alpha1_mean[t] ~ dnorm(0, 1 / 9) # fixed intercept differing by sex
    }
    
    for(g in 1:n_sites) {
      
      for(i in 1:M) {
        Sex[g, i] ~ dbern(psi.sex[g])
        Sex2[g, i] <- Sex[g, i] + 1
      }
      
      for(t in 1:2){
        # alpha1_mean[g, t] ~ dnorm(0, 1 / (5^2))T(0, ) ## half normal independent across sites and sexes
        sigma[g, t] <- pow(1 / (2*alpha1[g, t]), 0.5) # sd of half normal - derived parameter
        
        # log_alpha1 <- beta1
        # alpha1[g, t] <- alpha1_mean[g, t] + exp(log_alpha1[g, t])
        log(alpha1[g, t]) <- alpha1_mean[t] + beta1 * depth[g] # + beta2 * forest[g] # linear regression on home range hence density? - maybe should put regression on psi instead?
        # alpha1_mean[g, t] ~ dnorm(mu_a1, 1 / sd_a1 / sd_a1) # random intercept
        
        
      } # t
      
      
      for(i in 1:M) {
        alpha2[g, i] ~ dnorm(mu_a2, 1 / sd_a2 / sd_a2) # Trap behavior universal distribution across sites
      } # m
      
      psi[g] ~ dunif(0, 1) # prob of individual being in the population (for augmentation since N unknown)
      # logit(psi[g]) ~ dnorm(mu_psi, 1 / sd_psi / sd_psi) # consider drawing from a normal distribution across sites
      # mu_psi ~ dnorm(0, 0.01)
      # sd_psi ~ dunif(0, 10)
      psi.sex[g] ~ dunif(0, 1)
      
      
      alpha0[g] ~ dnorm(mu0, sigma0) 
      
      for(i in 1:M) {
        z[g, i] ~ dbern(psi[g])
        s[g, i] ~ dunif(xlim[g, 1], xlim[g, 2]) ##??
        
        for(j in 1:max_trap[g]) { 
          d[g,i,j] <- abs(s[g, i] - trap_locs[g, j])
          
          for(k in 1:K) {
            for(t in 1:2) {
              logit(p0[g, i, j, k, t]) <- alpha0[g] + (alpha2[g, i] * C[i, k, g])  # alpha2*C to rep. global behav. response
            } # t
          } # k
        } # j
      } # i    
      
      for(i in 1:M) {
        for (j in 1:max_trap[g]) {
          for (k in 1:K) {
            y[i, j, k, g] ~ dbern(p[g, i, j, k])
            p[g, i, j, k] <- z[g, i] * p0[g, i, j, k, Sex2[g, i]] * exp(- alpha1[g, Sex2[g, i]] * d[g, i,j] * d[g, i,j])
          } # i
        } # j
      } # k
      
      # Derived parameters
      # N[g] <- sum(z[g , ])
      # N[g] <- inprod(z[1:M_allsites], sitedummy[ , t]) ## see panel 9.2
      density[g] <- sum(z[g , ]) / (xlim[g, 2] - xlim[g, 1]) # divided distances by 100 so calculates turtles per 100 m of canal
      
    } # g
    
  })
  
######### Set data and MCMC Conditions ########

############ TEMP #################
# Fake covariate data for testing
forest_std <- rnorm(12, 0, 2)
depth_std <- rnorm(12, 0, 2)
#############################


jags_data_site <- list(y = EM_array, 
                       Sex = Sex, 
                       trap_locs = trap_locs, 
                       K=K, 
                       M=M, 
                       xlim=xlim, 
                       max_trap = max_trap, 
                       forest = forest_std,
                       depth = depth_std,
                       C = C, 
                       n_sites = G) #, n_ind = n_ind)
# "initial values for the observed data have to be specified as NA"
initsf <- function() {
  list(alpha0 = rnorm(n_sites, -2, 0.5), 
       alpha1 = matrix(abs(rnorm(n_sites * 2, 1, 2)), n_sites, 2),
       alpha2 = matrix(rnorm(n_sites * 2, 1, 2), n_sites, M),
       s = t(sst), 
       z = z, 
       psi = runif(n_sites), 
       psi.sex = runif(n_sites)) #, Sex = c(rep(NA, n_ind))) ## Error = "Invalid parameters for chain 1: non-numeric intial values supplied for variable(s) Sex"   #### ALPHA2????
}

parameters <- c("sigma", "density", "sigma_ind", "alpha2", "alpha0", "alpha1", "sigma", "beta1", "beta2") # "C", maybe C or a summary stat, might blow up if saving each activity center "s". "N", 

samples <- nimbleMCMC(
  code = scr_reg,
  constants = jags_data_site, ## provide the combined data & constants as constants
  inits = inits,
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




