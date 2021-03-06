######### Load Libraries #########
##

library(dplyr)
library(nimble)

######### Load Data from Previous script #########

testing <- TRUE
run_date <- "2019-12-19"

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


######### NIMBLE Model ##########


scr_zeros <- nimbleCode( {
  
  alpha2 ~ dnorm(0, pow(1.5, -2)) # Trap behavior universal distribution across sites
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
      s[g, i] ~ dunif(xlim[g, 1], xlim[g, 2])
      
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
      for(k in 1:K) {
        for(j in 1:max_trap[g]) {
          logit(p0[g, i, j, k]) <- alpha0[g, k]
        } # j
      } # k
      zeros[i, g] ~ dbern(1 - prod(p[i, g, 1:max_trap[g], 1:K] * z[g, i]))
    } # i
    
    # Derived parameters
    N[g] <- sum(z[g , 1:M[g]])
    density[g] <- sum(z[g , 1:M[g]]) / (xlim[g, 2] - xlim[g, 1]) # divided distances by 100 so calculates turtles per 100 m of canal
    
    site_zeros[g] ~ dnorm(density[g] - (beta_0 + beta_1 * forest[g] + beta_2 * depth[g]), pow(3, -2))
  } # g
  
})


scr_reg <- nimbleCode( {
  
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
    psi[g] ~ dunif(0, 1) # prob of individual being in the population
    psi.sex[g] ~ dunif(0, 1) 
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
          logit(p0[g, i, j, k]) <- alpha0[g, k] + alpha2 * C[i, k, g]
          y[i, j, k, g] ~ dbern(p[g, i, j, k])
          p[g, i, j, k] <- z[g, i] * p0[g, i, j, k] * exp(-1 * alpha1[g, Sex2[g, i]] * d[g, i, j] * d[g, i, j])
        } # i
      } # j
    } # k
    
    # Derived parameters
    N[g] <- sum(z[g , ])
    density[g] <- sum(z[g , ]) / (xlim[g, 2] - xlim[g, 1]) # divided distances by 100 so calculates turtles per 100 m of canal
    
    site_zeros[g] ~ dnorm(density[g] - (beta_0 + beta_1 * forest[g] + beta_2 * depth[g]), pow(3, -2))
  } # g
  
})

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

#-------- test regular model ------------

# need to expand recaptures to include zeros for augmented individuals
recaptured2 <- base::array(0, dim = c(max(M), n_days, n_sites))
for(i in 1:dim(recaptured)[1]) {
  for(j in 1:dim(recaptured)[2]) {
    for(k in 1:dim(recaptured)[3]) {
      recaptured2[i,j,k] <- recaptured[i,j,k]
    }
  }
}

# need to expand observations to include zeros for augmented individuals
y2 <- base::array(NA, dim = c(max(M), max(n_traps_site$max_traps), n_days, n_sites))
for(i in 1:dim(EM_array)[1]) {
  for(j in 1:dim(EM_array)[2]) {
    for(k in 1:dim(EM_array)[3]) {
      for(g in 1:dim(EM_array)[4]) {
        y2[i,j,k,g] <- EM_array[i,j,k,g]
      }
    }
  }
}

jags_data_site <- list(y = y2, 
                       Sex = sex, 
                       trap_locs = trap_locs, 
                       K=n_days, 
                       M=M, 
                       xlim=xlim, 
                       max_trap = n_traps_site$max_trap, 
                       forest = forest_std,
                       depth = depth_std,
                       site_zeros = rep(0, n_sites),
                       C = recaptured2, 
                       zeros = matrix(0, max(M), n_sites),
                       n_sites = n_sites,
                       n0 = n_ind_site$n)
start_reg <- Sys.time()
samples <- nimbleMCMC(
  code = scr_reg,
  constants = jags_data_site, ## provide the combined data & constants as constants
  inits = initsf,
  monitors = parameters,
  niter = 200,
  nburnin = 100,
  thin = 1,
  nchains = 1)
end_reg <- Sys.time()

#--------- test model with zeros trick ------------
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
                       n0 = n_ind_site$n)
start_zeros <- Sys.time()
samples <- nimbleMCMC(
  code = scr_zeros,
  constants = jags_data_site, ## provide the combined data & constants as constants
  inits = initsf,
  monitors = parameters,
  niter = 200,
  nburnin = 100,
  thin = 1,
  nchains = 1)
end_zeros <- Sys.time()

(time_reg <- end_reg - start_reg)
(time_zeros <- end_zeros - start_zeros)

as.numeric(time_reg - time_zeros) / as.numeric(time_reg)



plot(samples[ , "alpha2"], type = "l", xlab = "iteration",
     ylab = expression(alpha))

plot(samples[ , "density[1]"], type = "l", xlab = "iteration",
     ylab = expression(alpha))

plot(samples[ , "density[3]"], type = "l", xlab = "iteration",
     ylab = expression(alpha))

plot(samples[ , "beta_1"], type = "l", xlab = "iteration",
     ylab = expression(alpha))
plot(samples[ , "beta_2"], type = "l", xlab = "iteration",
     ylab = expression(alpha))

plot(samples[ , "sigma[1, 1]"], type = "l", xlab = "iteration",
     ylab = expression(alpha))

if(run_par) {
  
  nc <- 3
  ni <- 200
  na <- 100
  nt <- 1
  
  cl <- makeCluster(nc)                        # Request # cores
  clusterExport(cl, c("scr_reg", "jags_data_site", "initsf", "parameters", "z", "sst", "Sex", "ni", "na", "nt", "K", "C", "M", "G", "n_sites")) # Make these available
  clusterSetRNGStream(cl = cl, 54354354)
  
  system.time({ 
    out <- clusterEvalQ(cl, {
      library(nimble)
      samples <- nimbleMCMC(
        code = scr_reg,
        constants = jags_data_site, ## provide the combined data & constants as constants
        inits = initsf,
        monitors = parameters,
        niter = 200,
        nburnin = 100,
        thin = 1,
        nchains = 1)
      return(samples)
    })
  }) #
  
  stopCluster(cl)
}

if(!dir.exists("Results/nimble")) dir.create("Results/nimble", recursive = TRUE)
saveRDS(out, "Results/nimble/all_site_test_", run_date, ".rds")
saveRDS(out, paste0("Results/nimble/all_site_final_", run_date, ".rds"))

library(bayesplot)
post_mat <- as.matrix(out)

color_scheme_set("mix-blue-pink")

p <- mcmc_trace(post_mat,  regex_pars = c("N"), n_warmup = 10) #, facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)
p <- mcmc_trace(post_mat,  regex_pars = c("sigma"), n_warmup = 10) #, facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)
p <- mcmc_trace(post_mat,  regex_pars = c("alpha_1_sex"), n_warmup = 10) #, facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)

p <- mcmc_trace(post_mat,  regex_pars = c("alpha2"), n_warmup = 10) #, facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)




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





