########## Load Libraries ##########
library(coda)
library(rjags)
# library(devtools)
# install_github("https://github.com/stan-dev/bayesplot" # if want latest development version of bayesplot)
library(bayesplot)

######### Load MCMC Object #########

out <- readRDS("Results/JAGS/all_site_reg_final.rds")

samples <- out$samples

# data <- bayesplot:::mcmc_trace_data(samples, regex_pars = c("density", "N", "alpha2", "alpha0", "alpha1", "mu_0", "sd_0", "mu_1", "sd_1", "beta_0", "beta_1", "beta_2", "beta_3", "psi_sex", "p_cap_day", "mu_psi", "sd_psi", "sigma_mean", "sigma", "p_cap_site"))


######### Check MCMC ########

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(samples, regex_pars = c("mu"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("psi.sex"))
p + facet_text(size = 15)

# p <- mcmc_trace(samples, regex_pars = c("alpha1"))
# p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("alpha2"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("alpha3"))
p + facet_text(size = 15)

# p <- mcmc_trace(samples, regex_pars = c("sigma"))
# p + facet_text(size = 15)

# p <- mcmc_dens(samples, regex_pars = c("sigma"))
# p + facet_text(size = 15) + lim(0, 10)

p <- mcmc_trace(samples, regex_pars = c("beta"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("sd"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("density"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("N"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("mu_psi_site"))
p + facet_text(size = 15)

# p <- mcmc_trace(samples, regex_pars = c("p_cap"))
# p + facet_text(size = 15)

effectiveSize(samples)

gelman.diag(samples)

gelman.diag(samples[ , c("density[1]", "density[2]", "density[3]", "density[4]", "density[5]", "density[6]", "density[7]", "density[8]", "density[9]", "density[10]", "density[11]", "density[12]")])

# Detection intercept 
# p <- mcmc_trace(samples, regex_pars = c("alpha0"))
# p + facet_text(size = 15)

# p <- mcmc_trace(samples, regex_pars = c("alpha3"))
# p + facet_text(size = 15)

# p <- mcmc_trace(samples, regex_pars = c("p_cap"))
# p + facet_text(size = 15)

library(rstan)
sample_array <- as.array(samples)
sample_array <- aperm(sample_array, c(1, 3, 2))
mon <- monitor(sample_array, warmup = 0)
rstan:::print.simsummary(mon)
max(mon$Rhat)
min(mon$Bulk_ESS)
min(mon$Tail_ESS)

