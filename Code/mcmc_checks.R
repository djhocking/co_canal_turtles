########## Load Libraries ##########
library(coda)
library(rjags)
# library(devtools)
# install_github("https://github.com/stan-dev/bayesplot" # if want latest development version of bayesplot)
library(bayesplot)

######### Load MCMC Object #########

samples <- readRDS("Results/JAGS/all_site_reg.rds")

######### Check MCMC ########

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(samples, regex_pars = c("mu"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("alpha1"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("alpha2"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("alpha3"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("sigma"))
p + facet_text(size = 15)

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

p <- mcmc_trace(samples, regex_pars = c("p_cap"))
p + facet_text(size = 15)

# Detection intercept 
p <- mcmc_trace(samples, regex_pars = c("alpha0"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("alpha3"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("p_cap"))
p + facet_text(size = 15)


