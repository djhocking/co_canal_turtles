########## Load Libraries ##########
library(coda)
library(rjags)
# library(nimble)
# library(devtools)
# install_github("https://github.com/stan-dev/bayesplot" # if want latest development version of bayesplot)
library(ggplot2)
library(bayesplot)
library(dplyr)
# library(tibble)

######### Load MCMC Object #########

# out <- readRDS("Results/JAGS/all_sites_reg_2020-05-01.rds")
Species <- "PRUB"
run_date <- "2020-05-21"
out <- readRDS(paste0("Results/JAGS/", Species, "_all_sites_reg_", run_date, ".rds"))

samples <- out$samples

rm(out)

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

p <- mcmc_trace(samples, regex_pars = c("sigma"))
p + facet_text(size = 15)

# 50% kernal density home range (+/- 0.68 SD)

# mat <- t(as.matrix(samples))
# df <- data.frame(mat, stringsAsFactors = FALSE)
# names(df) <- dimnames(mat)[[1]]
# df <- df %>%
#   # dplyr::select(sigma) %>%
#   dplyr::mutate(hr50 = 2 * sigma * 0.68,
#          hr95 = 2 * sigma * 2)

# ggplot(data = df, aes(hr50)) + geom_histogram()
# p <- mcmc_dens(samples, regex_pars = c("home_50")) # + panel_cols(color = "gray20", fill = "gray30")
# p + facet_text(size = 15)

p <- mcmc_dens(samples, regex_pars = c("home_"))
p + facet_text(size = 15)

p <- mcmc_dens_chains(samples, regex_pars = c("home_"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("alpha1"))
p + facet_text(size = 15)


p <- mcmc_trace(samples, regex_pars = c("alpha2"))
p + facet_text(size = 15)

p <- mcmc_dens(samples, regex_pars = c("alpha1"))
p + facet_text(size = 15) #+ lim(0, 10)

p <- mcmc_trace(samples, regex_pars = c("beta"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("sd"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("density"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("density_ha"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("N"))
p + facet_text(size = 15)

p <- mcmc_hist(samples, regex_pars = c("N"))
p + facet_text(size = 15)
# p + geom_vline(data = n_ind_site, aes(xintercept = n), color = "red") + facet_wrap(~site_num) # doesn't work
# 
# EM_array[ , , , 1]
# EM_array[ , , 3, 1]
# EM_array[ , , 4, 1]

# p <- mcmc_trace(samples, regex_pars = c("mu_psi_site"))
# p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("p_cap_site"))
p + facet_text(size = 15)

p <- mcmc_trace(samples, regex_pars = c("p_cap_day"))
p + facet_text(size = 15)

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
nc <- 8
sample_array <- as.array(samples)
sample_array <- aperm(sample_array, c(1, 3, 2))
mon <- rstan::monitor(sample_array, warmup = 0)
rstan:::print.simsummary(mon)
max(mon$Rhat)
min(mon$Bulk_ESS)
min(mon$Tail_ESS)

# for nimble
sample_array <- base::array(as.numeric(unlist(samples)), dim=c(nrow(samples[[1]]), length(samples), ncol(samples[[1]])), dimnames = list(NULL, NULL, colnames(samples[[1]])))
mon <- rstan::monitor(sample_array, warmup = 0)
rstan:::print.simsummary(mon)
max(mon$Rhat)
min(mon$Bulk_ESS)
min(mon$Tail_ESS)

