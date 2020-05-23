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

nowt <- function(x = NULL) x

######### Load MCMC Object #########

# out <- readRDS("Results/JAGS/all_sites_reg_2020-05-01.rds")
Species <- "CSER"
run_date <- "2020-05-21"
out <- readRDS(paste0("Results/JAGS/", Species, "_all_sites_reg_", run_date, ".rds"))

samples <- out$samples
rm(out)

#----- Check that N >= num unique caps -----
samples_df<- samples %>%
  as.matrix(.) %>%
  data.frame(.)
# names(samples_df) <- colnames

samples_df %>%
  dplyr::select(starts_with("N.")) %>%
  dplyr::summarise_all(list(min = min, 
                            # q25 = quantile(., 0.25), 
                            median = median, 
                            # q75 = quantile(., 0.75), 
                            max = max,
                            mean = mean, 
                            sd = sd))
  
samples_N <- samples_df %>%
  dplyr::select(starts_with("N."))
summary_N <- samples_N %>%
  purrr::map_df(.f = ~ broom::tidy(summary(.x)), .id = "variable") %>%
  # purrr::map_df(.f = ~ quantile(.x, probs = c(0, 0.05, 0.5, 0.95, 1)), .id = "variable") %>%
  dplyr::mutate(n = n_ind_site$n,
                check = minimum >= n) %>%
  nowt()

#----- Get proportion of total individuals captured over all sample periods and traps -----
samples_prop_cap <- samples_N
samples_prop_cap[ , ] <- NA_real_
for(i in 1:ncol(samples_N)) {
  samples_prop_cap[ , i] <- n_ind_site$n[i] / samples_N[ , i]
}
summary(samples_prop_cap)

samples_cap_long <- samples_prop_cap %>%
  tidyr::pivot_longer(cols = starts_with("N"), names_to = "site", values_to = "p")

ggplot(samples_cap_long, aes(p)) + geom_histogram(bins = 50, aes(y = ..density..)) + coord_cartesian(ylim = c(0, 20), xlim = c(0, 0.5)) + facet_wrap(~site)

# Change Order of Facets
# temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))
# or
# arrange(size_num) %>% # sort
#   mutate_at(vars(size), funs(factor(., levels=unique(.)))) %>% # convert to factor


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
# p + vline_at(n_ind_site$n)
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

