# co_canal_turtles

## Workflow

1. Set the species and testing in `Code/config.yml`
2. Run `Code/prep_data_species.R` to organize data, augment, set starting values
3. Run `Code/run_jags.R` to implement JAGS in parallel on 8 cores
4. Use `Code/mcmc_check.R` to visualize posteriors and check ESS and $\hat{R}$
5. The `mcmc_Analysis.Rmd` file can be used to create figures and tables (refine & turn into supplemental materials?)