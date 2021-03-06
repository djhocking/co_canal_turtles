
# sex ratios vary by site
scr_zeros <- nimbleCode( {
  
  alpha2 ~ dnorm(0, pow(1.5, -2)) # Trap behavior universal distribution across sites
  mu_0 ~ dnorm(0, pow(1.5, -2))
  sd_0 ~ dunif(0, 3)
  # alpha_1_sex ~ dnorm(0, pow(1.5, -2))
  beta_0 ~ dnorm(0, pow(10, -2))
  beta_1 ~ dnorm(0, pow(10, -2))
  beta_2 ~ dnorm(0, pow(10, -2))
  # alpha_1_int ~ dnorm(0, pow(3, -2))
  
  for(g in 1:n_sites) {
    psi[g] ~ dbeta(1, 1) # prob of individual being in the population
    psi.sex[g] ~ dbeta(1, 1) 
    
    for(t in 1:2) {
      sigma[t] <- pow(1 / (2 * alpha1[t]), 0.5) # sd of half normal
      # log(alpha1[t]) <- alpha_1_int + alpha_1_sex * (t-1) # affect of being female on home range
      alpha1[t] ~ dunif(0, 20)
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