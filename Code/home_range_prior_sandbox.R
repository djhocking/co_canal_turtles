library(magrittr)

set.seed(2439761)
alpha1 <- rnorm(1000000, 1.5, 2) %>%
  .[. >= 0.02] # truncates 95% home rante at 1km if 0.08, 0.02 truncates at 2km

# alpha1 <- runif(10000, 0.01, 2)

sigma <- (1 / (2 * alpha1)) ^ 0.5

home_range <- 4 * sigma * 100


summary(alpha1)
# plot(density(alpha1))
summary(sigma)
summary(home_range)
sd(home_range)
# hist(home_range)
plot(density(home_range), xlim = c(0, 1000)) # prior on home range


# if typical home range size from literature is ~2 ha, then in a 20 m width section of the canal, linear home range would be 1000 m (seems large). But less than 100m moved per day, so we wouldn't expect more than 400 m for ours?