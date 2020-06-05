
# Assume 25 m canal width

# If 95% home range = 6 ha 
sig <- 60000 / 25 / 4 / 100 # 25 m, 4 x sd for 95%, scale by 100 m length
sig

# If 50% home range = 1 ha
sig <- 10000 / 25 / 2 / 100
sig

# If linear home range = 500 m
sig <- 500 / 4 / 100
sig

# If linear home range = 250 m
sig <- 250 / 4 / 100
sig

# range of values from literature 0.6 - 6 for sigma, but most weight maybe around ~2-3. Sigma will be smaller than longer term home range estimates because of 4 day limits on movement even if they are otherwise randomly within the home range.  

# parameterized by mean (m) and standard deviation (sd)
sig_mean <- 3
sig_sd <- 1.5

sh <- sig_mean^2 / sig_sd^2
ra <- sig_mean / sig_sd^2

x <- seq(0, 10, 0.05)
plot(x = x, dgamma(x = x, shape = sh, rate = ra), type = "l")

plot(x = x, dgamma(x = x, shape = 2, rate = 0.5), type = "l")

# rate = 4 and shape = 1.3 seems reasonable based on the literature for priors for SODO and hopefully with so little weight near zero that it prevents the bimodal posterior. 

# tried 4, 1.3 and still got a small second mode between 0 & 1

pgamma(1, 4, rate = 1.3)
pgamma(1, 4, rate = 1)
pgamma(1, 4, rate = 0.75)
pgamma(1, 5, rate = 1.25)
pgamma(1, 6, rate = 1.5)

pgamma(2, 4, rate = 1.3)
pgamma(2, 4, rate = 1)
pgamma(2, 4, rate = 0.75)
pgamma(2, 5, rate = 1.25)
pgamma(2, 6, rate = 1.5)

1 - pgamma(2, 4, rate = 1) - pgamma(6, 4, rate = 1, lower.tail = FALSE) # 70% b.w 2-6
1 - pgamma(2, 5, rate = 1.25) - pgamma(6, 5, rate = 1.25, lower.tail = FALSE) # 80% bw 2-6

plot(x = x, dgamma(x = x, shape = 4, rate = 1.3), type = "l")
lines(x = x, dgamma(x = x, shape = 4, rate = 1), col = "red")
lines(x = x, dgamma(x = x, shape = 4, rate = 0.75), col = "blue")
lines(x = x, dgamma(x = x, shape = 5, rate = 1.25), col = "green")
lines(x = x, dgamma(x = x, shape = 6, rate = 1.5), col = "purple")

