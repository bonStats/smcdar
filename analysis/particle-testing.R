num_particles <- 10
len <- 2
rv <- matrix(rnorm(num_particles * len), nrow = num_particles, ncol = len)
prts <- particles(beta = rv)

weights(prts, log = F) <- runif(10)

sum( exp(weights(prts, log = T) ) )

w <- weights(prts, log = F)
