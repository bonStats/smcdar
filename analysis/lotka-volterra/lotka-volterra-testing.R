## matlab
dat <- read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test1/data.csv", header = F, col.names = c("y1","y2","t"))

gen1 <- read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test1/generator.csv", header = F, col.names = c("i","j","x"))
theta1 <-  as.matrix(read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test1/theta.csv",header = F))
max_vals1 <-  as.matrix(read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test1/max_vals.csv",header = F))
ll1_m <-  as.matrix(read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test1/log_like.csv",header = F))

gen2 <- read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test2/generator.csv",header = F, col.names = c("i","j","x"))
theta2 <-  as.matrix(read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test2/theta.csv",header = F))
max_vals2 <-  as.matrix(read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test2/max_vals.csv",header = F))
ll2_m <-  as.matrix(read.csv("~/sdrive/Research/smc-da/lotka_volterra_matlab/test2/log_like.csv",header = F))


# check
G1_r <- generator_matrix_lotka_volterra_ctmc(theta1[1,], y1_max = max_vals1[1], y2_max = max_vals1[2])

G1_m <- Matrix::sparseMatrix(i = gen1$i,
                     j =  gen1$j,
                     x = gen1$x,
                     giveCsparse = FALSE)

sum( abs(G1_r - G1_m) > 1e-8)

G2_r <- generator_matrix_lotka_volterra_ctmc(theta2[1,], y1_max = max_vals2[1], y2_max = max_vals2[2])

G2_m <- Matrix::sparseMatrix(i = gen2$i,
                             j =  gen2$j,
                             x = gen2$x,
                             giveCsparse = FALSE)

sum( abs(G2_r - G2_m) > 1e-8)


ll1_r1 <- log_lhood_lotka_volterra_ctmc(theta1, y1 = dat$y1, y2 = dat$y2, times = dat$t,   y1_max = max_vals1[1], y2_max = max_vals1[2], package = "expoRkit")
ll1_r2 <- log_lhood_lotka_volterra_ctmc(theta1, y1 = dat$y1, y2 = dat$y2, times = dat$t,   y1_max = max_vals1[1], y2_max = max_vals1[2], package = "expm")

ll1_r - ll1_m

v <- rep(0, nrow(G1_r))
v[123] <- 1
e1 <- expm::expAtv(G1_r, v = v, t = 1)
e2 <- expoRkit::expv(G1_r, v, 1)

system.time(
replicate(100, expm::expAtv(G1_r, v = v, t = 1))
)

system.time(
  replicate(100, expoRkit::expv(G1_r, v, 1))
)

## change to unit testing

library(Matrix)
library(expm)
theta <- runif(3)
y1_min <- 0
y2_min <- 0
y1_max <- 10
y2_max <- 10

G <- generator_matrix_lotka_volterra_ctmc(theta = theta, y1_min = y1_min, y1_max = y1_max, y2_min = y2_min, y2_max = y2_max)

# check probs
f_init_s <- 11
f_final_s <- 10
time <- 2
P <- expm(time * G, method = "R_Eigen")

rowSums(P)

P[1:3,1:3]

v <- rep(0, times = nrow(P))
v[f_final_s] <- 1

exp_tG_v1 <- expAtv(A = G, t = time, v = v)
exp_tG_v2 <- P %*% v

max(abs(exp_tG_v1$eAtv - exp_tG_v2))

smcdar:::prob_lotka_volterra_ctmc(flat_init_state = f_init_s, flat_final_state = f_final_s,elapsed_time = time, generator_matrix = G)

exp_tG_v1$eAtv[f_init_s]

P[f_init_s,f_final_s]


#### LL ####

y1 <- c(4,5,6,2)
y2 <- c(3,4,6,1)
times <- 1:4

log_lhood_lotka_volterra_ctmc(theta = theta, y1 = y1, y2 = y2, times = times, y1_max = y1_max, y2_max = y2_max)
