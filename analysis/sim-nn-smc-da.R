.libPaths(Sys.getenv("R_LIB_USER"))
outdir <- Sys.getenv("OUT_DIRECTORY")

devtools::install_github("bonStats/smcdar")

#### load packages ####

library(smcdar)
library(dplyr)
#library(tidyr)
#library(ggplot2)
library(future)
#library(future.apply)
library(furrr)

#### Settings ####

## Priors, likelihoods

# prior: beta ~ Normal(0,1)
log_prior <- function(beta){

  sum( dnorm(beta, mean = rep(0, times = 5), sd = 5, log = T) )

}

draw_prior <- function(n){
  matrix(
    rnorm(n * 5, mean = rep(0, times = 5), sd = 5), ncol = 5, byrow = T
  )
}

log_like <- function(beta, X, y){
  Sys.sleep(0.1)
  sum( dnorm(y, mean = X %*% beta, sd = 1, log = T) )

}

log_like_approx <- function(beta, X, y, bias){

  sum( dnorm(y, mean = X %*% (beta + bias), sd = 1.5, log = T) )

}

nn_optimise_pre_approx_llhood_transformation <- function(particles, b_s_start, loglike, temp){

  #log_theta is matrix
  # should only send values of log_like that have been cached already
  true_ll <- apply(particles, MARGIN = 1, loglike, type = "full_likelihood", temp = temp)

  len <- ncol(particles)

  if(missing(b_s_start)){

    b_s_start <- rep(0, times = len)

  }

  topt <- function(b_s){

    b_mat <- matrix(b_s[1:len], ncol = len, nrow = nrow(particles), byrow = T)
    s_mat <- matrix(b_s[(len+1):(2*len)], ncol = len, nrow = nrow(particles), byrow = T)

    approx_ll <- apply((particles - b_mat)/exp(s_mat), MARGIN = 1, loglike, type = "approx_likelihood", temp = temp)

    keep_ll <- is.finite(approx_ll)

    # KL discrepancy
    sum( # min{ p * (log ~p - log ~q) }
      weights(particles)[keep_ll] * ( true_ll[keep_ll] - approx_ll[keep_ll] )
    )

  }

  optim(par = b_s_start, fn = topt, control = list(maxit = 100), method = "BFGS")

}

nn_posterior <- function(y, X, sigma, tau){
  p <- ncol(X)
  n <- length(y)
  XtX <- crossprod(X)
  Xty <- crossprod(X,y)
  mu <- solve(XtX + diag(1/(tau^2), nrow = p), Xty)
  V <- solve(XtX + diag(1/(tau^2), nrow = p)) * (sigma^2)

  log_z_lhood <- -(n/2) * ( log(2) + log(pi) + 2 * log(sigma) )
  log_z_prior <- -(p/2) * ( log(2) + log(pi) + 2 * log(tau) )
  log_z_post <- -(p/2) * ( (log(2) + log(pi)) + log(det(V)) )

  return(list(
    mu = mu,
    var = V,
    log_z =  log_z_prior + log_z_lhood - log_z_post
  ))

}

# true_theta = (birth rate prey, death rate prey/ birth rate pred, death rate pred)
g_pars <- list(N = 100, N_approx = 100, true_beta = c(0, 0.5, -1.5, 1.5, -3),
               log_prior = log_prior, log_like = log_like, log_like_approx = log_like_approx,
               draw_prior = draw_prior, optimise_pre_approx_llhood_transformation = nn_optimise_pre_approx_llhood_transformation
)

sim_settings <- rep(list(list(f_pars = NULL, g_pars = g_pars)), 2)

sim_settings[[1]]$f_pars <- list(
  num_p = 100,
  step_scale_set = c(0.05, 0.1, 0.2, 0.4, 0.6),
  b_s_start = rep(0, 10),
  refresh_ejd_threshold = 1,
  approx_ll_bias = 1
)

sim_settings[[2]]$f_pars <- list(
  num_p = 200,
  step_scale_set =  c(0.05, 0.1, 0.2, 0.4, 0.6),
  b_s_start = rep(0, 10),
  refresh_ejd_threshold = 1,
  approx_ll_bias = 1
)

simulate_regr <- function(N, beta){

  X <- matrix(rnorm(length(beta)*N), nrow = N)
  y <- X %*% beta + rnorm(N)

  return(
    list(y = y, X = X, beta = beta)
  )

}

####

#### Code for sim ####

run_sim <- function(ss, verbose = F){

  ## simulate data

  sim <- simulate_regr(N = ss$g_pars$N, beta = ss$g_pars$true_beta)

  source("analysis/nn-smc-da.R")

  log_like <- function(beta) ss$g_pars$log_like(beta, X = sim$X, y = sim$y)
  log_like_approx <- function(beta) ss$g_pars$log_like_approx(beta, X = sim$X[1:ss$g_pars$N_approx,],
                                                              y = sim$y[1:ss$g_pars$N_approx],
                                                              bias = ss$f_pars$approx_ll_bias)

  ## run
  cat("smc da")
  smc_da <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    b_s_start = b_s_start,
    refresh_ejd_threshold = refresh_ejd_threshold,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like,
    log_like_approx = log_like_approx,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  ss$g_pars$optimise_pre_approx_llhood_transformation,
    verbose = verbose
  )
  )
  cat("smc full")
  smc_full <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    b_s_start = b_s_start,
    refresh_ejd_threshold = refresh_ejd_threshold,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like,
    log_like_approx = log_like_approx,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  ss$g_pars$optimise_pre_approx_llhood_transformation,
    verbose = verbose
  )
   )

  cat("smc approx")
  smc_approx <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = T,
    num_p = num_p,
    step_scale_set = step_scale_set,
    b_s_start = b_s_start,
    refresh_ejd_threshold = refresh_ejd_threshold,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like,
    log_like_approx = log_like_approx,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  ss$g_pars$optimise_pre_approx_llhood_transformation,
    verbose = verbose
  )
  )

  return(
    list(sim_settings = ss,
         sim_data = sim,
         smc_da = smc_da,
         smc_full = smc_full,
         smc_approx = smc_approx
    )
  )

}

####

res <- run_sim(ss = sim_settings[[1]], verbose = T)

#res <- future_map(1:100, .f = function(x) , .progress = T)
