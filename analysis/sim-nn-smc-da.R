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
library(mvtnorm)

#### Settings ####

## Priors, likelihoods

# prior: beta ~ Normal(0,1)
log_prior <- function(beta){

  sum( dnorm(beta, mean = rep(0, times = 5), sd = 2, log = T) )

}

draw_prior <- function(n){
  matrix(
    rnorm(n * 5, mean = rep(0, times = 5), sd = 2), ncol = 5, byrow = T
  )
}

log_like <- function(beta, X, y){
  Sys.sleep(0.1)
  sum( dnorm(y, mean = X %*% beta, sd = 1, log = T) )

}

log_like_approx <- function(beta, X, y, bias){

  sum( dnorm(y, mean = X %*% ( (beta + bias) / 0.5) , sd = 1, log = T) )

}

nn_optimise_pre_approx_llhood_transformation <- function(particles, b_s_start, loglike, temp){

  #log_theta is matrix
  # should only send values of log_like that have been cached already
  true_ll <- apply(particles, MARGIN = 1, loglike, type = "full_likelihood", temp = temp)

  len <- ncol(particles)

  if(missing(b_s_start) | is.null(b_s_start)){

    b_s_start <- rep(0, times = len * 2)

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

  optim(par = b_s_start, fn = topt, control = list(maxit = 5), method = "BFGS")

}

best_step_scale_ejd_v_time <- function(step_scale, dist, comptime){

  tibble(step_scale = step_scale, dist = dist, comptime = comptime) %>%
    group_by(step_scale) %>%
    summarise(objective = mean(dist/comptime), mean_dist = mean(dist)) %>%
    {list(
      step_scale = .$step_scale[which.max(.$objective)],
      dist = .$mean_dist[which.max(.$objective)]
    )}

}

best_step_scale <- function(eta, dist, D, rho, max_T = 10, surrogate_acceptance, surrogate_cost, full_cost, model = "empirical", da = T){

  find_min_iter <- switch(model, # nned to make more robust to situations where only 1 is accepted.
                          empirical = time_steps_to_min_quantile_dist_emp,
                          normal = time_steps_to_min_quantile_dist_normal,
                          gamma = time_steps_to_min_quantile_dist_gamma,
                          bootstrap = time_steps_to_min_quantile_dist_bootstrap
  )

  unq_eta <- sort(unique(eta), decreasing = T)
  min_T_res <- vector(mode = "list", length = length(unq_eta))

  saccept_tb <- tibble(accept = surrogate_acceptance, eta = eta) %>%
    group_by(eta) %>%
    summarise(surrogate_accept_rate = mean(accept)) %>%
    arrange(-eta)
  # to be same order as unq_eta

  surrogate_acceptance_rate <- saccept_tb$surrogate_accept_rate

  for(i in 1:length(unq_eta)){

    ue <- unq_eta[i]
    min_T_res[[i]] <- find_min_iter(dist = dist[eta == ue], D = D, rho = rho, max_T = max_T)

  }

  which_eta_consider <- sapply(min_T_res, getElement, name = "sufficient_iter")

  if(all(!which_eta_consider)){
    warning("No MH tuning parameters have sufficient iterations to reach target median.")
    which_eta_consider[] <-  T
    # consider all but last and find cheapest after 10 iterations,
    # consider making max_T the maximum iterations possible.
  }

  min_T <- sapply(min_T_res, getElement, name = "iter")
  min_T[!which_eta_consider] <- Inf

  if(da){
    which_mintotal_cost <- min_da_mh_cost(min_T, surrogate_acceptance_rate, surrogate_cost, full_cost)
  } else {
    which_mintotal_cost <- min_mh_cost(min_T)
  }

  return(list(step_scale = eta[which_mintotal_cost], expected_iter = min_T[which_mintotal_cost]))

}

nn_posterior <- function(y, X, sigma, tau){
  p <- ncol(X)
  n <- length(y)
  XtX <- crossprod(X)
  Xty <- crossprod(X,y)
  mu <- solve(XtX + diag(1/(tau^2), nrow = p), Xty)
  V <- solve(XtX + diag(1/(tau^2), nrow = p)) * (sigma^2)

  # log_z_lhood <- -( n/2 ) * ( log(2) + log(pi) + 2 * log(sigma) )
  # log_z_prior <- -( p/2 ) * ( log(2) + log(pi) + 2 * log(tau) )
  # log_z_post <-  -( p/2 ) * ( log(2) + log(pi) ) - ( log(det(V))/2 )

  log_z <-
    sum( dnorm(x = rep(0, times = p), mean = rep(0, times = p), sd = tau, log = T) ) +
    sum( dnorm(x = y, mean = rep(0, times = n), sd = sigma, log = T) ) -
    sum( dmvnorm(x = rep(0, times = p), mean = mu, sigma = V, log = T) )

  return(list(
    mu = mu,
    var = V,
    log_z = log_z # log_z_prior + log_z_lhood - log_z_post
  ))

}

# true_theta = (birth rate prey, death rate prey/ birth rate pred, death rate pred)
g_pars <- list(N = 100, N_approx = 100, true_beta = c(0, 0.5, -1.5, 1.5, -3),
               log_prior = log_prior, log_like = log_like, log_like_approx = log_like_approx,
               draw_prior = draw_prior, optimise_pre_approx_llhood_transformation = nn_optimise_pre_approx_llhood_transformation,
               best_step_scale = best_step_scale
)

sim_settings <- rep(list(list(f_pars = NULL, g_pars = g_pars)), 2)

sim_settings[[1]]$f_pars <- list(
  num_p = 100,
  step_scale_set = c(0.25, 0.4, 0.55, 0.7),
  b_s_start = rep(0, 10),
  approx_ll_bias = 1,
  bss_model = "normal",
  bss_D = 1,
  bss_rho = 0.5
)

sim_settings[[2]]$f_pars <- list(
  num_p = 200,
  step_scale_set =  c(0.25, 0.4, 0.55, 0.7),
  b_s_start = rep(0, 10),
  approx_ll_bias = 1,
  bss_model = "normal",
  bss_D = 1,
  bss_rho = 0.5
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

# ss = sim_settings[[1]]; verbose = T

run_sim <- function(ss, verbose = F){

  ## simulate data

  sim <- simulate_regr(N = ss$g_pars$N, beta = ss$g_pars$true_beta)

  source("analysis/nn-smc-da.R")

  log_like <- function(beta) ss$g_pars$log_like(beta, X = sim$X, y = sim$y)
  log_like_approx <- function(beta) ss$g_pars$log_like_approx(beta, X = sim$X[1:ss$g_pars$N_approx,],
                                                              y = sim$y[1:ss$g_pars$N_approx],
                                                              bias = ss$f_pars$approx_ll_bias)

  best_step_scale_f <- function(eta, dist, surrogate_acceptance, surrogate_cost, full_cost, da = T){
    best_step_scale(eta = eta, dist = dist, D = ss$f_pars$bss_D, rho = ss$f_pars$bss_rho, max_T = 10,
                    surrogate_acceptance = surrogate_acceptance, surrogate_cost = surrogate_cost, full_cost = full_cost,
                    model = ss$f_pars$bss_model, da = da)
  }

  ## run
  cat("smc da")
  smc_da <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    b_s_start = b_s_start,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like,
    log_like_approx = log_like_approx,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  ss$g_pars$optimise_pre_approx_llhood_transformation,
    find_best_step_scale = best_step_scale_f,
    verbose = verbose
  )
  )
  cat("smc full")
  smc_full <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    b_s_start = b_s_start,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like,
    log_like_approx = log_like_approx,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  NULL,
    find_best_step_scale = best_step_scale_f,
    verbose = verbose
  )
   )

  cat("smc approx")
  smc_approx <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = T,
    num_p = num_p,
    step_scale_set = step_scale_set,
    b_s_start = b_s_start,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like,
    log_like_approx = log_like_approx,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  NULL,
    find_best_step_scale = best_step_scale_f,
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

res <- run_sim(ss = sim_settings[[2]], verbose = T)

#res <- future_map(1:100, .f = function(x) , .progress = T)

da_opt_approx_ll_hood <- function(b, res, i){

  # pre-trans is identity if not turned on.
  res$smc_da$log_ann_post_ctmc_da(
    x = b, temp = 1,
    lh_trans = res$smc_da$iter_summary[[i]]$approx_llh_pre_trans,
    type = "approx_likelihood"
    )

}

ll_hood <- function(b, res, i){

  # pre-trans is identity if not turned on.
  res$smc_da$log_ann_post_ctmc_da(
    x = b, temp = 1,
    type = "full_likelihood"
  )

}

mean_beta <- colMeans(res$smc_da$particles)

approx_ll <- function(x) sapply(x, function(x) da_opt_approx_ll_hood(b = replace(mean_beta, list = 1, x), res = res, i = 1))
full_ll <- function(x) sapply(x, function(x) ll_hood(b = replace(mean_beta, list = 1, x), res = res, i = 1))

x <- seq(-10,10, length.out = 200)

plot(x, exp(approx_ll(x)), type = "l")
lines(x, exp(full_ll(x)), col = "red")
