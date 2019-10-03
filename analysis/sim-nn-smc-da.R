.libPaths(Sys.getenv("R_LIB_USER"))
outdir <- Sys.getenv("OUT_DIRECTORY")

devtools::install_github("bonStats/smcdar", lib = Sys.getenv("R_LIB_USER"))

#### load packages ####

library(smcdar)
library(dplyr)
library(parallel)

#### Settings ####

## Sim

simulate_regr <- function(N, beta, error_dist = rnorm, ...){

  X <- matrix(rnorm(length(beta)*N), nrow = N)
  y <- X %*% beta + error_dist(n = N, ...)

  return(
    list(y = y, X = X, beta = beta)
  )

}

simulate_regr_norm  <- function(N, beta) simulate_regr(N, beta, error_dist = rnorm, sd = 0.5)
simulate_regr_tdist  <- function(N, beta) simulate_regr(N, beta, error_dist = rt, df = 3)

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

log_like_norm <- function(beta, X, y, sleep_time = 0.1){
  if(sleep_time > 0) Sys.sleep(sleep_time)
  sum( dnorm(y, mean = X %*% beta, sd = 0.5, log = T) )
}

log_like_tdist <- function(beta, X, y, sleep_time = 0.1){
  if(sleep_time > 0) Sys.sleep(sleep_time)
  sum( dt(y - X %*% beta, df = 3, log = T) )
}

log_like_approx <- function(beta, X, y, bias_mean, bias_scale = 1){

  -50 + sum( dnorm(y, mean = X %*% ( bias_scale * beta + bias_mean), sd = 1, log = T) )

}
# need to change Dlog_like_approx_Dbeta inorder to update this
# log_like_approx <- function(beta, X, y, bias_mean, bias_scale = 1){
#
#   1 + sum( dt(y -  X %*% ( bias_scale * beta + bias_mean), df = 2, log = T) )
#
# }

Dlog_like_approx_Dbeta <- function(beta, X, y, bias_mean, bias_scale = 1){
  # sd = 1
  mu = X %*% ( bias_scale * beta + bias_mean)

  as.vector( bias_scale * crossprod(X, y - mu) )

}

particle_optimise_pre_approx_llhood_transformation <- function(particles, b_s_start, loglike, temp, max_iter = 20, ...){

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

  optim(par = b_s_start, fn = topt, control = list(maxit = max_iter), method = "BFGS")

}

visited_optimise_pre_approx_llhood_transformation <- function(b_s_start, loglike, temp, max_iter = 20, max_strata_size = 25, ...){

  # should only send values of log_like that have been cached already
  unique_x_val <- get("unique_x", envir = environment(loglike))
  true_ll_val <- apply(unique_x_val, MARGIN = 1, loglike, type = "full_likelihood", temp = 1)

  tlv_quant <- quantile(true_ll_val, probs = 0:10/10)
  tlv_quant <- cbind(lower = tlv_quant[-11], upper = tlv_quant[-1])

  strata_size <- min(max_strata_size, floor(length(true_ll_val)/10) )

  true_ll_sample_ind <- c(apply(tlv_quant, MARGIN = 1, function(l_u){
    sample(which(l_u[1] <= true_ll_val & true_ll_val <= l_u[2]),
           size = strata_size, replace = F)
  }))

  xval <- unique_x_val[true_ll_sample_ind,]
  true_ll <- true_ll_val[true_ll_sample_ind]

  len <- ncol(xval)

  if(missing(b_s_start) | is.null(b_s_start)){

    b_s_start <- rep(0, times = len * 2)

  }

  topt <- function(b_s){

    b_mat <- matrix(b_s[1:len], ncol = len, nrow = nrow(xval), byrow = T)
    s_mat <- matrix(b_s[(len+1):(2*len)], ncol = len, nrow = nrow(xval), byrow = T)

    approx_ll <- apply((xval - b_mat)/exp(s_mat), MARGIN = 1, loglike, type = "approx_likelihood", temp = 1)

    keep_ll <- is.finite(approx_ll)

    dm_true_ll <- true_ll[keep_ll] - mean(true_ll[keep_ll])
    dm_approx_ll <- approx_ll[keep_ll] - mean(approx_ll[keep_ll])

    # Least-squares
    sum(
      ( dm_true_ll - dm_approx_ll )^2
    )

  }

  optim(par = b_s_start, fn = topt, control = list(maxit = max_iter), method = "BFGS")

}

visit_gr_optimise_pre_approx_llhood_transformation <- function(par_start, loglike, D_approx_log_like, temp, max_iter = 40, max_strata_size = 20, ...){

  # should only send values of log_like that have been cached already
  unique_x_val <- get("unique_x", envir = environment(loglike))
  true_ll_val <- apply(unique_x_val, MARGIN = 1, loglike, type = "full_likelihood", temp = 1)

  tlv_quant <- quantile(true_ll_val, probs = 0:10/10)
  tlv_quant <- cbind(lower = tlv_quant[-11], upper = tlv_quant[-1])

  strata_size <- min(max_strata_size, floor(length(true_ll_val)/10) )

  true_ll_sample_ind <- c(apply(tlv_quant, MARGIN = 1, function(l_u){
    sample(which(l_u[1] <= true_ll_val & true_ll_val <= l_u[2]),
           size = strata_size, replace = F)
  }))

  xval <- unique_x_val[true_ll_sample_ind,]
  true_ll <- true_ll_val[true_ll_sample_ind]

  len <- ncol(xval)

  if(missing(par_start) | is.null(par_start)){

    par_start <- rep(0, times = len * 2 + 1)

  }


  gropt <- function(vr, penalty){
    # xval, loglike, len
    mu <- vr[1]
    b <- vr[2:(len+1)]
    s <- vr[(len+2):(2*len + 1)]


    beta_trf <- apply(xval, MARGIN = 1, function(x) exp(-s) * (x - b) )
    ll_beta_trf <- apply(beta_trf, MARGIN = 2, loglike, type = "approx_likelihood", temp = 1)

    Dbeta_db <- - matrix(exp(-s), ncol = ncol(beta_trf), nrow = nrow(beta_trf))
    Dbeta_ds <- - beta_trf
    Dll_dbeta <- apply(beta_trf, MARGIN = 2, D_approx_log_like)
    Dr_dll <- -2 * (true_ll - ll_beta_trf - mu)
    Dr_dmu <- -2 * (true_ll - ll_beta_trf - mu)

    Dpen_b <- 2 * penalty * b
    Dpen_s <-  - 2 * penalty * (exp(-s) - 1) * exp(-s)

    c(
      sum(Dr_dmu),
      (Dll_dbeta * Dbeta_db) %*% Dr_dll + Dpen_b,
      (Dll_dbeta * Dbeta_ds) %*% Dr_dll + Dpen_s
    )

  }

  topt <- function(vr, penalty){

    mu <- vr[1]
    b <- vr[2:(len+1)]
    s <- vr[(len+2):(2*len + 1)]

    approx_ll <- apply(xval, MARGIN = 1, function(x) loglike( exp(-s) * (x - b), type = "approx_likelihood", temp = 1))
    # Least-squares
    sum( ( true_ll - approx_ll - mu )^2 ) +
      penalty * ( sum(b^2) + sum((exp(-s) - 1)^2) )

  }

  optim(par = par_start, fn = topt, gr = gropt, control = list(maxit = max_iter), method = "BFGS", penalty = 10)

}

visit_nls_optimise_pre_approx_llhood_transformation <- function(par_start, penalty, loglike, D_approx_log_like, temp, max_iter = 50, max_strata_size = 25, add_noise = T, ...){

  # should only send values of log_like that have been cached already
  unique_x_val <- get("unique_x", envir = environment(loglike))
  true_ll_val <- apply(unique_x_val, MARGIN = 1, loglike, type = "full_likelihood", temp = 1)

  tlv_quant <- quantile(true_ll_val, probs = 0:10/10)
  tlv_quant <- cbind(lower = tlv_quant[-11], upper = tlv_quant[-1])

  strata_size <- min(max_strata_size, floor(length(true_ll_val)/10) )

  true_ll_sample_ind <- c(apply(tlv_quant, MARGIN = 1, function(l_u){
    sample(which(l_u[1] <= true_ll_val & true_ll_val <= l_u[2]),
           size = strata_size, replace = F)
  }))

  xval <- unique_x_val[true_ll_sample_ind,]
  true_ll <- true_ll_val[true_ll_sample_ind]

  # add noise to make optimisation more realistic
  if(add_noise) true_ll <- true_ll + rnorm(n = length(true_ll), sd = sd(true_ll)/50)

  len <- ncol(xval)

  if(missing(par_start) | is.null(par_start)){

    par_start <- rep(0, times = len * 2 + 1)

  }

  approx_log_like <- function(b, s, ...){
    # value
    xv <- cbind(...)

    xv_trf <- t(apply(xv, MARGIN = 1, function(x) exp(-s) * (x - b) ))
    ll_val <-  apply(xv_trf, MARGIN = 1, loglike, type = "approx_likelihood", temp = 1)

    Dbeta_db <- - matrix(exp(-s), ncol = ncol(xv_trf), nrow = nrow(xv_trf))
    Dbeta_ds <- - xv_trf
    Dll_dbeta <- t(apply(xv_trf, MARGIN = 1, D_approx_log_like))

    # Somehow broke the gradient, but not neccesary for nls
    # grad <- cbind(
    #   (Dll_dbeta * Dbeta_db),
    #   (Dll_dbeta * Dbeta_ds)
    # )
    #
    # attr(ll_val, "gradient") <- grad

    return(ll_val)

  }


  dat <- rbind(
    data.frame(y = true_ll, xval, ll = TRUE, w = 1),
    data.frame(y = 0, matrix(0, nrow = len*2, ncol = ncol(xval)), ll = FALSE, w = sqrt(penalty))
  )

  nls_ridge <- nls(formula =
                     y ~ ll * ( approx_log_like(b, s, X1, X2, X3, X4, X5) + mu ) + # least squares
                     I(!ll) * ( identity(b) + identity(s) ), #ridge
                   data = dat, start = list(b = par_start[1:5], s = par_start[6:10], mu = par_start[11]),
                   weights = w,
                   control = list(tol = 1e-02, maxiter = max_iter)
  )

  # mu first
  list(par = coef(nls_ridge))

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

g_pars1 <- list(N = 100, N_approx = 100, true_beta = c(0, 0.5, -1.5, 1.5, -3),
                log_prior = log_prior, log_like = log_like_norm, sim_func = simulate_regr_norm,
                log_like_approx = log_like_approx,
                Dlog_like_approx_Dbeta = Dlog_like_approx_Dbeta,
                draw_prior = draw_prior, optimise_pre_approx_llhood_transformation = visit_nls_optimise_pre_approx_llhood_transformation,
                best_step_scale = best_step_scale,
                save_post_interface = F
)

g_pars2 <- within(g_pars1,
                  {
                    log_like = log_like_tdist
                    sim_func = simulate_regr_tdist
                  }
)

sim_settings <- c(rep(list(list(f_pars = NULL, g_pars = g_pars1)), 2),
                  rep(list(list(f_pars = NULL, g_pars = g_pars2)), 2)
)

sim_settings[[1]]$f_pars <- sim_settings[[3]]$f_pars <- list(
  num_p = 250,
  step_scale_set = c(0.25, 0.4, 0.55, 0.7),
  par_start = rep(0, 11),
  approx_ll_bias_mean = 1,
  approx_ll_bias_scale = exp(0.1),
  bss_model = "normal",
  bss_D = 1,
  bss_rho = 0.5,
  ll_tune_shrinage_penalty = 5
)

sim_settings[[2]]$f_pars <- sim_settings[[4]]$f_pars <- list(
  num_p = 500,
  step_scale_set =  c(0.1, 0.25, 0.4, 0.55, 0.7, 0.85),
  par_start = rep(0, 11),
  approx_ll_bias_mean = 1,
  approx_ll_bias_scale = exp(0.1),
  bss_model = "normal",
  bss_D = 1,
  bss_rho = 0.5,
  ll_tune_shrinage_penalty = 5
)

####

#### Code for sim ####

# ss = sim_settings[[1]]; verbose = T

run_sim <- function(ss, verbose = F){

  ## simulate data

  sim <- ss$g_pars$sim_func(N = ss$g_pars$N, beta = ss$g_pars$true_beta)

  source("nn-smc-da-server.R")

  log_like_f <- function(beta) ss$g_pars$log_like(beta, X = sim$X, y = sim$y)
  log_like_approx_f <- function(beta) ss$g_pars$log_like_approx(beta, X = sim$X[1:ss$g_pars$N_approx,],
                                                                y = sim$y[1:ss$g_pars$N_approx],
                                                                bias_mean = ss$f_pars$approx_ll_bias_mean,
                                                                bias_scale = ss$f_pars$approx_ll_bias_scale)

  Dlog_like_approx_Dbeta_f <- function(beta) ss$g_pars$Dlog_like_approx_Dbeta(beta, X = sim$X[1:ss$g_pars$N_approx,],
                                                                              y = sim$y[1:ss$g_pars$N_approx],
                                                                              bias_mean = ss$f_pars$approx_ll_bias_mean,
                                                                              bias_scale = ss$f_pars$approx_ll_bias_scale)

  best_step_scale_f <- function(eta, dist, surrogate_acceptance, surrogate_cost, full_cost, da = T){
    best_step_scale(eta = eta, dist = dist, D = ss$f_pars$bss_D, rho = ss$f_pars$bss_rho, max_T = 10,
                    surrogate_acceptance = surrogate_acceptance, surrogate_cost = surrogate_cost, full_cost = full_cost,
                    model = ss$f_pars$bss_model, da = da)
  }


  optimise_pre_approx_llhood_transformation_f <- function(par_start, loglike, D_approx_log_like, temp, ...){
    ss$g_pars$optimise_pre_approx_llhood_transformation(par_start = par_start, penalty = ss$f_pars$ll_tune_shrinage_penalty,
                                                        loglike = loglike, D_approx_log_like = D_approx_log_like, temp = temp,
                                                        max_iter = 50, max_strata_size = 25, add_noise = F)
  }
  ## run
  if(verbose) cat("smc da\n")
  smc_da <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = par_start,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    Dlog_like_approx_Dbeta = Dlog_like_approx_Dbeta_f,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  optimise_pre_approx_llhood_transformation_f,
    find_best_step_scale = best_step_scale_f,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )
  if(verbose) cat("smc full\n")
  smc_full <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    Dlog_like_approx_Dbeta = NULL,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  NULL,
    find_best_step_scale = best_step_scale_f,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  if(verbose) cat("smc approx\n")
  smc_approx <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = T,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    Dlog_like_approx_Dbeta = NULL,
    draw_prior = ss$g_pars$draw_prior,
    optimise_pre_approx_llhood_transformation =  NULL,
    find_best_step_scale = best_step_scale_f,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  if(verbose) cat("smc approx + da\n")
  smc_approx_then_da <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F, start_from_approx = T,
    start_from_approx_fit = smc_approx,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = par_start,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    Dlog_like_approx_Dbeta = Dlog_like_approx_Dbeta_f,
    draw_prior = NULL,
    optimise_pre_approx_llhood_transformation =  optimise_pre_approx_llhood_transformation_f,
    find_best_step_scale = best_step_scale_f,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  cat("done")

  return(
    list(sim_settings = ss,
         sim_data = sim,
         smc_da = smc_da,
         smc_full = smc_full,
         smc_approx = smc_approx,
         smc_approx_then_da = smc_approx_then_da
    )
  )

}

####

cat("parallel::detectCores():", parallel::detectCores(),"\n")

for(ssi in 2:length(sim_settings)){

  cat("*SIM",ssi,"*\n")
  res <- mclapply(1:100, function(...)
  {
    tryCatch(run_sim(ss = sim_settings[[ssi]]),
             error=function(e) {
               cat(paste(e))
               paste(e)
             })
  }, mc.cores = 15
  )
  saveRDS(res, file.path(outdir, paste0("sim_res_",ssi,"_20190926.rds")))
  rm(res)
}
