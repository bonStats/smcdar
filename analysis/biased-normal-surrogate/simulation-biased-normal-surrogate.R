#### load packages ####
library(smcdar)
library(dplyr)
library(glmnet)

#### Settings ####

# Select settings to use. See variable `sim_settings` below.
job_n <- 3

# How long should full likelihood (artifically) sleep for?
# (This creates an expensive likelihood to evaluate)
sleep_ll <- 10000

# the surrogate likelihood is currently set at sleep = 0.01 units,
# so the computation ratio (rho) is:
sleep_ll/0.01

# Total jumping distance required in mutation steps
set_bss_D <- qchisq(0.2, 5)

## Simulate data

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
  b <- matrix(
    rnorm(n * 5, mean = rep(0, times = 5), sd = 2), ncol = 5, byrow = T
  )

  list(beta = b)
}

# Likelihood: Normal
log_like_norm <- function(beta, X, y, sleep_time = 10, real_sleep = F){

  ll <- sum( dnorm(y, mean = X %*% beta, sd = 0.5, log = T) )

  if(sleep_time > 0){
    if(real_sleep) Sys.sleep(sleep_time)
    else attr(ll, "comptime") <- sleep_time
  }

  return(ll)

}

# Likelihood: Student-t
log_like_tdist <- function(beta, X, y, sleep_time = 10, real_sleep = F){

  ll <- sum( dt(y - X %*% beta, df = 3, log = T) )

  if(sleep_time > 0){
    if(real_sleep) Sys.sleep(sleep_time)
    else attr(ll, "comptime") <- sleep_time
  }

  return(ll)
}

# Surrogate Likelihood: Normal
log_like_approx <- function(beta, X, y, bias_mean, bias_scale = 1, sleep_time = 1, real_sleep = F, component_wise = F, calc_comps = T){

  stopifnot(is.logical(calc_comps))
  if(length(calc_comps) == 1) calc_comps <- rep(calc_comps, length(y))

  ll_comp <- rep(0, times = length(y))
  ll_comp[calc_comps] <- -( 50/length(y) ) + dnorm(y[calc_comps], mean = X[calc_comps,] %*% ( bias_scale * beta + bias_mean), sd = 1, log = T)

  if(component_wise) ll <- ll_comp
  else ll <- sum(ll_comp)

  if(sleep_time > 0){
    if(real_sleep) Sys.sleep( sleep_time * mean(calc_comps) )
    else attr(ll, "comptime") <- sleep_time * mean(calc_comps)
  }

  return(ll)
}

## Surrogate likelihood calibration function

optimise_approx_llhood_hybrid <- function(particles, par_start, initial_par_start, penalty, loglike, log_like_approx_comps, max_iter = 50, ...){

  res <- list(par = NULL, trans = identity, weights = NULL)

  xval <- particles[]
  true_ll <- loglike(particles, temp = 1, type = "full_likelihood", comp_time = F)

  len <- ncol(xval)

  if(missing(par_start) | is.null(par_start)){

    par_start <- rep(0, times = len + 1)

  }

  approx_log_like <- function(b, ...){
    # x value
    xv <- cbind(...)

    xv_trf <- t(apply(xv, MARGIN = 1, function(x) (x - b) ))
    ll_val <- loglike(xv_trf, temp = 1, type = "approx_likelihood", comp_time = F)

    return(ll_val)

  }

  dat <- data.frame(y = true_ll, xval)

  nls_ridge <- tryCatch(
    {nls(formula = y ~ approx_log_like(b, X1, X2, X3, X4, X5) + mu, #ridge
         data = dat, start = list(b = par_start[1:5], mu = par_start[6]),
         control = list(maxiter = max_iter))},
    error = function(e) paste(e)
  )

  if(class(nls_ridge) == "nls") {

    res$par <- coef(nls_ridge)
    res$trans <- function(x){ x - res$par[1:len] }

  } else {

    warning(paste("nls did not converge in", max_iter, "iterations"))
    res$par <- par_start * 0.5 + initial_par_start * 0.5

  }

  tr_approx_comps <- t(apply(xval, MARGIN = 1, FUN = function(x) log_like_approx_comps(res$trans(x), calc_comps = T)))

  # shrink to weights to 1
  glm_fit <- glmnet::cv.glmnet(x = tr_approx_comps, y = true_ll - rowSums(tr_approx_comps), nfolds = 20, alpha = 1) # alpha = 1 <=> lasso

  res$weights <- 1 + as.vector( coef(glm_fit, s = c("lambda.1se", "lambda.min")[1]) )[-1] # no intercept

  return(res)

}

## Function to choose best step scale

best_step_scale <- function(eta, dist, prob_accept, D, rho, max_T = 10, surrogate_expected_acceptance, surrogate_cost, full_cost, model = "empirical", da = T){
  # surrogate_expected_acceptance is a probability
  find_min_iter <- switch(model, # nned to make more robust to situations where only 1 is accepted.
                          empirical = time_steps_to_min_quantile_dist_emp,
                          normal = time_steps_to_min_quantile_dist_normal,
                          gamma = time_steps_to_min_quantile_dist_gamma,
                          bootstrap = time_steps_to_min_quantile_dist_bootstrap,
                          median = time_steps_to_min_quantile_dist_median
  )

  unq_eta <- sort(unique(eta), decreasing = T)
  min_T_res <- vector(mode = "list", length = length(unq_eta))

  for(i in 1:length(unq_eta)){

    ue <- unq_eta[i]
    min_T_res[[i]] <- find_min_iter(dist = dist[eta == ue], prob_accept = prob_accept[eta == ue],  D = D, rho = rho, max_T = max_T)

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
    saccept_tb <- tibble(accept = surrogate_expected_acceptance, eta = eta) %>%
      group_by(eta) %>%
      summarise(surrogate_accept_rate = mean(accept), .groups = "drop") %>%
      arrange(-eta)
    # to be same order as unq_eta

    surrogate_acceptance_rate <- saccept_tb$surrogate_accept_rate

    which_mintotal_cost <- min_da_mh_cost(min_T, surrogate_acceptance_rate, surrogate_cost, full_cost)
  } else {
    which_mintotal_cost <- min_mh_cost(min_T)
  }

  return(list(step_scale = unq_eta[which_mintotal_cost], expected_iter = min_T[which_mintotal_cost]))

}

## Simulation settings

# Specify the possible simulation settings (`job_n` specifies which `sim_setting` is used.)

g_pars1 <- list(N = 100, N_approx = 100, true_beta = c(0, 0.5, -1.5, 1.5, -3),
                log_prior = log_prior, log_like = log_like_norm, sim_func = simulate_regr_norm,
                log_like_approx = log_like_approx,
                draw_prior = draw_prior,
                best_step_scale = best_step_scale,
                save_post_interface = F,
                max_mh_steps = 100
)

g_pars2 <- within(g_pars1,{
                    log_like = log_like_tdist
                    sim_func = simulate_regr_tdist
                  })

sim_settings <- c(rep(list(list(f_pars = NULL, g_pars = g_pars1)), 2),
                  rep(list(list(f_pars = NULL, g_pars = g_pars2)), 2)
)

sim_settings <- c(sim_settings, sim_settings)


sim_settings[[1]]$f_pars <- sim_settings[[3]]$f_pars <- list(
  num_p = 2000,
  step_scale_set = c(0.1, seq(from = 0.25, to = 3.25, by = 0.5)),
  par_start = rep(0, 11),
  approx_ll_bias_mean = 0.25,
  approx_ll_bias_scale = exp(0.1),
  sleep_ll = sleep_ll,
  sleep_ll_approx = 0.01,
  bss_model_std = "median",
  bss_model = "median",
  bss_D = set_bss_D,
  bss_rho = 0.5,
  max_anneal_temp = 0.1,
  calibration_method = optimise_approx_llhood_hybrid,
  ll_tune_shrinage_penalty = 15
)


sim_settings[[2]]$f_pars <- sim_settings[[4]]$f_pars <- list(
  num_p = 2000,
  step_scale_set =  c(0.1, seq(from = 0.25, to = 3.25, by = 0.5)),
  par_start = rep(0, 11),
  approx_ll_bias_mean = 0.25,
  approx_ll_bias_scale = exp(0.1),
  sleep_ll = sleep_ll,
  sleep_ll_approx = 0.01,
  bss_model_std = "median",
  bss_model = "gamma",
  bss_D = set_bss_D,
  bss_rho = 0.5,
  max_anneal_temp = 0.1,
  calibration_method = optimise_approx_llhood_hybrid,
  ll_tune_shrinage_penalty = 15
)


sim_settings[[5]]$f_pars <- sim_settings[[7]]$f_pars <- list(
  num_p = 2000,
  step_scale_set = c(0.1, seq(from = 0.25, to = 3.25, by = 0.5)),
  par_start = rep(0, 11),
  approx_ll_bias_mean = 0.25,
  approx_ll_bias_scale = exp(0.1),
  sleep_ll = sleep_ll,
  sleep_ll_approx = 0.01,
  bss_model_std = "median",
  bss_model = "bootstrap",
  bss_D = set_bss_D,
  bss_rho = 0.5,
  max_anneal_temp = 0.1,
  calibration_method = optimise_approx_llhood_hybrid,
  ll_tune_shrinage_penalty = 15
)

# not in use
sim_settings[[6]] <- sim_settings[[8]] <- NULL

####

#### Code for sim ####

# To test inner code use:
# ss = sim_settings[[job_n]]; verbose = T

# use function to run set of algorithms multiple times

run_sim <- function(ss, verbose = F){

  stopifnot(is.list(ss))

  # simulate data
  sim <- ss$g_pars$sim_func(N = ss$g_pars$N, beta = ss$g_pars$true_beta)

  # log likelihood function
  log_like_f <- function(beta) ss$g_pars$log_like(beta, X = sim$X, y = sim$y, sleep_time = ss$f_pars$sleep_ll)
  # surrogate log likelihood function
  log_like_approx_f <- function(beta, ...) ss$g_pars$log_like_approx(beta, X = sim$X[1:ss$g_pars$N_approx,],
                                                                y = sim$y[1:ss$g_pars$N_approx],
                                                                bias_mean = ss$f_pars$approx_ll_bias_mean,
                                                                bias_scale = ss$f_pars$approx_ll_bias_scale,
                                                                sleep_time = ss$f_pars$sleep_ll_approx,
                                                                ...)
  # component-wise surrogate log likelihood function
  log_like_approx_comps <- function(beta, calc_comps) log_like_approx_f(beta = beta, calc_comps = calc_comps, component_wise = T)

  # choose best step scale function (all SMC algorithms except standard)
  best_step_scale_f <- function(eta, dist, prob_accept, adjust_D, surrogate_expected_acceptance, surrogate_cost, full_cost, da = T){
    best_step_scale(eta = eta, dist = dist, prob_accept = prob_accept, D = ss$f_pars$bss_D - adjust_D, rho = ss$f_pars$bss_rho, max_T = 100,
                    surrogate_expected_acceptance = surrogate_expected_acceptance, surrogate_cost = surrogate_cost, full_cost = full_cost,
                    model = ss$f_pars$bss_model, da = da)
  }

  # choose best step scale function (for standard SMC)
  best_step_scale_f_std <- function(eta, dist, prob_accept, adjust_D, surrogate_expected_acceptance, surrogate_cost, full_cost, da = T){
    best_step_scale(eta = eta, dist = dist, prob_accept = prob_accept, D = ss$f_pars$bss_D - adjust_D, rho = ss$f_pars$bss_rho, max_T = 100,
                    surrogate_expected_acceptance = surrogate_expected_acceptance, surrogate_cost = surrogate_cost, full_cost = full_cost,
                    model = ss$f_pars$bss_model_std, da = da)
  }

  ## run
  if(verbose) cat("smc da\n")
  smc_da <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = par_start,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    log_like_approx_comps = log_like_approx_comps,
    draw_prior = ss$g_pars$draw_prior,
    calibrate_approx_likelihood = calibration_method,
    find_best_step_scale = best_step_scale_f,
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  if(verbose) cat("smc da - no T\n")
  smc_da_no_trans <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = par_start,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = ss$g_pars$draw_prior,
    calibrate_approx_likelihood =  NULL,
    find_best_step_scale = best_step_scale_f,
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  if(verbose) cat("smc standard\n")
  smc_standard <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = ss$g_pars$draw_prior,
    calibrate_approx_likelihood = NULL,
    find_best_step_scale = best_step_scale_f_std,
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  # find best step scale, test +-(0.5) above it with standard SMC algorithm
  avg_step_scale <- sapply(smc_standard$iter_summary, getElement, name = "step_scale") %>% mean()

  step_scale_fixed <- function(..., offset){
    return(list(step_scale = max(0.1, avg_step_scale + offset), expected_iter = NA))
  }

  if(verbose) cat("smc standard avg+0.5 (2) \n")
  smc_standard_fixed_step_p5 <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = ss$g_pars$draw_prior,
    calibrate_approx_likelihood = NULL,
    find_best_step_scale = function(...) step_scale_fixed(..., offset = 0.5),
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  # Don't save, just want timings
  smc_standard_fixed_step_p5$particle_list <- NULL
  smc_standard_fixed_step_p5$iter_summary <- NULL

  if(verbose) cat("smc standard avg-0.5 (3) \n")
  smc_standard_fixed_step_m5 <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = ss$g_pars$draw_prior,
    calibrate_approx_likelihood = NULL,
    find_best_step_scale = function(...) step_scale_fixed(..., offset = -0.5),
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  # Don't save, just want timings
  smc_standard_fixed_step_m5$particle_list <- NULL
  smc_standard_fixed_step_m5$iter_summary <- NULL

  if(verbose) cat("smc approx\n")
  smc_approx <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = T,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = ss$g_pars$draw_prior,
    calibrate_approx_likelihood =  NULL,
    find_best_step_scale = best_step_scale_f,
    max_anneal_temp = max_anneal_temp,
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  if(verbose) cat("smc + sfa\n")
  smc_sfa <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = F, start_from_approx = T,
    start_from_approx_fit = smc_approx,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = par_start,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = NULL,
    calibrate_approx_likelihood = NULL,
    find_best_step_scale = best_step_scale_f,
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  if(verbose) cat("smc + sfa + da\n")
  smc_sfa_da <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F, start_from_approx = T,
    start_from_approx_fit = smc_approx,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = par_start,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    log_like_approx_comps = log_like_approx_comps,
    draw_prior = NULL,
    calibrate_approx_likelihood = calibration_method,
    find_best_step_scale = best_step_scale_f,
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  if(verbose) cat("smc + sfa + da no trans\n")
  smc_sfa_da_no_trans <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F, start_from_approx = T,
    start_from_approx_fit = smc_approx,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = par_start,
    refresh_ejd_threshold = bss_D,
    log_prior = ss$g_pars$log_prior,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = NULL,
    calibrate_approx_likelihood = NULL,
    find_best_step_scale = best_step_scale_f,
    max_mh_steps = ss$g_pars$max_mh_steps,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose
  )
  )

  cat("done\n")

  return(
    list(sim_settings = ss,
         sim_data = sim,
         smc_da = smc_da,
         smc_da_no_trans = smc_da_no_trans,
         smc_sfa = smc_sfa,
         smc_sfa_da = smc_sfa_da,
         smc_sfa_da_no_trans = smc_sfa_da_no_trans,
         smc_standard = smc_standard,
         smc_standard_fixed_step_p5 = smc_standard_fixed_step_p5,
         smc_standard_fixed_step_m5 = smc_standard_fixed_step_m5,
         smc_approx = smc_approx
    )
  )

}



####

try_catch_all <- function(expr) {
  warn <- err <- NULL
  val <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=val, warning=warn, error=err)
}

smc_result <- try_catch_all(
    run_sim(ss = sim_settings[[job_n]], verbose = T)
    )

# The following warning message is expected for the SMC algorithms which don't tune
# the surrogate likelihood (because the surrogate likelihood is biased in this
# example. It can be ignored in SMC algorithms that DO use surrogate likelihood calibration
# if it occurs infrequently.
#   In mh_da_step_bglr(...)
#     Positive slope for step scale versus full acceptance ratio:
#     May indicate bad surrogate or too few first stage acceptances.
# It refers to the procedure which estimates the second stage acceptance probabilities in DA-MH
# when the particular particle proposal is rejected in the first stage.
