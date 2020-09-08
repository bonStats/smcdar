#### load packages ####

library(smcdar)
library(dplyr)
library(parallel)
library(arfima)
library(TSA)
library(matrixStats)
library(mvtnorm)
library(glmnet)

#### Settings ####

# number of parallel cores to use for mutation step
par_cores <- 12

# load data (if not simulating)
whittle_data <- readRDS("analysis/whittle-likelihood/whittle_sim_data_8001.rds")

## Simulate data (if required, default data is in directory)

simulate_arfima <- function(N, phi, dfrac, theta, sigma2, ...){
  # phi = AR parameters
  # theta = MA paramters (Box-Jenkins notation)
  # dfrac = fractional differencing (long memory)

  y <- arfima::arfima.sim(
    n = N,
    model = list(phi = phi,
                 dfrac = dfrac,
                 theta = theta),
    sigma2 = sigma2
  )

  return(
    list(y = y, phi = phi, dfrac = dfrac, theta = theta, sigma2 = sigma2)
  )

}

## Priors

# prior:
# pphi ~ [-1,1] ( phi = arfima::PacfToAR(pphi) )
# ==> t_pphi ~ cosh(p)^{-2}
# ptheta ~ [-1,1]
# ==> ptheta ~ cosh(t)^{-2} ( theta = arfima::PacfToAR(ptheta) )
# t_dfrac ~ N(0,1) ( dfrac = 0.5 * tanh(t_dfrac) )
# log_sigma2 ~ N(0,1)
log_prior <- function(t_pphi, t_dfrac, t_ptheta, log_sigma2){

  -2 * sum(log(cosh(t_pphi))) +
    -2 * sum(log(cosh(t_ptheta))) +
    dnorm(t_dfrac, mean = 0, sd = 1, log = T) +
    dnorm(log_sigma2, mean = 0, sd = 1, log = T)

}

draw_prior <- function(n, len_phi, len_theta){

  list(
    t_pphi = matrix(atanh(runif(n = len_phi*n, min = -1, max = 1)), nrow = n),
    t_dfrac = matrix(rnorm(n = n, mean = 0, sd = 1), nrow = n),
    t_ptheta = matrix(atanh(runif(n = len_theta*n, min = -1, max = 1)), nrow = n),
    log_sigma2 = matrix(rnorm(n = n, mean = 0, sd = 1), nrow = n)
  )

}

## Functions for full likelihood

# redefine lARFIMA  from arfima package to handle variance (sigma)
lARFIMA2 <- function(z, phi = numeric(0), theta = numeric(0), dfrac = numeric(0), sigma2 = 1){

  EPS <- .Machine$double.eps

  if (is.null(phi) || any(is.na(phi))) phi <- numeric(0)
  if (is.null(theta) || any(is.na(theta))) theta <- numeric(0)

  n <- length(z)

  r <- tryCatch(
    arfima::tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, sigma2 = sigma2,
                        maxlag = n - 1),
    error = function(e) paste(e)
  )

  if(is.numeric(r)){
    r <- as.double(r) #not sure if neccesary
  } else {
    warning("In arfima::tacvfARFIMA: ", paste(r))
    return(-1e+10)
  }

  r0 <- r[1]
  r <- r/r[1]
  if (r[1] < EPS)  stop("error: r[1] the variance is <= 0")
  if (length(r) != length(z)) stop("error: r and z have unequal lengths")
  out <- .C("trenchQFR", as.double(r), as.integer(length(r)),
            as.double(z), as.integer(length(z)), EPS,
            tr = array(0, dim = c(1, 2)), fault = as.integer(1), PACKAGE = "ltsa")
  fault <- out$fault
  if (fault != 0){
    warning(c("Singular Matrix",
              "Input r[0] is not equal to 1.",
              "The length(r) is not equal to the length(z))")[fault]
    )
    return(-1e+10)
  }
  norm_quad_prod <- (out$tr)[1]
  norm_log_det <- (out$tr)[2]

  quad_prod <- norm_quad_prod / r0
  log_det <- norm_log_det + n * log(r0)

  #logl <- - 0.5 * (log_det  + quad_prod)
  logl <- - 0.5 * (n * log(2 * pi) + log_det  + quad_prod)
  # Section 2.3. http://dx.doi.org/10.18637/jss.v023.i05
  return(logl)
}

log_like_arfima <- function(y, t_pphi, t_dfrac, t_ptheta, log_sigma2){

  phi <- arfima::PacfToAR(tanh(t_pphi))
  theta <- arfima::PacfToAR(tanh(t_ptheta))
  dfrac <- 0.5 * tanh(t_dfrac)
  sigma2 <- exp(log_sigma2)

  # log-like
  ll <- lARFIMA2( #arfima::lARFIMA() doesn't have sigma term built in
    z = y,
    phi = phi,
    theta = theta,
    dfrac = dfrac,
    sigma2 = sigma2
  )

  return(ll)

}

## Functions for surrogate likelihood

spectral_farima <- function(phi, dfrac, theta, log_sigma2, spec_den_matrix, w){

  p <- length(phi);
  q <- length(theta);

  if(missing(spec_den_matrix)){
    spec_den_matrix <- outer(1:max(p,q), w, FUN = function(j,w) exp(-1i*w*j) )
  }

  if (p >= 1){
    denom <- 1 - as.complex(phi %*% spec_den_matrix[1:p,])
  } else {
    denom <- rep(1, ncol(spec_den_matrix));
  }

  if (q >= 1){
    numer <- 1 - as.complex(theta %*% spec_den_matrix[1:q,])
  } else {
    numer <- rep(1, ncol(spec_den_matrix));
  }

  log_f <- (-2*dfrac) * log(abs(1 - exp(-1i*w))) + 2 * (log(abs(numer)) - log(abs(denom)) ) + log_sigma2 - log(2*pi)

  return(exp(log_f));

}

log_like_whittle <- function(t_pphi, t_dfrac, t_ptheta, log_sigma2, pgram, w, component_wise = F, calc_comps = T, ...){
  # For '...' argument use: "spec_den_matrix" OR "w"
  phi <- arfima::PacfToAR(tanh(t_pphi))
  theta <- arfima::PacfToAR(tanh(t_ptheta))
  dfrac <- 0.5 * tanh(t_dfrac)

  stopifnot(is.logical(calc_comps))
  if(length(calc_comps) == 1) calc_comps <- rep(calc_comps, length(w))

  sp <- spectral_farima(phi = phi,
                        dfrac = dfrac,
                        theta = theta,
                        log_sigma2 = log_sigma2,
                        w = w[calc_comps],
                        ...
  )

  # whittle log-like
  ll_comp <- rep(NA_real_, length(w))
  ll_comp[calc_comps] <- log(sp) + pgram[calc_comps]/sp
  ll_comp[!calc_comps] <- 0

  if(component_wise) ll <- -2 * ll_comp
  else ll <- -2 * sum(ll_comp)

  return(ll)

}

calc_pgram <-function(y){

  pgram <- spec.pgram(y, plot = F)
  n <- length(y)

  list(
    w_freq = 2*pi*pgram$freq,
    p_spec = pgram$spec/(2*pi)
  )

}

optimise_approx_llhood_weights_lasso <- function(particles, loglike, log_like_approx_comps, ...){

  res <- list()

  true_ll <- loglike(particles, temp = 1, type = "full_likelihood", comp_time = F)

  xval <- particles[]
  approx_ll_comps <- t(apply(xval, MARGIN = 1, FUN = function(x) log_like_approx_comps(x, calc_comps = T)))

  # shrink to weights to 1
  glm_fit <- glmnet::cv.glmnet(x = approx_ll_comps, y = true_ll - rowSums(approx_ll_comps), nfolds = 20, alpha = 1) # alpha = 1 <=> lasso
  res$weights <- 1 + as.vector( coef(glm_fit, s = c("lambda.1se", "lambda.min")[1]) )[-1] # no intercept

  res$trans <- identity

  return(res)

}


best_step_scale <- function(eta, dist, prob_accept, D, rho, max_T = 10, surrogate_expected_acceptance, surrogate_cost, full_cost, model = "median", da = T){
  # surrogate_expected_acceptance is a probability
  find_min_iter <- switch(model, # need to make more robust to situations where only 1 is accepted.
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
    min_T_res[[i]] <- find_min_iter(
      dist = dist[eta == ue],
      prob_accept = prob_accept[eta == ue],
      D = D, rho = rho, max_T = max_T)

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

## Data and algorithm settings

g_pars <- list(N = length(whittle_data$y), # or specify: e.g. 8001
                true_phi = c(0.45,0.1), true_dfrac = 0.4,
                true_theta = -0.4, true_sigma2 = 1,
                log_prior = log_prior, log_like = log_like_arfima,
                sim_func = function(...) whittle_data, # or use: simulate_arfima,
                log_like_approx = log_like_whittle,
                draw_prior = draw_prior,
                best_step_scale = best_step_scale,
                use_robust_cov = T,
                save_post_interface = F
)

g_pars$len_post_par <- with(g_pars,
                             length(true_phi) + length(true_dfrac) +
                               length(true_theta) + length(true_sigma2)
)

g_pars$grp_post_par <- with(g_pars,
                             c(
                               rep(1, length(true_phi)),
                               rep(2, length(true_dfrac)),
                               rep(3, length(true_theta)),
                               rep(4, length(true_sigma2))
                             )
)


f_pars <- list(
  num_p = 5000,
  step_scale_set = seq(0.1, 2, by = 0.2),
  par_start = rep(0, 2 * g_pars$len_post_par + 1),
  bss_model_noda = "median",
  bss_model_da = "median",
  whittle_thin = 1, # 5, 10 didn't seem to work well
  smooth_pgram = F, # smoothing didn't seem to work well
  bss_D = 3,
  bss_rho = 0.5,
  sfa_max_anneal_temp = 0.01,
  calibration_method = optimise_approx_llhood_weights_lasso,
  cores = par_cores
)

sim_settings <- list()

sim_settings[[1]] <- list(f_pars = f_pars, g_pars = g_pars)

####

#### Code for sim ####

# To test inner code use:
# ss = sim_settings[[1]]; verbose = T

# use function to run set of algorithms multiple times

run_sim <- function(ss, verbose = F){

  # list of pars to vector
  pars_to_vec <- function(phi, dfrac, theta, sigma2){
    c(phi, dfrac, theta, sigma2)
  }

  # par vector to list
  pars_to_list <-function(x){

    setNames(split(x, f = ss$g_pars$grp_post_par),  c("phi", "dfrac", "theta", "sigma2"))

  }

  # simulate data (if necessary)
  sim <- ss$g_pars$sim_func(N = ss$g_pars$N,
                            phi = ss$g_pars$true_phi,
                            dfrac = ss$g_pars$true_dfrac,
                            theta = ss$g_pars$true_theta,
                            sigma2 = ss$g_pars$true_sigma2
  )

  # calculate periodogram
  pgram <- calc_pgram(sim$y)

  # smooth if wanted
  if(ss$f_pars$smooth_pgram){
    pgram$p_spec <- smooth.spline(x = pgram$w_freq, y = pgram$p_spec)$y
  }

  draw_prior_f <- function(n){
    ss$g_pars$draw_prior(n = n,
                         len_phi = length(ss$g_pars$true_phi),
                         len_theta = length(ss$g_pars$true_theta)
    )
  }

  # prior
  log_prior_f <- function(x){
    pars <- pars_to_list(x) # should already be transformed vars
    ss$g_pars$log_prior(t_pphi = pars$phi,
                        t_dfrac = pars$dfrac,
                        t_ptheta = pars$theta,
                        log_sigma2 = pars$sigma2)
  }

  # full log-likelihood
  log_like_f <- function(x){
    pars <- pars_to_list(x) # should already be transformed vars
    ss$g_pars$log_like(y = sim$y,
                       t_pphi = pars$phi,
                       t_dfrac = pars$dfrac,
                       t_ptheta = pars$theta,
                       log_sigma2 = pars$sigma2)
  }

  # surrogate log-likelihood
  log_like_approx_f <- function(x){
    pars <- pars_to_list(x) # should already be transformed vars
    thin <- seq(from = 1, to = length(pgram$w_freq), by = ss$f_pars$whittle_thin)
    ss$g_pars$log_like_approx(
      t_pphi = pars$phi,
      t_dfrac = pars$dfrac,
      t_ptheta = pars$theta,
      log_sigma2 = pars$sigma2,
      pgram = pgram$p_spec[thin],
      w = pgram$w_freq[thin],
      component_wise = F
    )
  }

  # component-wise surrogate log-likelihood
  log_like_approx_c <- function(x, calc_comps){
    pars <- pars_to_list(x) # should already be transformed vars
    thin <- seq(from = 1, to = length(pgram$w_freq), by = 1) # if thinning check makes sense
    ss$g_pars$log_like_approx(
      t_pphi = pars$phi,
      t_dfrac = pars$dfrac,
      t_ptheta = pars$theta,
      log_sigma2 = pars$sigma2,
      pgram = pgram$p_spec[thin],
      w = pgram$w_freq[thin],
      component_wise = T,
      calc_comps = calc_comps
    )
  }

  # choose step scale when no DA
  best_step_scale_f_noda <- function(eta, dist, prob_accept, adjust_D, surrogate_expected_acceptance, surrogate_cost, full_cost, da = F){
    best_step_scale(eta = eta, dist = dist, prob_accept = prob_accept, D = ss$f_pars$bss_D - adjust_D, rho = ss$f_pars$bss_rho, max_T = 200,
                    surrogate_expected_acceptance = surrogate_expected_acceptance, surrogate_cost = surrogate_cost, full_cost = full_cost,
                    model = ss$f_pars$bss_model_noda, da = da)
  }

  # choose step scale with DA
  best_step_scale_f_da <- function(eta, dist, prob_accept, adjust_D, surrogate_expected_acceptance, surrogate_cost, full_cost, da = T){
    best_step_scale(eta = eta, dist = dist, prob_accept = prob_accept, D = ss$f_pars$bss_D - adjust_D, rho = ss$f_pars$bss_rho, max_T = 200,
                    surrogate_expected_acceptance = surrogate_expected_acceptance, surrogate_cost = surrogate_cost, full_cost = full_cost,
                    model = ss$f_pars$bss_model_da, da = da)
  }

 if(verbose) cat("smc approx\n")
  smc_approx <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = T,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = log_prior_f,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = draw_prior_f,
    calibrate_approx_likelihood = NULL,
    find_best_step_scale = best_step_scale_f_noda,
    max_anneal_temp = sfa_max_anneal_temp,
    use_robust_cov = ss$g_pars$use_robust_cov,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose,
    cores = ss$f_pars$cores
  )
  )


  if(verbose) cat("smc + sfa + da\n")
  smc_sfa_da <- with(ss$f_pars, run_smc_da(
    use_da = T, use_approx = F, start_from_approx = T,
    start_from_approx_fit = smc_approx,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = par_start,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = log_prior_f,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    log_like_approx_comps = log_like_approx_c,
    draw_prior = NULL,
    calibrate_approx_likelihood = calibration_method,
    find_best_step_scale = best_step_scale_f_da,
    use_robust_cov = ss$g_pars$use_robust_cov,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose,
    cores = ss$f_pars$cores
  )
  )


  if(verbose) cat("smc standard\n")
  smc_standard <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = F,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = log_prior_f,
    log_like = log_like_f,
    log_like_approx = NULL,
    draw_prior = draw_prior_f,
    calibrate_approx_likelihood = NULL,
    find_best_step_scale = best_step_scale_f_noda,
    use_robust_cov = ss$g_pars$use_robust_cov,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose,
    cores = ss$f_pars$cores
  )
  )


  if(verbose) cat("smc approx full\n")
  smc_approx_full <- with(ss$f_pars, run_smc_da(
    use_da = F, use_approx = T,
    num_p = num_p,
    step_scale_set = step_scale_set,
    par_start = NULL,
    refresh_ejd_threshold = ss$f_pars$bss_D,
    log_prior = log_prior_f,
    log_like = log_like_f,
    log_like_approx = log_like_approx_f,
    draw_prior = draw_prior_f,
    calibrate_approx_likelihood = NULL,
    find_best_step_scale = best_step_scale_f_noda,
    max_anneal_temp = 1,
    use_robust_cov = ss$g_pars$use_robust_cov,
    save_post_interface = ss$g_pars$save_post_interface,
    verbose = verbose,
    cores = ss$f_pars$cores
  )
  )


  return(
    list(sim_settings = ss,
         sim_data = sim,
         smc_sfa_da = smc_sfa_da,
         smc_standard = smc_standard,
         smc_approx = smc_approx,
         smc_approx_full = smc_approx_full
    )
  )

}

####

# Can use try_catch_all to catch warnings from multiple runs:

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


smc_results <- try_catch_all( run_sim(ss = sim_settings[[1]], verbose = T) )


