#' Delayed-acceptance MH for symmetric kernel.
#'
#' Uses http://aimsciences.org//article/doi/10.3934/fods.2019005
#'
#' Assumes symmetric or independent kernel.
#'
#' @param new_particles Locations for proposed new particles.
#' @param old_particles Locations for current/old particles.
#' @param var Empirical variance for Mahalanobis distance calculation.
#' @param temp Current SMC temperature.
#' @param loglike Generic log-likelihood interface. Vectorised!
#' @param pre_trans Pre-transformation for approximate log-likelihood.
#' @param time_on Logical. Should computation time be recorded?
#' @param approx_ll_weights Weights for approximate likelihood.
#' @param mvn_step_scale Step scale used for to generate each proposed particle.
#' @param reg_warnings Print regression warnings.
#'
#' @return List.
#' @export
mh_da_step_bglr <- function(new_particles, old_particles, var, temp, loglike, pre_trans = identity, time_on = T, approx_ll_weights = NULL, mvn_step_scale = NULL, reg_warnings = T){

  c_const <- 0.01 # how to choose?
  d_const <- 2
  log_b_const <- ( 1 / ( d_const - 1) ) * log(c_const)

  #trans_new_particles <- t(papply(new_particles, pre_trans))
  #trans_old_particles <- t(papply(old_particles, pre_trans))

  approx_trans_post_new_particles <- loglike(new_particles, lh_trans = pre_trans, type = "approx_posterior", temp = temp, comp_time = T & time_on, weights = approx_ll_weights)
  approx_trans_post_old_particles <- loglike(old_particles, lh_trans = pre_trans, type = "approx_posterior", temp = temp, comp_time = T & time_on, weights = approx_ll_weights)

  # approximate likelihood threshold (surrogate model)
  log_rho_tilde_1 <- approx_trans_post_new_particles - approx_trans_post_old_particles

  # conservative version
  log_rho_1 <- pmin( -log_b_const, pmax( log_b_const, log_rho_tilde_1 ) )

  rho_1 <- exp(log_rho_1)

  accept_1 <- rho_1 > stats::runif(num_particles(new_particles))

  full_post_new_particles <- loglike(new_particles[accept_1,,drop = F], type = "full_likelihood", temp = temp, comp_time = T & time_on)
  full_post_old_particles <- loglike(old_particles[accept_1,,drop = F],  type = "full_likelihood", temp = temp, comp_time = F) # should be memoised.
  approx_trans_llh_new_particles <- loglike(new_particles[accept_1,,drop = F], lh_trans = pre_trans, type = "approx_likelihood", temp = temp, comp_time = T & time_on, weights = approx_ll_weights)
  approx_trans_llh_old_particles <- loglike(old_particles[accept_1,,drop = F], lh_trans = pre_trans, type = "approx_likelihood", temp = temp, comp_time = T & time_on, weights = approx_ll_weights)

  log_rho_tilde_2 <-
    ( full_post_new_particles - full_post_old_particles ) -
    ( approx_trans_llh_new_particles - approx_trans_llh_old_particles)

  log_rho_2 <- log_rho_tilde_1[accept_1] + log_rho_tilde_2 - log_rho_1[accept_1]
  rho_2 <- exp(log_rho_2) # only for accepted particles

  accept_2 <- rho_2 > stats::runif(sum(accept_1))

  accept <- accept_1
  accept[accept_1] <- accept_2

  maha_dist <- papply(new_particles - old_particles, function(x){t(x) %*% solve(var, x)})


  est_true_log_mh <- rep(NA_real_, length = length(accept))
  est_true_log_mh[accept_1] <- full_post_new_particles - full_post_old_particles

  # estimate second stage acceptance probability:
  if(is.null(mvn_step_scale)){
  # linear reg:
    log_mh <- data.frame(true_mh = est_true_log_mh[accept_1],
                     surr_mh = log_rho_tilde_1[accept_1])

    est_ts_md <- stats::lm(true_mh ~ surr_mh, data = log_mh)

    est_true_log_mh[!accept_1] <- predict(est_ts_md, data.frame(surr_mh = log_rho_tilde_1[!accept_1]))

    if(reg_warnings & coef(est_ts_md)["surr_mh"] < 0) warning("Negative slope for surrogate versus full acceptance ratio: May indicate bad surrogate or too few first stage acceptances.")

  } else {
    log_mh <- data.frame(true_mh = est_true_log_mh[accept_1],
                         surr_mh = log_rho_tilde_1[accept_1],
                         eta = mvn_step_scale[accept_1])

    est_ts_md <- stats::lm(true_mh ~ surr_mh + eta, data = log_mh)

    if(reg_warnings & coef(est_ts_md)["surr_mh"] < 0){
      warning("Negative slope for surrogate versus full acceptance ratio: May indicate bad surrogate or too few first stage acceptances.")
      est_ts_md <- stats::lm(true_mh ~ eta, data = log_mh)
    }

    if(reg_warnings & coef(est_ts_md)["eta"] > 0) warning("Positive slope for step scale versus full acceptance ratio: May indicate bad surrogate or too few first stage acceptances.")

    est_true_log_mh[!accept_1] <- predict(est_ts_md, data.frame(surr_mh = log_rho_tilde_1[!accept_1], eta = mvn_step_scale[!accept_1]))

  }

  expected_pre_accept <- pmin(rho_1, 1)
  est_prob_accept <- pmin(exp(est_true_log_mh), 1)

  if(time_on){

    comp_time <- attr(approx_trans_post_new_particles, "comptime") + attr(approx_trans_post_old_particles, "comptime")

    comp_time[accept_1] <- comp_time[accept_1] +
      attr(full_post_new_particles, "comptime") +
      attr(approx_trans_llh_new_particles, "comptime")

    avg_full_like_cost <- mean( attr(full_post_new_particles, "comptime") )
    avg_surr_like_cost <-
      mean(attr(approx_trans_post_new_particles, "comptime")) +
      mean(attr(approx_trans_post_old_particles, "comptime")) +
      mean(attr(approx_trans_llh_new_particles, "comptime"))

    extra_artifical_time <-
      sum(attr(full_post_new_particles, "artiftime")) +
      sum(attr(approx_trans_post_new_particles, "artiftime")) +
      sum(attr(approx_trans_post_old_particles, "artiftime")) +
      sum(attr(approx_trans_llh_new_particles, "artiftime"))


  } else {
    comp_time <- NA
    avg_full_like_cost <- NA
    avg_surr_like_cost <- NA
    extra_artifical_time <- 0
  }

  proposal_dist <- sqrt( maha_dist )

  return(list(
    expected_pre_accept = expected_pre_accept,
    pre_accept = accept_1,
    accept = accept,
    est_prob_accept = est_prob_accept,
    proposal_dist = proposal_dist,
    actual_dist = proposal_dist * accept,
    est_expected_dist = proposal_dist * est_prob_accept,
    comp_time = comp_time,
    avg_full_like_cost = avg_full_like_cost,
    avg_surr_like_cost = avg_surr_like_cost,
    extra_artifical_time = extra_artifical_time
  ))

}
