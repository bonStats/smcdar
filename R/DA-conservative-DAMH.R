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
#' @param loglike Generic log-likelihood interface.
#' @param pre_trans Pre-transformation for approximate log-likelihood.
#' @param time_on Logical. Should computation time be recorded?
#'
#' @return List.
#' @export
mh_da_step_bglr <- function(new_particles, old_particles, var, temp, loglike, pre_trans = identity, time_on = T){

  c_const <- 0.01 # how to choose?
  d_const <- 2
  log_b_const <- ( 1 / ( d_const - 1) ) * log(c_const)

  #trans_new_particles <- t(papply(new_particles, pre_trans))
  #trans_old_particles <- t(papply(old_particles, pre_trans))

  approx_trans_post_new_particles <- papply(new_particles, fun = loglike, lh_trans = pre_trans, type = "approx_posterior", temp = temp, comp_time = T & time_on)
  approx_trans_post_old_particles <- papply(old_particles, fun = loglike, lh_trans = pre_trans, type = "approx_posterior", temp = temp, comp_time = T & time_on)

  # approximate likelihood threshold (surrogate model)
  log_rho_tilde_1 <- approx_trans_post_new_particles - approx_trans_post_old_particles

  # conservative version
  log_rho_1 <- pmin( -log_b_const, pmax( log_b_const, log_rho_tilde_1 ) )

  rho_1 <- exp(log_rho_1)

  accept_1 <- rho_1 > stats::runif(num_particles(new_particles))

  full_post_new_particles <- papply(new_particles[accept_1,,drop = F], fun = loglike, type = "full_likelihood", temp = temp, comp_time = T & time_on)
  full_post_old_particles <- papply(old_particles[accept_1,,drop = F], fun = loglike, type = "full_likelihood", temp = temp, comp_time = F) # should be memoised.
  approx_trans_llh_new_particles <- papply(new_particles[accept_1,,drop = F], fun = loglike, lh_trans = pre_trans, type = "approx_likelihood", temp = temp, comp_time = T & time_on)
  approx_trans_llh_old_particles <- papply(old_particles[accept_1,,drop = F], fun = loglike, lh_trans = pre_trans, type = "approx_likelihood", temp = temp, comp_time = T & time_on)

  log_rho_tilde_2 <-
    ( full_post_new_particles - full_post_old_particles ) -
    ( approx_trans_llh_new_particles - approx_trans_llh_old_particles)

  rho_2 <- exp(log_rho_tilde_1[accept_1] + log_rho_tilde_2 - log_rho_1[accept_1]) # only for accepted particles

  accept_2 <- rho_2 > stats::runif(sum(accept_1))

  accept <- accept_1
  accept[accept_1] <- accept_2

  maha_dist <- papply(new_particles - old_particles, function(x){t(x) %*% solve(var, x)})


  expected_pre_accept <- pmin(rho_1, 1)
  est_prob_accept <- expected_pre_accept
  accept_2_prob <- pmin(rho_2, 1)
  est_prob_accept[accept_1] <- est_prob_accept[accept_1] * accept_2_prob

  est_prob_accept[!accept_1] <- est_prob_accept[!accept_1] * mean(accept_2_prob)

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
