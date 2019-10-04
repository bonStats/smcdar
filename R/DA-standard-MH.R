#' MH for symmetric kernel.
#'
#' Assumes symmetric or independent kernel.
#'
#' @param new_particles Locations for proposed new particles.
#' @param old_particles Locations for current/old particles.
#' @param var Empirical variance for Mahalanobis distance calculation.
#' @param temp Current SMC temperature.
#' @param loglike Generic log-likelihood interface.
#' @param type Which log-likelihood to use. approx or full
#' @param time_on Logical. Should computation time be recorded?
#'
#' @return List.
#' @export
mh_step <- function(new_particles, old_particles, var, temp, loglike, type, time_on = T){
  # using MVN (symmetric kernel)

  new_loglike_type <- papply(new_particles, fun = loglike, temp = temp, type = type, comp_time = T & time_on)
  old_loglike_type <- papply(old_particles, fun = loglike, temp = temp, type = type, comp_time = F)

  accept <- exp(
    new_loglike_type - old_loglike_type
  ) > stats::runif(num_particles(new_particles))

  maha_dist <- papply(new_particles - old_particles, function(x){t(x) %*% solve(var, x)})

  proposal_dist <- sqrt( maha_dist )

  if(time_on){
    avg_full_like_cost <- mean(attr(new_loglike_type, "comptime"))
    comp_time <- attr(new_loglike_type, "comptime")
    extra_artifical_time <- sum(attr(new_loglike_type, "artiftime"))

  } else {
    avg_full_like_cost <- NA
    comp_time <- NA
    extra_artifical_time <- 0
  }

  return(list(
    expected_pre_accept = NA,
    pre_accept = NA,
    accept = accept,
    proposal_dist = proposal_dist,
    actual_dist = proposal_dist * accept,
    comp_time = comp_time,
    avg_full_like_cost = avg_full_like_cost,
    avg_surr_like_cost = NA,
    extra_artifical_time = extra_artifical_time
  ))

}

