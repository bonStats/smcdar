#' Variance estimator using ancestory
#'
#' @param particles Particles to calculate variance estimate for.
#' @param log_z log of normalising constant estimated by algorithm.
#' @param num_iter Number of SMC iterations.
#'
#' @return estimate of MC variance.
#' @export
#'
eve_var_est <- function(particles, log_z, num_iter){

  N <- num_particles(particles)
  unq_eve <- length(unique(eve(particles)))


  Vn <- 1 - exp(num_iter * (log(N) - log(N-1)) - 2 * log(N) + log(2 * unq_eve))

  exp(2 * log_z) * Vn

}
