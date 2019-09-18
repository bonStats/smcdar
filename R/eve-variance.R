#' Variance estimator using ancestory
#'
#' See: Supplementary material for 'Variance estimation in the particle filter'
#'
#' @param particles Particles to calculate variance estimate for.
#' @param log_z log of normalising constant estimated by algorithm.
#' @param num_iter Number of SMC iterations.
#'
#' @return estimate of MC variance of normalising constant.
#' @export
#'
eve_var_est <- function(particles, log_z, num_iter){

  N <- num_particles(particles)

  Z_sq <- exp(2 * log_z)

  S_0_n <- table(factor(eve(particles), levels = 1:N)) / N

  m_star <- exp(num_iter * ( log(N) - log(N-1) )) * ( Z_sq - sum( S_0_n^2 ) )

  Z_sq - m_star

}
