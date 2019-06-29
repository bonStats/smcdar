#' Jitter particles object with MVN kernel.
#'
#' \code{mvn_jitter} returns particles that have been jittered.
#'
#' @param particles Particles to be jittered.
#' @param step_scale How much to scale the standard deviation by.
#' @param var Variance matrix.
#' @param prec Precision matrix.
#' @param var_chol Cholesky decomposition of variance matrix (upper triangular).
#' @param prec_chol Cholesky decomposition of precision matrix (upper triangular).
#'
#' @export
#'
mvn_jitter <- function(particles, step_scale = 1, var = NULL, prec = NULL, var_chol = NULL, prec_chol = NULL){

  stopifnot( is_particles_obj(particles) )

  # only one to be specified...
  stopifnot( (is.null(var) + is.null(prec) + is.null(var_chol) + is.null(prec_chol) ) == 3)

  stopifnot(length(step_scale) == 1 | length(step_scale) == num_particles(particles))

  if(!is.null(var_chol))  stopifnot( is_upper_triangular_matrix(var_chol) )

  if(!is.null(prec_chol)) stopifnot( is_upper_triangular_matrix(prec_chol) )

  n <- num_particles(particles)
  k <- dim(particles)[2]

  z <- stats::rnorm(n * k)

  if( !is.null(var) |  !is.null(var_chol) ){

    if(is.null(var_chol)) var_chol <- chol(var)

    stopifnot( k == nrow(var_chol) )

    jit <- step_scale * matrix(z, n, k) %*% var_chol

  } else if( !is.null(prec) | !is.null(prec_chol) ){

    if(is.null(prec_chol)) prec_chol <- chol(prec)

    stopifnot( k == nrow(prec_chol) )

    jit <- step_scale * t( backsolve(prec_chol, matrix(z, k, n)) )

  }

  new_prts <- particles
  new_prts[] <- new_prts[] + jit

  return(new_prts)

}
