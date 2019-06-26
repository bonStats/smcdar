#' Generic stratified resample
#'
#' @param x Object to resample elements of.
#' @param num_strata Number of strata.
#'
#' @return Resampled object
#'
#' @export
resample_stratified <- function(x, num_strata){

  UseMethod("resample_stratified")

}

#' @export
resample_stratified.numeric <- function(x, num_strata){
  # may not be numerically stable see Appendix A: Parallel resampling in the particle filter
  stopifnot( all(x >= 0), all(is.finite(x)), num_strata > 0 )

  N <- length(x)
  M <- as.integer(num_strata)

  x_order <- order(x)
  w <- cumsum(x[x_order])
  w <- w / w[N]

  k <- sample.int(n = M, replace = T)
  u <- ( k - 1 + runif(N) ) / M

  out <- rowSums(outer(u, w, ">=")) + 1

  return(out[x_order])

}

#' @export
resample_stratified.particles <- function(x, num_strata){

  sample_index <- resample_stratified.numeric(weights(x), num_strata = num_strata)

  select_particles(x, index = sample_index, reweight = T)

}
