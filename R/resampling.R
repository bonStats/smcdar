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
  stopifnot( all(x >= 0),
             all(is.finite(x)),
             num_strata > 0,
             length(x) >= num_strata
             )

  N <- length(x)
  M <- as.integer(num_strata)

  x_order <- order(x)
  w <- cumsum(x[x_order])
  w <- w / tail(w,1)

  k <- c( rep(1:M, times = N %/% M), sample(x = 1:M, size = N %% M, replace = T) )
  u <- ( k - 1 + runif(N) ) / M

  smpl <- rowSums(outer(u, w, ">=")) +1

  return(x_order[smpl])

}

#' @export
resample_stratified.particles <- function(x, num_strata){

  sample_index <- resample_stratified.numeric(weights(x), num_strata = num_strata)

  list(
    particles = select_particles(x, index = sample_index, reweight = T),
    index = sample_index
    )

}

#' Generic multinomial resample
#'
#' @param x Object to resample elements of.
#'
#' @return Resampled object
#'
#' @export
resample_multinomial <- function(x){

  UseMethod("resample_multinomial")

}

#' @export
resample_multinomial.numeric <- function(x){

  stopifnot( all(x >= 0), all(is.finite(x)) )

  sample.int(n = length(x), replace = T, prob = x)

}

#' @export
resample_multinomial.particles <- function(x){

  sample_index <- resample_multinomial.numeric(weights(x))

  list(
    particles = select_particles(x, index = sample_index, reweight = T),
    index = sample_index
  )

}

