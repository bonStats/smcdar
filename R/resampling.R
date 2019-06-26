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

# resample_stratified.numeric <- function(x, num_strata){
#   # may not be numerically stable see Appendix A: Parallel resampling in the particle filter
#   stopifnot( all(x >= 0), all(is.finite(x)), num_strata > 0 )
#
#   x_order <- order(x)
#
#   N <- length(x)
#   u <- runif(as.integer(num_strata))
#   w <- cumsum(x[x_order])
#   w <- w / w[N]
#
#   out <- numeric(length = N)
#
#   for(i in 1:N){
#     r <- N * w[i] / w[N]
#     k <- min(N, floor(r) + 1)
#     out[i] <- min(N, floor(r + u[k]) )
#   }
#
#   return(out[x_order])
#
# }
