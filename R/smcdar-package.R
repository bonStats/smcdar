#' \code{smcdar}: A package for delayed-acceptance SMC.
#'
#' smcdar...
#'
#' @section \code{smcdar} objects:
#' \itemize{
#'    \item \code{\link{particles}}
#' }
#' @section \code{smcdar} functions:
#' \itemize{
#'    \item \code{\link{mvn_jitter}}
#'    \item \code{\link{num_particles}}
#'    \item \code{\link{p_components}}
#'    \item \code{\link{p_dim}}
#'    \item \code{\link{p_shape}}
#' }
#'
#' @docType package
#' @name smcdar
#'
#' @importFrom stats weights
#'
## usethis namespace: start
#' @useDynLib smcdar, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
#' @importFrom utils tail
NULL
