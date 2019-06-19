#' Apply function to all particles
#'
#' Function must be formatted to take rows from particle matrix.
#'
#' @param particles Particle object.
#' @param fun Function to apply.
#' @param ...
#'
#' @return Function evaluated at each row.
#' @export
#'
papply <- function(particles, fun, ...){

  apply(particles, MARGIN = 1, FUN = fun, ...)

}
