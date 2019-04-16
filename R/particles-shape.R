#' Reshape components of particles to matrix.
#'
#' @param x Object to shape.
#' @param ... Extra arguments.
#'
#' @export
p_shape <- function(x, ...){
  UseMethod("p_shape", x)
}

p_shape.numeric <- function(x, ...){

  array(x, dim = c(length(x),1) )

}

p_shape.matrix <- function(x, ...){

  array(x, dim = dim(x) )

}

p_shape.array <- function(x, ...){

  new_dim <- c(dim(x)[1], prod(dim(x)[-1]))
  array(x, dim = new_dim )

}
