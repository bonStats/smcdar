#' Get the dimensions of the components of particles.
#'
#' @param x Object to check dimesion of.
#' @param ... Extra arguments.
#'
#' @export
p_dim <- function(x, ...){
  UseMethod("p_dim", x)
}

p_dim.numeric <- function(x, ...){

  length(x)

}

p_dim.matrix <- function(x, ...){

  dim(x)

}

p_dim.array <- function(x, ...){

  dim(x)

}

p_dim.list <- function(x, ...){

  lapply(x, FUN = function(s) p_dim(s, ...) )

}

p_dim.default <- function(x, ...){

  dim(x)

}
