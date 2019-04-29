#' Generator matrix for Lotka-Volterra continuous-time Markov chain
#'
#' @param theta Numeric vector of length 3.
#' @param y1_min Min for y1, integer, default is 0.
#' @param y1_max Max for y1, integer.
#' @param y2_min Min for y2, integer, default is 0.
#' @param y2_max Max for y2, integer.
#'
#' @return Sparse matrix
#'
#' @export
#'
generator_matrix_lotka_volterra_ctmc <- function(theta, y1_min = 0L, y1_max, y2_min = 0L, y2_max){

  stopifnot(length(theta) == 3)

  sparse_mat_list <- lotka_volterra_generator(theta = theta, y1_min = y1_min, y1_max = y1_max, y2_min = y2_min, y2_max = y2_max)

  Matrix::sparseMatrix(i = sparse_mat_list$rows,
                       j =  sparse_mat_list$cols,
                       x = sparse_mat_list$vals
  )

}

log_lhood_lotka_volterra_ctmc <- function(generator_matrix, y1, y2, times){






}
