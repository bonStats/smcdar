#' Robust covariance and correlation estimation from matrix
#'
#' Uses a clustering algorithm to estimate means, then estimate the covariance or correlation matrix.
#' Better for multimodal distributions. Uses \code{fpc::kmeansruns} to determine number of clusters.
#'
#' @param X A matrix, rows are observations.
#' @param max_ncluster Max number of clusters to consider.
#' @param ... Settings for \code{stats::cov.wt}.
#'
#' @return Covariance or correlation matrix.
#' @export
#'
robust_mean_cov <- function(X, max_ncluster = 5, ...){

  opt_cluster <- fpc::kmeansruns(data = X, krange = 1:max_ncluster, criterion = "ch")

  if(opt_cluster$bestk > 1){
    X <- X - opt_cluster$centers[opt_cluster$cluster,]
  }

  c(cov.wt(X, ...), nclusters = opt_cluster$bestk)

}
