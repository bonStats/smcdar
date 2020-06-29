#' Robust covariance and correlation estimation from matrix
#'
#' Uses a clustering algorithm to estimate means, then estimate the covariance or correlation matrix.
#' Better for multimodal distributions.
#'
#' @param X A matrix, rows are observations.
#' @param ncluster Number of clusters to use. If missing \code{ClusterR::Optimal_Clusters_KMeans} is used.
#' @param ... Settings for \code{stats::cov.wt}.
#'
#' @return Covariance or correlation matrix.
#' @export
#'
robust_mean_cov <- function(X, ncluster, ...){

  if(missing(ncluster)){
    opt_cluster <- ClusterR::Optimal_Clusters_KMeans(
      X,
      max_clusters = 5,
      criterion = "BIC",
      num_init = 10,
      plot_clusters = F)

    ncluster <- which.min(opt_cluster)
  }

  opt_kmeans <- kmeans(X, centers = ncluster)
  X_centered <- X - opt_kmeans$centers[opt_kmeans$cluster,]

  cov.wt(X_centered, ...)

}
