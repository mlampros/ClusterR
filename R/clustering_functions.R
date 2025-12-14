
utils::globalVariables(c("x", "y"))           # to avoid the following NOTE when package checking takes place --> plot_2d: no visible binding for global variables 'x', 'y'


#' tryCatch function to prevent armadillo errors
#'
#' @keywords internal

tryCatch_GMM <- function(data,
                         gaussian_comps,
                         dist_mode,
                         seed_mode,
                         km_iter,
                         em_iter,
                         verbose,
                         var_floor,
                         seed,
                         full_covariance_matrices) {

  Error = tryCatch(GMM_arma(data,
                            gaussian_comps,
                            dist_mode,
                            seed_mode,
                            km_iter,
                            em_iter,
                            verbose,
                            var_floor,
                            seed,
                            full_covariance_matrices),

                   error = function(e) e)

  if (inherits(Error, "error")) {

    return(list(Error = Error, warning = "probable causes of error: 'warning: gmm_diag::learn(): number of vectors is less than number of gaussians' OR 'warning: gmm_diag::learn(): EM algorithm failed'"))}

  else {

    return(Error)
  }
}


#' Gaussian Mixture Model clustering
#'
#' @param data matrix or data frame
#' @param gaussian_comps the number of gaussian mixture components
#' @param dist_mode the distance used during the seeding of initial means and k-means clustering. One of, \emph{eucl_dist}, \emph{maha_dist}.
#' @param seed_mode how the initial means are seeded prior to running k-means and/or EM algorithms. One of, \emph{static_subset}, \emph{random_subset}, \emph{static_spread}, \emph{random_spread}.
#' @param km_iter the number of iterations of the k-means algorithm
#' @param em_iter the number of iterations of the EM algorithm
#' @param verbose either TRUE or FALSE; enable or disable printing of progress during the k-means and EM algorithms
#' @param var_floor the variance floor (smallest allowed value) for the diagonal covariances
#' @param seed integer value for random number generator (RNG)
#' @param full_covariance_matrices a boolean. If FALSE "diagonal" covariance matrices (i.e. in each covariance matrix, all entries outside the main diagonal are assumed to be zero) otherwise "full" covariance matrices will be returned. Be aware in case of "full" covariance matrices a cube (3-dimensional) rather than a matrix for the output "covariance_matrices" value will be returned.
#' @return a list consisting of the centroids, covariance matrix ( where each row of the matrix represents a diagonal covariance matrix), weights and the log-likelihoods for each gaussian component. In case of Error it returns the error message and the possible causes.
#' @details
#' This function is an R implementation of the 'gmm_diag' class of the Armadillo library. The only exception is that user defined parameter settings are not supported, such as seed_mode = 'keep_existing'.
#' For probabilistic applications, better model parameters are typically learned with dist_mode set to maha_dist.
#' For vector quantisation applications, model parameters should be learned with dist_mode set to eucl_dist, and the number of EM iterations set to zero.
#' In general, a sufficient number of k-means and EM iterations is typically about 10.
#' The number of training samples should be much larger than the number of Gaussians.
#' Seeding the initial means with static_spread and random_spread can be much more time consuming than with static_subset and random_subset.
#' The k-means and EM algorithms will run faster on multi-core machines when OpenMP is enabled in your compiler (eg. -fopenmp in GCC)
#' @references
#' http://arma.sourceforge.net/docs.html
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = as.matrix(dietary_survey_IBS[, -ncol(dietary_survey_IBS)])
#'
#' dat = center_scale(dat)
#'
#' gmm = GMM(dat, 2, "maha_dist", "random_subset", 10, 10)

GMM = function(data,
               gaussian_comps = 1,
               dist_mode = 'eucl_dist',
               seed_mode = 'random_subset',
               km_iter = 10,
               em_iter = 5,
               verbose = FALSE,
               var_floor = 1e-10,
               seed = 1,
               full_covariance_matrices = FALSE) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (gaussian_comps < 1) stop('the number of gaussian mixture components should be greater than 0')
  if (!dist_mode %in% c('eucl_dist', 'maha_dist')) stop("available distance modes are 'eucl_dist' and 'maha_dist'")
  if (!seed_mode %in% c('static_subset','random_subset','static_spread','random_spread')) stop("available seed modes are 'static_subset','random_subset','static_spread' and 'random_spread'")
  if (km_iter < 0 ) stop('the km_iter parameter can not be negative')
  if (em_iter < 0 ) stop('the em_iter parameter can not be negative')
  if (!is.logical(verbose)) stop('the verbose parameter should be either TRUE or FALSE')
  if (var_floor < 0 ) stop('the var_floor parameter can not be negative')
  if (!inherits(full_covariance_matrices, 'logical')) stop('The full_covariance_matrices parameter must be a boolean!')

  flag_non_finite = check_NaN_Inf(data)

  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  res = tryCatch_GMM(data,
                     gaussian_comps,
                     dist_mode,
                     seed_mode,
                     km_iter,
                     em_iter,
                     verbose,
                     var_floor,
                     seed,
                     full_covariance_matrices)

  if ('Error' %in% names(res)) {

    return(res)

  } else {

    structure(list(call = match.call(),
                   centroids = res$centroids,
                   covariance_matrices = res$covariance_matrices,
                   weights = as.vector(res$weights),
                   Log_likelihood = res$Log_likelihood_raw),
              class = c("GMMCluster", 'Gaussian Mixture Models'))
  }
}

#' Prediction function for a Gaussian Mixture Model object
#'
#' @param data matrix or data frame
#' @param CENTROIDS matrix or data frame containing the centroids (means), stored as row vectors
#' @param COVARIANCE matrix or data frame (for diagonal covariance) or 3D array (for full covariance matrices)
#' @param WEIGHTS vector containing the weights
#' @return a list consisting of the log-likelihoods, cluster probabilities and cluster labels.
#' @author Lampros Mouselimis
#' @details
#' This function takes the centroids, covariance matrix and weights from a trained model and returns the log-likelihoods, cluster probabilities and cluster labels for new data.
#' The function handles both diagonal covariance matrices (2D matrix) and full covariance matrices (3D array/cube).
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = as.matrix(dietary_survey_IBS[, -ncol(dietary_survey_IBS)])
#'
#' dat = center_scale(dat)
#'
#' gmm = GMM(dat, 2, "maha_dist", "random_subset", 10, 10)
#'
#' # pr = predict_GMM(dat, gmm$centroids, gmm$covariance_matrices, gmm$weights)
#'
predict_GMM = function(data, CENTROIDS, COVARIANCE, WEIGHTS) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if ('data.frame' %in% class(CENTROIDS)) CENTROIDS = as.matrix(CENTROIDS)
  if (!inherits(CENTROIDS, 'matrix')) stop('CENTROIDS should be either a matrix or a data frame')
  if (!inherits(WEIGHTS, 'numeric') || !is.vector(WEIGHTS))
    stop('WEIGHTS should be a numeric vector')

  flag_non_finite = check_NaN_Inf(data)

  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  # Check if COVARIANCE is a 3D array (full covariance) or 2D matrix (diagonal covariance)
  is_full_covariance = length(dim(COVARIANCE)) == 3

  if (is_full_covariance) {
    # Full covariance matrices - COVARIANCE is a 3D array (cube)
    if (!is.array(COVARIANCE)) stop('COVARIANCE should be a 3D array for full covariance matrices')
    if (dim(COVARIANCE)[1] != ncol(data) || dim(COVARIANCE)[2] != ncol(data) || dim(COVARIANCE)[3] != length(WEIGHTS))
      stop('for full covariance: dim(COVARIANCE)[1] and dim(COVARIANCE)[2] should equal ncol(data), and dim(COVARIANCE)[3] should equal length(WEIGHTS)')
    if (ncol(data) != ncol(CENTROIDS) || length(WEIGHTS) != nrow(CENTROIDS))
      stop('the number of columns of the data and CENTROIDS should match and the number of rows of CENTROIDS should equal the length of the WEIGHTS vector')

    res = predict_MGausDPDF_full(data, CENTROIDS, COVARIANCE, WEIGHTS, eps = 1.0e-8)
  } else {
    # Diagonal covariance matrices - COVARIANCE is a 2D matrix
    if ('data.frame' %in% class(COVARIANCE)) COVARIANCE = as.matrix(COVARIANCE)
    if (!inherits(COVARIANCE, 'matrix')) stop('COVARIANCE should be either a matrix or a data frame for diagonal covariance')
    if (ncol(data) != ncol(CENTROIDS) || ncol(data) != ncol(COVARIANCE) || length(WEIGHTS) != nrow(CENTROIDS) || length(WEIGHTS) != nrow(COVARIANCE))
      stop('the number of columns of the data, CENTROIDS and COVARIANCE should match and the number of rows of the CENTROIDS AND COVARIANCE should be equal to the length of the WEIGHTS vector')

    res = predict_MGausDPDF(data, CENTROIDS, COVARIANCE, WEIGHTS, eps = 1.0e-8)
  }

  # I've added 1 to the output cluster labels to account for the difference in indexing between R and C++
  list(log_likelihood = res$Log_likelihood_raw,
       cluster_proba = res$cluster_proba,
       cluster_labels = as.vector(res$cluster_labels) + 1)
}


#' @rdname predict_GMM
#' @param object,newdata,... arguments for the `predict` generic
#' @export
predict.GMMCluster <- function(object, newdata, ...) {
  predict_GMM(newdata, object$centroids, object$covariance_matrices, object$weights)$cluster_labels
}

#' @export
print.GMMCluster <- function(x, ...) {
  cat("GMM Cluster\n",
      "Call:", deparse(x$call), "\n",
      "Data cols:", ncol(x$centroids), "\n",
      "Centroids:", nrow(x$centroids), "\n")
}

#' tryCatch function to prevent armadillo errors in GMM_arma_AIC_BIC
#'
#' @keywords internal
tryCatch_optimal_clust_GMM <- function(data, max_clusters, dist_mode, seed_mode, km_iter, em_iter, verbose, var_floor, criterion, seed, full_covariance_matrices) {

  Error = tryCatch(GMM_arma_AIC_BIC(data, max_clusters, dist_mode, seed_mode, km_iter, em_iter, verbose, var_floor, criterion, seed, full_covariance_matrices),

                   error = function(e) e)

  if (inherits(Error, "error")) {

    return(list(Error = Error, warning = "probable causes of error: 'warning: gmm_diag::learn(): number of vectors is less than number of gaussians' OR 'warning: gmm_diag::learn(): EM algorithm failed'"))}

  else {

    return(Error)
  }
}



#' Optimal number of Clusters for the gaussian mixture models
#'
#' @param data matrix or data frame
#' @param max_clusters either a numeric value, a contiguous or non-continguous numeric vector specifying the cluster search space
#' @param criterion one of 'AIC' or 'BIC'
#' @param dist_mode the distance used during the seeding of initial means and k-means clustering. One of, \emph{eucl_dist}, \emph{maha_dist}.
#' @param seed_mode how the initial means are seeded prior to running k-means and/or EM algorithms. One of, \emph{static_subset}, \emph{random_subset}, \emph{static_spread}, \emph{random_spread}.
#' @param km_iter the number of iterations of the k-means algorithm
#' @param em_iter the number of iterations of the EM algorithm
#' @param verbose either TRUE or FALSE; enable or disable printing of progress during the k-means and EM algorithms
#' @param var_floor the variance floor (smallest allowed value) for the diagonal covariances
#' @param plot_data either TRUE or FALSE indicating whether the results of the function should be plotted
#' @param seed integer value for random number generator (RNG)
#' @param full_covariance_matrices a boolean. If FALSE "diagonal" covariance matrices (i.e. in each covariance matrix, all entries outside the main diagonal are assumed to be zero) otherwise "full" covariance matrices will be used. Note: when using full covariance matrices, the AIC/BIC calculation accounts for the increased number of parameters.
#' @return a vector with either the AIC or BIC for each iteration. In case of Error it returns the error message and the possible causes.
#' @author Lampros Mouselimis
#' @details
#' \strong{AIC}  : the Akaike information criterion
#'
#' \strong{BIC}  : the Bayesian information criterion
#'
#' In case that the \emph{max_clusters} parameter is a contiguous or non-contiguous vector then plotting is disabled. Therefore, plotting is enabled only if the \emph{max_clusters} parameter is of length 1.
#'
#' When \emph{full_covariance_matrices} is TRUE, the AIC/BIC values will be different from when it is FALSE because full covariance matrices have more free parameters (k*(d + d*(d+1)/2)) compared to diagonal covariance matrices (k*2*d), where k is the number of clusters and d is the number of dimensions.
#'
#' @importFrom grDevices dev.cur
#' @importFrom grDevices dev.off
#' @importFrom graphics plot
#' @importFrom graphics axis
#' @importFrom graphics abline
#' @importFrom graphics text
#'
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' opt_gmm = Optimal_Clusters_GMM(dat, 10, criterion = "AIC", plot_data = FALSE)
#'
#'
#' #----------------------------
#' # non-contiguous search space
#' #----------------------------
#'
#' search_space = c(2,5)
#'
#' opt_gmm = Optimal_Clusters_GMM(dat, search_space, criterion = "AIC", plot_data = FALSE)
#'

Optimal_Clusters_GMM = function(data,
                                max_clusters,
                                criterion = "AIC",
                                dist_mode = 'eucl_dist',
                                seed_mode = 'random_subset',
                                km_iter = 10,
                                em_iter = 5,
                                verbose = FALSE,
                                var_floor = 1e-10,
                                plot_data = TRUE,
                                seed = 1,
                                full_covariance_matrices = FALSE) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!inherits(max_clusters, c('numeric', 'integer'))) stop('max_clusters should be a numeric or integer vector')
  if (length(max_clusters) == 1) {
    if (plot_data && max_clusters < 2) stop('if plot_data is TRUE the max_clusters parameter should be at least 2')
  }
  if (!criterion %in% c("AIC", "BIC")) stop("supported criteria are 'AIC' or 'BIC'")
  if (!dist_mode %in% c('eucl_dist', 'maha_dist')) stop("available distance modes are 'eucl_dist' and 'maha_dist'")
  if (!seed_mode %in% c('static_subset','random_subset','static_spread','random_spread'))
    stop("available seed modes are 'static_subset','random_subset','static_spread' and 'random_spread'")
  if (km_iter < 0 ) stop('the km_iter parameter can not be negative')
  if (em_iter < 0 ) stop('the em_iter parameter can not be negative')
  if (!is.logical(verbose)) stop('the verbose parameter should be either TRUE or FALSE')
  if (var_floor < 0 ) stop('the var_floor parameter can not be negative')
  if (!is.logical(full_covariance_matrices)) stop('The full_covariance_matrices parameter must be a boolean!')

  if (length(max_clusters) != 1) {
    plot_data = FALSE                       # set "plot_data" to FALSE if the "max_clusters" parameter is not of length 1
    if (nrow(data) < max(max_clusters) && verbose) {
      warning("the number of rows of the data should be larger than the maximum value of 'max_clusters'", call. = F)
      cat(" ", '\n')
    }
  }
  else {
    if (nrow(data) < max_clusters && verbose) {
      warning("the number of rows of the data should be larger than 'max_clusters'", call. = F)
      cat(" ", '\n')
    }
  }

  flag_non_finite = check_NaN_Inf(data)
  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  if (length(max_clusters) == 1) {
    pass_vector = 1:max_clusters}
  else {
    pass_vector = max_clusters
  }
  if (0 %in% pass_vector) {
    stop("The 'max_clusters' vector can not include a 0 value !", call. = F)
  }

  gmm = tryCatch_optimal_clust_GMM(data, pass_vector, dist_mode, seed_mode, km_iter, em_iter, verbose, var_floor, criterion, seed, full_covariance_matrices)

  if ('Error' %in% names(gmm)) {
    return(gmm)
  }
  else {
    if (plot_data) {
      if (grDevices::dev.cur() != 1) {
        grDevices::dev.off()                          # reset par()
      }

      vec_out = as.vector(gmm)
      tmp_VAL = as.vector(stats::na.omit(vec_out))

      if (length(which(is.na(vec_out))) > 0) {
        x_dis = (1:length(vec_out))[-which(is.na(vec_out))]
        y_dis = vec_out[-which(is.na(vec_out))]
      }
      else {
        x_dis = 1:length(vec_out)
        y_dis = vec_out
      }

      y_MAX = max(tmp_VAL)
      graphics::plot(x = x_dis, y = y_dis, type = 'l', xlab = 'clusters', ylab = criterion, col = 'blue', lty = 3, axes = FALSE)
      graphics::axis(1, at = seq(1, length(vec_out) , by = 1))
      graphics::axis(2, at = seq(round(min(tmp_VAL) - round(summary(y_MAX)[['Max.']]) / 10), y_MAX + round(summary(y_MAX)[['Max.']]) / 10, by = round((summary(tmp_VAL)['Max.'] - summary(tmp_VAL)['Min.']) / 5)), las = 1, cex.axis = 0.8)
      graphics::abline(h = seq(round(min(tmp_VAL) - round(summary(y_MAX)[['Max.']]) / 10), y_MAX + round(summary(y_MAX)[['Max.']]) / 10, by = round((summary(tmp_VAL)['Max.'] - summary(tmp_VAL)['Min.']) / 5)), v = seq(1, length(vec_out) , by = 1), col = "gray", lty = 3)
      graphics::text(x = 1:length(vec_out), y = vec_out, labels = round(vec_out, 1), cex = 0.8, font = 2)
    }

    res = as.vector(gmm)
    return(res)
  }
}



#' tryCatch function to prevent armadillo errors in KMEANS_arma
#'
#' @keywords internal

tryCatch_KMEANS_arma <- function(data, clusters, n_iter, verbose, seed_mode, CENTROIDS, seed) {

  Error = tryCatch(KMEANS_arma(data, clusters, n_iter, verbose, seed_mode, CENTROIDS, seed),

                   error = function(e) e)

  if (inherits(Error, "error")) {

    return(list(Error = Error, message = Error$message))}

  else if (sum(dim(Error)) == 0) {

    return("warning: kmeans(): number of vectors is less than number of means")}

  else {

    return(Error)
  }
}




#' k-means using the Armadillo library
#'
#' @param data matrix or data frame
#' @param clusters the number of clusters
#' @param n_iter the number of clustering iterations (about 10 is typically sufficient)
#' @param seed_mode how the initial centroids are seeded. One of, \emph{keep_existing}, \emph{static_subset}, \emph{random_subset}, \emph{static_spread}, \emph{random_spread}.
#' @param verbose either TRUE or FALSE, indicating whether progress is printed during clustering
#' @param CENTROIDS a matrix of initial cluster centroids. The rows of the CENTROIDS matrix should be equal to the number of clusters and the columns should be equal to the columns of the data. CENTROIDS should be used in combination with seed_mode 'keep_existing'.
#' @param seed integer value for random number generator (RNG)
#' @return the centroids as a matrix. In case of Error it returns the error message, whereas in case of an empty centroids-matrix it returns a warning-message.
#' @details
#' This function is an R implementation of the 'kmeans' class of the Armadillo library.
#' It is faster than the KMeans_rcpp function but it lacks some features. For more info see the details section of the KMeans_rcpp function.
#' The number of columns should be larger than the number of clusters or CENTROIDS.
#' If the clustering fails, the means matrix is reset and a bool set to false is returned.
#' The clustering will run faster on multi-core machines when OpenMP is enabled in your compiler (eg. -fopenmp in GCC)
#' @references
#' http://arma.sourceforge.net/docs.html
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' km = KMeans_arma(dat, clusters = 2, n_iter = 10, "random_subset")
#'
KMeans_arma = function(data,
                       clusters,
                       n_iter = 10,
                       seed_mode = "random_subset",
                       verbose = FALSE,
                       CENTROIDS = NULL,
                       seed = 1) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!is.numeric(clusters) || length(clusters) != 1 || clusters < 1) stop('clusters should be numeric and greater than 0')
  if (n_iter < 0) stop('the n_iter parameter can not be negative')
  if (!seed_mode %in% c('keep_existing','static_subset','random_subset','static_spread','random_spread'))
    stop("available seed modes are 'keep_existing','static_subset','random_subset','static_spread' and 'random_spread'")
  if ((seed_mode == 'keep_existing' && is.null(CENTROIDS)) || (seed_mode != 'keep_existing' && !is.null(CENTROIDS)))
    stop('the keep_existing seed_mode should be used when CENTROIDS is not NULL')
  if (!is.logical(verbose)) stop('the verbose parameter should be either TRUE or FALSE')
  if (!is.null(CENTROIDS) && (!inherits(CENTROIDS, 'matrix') || nrow(CENTROIDS) != clusters || ncol(CENTROIDS) != ncol(data)))
    stop('CENTROIDS should be a matrix with number of rows equal to the number of clusters and number of columns equal to the number of columns of the data')

  flag_non_finite = check_NaN_Inf(data)
  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  res = tryCatch_KMEANS_arma(data, clusters, n_iter, verbose, seed_mode, CENTROIDS, seed)

  if ('Error' %in% names(res) || is.character(res)) {

    return(res)

  } else {

    ## FIXME: this function currently returns centroids. It should probably
    ## return the same data structure as KMeans_cpp.
    ## return(structure(res, class = c("KMeansCluster", "k-means clustering")))
    return (structure(res, class = "k-means clustering"))
  }
}



#' k-means using RcppArmadillo
#'
#' @param data matrix or data frame
#' @param clusters the number of clusters
#' @param num_init number of times the algorithm will be run with different centroid seeds
#' @param max_iters the maximum number of clustering iterations
#' @param initializer the method of initialization. One of, \emph{optimal_init}, \emph{quantile_init}, \emph{kmeans++} and \emph{random}. See details for more information
#' @param fuzzy either TRUE or FALSE. If TRUE, then prediction probabilities will be calculated using the distance between observations and centroids
#' @param verbose either TRUE or FALSE, indicating whether progress is printed during clustering.
#' @param CENTROIDS a matrix of initial cluster centroids. The rows of the CENTROIDS matrix should be equal to the number of clusters and the columns should be equal to the columns of the data.
#' @param tol a float number. If, in case of an iteration (iteration > 1 and iteration < max_iters) 'tol' is greater than the squared norm of the centroids, then kmeans has converged
#' @param tol_optimal_init tolerance value for the 'optimal_init' initializer. The higher this value is, the far appart from each other the centroids are.
#' @param seed integer value for random number generator (RNG)
#' @return a list with the following attributes: clusters, fuzzy_clusters (if fuzzy = TRUE), centroids, total_SSE, best_initialization, WCSS_per_cluster, obs_per_cluster, between.SS_DIV_total.SS
#' @author Lampros Mouselimis
#' @details
#' This function has the following features in comparison to the KMeans_arma function:
#'
#' Besides optimal_init, quantile_init, random and kmeans++ initilizations one can specify the centroids using the CENTROIDS parameter.
#'
#' The running time and convergence of the algorithm can be adjusted using the num_init, max_iters and tol parameters.
#'
#' If num_init > 1 then KMeans_rcpp returns the attributes of the best initialization using as criterion the within-cluster-sum-of-squared-error.
#'
#'
#' ---------------initializers----------------------
#'
#' \strong{optimal_init}   : this initializer adds rows of the data incrementally, while checking that they do not already exist in the centroid-matrix [ experimental ]
#'
#' \strong{quantile_init}  : initialization of centroids by using the cummulative distance between observations and by removing potential duplicates [ experimental ]
#'
#' \strong{kmeans++}       : kmeans++ initialization. Reference : http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf AND http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work
#'
#' \strong{random}         : random selection of data rows as initial centroids
#'
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' km = KMeans_rcpp(dat, clusters = 2, num_init = 5, max_iters = 100, initializer = 'kmeans++')
#'
KMeans_rcpp = function(data,
                       clusters,
                       num_init = 1,
                       max_iters = 100,
                       initializer = 'kmeans++',
                       fuzzy = FALSE,
                       verbose = FALSE,
                       CENTROIDS = NULL,
                       tol = 1e-4,
                       tol_optimal_init = 0.3,
                       seed = 1) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!is.numeric(clusters) || length(clusters) != 1 || clusters < 1) stop('clusters should be numeric and greater than 0')
  if (num_init < 1) stop('the num_init parameter should be greater than 0')
  if (max_iters < 1) stop('the max_iters parameter should be greater than 0')
  if (!initializer %in% c('kmeans++', 'random', 'optimal_init', 'quantile_init')) stop("available initializer methods are 'kmeans++', 'random', 'optimal_init' and 'quantile_init'")
  if (!is.logical(fuzzy)) stop('the fuzzy parameter should be either TRUE or FALSE')
  if (!is.logical(verbose)) stop('the verbose parameter should be either TRUE or FALSE')
  if (!is.null(CENTROIDS) && (!inherits(CENTROIDS, 'matrix') || nrow(CENTROIDS) != clusters || ncol(CENTROIDS) != ncol(data)))
    stop('CENTROIDS should be a matrix with number of rows equal to the number of clusters and number of columns equal to the number of columns of the data')
  if (tol <= 0.0) stop('tol should be a float number greater than 0.0')
  if (tol_optimal_init <= 0.0) stop('tol_optimal_init should be a float number greater than 0.0')

  flag_non_finite = check_NaN_Inf(data)
  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  res = KMEANS_rcpp(data, clusters, num_init, max_iters, initializer, fuzzy, verbose, CENTROIDS, tol, eps = 1.0e-6, tol_optimal_init, seed)

  if (fuzzy) {

    return(structure(list(call = match.call(),
                          clusters = as.vector(res$clusters + 1),
                          fuzzy_clusters = res$fuzzy_clusters,
                          centroids = res$centers,
                          total_SSE = res$total_SSE,
                          best_initialization = res$best_initialization,
                          WCSS_per_cluster = res$WCSS_per_cluster,
                          obs_per_cluster = res$obs_per_cluster,
                          between.SS_DIV_total.SS = (res$total_SSE - sum(res$WCSS_per_cluster)) / res$total_SSE),
                     class = c("KMeansCluster", "k-means clustering")))

  } else {

    return(structure(list(call = match.call(),
                          clusters = as.vector(res$clusters + 1),
                          centroids = res$centers,
                          total_SSE = res$total_SSE,
                          best_initialization = res$best_initialization,
                          WCSS_per_cluster = res$WCSS_per_cluster,
                          obs_per_cluster = res$obs_per_cluster,
                          between.SS_DIV_total.SS = (res$total_SSE - sum(res$WCSS_per_cluster)) / res$total_SSE),
                     class = c("KMeansCluster", "k-means clustering")))
  }
}



#' Prediction function for the k-means
#'
#' @param data matrix or data frame
#' @param CENTROIDS a matrix of initial cluster centroids. The rows of the CENTROIDS matrix should be equal to the number of clusters and the columns should be equal to the columns of the data.
#' @param threads an integer specifying the number of cores to run in parallel
#' @param fuzzy either TRUE or FALSE. If TRUE, then probabilities for each cluster will be returned based on the distance between observations and centroids.
#' @return a vector (clusters)
#' @author Lampros Mouselimis
#' @details
#' This function takes the data and the output centroids and returns the clusters.
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' km = KMeans_rcpp(dat, clusters = 2, num_init = 5, max_iters = 100, initializer = 'kmeans++')
#'
#' pr = predict_KMeans(dat, km$centroids, threads = 1)
predict_KMeans = function(data, CENTROIDS, threads = 1, fuzzy = FALSE) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!is.matrix(CENTROIDS)) stop("CENTROIDS should be a matrix")
  if (ncol(data) != ncol(CENTROIDS))
    stop('the number of columns of the data should match the number of columns of the CENTROIDS ')
  if (threads < 1) stop('the number of threads should be greater or equal to 1')
  if (!is.logical(fuzzy)) stop('fuzzy should be either TRUE or FALSE')

  flag_non_finite = check_NaN_Inf(data)
  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values", call. = F)

  if (!is.null(class(CENTROIDS))) class(CENTROIDS) = NULL                                                  # set the class of the input 'CENTROIDS' to NULL otherwise the 'duplicated()' function might check it column-wise rather than row-wise
  flag_dups = duplicated(CENTROIDS)
  if (sum(flag_dups) > 0) stop("The 'CENTROIDS' input matrix includes duplicated rows!", call. = F)

  res = validate_centroids(data = data,
                           init_centroids = CENTROIDS,
                           threads = threads,
                           fuzzy = fuzzy,
                           eps = 1.0e-6)
  if (fuzzy) {
    return(res$fuzzy_probs)
  }
  else {
    return(as.vector(res$clusters) + 1)
  }
}

#' @rdname predict_KMeans
#' @param object,newdata,... arguments for the `predict` generic
#' @export
predict.KMeansCluster <- function(object, newdata, fuzzy = FALSE, threads = 1, ...) {

  out <- predict_KMeans(newdata,
                        CENTROIDS = object$centroids,
                        threads = threads,
                        fuzzy = fuzzy)
  return(out)
}

#' @export
print.KMeansCluster <- function(x, ...) {
  WSSE <- sum(x$WCSS_per_cluster)
  BSSE <- x$total_SSE - WSSE
  cat("KMeans Cluster\n",
      "Call:", deparse(x$call), "\n",
      "Data cols:", ncol(x$centroids), "\n",
      "Centroids:", nrow(x$centroids), "\n",
      "BSS/SS:", BSSE/x$total_SSE, "\n",
      "SS:", x$total_SSE, "=", WSSE, "(WSS) +", BSSE, "(BSS)\n")
}



#' Silhouette width based on pre-computed clusters
#'
#' @param data a matrix or a data frame
#' @param clusters a numeric vector which corresponds to the pre-computed clusters (see the example section for more details). The size of the 'clusters' vector must be equal to the number of rows of the input data
#' @return a list object where the first sublist is the 'silhouette summary', the second sublist is the 'silhouette matrix' and the third sublist is the 'global average silhouette' (based on the silhouette values of all observations)
#' @author Lampros Mouselimis
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#' dat = center_scale(dat)
#'
#' clusters = 2
#'
#' # compute k-means
#' km = KMeans_rcpp(dat, clusters = clusters, num_init = 5, max_iters = 100, initializer = 'kmeans++')
#'
#' # compute the silhouette width
#' silh_km = silhouette_of_clusters(data = dat, clusters = km$clusters)
#'
#' # silhouette summary
#' silh_summary = silh_km$silhouette_summary
#'
#' # silhouette matrix (including cluster & dissimilarity)
#' silh_mtrx = silh_km$silhouette_matrix
#'
#' # global average silhouette
#' glob_avg = silh_km$silhouette_global_average
#'

silhouette_of_clusters = function(data, clusters) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('the "data" parameter must be either a matrix or a data frame!')
  if (!inherits(clusters, c('numeric', 'integer'))) stop('the "clusters" parameter must be either of type numeric or of type integer!')
  if (nrow(data) != length(clusters)) stop("I expect that the number of observations of the 'clusters' parameter is equal to the number of rows of the input 'data'!")

  silh_lst = silhouette_clusters(data, clusters)
  colnames(silh_lst[["silhouette_matrix"]]) = c('cluster', 'intra_cluster_dissim', 'silhouette')
  return(silh_lst)
}


#' Optimal number of Clusters for Kmeans or Mini-Batch-Kmeans
#'
#' @param data matrix or data frame
#' @param max_clusters either a numeric value, a contiguous or non-continguous numeric vector specifying the cluster search space
#' @param criterion one of \emph{variance_explained}, \emph{WCSSE}, \emph{dissimilarity}, \emph{silhouette}, \emph{distortion_fK}, \emph{AIC}, \emph{BIC} and \emph{Adjusted_Rsquared}. See details for more information.
#' @param fK_threshold a float number used in the 'distortion_fK' criterion
#' @param num_init number of times the algorithm will be run with different centroid seeds
#' @param max_iters the maximum number of clustering iterations
#' @param initializer the method of initialization. One of, \emph{optimal_init}, \emph{quantile_init}, \emph{kmeans++} and \emph{random}. See details for more information
#' @param tol a float number. If, in case of an iteration (iteration > 1 and iteration < max_iters) 'tol' is greater than the squared norm of the centroids, then kmeans has converged
#' @param plot_clusters either TRUE or FALSE, indicating whether the results of the \emph{Optimal_Clusters_KMeans} function should be plotted
#' @param verbose either TRUE or FALSE, indicating whether progress is printed during clustering
#' @param tol_optimal_init tolerance value for the 'optimal_init' initializer. The higher this value is, the far appart from each other the centroids are.
#' @param seed integer value for random number generator (RNG)
#' @param mini_batch_params either NULL or a list of the following parameters : \emph{batch_size}, \emph{init_fraction}, \emph{early_stop_iter}. If not NULL then the optimal number of clusters will be found based on the Mini-Batch-Kmeans. See the details and examples sections for more information.
#' @return a vector with the results for the specified criterion. If plot_clusters is TRUE then it plots also the results.
#' @author Lampros Mouselimis
#' @details
#' ---------------criteria--------------------------
#'
#' \strong{variance_explained} : the sum of the within-cluster-sum-of-squares-of-all-clusters divided by the total sum of squares
#'
#' \strong{WCSSE}              : the sum of the within-cluster-sum-of-squares-of-all-clusters
#'
#' \strong{dissimilarity}      : the average intra-cluster-dissimilarity of all clusters (the distance metric defaults to euclidean)
#'
#' \strong{silhouette}         : the average silhouette width where first the average per cluster silhouette is computed and then the global average (the distance metric defaults to euclidean). To compute the silhouette width for each cluster separately see the 'silhouette_of_clusters()' function
#'
#' \strong{distortion_fK}      : this criterion is based on the following paper, 'Selection of K in K-means clustering' (https://www.ee.columbia.edu/~dpwe/papers/PhamDN05-kmeans.pdf)
#'
#' \strong{AIC}                : the Akaike information criterion
#'
#' \strong{BIC}                : the Bayesian information criterion
#'
#' \strong{Adjusted_Rsquared}  : the adjusted R^2 statistic
#'
#'
#' ---------------initializers----------------------
#'
#' \strong{optimal_init}   : this initializer adds rows of the data incrementally, while checking that they do not already exist in the centroid-matrix  [ experimental ]
#'
#' \strong{quantile_init}  : initialization of centroids by using the cummulative distance between observations and by removing potential duplicates   [ experimental ]
#'
#' \strong{kmeans++}       : kmeans++ initialization. Reference : http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf AND http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work
#'
#' \strong{random}         : random selection of data rows as initial centroids
#'
#'
#' If the \emph{mini_batch_params} parameter is not NULL then the optimal number of clusters will be found based on the Mini-batch-Kmeans algorithm, otherwise based on the Kmeans. The higher the \emph{init_fraction}
#' parameter is the more close the results between Mini-Batch-Kmeans and Kmeans will be.
#'
#' In case that the \emph{max_clusters} parameter is a contiguous or non-contiguous vector then plotting is disabled. Therefore, plotting is enabled only if the \emph{max_clusters} parameter is of length 1.
#' Moreover, the \emph{distortion_fK} criterion can't be computed if the \emph{max_clusters} parameter is a contiguous or non-continguous vector ( the \emph{distortion_fK} criterion requires consecutive clusters ).
#' The same applies also to the \emph{Adjusted_Rsquared} criterion which returns incorrect output.
#'
#' @references
#'
#' https://www.ee.columbia.edu/~dpwe/papers/PhamDN05-kmeans.pdf
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom graphics par
#' @importFrom graphics mtext
#' @importFrom stats na.omit
#'
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#'
#' #-------
#' # kmeans
#' #-------
#'
#' opt_km = Optimal_Clusters_KMeans(dat, max_clusters = 10, criterion = "distortion_fK",
#'
#'                                  plot_clusters = FALSE)
#'
#' #------------------
#' # mini-batch-kmeans
#' #------------------
#'
#'
#' params_mbkm = list(batch_size = 10, init_fraction = 0.3, early_stop_iter = 10)
#'
#' opt_mbkm = Optimal_Clusters_KMeans(dat, max_clusters = 10, criterion = "distortion_fK",
#'
#'                                    plot_clusters = FALSE, mini_batch_params = params_mbkm)
#'
#'
#' #----------------------------
#' # non-contiguous search space
#' #----------------------------
#'
#' search_space = c(2,5)
#'
#' opt_km = Optimal_Clusters_KMeans(dat, max_clusters = search_space,
#'
#'                                  criterion = "variance_explained",
#'
#'                                  plot_clusters = FALSE)
#'


Optimal_Clusters_KMeans = function(data,
                                   max_clusters,
                                   criterion = "variance_explained",
                                   fK_threshold = 0.85,
                                   num_init = 1,
                                   max_iters = 200,
                                   initializer = 'kmeans++',
                                   tol = 1e-4,
                                   plot_clusters = TRUE,
                                   verbose = FALSE,
                                   tol_optimal_init = 0.3,
                                   seed = 1,
                                   mini_batch_params = NULL) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!inherits(max_clusters, c('numeric', 'integer'))) stop('max_clusters should be a numeric or integer vector')
  if (length(max_clusters) == 1) {
    if (max_clusters < 1) {
      stop('In case that max_clusters is of length 1 it should be greater than 0')
    }
  }
  if (!criterion %in% c('variance_explained', 'WCSSE', 'dissimilarity', 'silhouette', 'distortion_fK', 'AIC', 'BIC', 'Adjusted_Rsquared'))
    stop("available criteria are 'variance_explained', 'WCSSE', 'dissimilarity', 'silhouette', 'distortion_fK', 'AIC', 'BIC' and 'Adjusted_Rsquared'")
  if (num_init < 1) stop('the num_init parameter should be greater than 0')
  if (max_iters < 1) stop('the max_iters parameter should be greater than 0')
  if (!initializer %in% c('kmeans++', 'random', 'optimal_init', 'quantile_init'))
    stop("available initializer methods are 'kmeans++', 'random', 'quantile_init' and 'optimal_init'")
  if (tol <= 0.0) stop('tol should be a float number greater than 0.0')
  if (!is.logical(plot_clusters)) stop('the plot_clusters parameter should be either TRUE or FALSE')
  if (!is.logical(verbose)) stop('the verbose parameter should be either TRUE or FALSE')
  if (tol_optimal_init <= 0.0) stop('tol_optimal_init should be a float number greater than 0.0')
  if (!is.null(mini_batch_params)) {
    if (!all(names(mini_batch_params) %in% c("batch_size", "init_fraction", "early_stop_iter"))) {
      stop("The 'mini_batch_params' parameter should be of type list and valid inputs to the 'mini_batch_params' are: 'batch_size', 'init_fraction' and 'early_stop_iter'!", call. = F)
    }
    if (criterion == "variance_explained") {
      stop("The 'variance_explained' criterion is not supported in case of mini-batch-kmeans (when 'mini_batch_params' is not NULL)!", call. = F)
    }
  }

  if (length(max_clusters) != 1) plot_clusters = FALSE                       # set "plot_clusters" to FALSE if the "max_clusters" parameter is not of length 1

  flag_non_finite = check_NaN_Inf(data)
  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  LEN_CLUST = ITER_CLUST = NA
  if (length(max_clusters) == 1) {
    LEN_CLUST = max_clusters
    ITER_CLUST = 1:max_clusters}
  else {
    LEN_CLUST = length(max_clusters)
    ITER_CLUST = max_clusters
  }

  vec_out = rep(NA, LEN_CLUST)
  if (verbose) { cat("", '\n'); pb = utils::txtProgressBar(min = 1, max = LEN_CLUST, style = 3); cat("", '\n') }

  COUNT = 1
  for (i in ITER_CLUST) {

    if (is.null(mini_batch_params)) {
      km = KMEANS_rcpp(data, i, num_init, max_iters, initializer, FALSE, FALSE, NULL, tol, 1.0e-6, tol_optimal_init, seed)
    }

    else {
      km = MiniBatchKmeans(data, i, mini_batch_params[["batch_size"]], num_init, max_iters, mini_batch_params[["init_fraction"]],
                           initializer, mini_batch_params[["early_stop_iter"]], FALSE, NULL, tol, tol_optimal_init, seed)

      tmp_cent = km$centroids
      km["centroids"] = NULL
      km[["centers"]] = tmp_cent                                           # rename the mini-batch-kmeans centroids-name to match the one of the kmeans algorithm

      if (criterion %in% c("dissimilarity", "silhouette", "BIC")) {        # in these cases call also the 'predict_MBatchKMeans' function to receive the clusters

        km_preds = predict_MBatchKMeans(data, tmp_cent, FALSE)
        km[["clusters"]] = as.vector(km_preds)
      }
    }

    if (criterion == "variance_explained") {
      vec_out[COUNT] = sum(stats::na.omit(as.vector(km$WCSS_per_cluster))) / km$total_SSE
    }

    if (criterion == "WCSSE") {
      vec_out[COUNT] = sum(stats::na.omit(as.vector(km$WCSS_per_cluster)))
    }

    if (criterion == "dissimilarity") {
      eval_km = evaluation_rcpp(data, as.vector(km$clusters), FALSE)
      tmp_dis = mean(stats::na.omit(unlist(lapply(eval_km$INTRA_cluster_dissimilarity, mean))))
      vec_out[COUNT] = tmp_dis
    }

    if (criterion == "silhouette") {

      if (i == 1) {
        vec_out[COUNT] = 0.0
      }
      else {
        silh_out = silhouette_of_clusters(data = data, clusters = as.vector(km$clusters))
        vec_out[COUNT] = silh_out$silhouette_global_average          # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html#sklearn.metrics.silhouette_score
      }
    }

    if (criterion == "distortion_fK") {
      vec_out[COUNT] = sum(stats::na.omit(as.vector(km$WCSS_per_cluster)))
    }

    if (criterion == "AIC") {                             # http://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
      m = ncol(km$centers)
      k = nrow(km$centers)
      D = sum(stats::na.omit(km$WCSS_per_cluster))
      vec_out[COUNT] = D + 2.0 * m * k
    }

    if (criterion == "BIC") {                             # http://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r
      m = ncol(km$centers)
      k = nrow(km$centers)
      n = length(km$clusters)
      D = sum(stats::na.omit(km$WCSS_per_cluster))
      vec_out[COUNT] = D + log(n) * m * k
    }

    if (criterion == 'Adjusted_Rsquared') {
      vec_out[COUNT] = sum(stats::na.omit(km$WCSS_per_cluster))
    }

    if (verbose) { utils::setTxtProgressBar(pb, COUNT) }
    COUNT = COUNT + 1
  }

  if (verbose) { close(pb); cat("", '\n') }

  if (criterion == 'Adjusted_Rsquared') {                 # http://sherrytowers.com/2013/10/24/k-means-clustering/

    if (length(max_clusters) != 1) {
      vec_out = "The 'Adjusted_Rsquared' criterion doesn't return the correct output if the 'max_clusters' parameter is greater than 1"
    }
    else {
      vec_out = 1.0 - (vec_out * (nrow(data) - 1)) / (vec_out[1] * (nrow(data) - ITER_CLUST))
    }
  }

  if (criterion %in% c('variance_explained', 'WCSSE', 'dissimilarity', 'silhouette', 'AIC', 'BIC', 'Adjusted_Rsquared')) {

    if (plot_clusters) {
      tmp_VAL = as.vector(stats::na.omit(vec_out))

      if (length(which(is.na(vec_out))) > 0) {
        x_dis = (1:length(vec_out))[-which(is.na(vec_out))]
        y_dis = vec_out[-which(is.na(vec_out))]
      }
      else {
        x_dis = 1:length(vec_out)
        y_dis = vec_out
      }

      y_MAX = max(tmp_VAL)
      graphics::plot(x = x_dis, y = y_dis, type = 'l', xlab = 'clusters', ylab = criterion, col = 'blue', lty = 3, axes = FALSE)
      graphics::axis(1, at = seq(1, length(vec_out) , by = 1))

      if (criterion == 'silhouette') {
        graphics::axis(2, at = seq(0, y_MAX + 0.05, by = 0.05 ), las = 1, cex.axis = 0.8)
        graphics::abline(h = seq(0.0, max(as.vector(stats::na.omit(vec_out))), 0.05), v = seq(1, length(vec_out) , by = 1), col = "gray", lty = 3)
      }
      else {
        tmp_summary = round(summary(y_MAX)[['Max.']])
        out_max_summary = ifelse(tmp_summary == 0, 1, tmp_summary)
        graphics::axis(2, at = seq(0, y_MAX + out_max_summary / 10, by = out_max_summary / 10), las = 1, cex.axis = 0.8)
        graphics::abline(h = seq(0.0, max(as.vector(stats::na.omit(vec_out))), out_max_summary / 10), v = seq(1, length(vec_out) , by = 1), col = "gray", lty = 3)
      }

      if (criterion %in% c("variance_explained", "Adjusted_Rsquared", "dissimilarity", "silhouette")) {
        graphics::text(x = 1:length(vec_out), y = vec_out, labels = round(vec_out, 2), cex = 0.8, font = 2)
      }
      else {
        graphics::text(x = 1:length(vec_out), y = vec_out, labels = round(vec_out, 1), cex = 0.8, font = 2)
      }
    }
  }

  else {                                                              # "distortion_fK"
    if (length(max_clusters) != 1) {
      fK_vec = "The 'distortion_fK' criterion can not be computed if the length of the 'max_clusters' parameter is greater than 1. See the details for more information!"
    }
    else {
      f_K = opt_clust_fK(vec_out, ncol(data), fK_threshold)
      fK_vec = as.vector(f_K$fK_evaluation)
    }

    if (plot_clusters) {
      if (length(which(is.na(fK_vec))) > 0) {
        x_fk = (1:length(fK_vec))[-which(is.na(fK_vec))]
        y_fk = fK_vec[-which(is.na(fK_vec))]
      }
      else {
        x_fk = 1:length(fK_vec)
        y_fk = fK_vec
      }

      graphics::par(oma = c(0, 2, 0, 0))
      graphics::plot(y_fk, type = 'l', xlab = 'clusters', ylab = 'f(K)', col = 'green', axes = FALSE)
      graphics::axis(1, at = x_fk)
      graphics::axis(2, at = seq(0, max(y_fk) + 0.1, by = round(summary(y_fk)[['Max.']]) / 10), las = 1, cex.axis = 0.8)
      graphics::abline(h = seq(0.0, max(y_fk), round(summary(y_fk)[['Max.']]) / 10), v = seq(1, length(y_fk) , by = 1), col = "gray", lty = 3)
      graphics::abline(h = fK_threshold, col = 'blue', lty = 3)
      graphics::mtext("threshold", side = 2, line = 2, at = fK_threshold, las = 1, cex = 0.9)
      graphics::text(x = x_fk, y = y_fk, labels = round(y_fk,2), cex = 0.8, font = 2)
    }
  }

  if (criterion %in% c('variance_explained', 'WCSSE', 'dissimilarity', 'silhouette', 'AIC', 'BIC', 'Adjusted_Rsquared')) {
    return(vec_out)
  } else {
    return(fK_vec)                                # "distortion_fK"
  }
}




#' Mini-batch-k-means using RcppArmadillo
#'
#' @param data matrix or data frame
#' @param clusters the number of clusters
#' @param batch_size the size of the mini batches
#' @param num_init number of times the algorithm will be run with different centroid seeds
#' @param max_iters the maximum number of clustering iterations
#' @param init_fraction percentage of data to use for the initialization centroids (applies if initializer is \emph{kmeans++} or \emph{optimal_init}). Should be a float number between 0.0 and 1.0.
#' @param initializer the method of initialization. One of, \emph{optimal_init}, \emph{quantile_init}, \emph{kmeans++} and \emph{random}. See details for more information
#' @param early_stop_iter continue that many iterations after calculation of the best within-cluster-sum-of-squared-error
#' @param verbose either TRUE or FALSE, indicating whether progress is printed during clustering
#' @param CENTROIDS a matrix of initial cluster centroids. The rows of the CENTROIDS matrix should be equal to the number of clusters and the columns should be equal to the columns of the data
#' @param tol a float number. If, in case of an iteration (iteration > 1 and iteration < max_iters) 'tol' is greater than the squared norm of the centroids, then kmeans has converged
#' @param tol_optimal_init tolerance value for the 'optimal_init' initializer. The higher this value is, the far appart from each other the centroids are.
#' @param seed integer value for random number generator (RNG)
#' @return a list with the following attributes: centroids, WCSS_per_cluster, best_initialization, iters_per_initialization
#' @author Lampros Mouselimis
#' @details
#' This function performs k-means clustering using mini batches.
#'
#' ---------------initializers----------------------
#'
#' \strong{optimal_init}   : this initializer adds rows of the data incrementally, while checking that they do not already exist in the centroid-matrix   [ experimental ]
#'
#' \strong{quantile_init}  : initialization of centroids by using the cummulative distance between observations and by removing potential duplicates  [ experimental ]
#'
#' \strong{kmeans++}       : kmeans++ initialization. Reference : http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf AND http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work
#'
#' \strong{random}         : random selection of data rows as initial centroids
#'
#' @references
#' http://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf, https://github.com/siddharth-agrawal/Mini-Batch-K-Means
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' MbatchKm = MiniBatchKmeans(dat, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10)
#'


MiniBatchKmeans = function(data,
                           clusters,
                           batch_size = 10,
                           num_init = 1,
                           max_iters = 100,
                           init_fraction = 1.0,
                           initializer = 'kmeans++',
                           early_stop_iter = 10,
                           verbose = FALSE,
                           CENTROIDS = NULL,
                           tol = 1e-4,
                           tol_optimal_init = 0.3,
                           seed = 1) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!is.numeric(clusters) || length(clusters) != 1 || clusters < 1) stop('clusters should be numeric and greater than 0')
  if (batch_size < 1) stop('batch_size should be greater than 0')
  if (num_init < 1) stop('the num_init parameter should be greater than 0')
  if (max_iters < 1) stop('the max_iters parameter should be greater than 0')
  if (init_fraction <= 0.0 || init_fraction > 1.0) stop('init_fraction is a float greater than 0 and less or equal to 1.0')
  if (!initializer %in% c('kmeans++', 'random', 'optimal_init', 'quantile_init')) stop("available initializer methods are 'kmeans++', 'random', 'optimal_init' and 'quantile_init'")
  if (early_stop_iter < 1) stop('early_stop_iter should be greater than 0')
  if (!is.logical(verbose)) stop('the verbose parameter should be either TRUE or FALSE')
  if (!is.null(CENTROIDS) && (!inherits(CENTROIDS, 'matrix') || nrow(CENTROIDS) != clusters || ncol(CENTROIDS) != ncol(data)))
    stop('CENTROIDS should be a matrix with number of rows equal to the number of clusters and number of columns equal to the number of columns of the data')
  if (tol <= 0.0) stop('tol should be a float number greater than 0.0')
  if (tol_optimal_init <= 0.0) stop('tol_optimal_init should be a float number greater than 0.0')

  flag_non_finite = check_NaN_Inf(data)
  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  res = mini_batch_kmeans(data, clusters, batch_size, max_iters, num_init, init_fraction, initializer, early_stop_iter, verbose, CENTROIDS, tol, tol_optimal_init, seed)

  structure(res, class = c("MBatchKMeans", "k-means clustering"))
}




#' Prediction function for Mini-Batch-k-means
#'
#' @param data matrix or data frame
#' @param CENTROIDS a matrix of initial cluster centroids. The rows of the CENTROIDS matrix should be equal to the number of clusters and the columns should equal the columns of the data.
#' @param fuzzy either TRUE or FALSE. If TRUE then prediction probabilities will be calculated using the distance between observations and centroids.
#' @param updated_output either TRUE or FALSE. If TRUE then the 'predict_MBatchKMeans' function will follow the same output object behaviour as the 'predict_KMeans' function (if fuzzy is TRUE it will return probabilities otherwise it will return the hard clusters). This parameter will be removed in version 1.4.0 because this will become the default output format.
#' @return if fuzzy = TRUE the function returns a list with two attributes: a vector with the clusters and a matrix with cluster probabilities. Otherwise, it returns a vector with the clusters.
#' @author Lampros Mouselimis
#' @importFrom lifecycle deprecate_warn is_present
#' @details
#' This function takes the data and the output centroids and returns the clusters.
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' MbatchKm = MiniBatchKmeans(dat, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10)
#'
#' pr = predict_MBatchKMeans(dat, MbatchKm$centroids, fuzzy = FALSE)
#'


predict_MBatchKMeans = function(data, CENTROIDS, fuzzy = FALSE, updated_output = FALSE) {

  if (lifecycle::is_present(fuzzy)) {

    lifecycle::deprecate_warn(
      when = "1.3.0",
      what = "predict_MBatchKMeans()",
      details = "Beginning from version 1.4.0, if the fuzzy parameter is TRUE the function 'predict_MBatchKMeans' will return only the probabilities, whereas currently it also returns the hard clusters"
    )
  }

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!inherits(CENTROIDS, 'matrix')) stop('CENTROIDS should be a matrix')
  if (!(ncol(data) == ncol(CENTROIDS)))
    stop('the number of columns of the data should match the number of columns of the CENTROIDS ')
  if (!is.logical(fuzzy)) stop('fuzzy should be either TRUE or FALSE')
  if (!is.logical(updated_output)) stop('updated_output should be either TRUE or FALSE')

  flag_non_finite = check_NaN_Inf(data)
  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  res = Predict_mini_batch_kmeans(data, CENTROIDS, fuzzy, eps = 1.0e-6)

  if (fuzzy) {
    if (updated_output) {
      return(res$fuzzy_clusters)
    }
    else {
      return(structure(list(clusters = as.vector(res$clusters + 1), fuzzy_clusters = res$fuzzy_clusters), class = "k-means clustering"))
    }
  }
  else {
    return(as.vector(res$clusters + 1))
  }
}


#' @rdname predict_MBatchKMeans
#' @param object,newdata,... arguments for the `predict` generic
#' @export
predict.MBatchKMeans <- function(object, newdata, fuzzy = FALSE, ...) {

  out <- predict_MBatchKMeans(newdata,
                              CENTROIDS = object$centroids,
                              fuzzy = fuzzy)
  if (fuzzy) {
    return(out$fuzzy_clusters)
  }
  else {
    return(out)
  }
}


#' Compute the cost and clusters based on an input dissimilarity matrix and medoids
#'
#' @param data a dissimilarity matrix, where the main diagonal equals 0.0 and the number of rows equals the number of columns
#' @param medoids a vector of output medoids of the 'Cluster_Medoids', 'Clara_Medoids' or any other 'partition around medoids' function
#' @return a list object that includes the cost and the clusters
#' @author Lampros Mouselimis
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#' dat = center_scale(dat)
#'
#' cm = Cluster_Medoids(dat, clusters = 3, distance_metric = 'euclidean', swap_phase = TRUE)
#' res = cost_clusters_from_dissim_medoids(data = cm$dissimilarity_matrix, medoids = cm$medoid_indices)
#'
#' # cm$best_dissimilarity == res$cost
#' # table(cm$clusters, res$clusters)

cost_clusters_from_dissim_medoids = function(data, medoids) {
  if (nrow(data) != ncol(data)) stop("I expect that the input data object is a dissimilarity matrix where the number of rows equals the number of columns!")
  if (!inherits(medoids, c('numeric', 'integer'))) stop("I expect the 'medoids' input object to be a numeric vector!")
  res_clust = cost_clusters_from_dis_meds(dissim_mat = data,
                                          medoids = medoids - 1)                     # adjust the medoids indexing to the Cpp-indexing by subtracting 1
  return(list(cost = res_clust$cost, clusters = as.vector(res_clust$clusters)))
}


#' Partitioning around medoids
#'
#' @param data matrix or data frame. The data parameter can be also a dissimilarity matrix, where the main diagonal equals 0.0 and the number of rows equals the number of columns
#' @param clusters the number of clusters
#' @param distance_metric a string specifying the distance method. One of,  \emph{euclidean},  \emph{manhattan},  \emph{chebyshev},  \emph{canberra},  \emph{braycurtis},  \emph{pearson_correlation},  \emph{simple_matching_coefficient},  \emph{minkowski},  \emph{hamming},  \emph{jaccard_coefficient},  \emph{Rao_coefficient},  \emph{mahalanobis}, \emph{cosine}
#' @param minkowski_p a numeric value specifying the minkowski parameter in case that distance_metric = "minkowski"
#' @param threads an integer specifying the number of cores to run in parallel
#' @param swap_phase either TRUE or FALSE. If TRUE then both phases ('build' and 'swap') will take place. The 'swap_phase' is considered more computationally intensive.
#' @param fuzzy either TRUE or FALSE. If TRUE, then probabilities for each cluster will be returned based on the distance between observations and medoids
#' @param verbose either TRUE or FALSE, indicating whether progress is printed during clustering
#' @param seed `r lifecycle::badge("deprecated")` `seed` (integer value for random number generator (RNG)) is no longer supported and will be removed in version 1.4.0
#' @return a list with the following attributes: medoids, medoid_indices, best_dissimilarity, dissimilarity_matrix, clusters, fuzzy_probs (if fuzzy = TRUE), silhouette_matrix, clustering_stats
#' @author Lampros Mouselimis
#' @details
#' Due to the fact that I didn't have access to the book 'Finding Groups in Data, Kaufman and Rousseeuw, 1990' (which includes the exact algorithm) I implemented the 'Cluster_Medoids' function based on the paper 'Clustering in an Object-Oriented Environment' (see 'References').
#' Therefore, the 'Cluster_Medoids' function is an approximate implementation and not an exact one. Furthermore, in comparison to k-means clustering, the function 'Cluster_Medoids' is more robust, because it minimizes the sum of unsquared dissimilarities. Moreover, it doesn't need initial guesses for the cluster centers.
#' @references
#' Anja Struyf, Mia Hubert, Peter J. Rousseeuw, (Feb. 1997), Clustering in an Object-Oriented Environment, Journal of Statistical Software, Vol 1, Issue 4
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' cm = Cluster_Medoids(dat, clusters = 3, distance_metric = 'euclidean', swap_phase = TRUE)
#'

Cluster_Medoids = function(data, clusters, distance_metric = 'euclidean', minkowski_p = 1.0, threads = 1, swap_phase = TRUE, fuzzy = FALSE, verbose = FALSE, seed = 1) {

  if (lifecycle::is_present(seed)) {

    lifecycle::deprecate_warn(
      when = "1.2.6",
      what = "Cluster_Medoids(seed)",
      details = "The 'seed' parameter will be removed in version 1.4.0"
    )
  }

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame or a dissimilarity matrix with equal number of rows and columns and a diagonal equal to 0.0')
  if (!is.numeric(clusters) || length(clusters) != 1 || clusters < 1) stop('clusters should be numeric and greater than 0')
  if (!distance_metric %in% c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis", "pearson_correlation", "simple_matching_coefficient", "minkowski",
                              "hamming", "jaccard_coefficient", "Rao_coefficient", "mahalanobis", "cosine"))
    stop("the distance_metric should be one of 'euclidean', 'manhattan', 'chebyshev', 'canberra', 'braycurtis', 'pearson_correlation', 'simple_matching_coefficient',
         'minkowski', 'hamming', 'jaccard_coefficient', 'Rao_coefficient', 'mahalanobis', 'cosine'")
  if (distance_metric == 'minkowski' && minkowski_p == 0.0) stop('if distance metric is minkowski then the minkowski_p should be either a positive or a negative number but not 0.0')
  if (threads < 1) stop('threads should be an integer greater than 0')
  if (!is.logical(swap_phase)) stop('swap_phase should be either TRUE or FALSE')
  if (!is.logical(fuzzy)) stop('fuzzy should be either TRUE or FALSE')
  if (!is.logical(verbose)) stop('verbose should be either TRUE or FALSE')

  flag_non_finite = check_NaN_Inf(data)

  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  medoids_mat = ClusterMedoids(data, clusters, distance_metric, minkowski_p, threads, verbose, swap_phase, fuzzy, seed)

  if (clusters > 1) {

    dsm = data.frame(medoids_mat$silhouette_matrix)

    colnames(dsm) = c('clusters', 'neighbor_clusters', 'intra_clust_dissim', 'outer_clust_dissim', 'silhouette_widths', 'diameter', 'separation')

    cs = data.frame(medoids_mat$clustering_stats)

    colnames(cs) = c('clusters', 'number_obs', 'max_dissimilarity', 'average_dissimilarity', 'diameter', 'separation')

  } else {

    dsm = NULL

    cs = NULL
  }

  if (medoids_mat$flag_dissim_mat) {

    tmp_rows = as.vector(medoids_mat$medoids) + 1

  } else {

    tmp_rows = data[as.vector(medoids_mat$medoids) + 1, ]
  }

  return(structure(list(call = match.call(),
                        medoids = tmp_rows,
                        medoid_indices = as.vector(medoids_mat$medoids) + 1,
                        best_dissimilarity = medoids_mat$cost,
                        dissimilarity_matrix = medoids_mat$dissimilarity_matrix,
                        clusters = as.vector(medoids_mat$clusters) + 1,
                        silhouette_matrix = dsm,
                        fuzzy_probs = medoids_mat$fuzzy_probs,
                        clustering_stats = cs,
                        distance_metric = distance_metric),
                   class = c("MedoidsCluster", "cluster medoids silhouette")))
}



#' Clustering large applications
#'
#' @param data matrix or data frame
#' @param clusters the number of clusters
#' @param samples number of samples to draw from the data set
#' @param sample_size fraction of data to draw in each sample iteration. It should be a float number greater than 0.0 and less or equal to 1.0
#' @param distance_metric a string specifying the distance method. One of,  \emph{euclidean},  \emph{manhattan},  \emph{chebyshev},  \emph{canberra},  \emph{braycurtis},  \emph{pearson_correlation},  \emph{simple_matching_coefficient},  \emph{minkowski},  \emph{hamming},  \emph{jaccard_coefficient},  \emph{Rao_coefficient},  \emph{mahalanobis}, \emph{cosine}
#' @param minkowski_p a numeric value specifying the minkowski parameter in case that distance_metric = "minkowski"
#' @param threads an integer specifying the number of cores to run in parallel. Openmp will be utilized to parallelize the number of the different sample draws
#' @param swap_phase either TRUE or FALSE. If TRUE then both phases ('build' and 'swap') will take place. The 'swap_phase' is considered more computationally intensive.
#' @param fuzzy either TRUE or FALSE. If TRUE, then probabilities for each cluster will be returned based on the distance between observations and medoids
#' @param verbose either TRUE or FALSE, indicating whether progress is printed during clustering
#' @param seed integer value for random number generator (RNG)
#' @return a list with the following attributes : medoids, medoid_indices, sample_indices, best_dissimilarity, clusters, fuzzy_probs (if fuzzy = TRUE), clustering_stats, dissimilarity_matrix, silhouette_matrix
#' @author Lampros Mouselimis
#' @details
#' The Clara_Medoids function is implemented in the same way as the 'clara' (clustering large applications) algorithm (Kaufman and Rousseeuw(1990)). In the 'Clara_Medoids'
#' the 'Cluster_Medoids' function will be applied to each sample draw.
#' @references
#' Anja Struyf, Mia Hubert, Peter J. Rousseeuw, (Feb. 1997), Clustering in an Object-Oriented Environment, Journal of Statistical Software, Vol 1, Issue 4
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' clm = Clara_Medoids(dat, clusters = 3, samples = 5, sample_size = 0.2, swap_phase = TRUE)
#'
Clara_Medoids = function(data, clusters, samples, sample_size, distance_metric = "euclidean", minkowski_p = 1.0, threads = 1, swap_phase = TRUE, fuzzy = FALSE, verbose = FALSE, seed = 1) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!is.numeric(clusters) || length(clusters) != 1 || clusters < 1) stop('clusters should be numeric and greater than 0')
  if (!is.numeric(samples) || length(samples) != 1 || samples < 1) stop("samples should be a numeric value greater than 0")
  if (!is.numeric(sample_size) || sample_size <= 0.0 || sample_size > 1.0 ) stop("sample_size should be a numeric value greater than 0.0 and less than or equal to 1.0")
  if (!distance_metric %in% c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis", "pearson_correlation", "simple_matching_coefficient", "minkowski",
                              "hamming", "jaccard_coefficient", "Rao_coefficient", "mahalanobis", "cosine"))
    stop("the distance_metric should be one of 'euclidean', 'manhattan', 'chebyshev', 'canberra', 'braycurtis', 'pearson_correlation', 'simple_matching_coefficient',
         'minkowski', 'hamming', 'jaccard_coefficient', 'Rao_coefficient', 'mahalanobis', 'cosine'")
  if (distance_metric == 'minkowski' && minkowski_p == 0.0) stop('if distance metric is minkowski then the minkowski_p should be either a positive or a negative number but not 0.0')
  if (threads < 1) stop('threads should be an integer greater than 0')
  if (!is.logical(verbose)) stop('verbose should be either TRUE or FALSE')
  if (!is.logical(swap_phase)) stop('swap_phase should be either TRUE or FALSE')
  if (!is.logical(fuzzy)) stop('fuzzy should be either TRUE or FALSE')

  flag_non_finite = check_NaN_Inf(data)

  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  medoids_mat = ClaraMedoids(data, clusters, distance_metric, samples, sample_size, minkowski_p, threads, verbose, swap_phase, fuzzy, seed)

  if (clusters > 1) {

    dsm = data.frame(medoids_mat$silhouette_matrix)

    colnames(dsm) = c('clusters', 'neighbor_clusters', 'intra_clust_dissim', 'outer_clust_dissim', 'silhouette_widths', 'diameter', 'separation')

  } else {

    dsm = NULL
  }

  cs = data.frame(medoids_mat$clustering_stats)

  colnames(cs) = c('clusters', 'number_obs', 'max_dissimilarity', 'average_dissimilarity', 'isolation')

  cs$clusters = cs$clusters + 1

  return(structure(list(call = match.call(),
                        medoids = medoids_mat$medoids,
                        medoid_indices = as.vector(medoids_mat$medoid_indices) + 1,
                        sample_indices = as.vector(medoids_mat$sample_indices) + 1,
                        best_dissimilarity = medoids_mat$best_dissimilarity,
                        clusters = as.vector(medoids_mat$clusters) + 1,
                        silhouette_matrix = dsm,
                        fuzzy_probs = medoids_mat$fuzzy_probs,
                        clustering_stats = cs,
                        dissimilarity_matrix = medoids_mat$dissimilarity_matrix,
                        distance_metric = distance_metric),
                   class = c("MedoidsCluster", "cluster medoids silhouette")))
}



#' Predictions for the Medoid functions
#'
#' @param data matrix or data frame
#' @param MEDOIDS a matrix of initial cluster medoids (data observations). The rows of the MEDOIDS matrix should be equal to the number of clusters and the columns of the MEDOIDS matrix should be equal to the columns of the data.
#' @param distance_metric a string specifying the distance method. One of,  \emph{euclidean},  \emph{manhattan},  \emph{chebyshev},  \emph{canberra},  \emph{braycurtis},  \emph{pearson_correlation},  \emph{simple_matching_coefficient},  \emph{minkowski},  \emph{hamming},  \emph{jaccard_coefficient},  \emph{Rao_coefficient},  \emph{mahalanobis}, \emph{cosine}
#' @param fuzzy either TRUE or FALSE. If TRUE, then probabilities for each cluster will be returned based on the distance between observations and medoids.
#' @param minkowski_p a numeric value specifying the minkowski parameter in case that distance_metric = "minkowski"
#' @param threads an integer specifying the number of cores to run in parallel. Openmp will be utilized to parallelize the number of initializations (num_init)
#' @return a list with the following attributes will be returned : clusters, fuzzy_clusters (if fuzzy = TRUE), dissimilarity.
#' @author Lampros Mouselimis
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat)
#'
#' cm = Cluster_Medoids(dat, clusters = 3, distance_metric = 'euclidean', swap_phase = TRUE)
#'
#' pm = predict_Medoids(dat, MEDOIDS = cm$medoids, 'euclidean', fuzzy = TRUE)
predict_Medoids = function(data, MEDOIDS = NULL, distance_metric = 'euclidean', fuzzy = FALSE, minkowski_p = 1.0, threads = 1) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if ('data.frame' %in% class(MEDOIDS)) MEDOIDS = as.matrix(MEDOIDS)
  if (is.null(MEDOIDS)) stop('the MEDOIDS should be a non-empty matrix or data frame')
  if (ncol(MEDOIDS) != ncol(data)) stop('the MEDOIDS columns should be equal to the number of columns of the data')
  if (!distance_metric %in% c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis", "pearson_correlation", "simple_matching_coefficient", "minkowski",
                              "hamming", "jaccard_coefficient", "Rao_coefficient", "mahalanobis", "cosine"))
    stop("the distance_metric should be one of 'euclidean', 'manhattan', 'chebyshev', 'canberra', 'braycurtis', 'pearson_correlation', 'simple_matching_coefficient',
         'minkowski', 'hamming', 'jaccard_coefficient', 'Rao_coefficient', 'mahalanobis', 'cosine'")
  if (!is.logical(fuzzy)) stop('fuzzy should be either TRUE or FALSE')
  if (distance_metric == 'minkowski' && minkowski_p == 0.0) stop('if distance metric is minkowski then the minkowski_p should be either a positive or a negative number but not 0.0')
  if (threads < 1) stop('threads should be an integer greater than 0')

  flag_non_finite = check_NaN_Inf(data)

  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  res = predict_medoids(data, distance_metric, MEDOIDS, minkowski_p, threads, fuzzy, 1.0e-6)

  structure(
    list(call = match.call(),
         clusters = as.vector(res$clusters) + 1,
         fuzzy_clusters = res$fuzzy_clusters,
         dissimilarity = res$dissimilarity,
         distance_metric = distance_metric),
    class = c("MedoidsCluster", "cluster medoids silhouette"))
}


#' @rdname predict_Medoids
#' @param object,newdata,... arguments for the `predict` generic
#' @export
predict.MedoidsCluster <- function(object, newdata, fuzzy = FALSE, threads = 1, ...) {
  out <- predict_Medoids(newdata,
                         MEDOIDS = object$medoids,
                         distance_metric = object$distance_metric,
                         fuzzy = fuzzy,
                         threads = threads,
                         ...)
  if (fuzzy) {
    return(out$fuzzy_clusters)
  }
  else {
    return(out$clusters)
  }
}

#' @export
print.MedoidsCluster <- function(x, ...) {
  stats <- as.data.frame(t(x$clustering_stats[-1]))
  for (nm in names(stats))
    stats[[nm]] <- format(stats[[nm]], drop0trailing = T)
  stats <- do.call(rbind, list(index = x$medoid_indices, stats))
  rownames(stats) <- paste("  ", rownames(stats))
  colnames(stats) <- x$clustering_stats$clusters
  cat("Medoids Cluster\n",
      "Call:", deparse(x$call), "\n",
      "Data cols:", ncol(x$medoids), "\n",
      "Medoids:", nrow(x$medoids), "\n",
      "Cluster stats:\n")
  print(stats)
}

#' Interactive function for consecutive plots ( using dissimilarities or the silhouette widths of the observations )
#'
#' @keywords internal
function_interactive = function(evaluation_objects, max_clusters, silhouette = FALSE) {

  cat(" ", '\n')

  x = readline("Based on the plot give the number of clusters (greater than 1) that you consider optimal? ")

  x = as.numeric(unlist(strsplit(x, ",")))

  if (x < 2 || x > max_clusters) { stop(paste0("The number of clusters should be at least 2 and at most ", max_clusters)) }      # silhouette and dissimilarity need at least two clusters

  sil_dis_obj = evaluation_objects[[x]]

  Silhouette_Dissimilarity_Plot(sil_dis_obj, silhouette)
}



#' Optimal number of Clusters for the partitioning around Medoids functions
#'
#' @param data matrix or data.frame. If both clara_samples and clara_sample_size equal 0, then the data parameter can be also a dissimilarity matrix, where the main diagonal equals 0.0 and the number of rows equals the number of columns
#' @param max_clusters either a numeric value, a contiguous or non-continguous numeric vector specifying the cluster search space
#' @param distance_metric a string specifying the distance method. One of,  \emph{euclidean},  \emph{manhattan},  \emph{chebyshev},  \emph{canberra},  \emph{braycurtis},  \emph{pearson_correlation},  \emph{simple_matching_coefficient},  \emph{minkowski},  \emph{hamming},  \emph{jaccard_coefficient},  \emph{Rao_coefficient},  \emph{mahalanobis}, \emph{cosine}
#' @param criterion one of 'dissimilarity' or 'silhouette'
#' @param clara_samples number of samples to draw from the data set in case of clustering large applications (clara)
#' @param clara_sample_size fraction of data to draw in each sample iteration in case of clustering large applications (clara). It should be a float number greater than 0.0 and less or equal to 1.0
#' @param minkowski_p a numeric value specifying the minkowski parameter in case that distance_metric = "minkowski"
#' @param swap_phase either TRUE or FALSE. If TRUE then both phases ('build' and 'swap') will take place. The 'swap_phase' is considered more computationally intensive.
#' @param threads an integer specifying the number of cores to run in parallel. Openmp will be utilized to parallelize the number of sample draws
#' @param verbose either TRUE or FALSE, indicating whether progress is printed during clustering
#' @param plot_clusters TRUE or FALSE, indicating whether the iterative results should be plotted. See the details section for more information
#' @param seed integer value for random number generator (RNG)
#' @return a list of length equal to the max_clusters parameter (the first sublist equals NULL, as dissimilarities and silhouette widths can be calculated if the number of clusters > 1). If plot_clusters is TRUE then the function plots also the results.
#' @author Lampros Mouselimis
#' @details
#' In case of plot_clusters = TRUE, the first plot will be either a plot of dissimilarities or both dissimilarities and silhouette widths giving an indication of the optimal number
#' of the clusters. Then, the user will be asked to give an optimal value for the number of the clusters and after that the second plot will appear with either the dissimilarities or the
#' silhouette widths belonging to each cluster.
#'
#' In case that the \emph{max_clusters} parameter is a contiguous or non-contiguous vector then plotting is disabled. Therefore, plotting is enabled only if the \emph{max_clusters} parameter is of length 1.
#'
#' @importFrom graphics legend
#' @importFrom graphics lines
#'
#' @export
#' @examples
#'
#' \dontrun{
#' data(soybean)
#'
#' dat = soybean[, -ncol(soybean)]
#'
#' opt_md = Optimal_Clusters_Medoids(dat, 10, 'jaccard_coefficient', plot_clusters = FALSE)
#'
#'
#' #----------------------------
#' # non-contiguous search space
#' #----------------------------
#'
#' search_space = c(2,5)
#'
#' opt_md = Optimal_Clusters_Medoids(dat, search_space, 'jaccard_coefficient', plot_clusters = FALSE)
#'
#' }


Optimal_Clusters_Medoids = function(data, max_clusters, distance_metric, criterion = "dissimilarity", clara_samples = 0, clara_sample_size = 0.0,

                                    minkowski_p = 1.0, swap_phase = TRUE, threads = 1, verbose = FALSE, plot_clusters = TRUE, seed = 1) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!inherits(max_clusters, c('numeric', 'integer'))) stop('max_clusters should be a numeric or integer vector')
  if (length(max_clusters) == 1) {
    if (max_clusters < 1) {
      stop('In case that max_clusters is of length 1 it should be greater than 0')
    }
  }
  if (!distance_metric %in% c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis", "pearson_correlation", "simple_matching_coefficient", "minkowski",
                              "hamming", "jaccard_coefficient", "Rao_coefficient", "mahalanobis", "cosine"))
    stop("the distance_metric should be one of 'euclidean', 'manhattan', 'chebyshev', 'canberra', 'braycurtis', 'pearson_correlation', 'simple_matching_coefficient',
         'minkowski', 'hamming', 'jaccard_coefficient', 'Rao_coefficient', 'mahalanobis', 'cosine'")
  if (!criterion %in% c("silhouette","dissimilarity")) stop("supported criteria are 'silhouette' and 'dissimilarity'")
  if (distance_metric == 'minkowski' && minkowski_p == 0.0) stop('if distance metric is minkowski then the minkowski_p should be either a positive or a negative number but not 0.0')
  if (!is.logical(swap_phase)) stop('swap_phase should be either TRUE or FALSE')
  if (threads < 1) stop('threads should be an integer greater than 0')
  if (!is.logical(verbose)) stop('verbose should be either TRUE or FALSE')
  if (!is.logical(plot_clusters)) stop('plot_clusters should be either TRUE or FALSE')
  if (clara_samples != 0 && (!is.numeric(clara_samples) || length(clara_samples) != 1 || clara_samples < 1))
    stop("clara_samples should be a numeric value greater than 0")
  if (clara_sample_size != 0.0 && (!is.numeric(clara_sample_size) || clara_sample_size < 0.0 || clara_sample_size > 1.0))
    stop("clara_sample_size should be a numeric value greater than 0.0 and less than or equal to 1.0")
  if ((clara_samples > 0 && clara_sample_size == 0.0) || (clara_samples == 0 && clara_sample_size > 0.0))
    stop("to run clustering for large applications (clara) both 'clara_samples' and 'clara_sample_size' should be greater than 0")
  if (clara_samples > 0 && clara_sample_size > 0.0 && sum(diag(data)) == 0.0 && nrow(data) == ncol(data))
    stop("a dissimilarity matrix is only allowed for the 'Cluster_Medoids' function")

  if (length(max_clusters) != 1) plot_clusters = FALSE                       # set "plot_clusters" to FALSE if the "max_clusters" parameter is not of length 1

  flag_non_finite = check_NaN_Inf(data)

  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  if (length(max_clusters) == 1) {
    pass_vector = 1:max_clusters}
  else {
    pass_vector = max_clusters
  }
  if (0 %in% pass_vector) {
    stop("The 'max_clusters' vector can not include a 0 value !", call. = F)
  }

  inter_bool = ifelse(criterion == "silhouette", TRUE, FALSE)

  if (clara_samples > 0 && clara_sample_size > 0.0) {

    opt_cl =  OptClust(data, pass_vector, distance_metric, TRUE, clara_samples, clara_sample_size, minkowski_p, criterion, threads, swap_phase, verbose, seed)        # Clara_Medoids
  }

  else {

    opt_cl =  OptClust(data, pass_vector, distance_metric, FALSE, clara_samples, clara_sample_size, minkowski_p, criterion, threads, swap_phase, verbose, seed)        # Cluster_Medoids
  }

  if (plot_clusters) {

    if (grDevices::dev.cur() != 1) {

      grDevices::dev.off()                          # reset par()
    }

    if (criterion == "dissimilarity") {

      tmp_dis = rep(NA, max_clusters)

      for (i in 2:max_clusters) { tmp_dis[i] = opt_cl[[i]]$avg_intra_clust_dissimilarity }

      SUM_dis = sum(stats::na.omit(tmp_dis))

      for (i in 2:max_clusters) { tmp_dis[i] = tmp_dis[i] / SUM_dis }           # the dissimilarities are divided by the sum, so that they are in the range 0 to 1

      tmp_VAL = as.vector(stats::na.omit(tmp_dis))

      if (length(which(is.na(tmp_dis))) > 0) {

        x_dis = (1:length(tmp_dis))[-which(is.na(tmp_dis))]

        y_dis = tmp_dis[-which(is.na(tmp_dis))]}

      else {

        x_dis = 1:length(tmp_dis)

        y_dis = tmp_dis
      }

      graphics::plot(x = x_dis, y = y_dis, type = 'l', xlab = 'clusters', ylab = criterion, col = 'red', xaxp = c(1, 10, 10), axes = FALSE)

      graphics::axis(1, at = seq(1, length(tmp_dis) , by = 1))

      graphics::axis(2, at = round(seq(0, max(tmp_VAL) + max(tmp_VAL)/10 , by = 0.01), 2))

      graphics::abline(h = seq(0.0, max(tmp_VAL) + max(tmp_VAL)/10, 0.1), v = seq(1, length(tmp_dis) , by = 1), col = "gray", lty = 3)

      graphics::text(x = x_dis, y = y_dis, labels = round(y_dis, 3), cex = 0.8)

      graphics::legend("topright", legend = 'avg. dissimilarity', col = "red", lty = 1, text.font = 1)


      # consecutive plots [ the 'tryCatch' function is necessary otherwise it raises an error of the type : "Error in plot.new() : figure margins too large" ]

      try_c = tryCatch(function_interactive(opt_cl, max_clusters, inter_bool), error = function(e) e)

      if (inherits(try_c, "error")) {

        warning("The plot can not be created for the specified number of clusters. This means the output data do not fit in the figure (plot) margins.", call. = F)
      }

      else {

        try_c
      }
    }

    if (criterion == "silhouette") {

      tmp_dis = rep(NA, max_clusters)

      for (i in 2:max_clusters) { tmp_dis[i] = opt_cl[[i]]$avg_intra_clust_dissimilarity }

      SUM_dis = sum(stats::na.omit(tmp_dis))

      for (i in 2:max_clusters) { tmp_dis[i] = tmp_dis[i] / SUM_dis }            # the dissimilarities are divided by the sum, so that they are in the range 0 to 1

      if (length(which(is.na(tmp_dis))) > 0) {

        x_dis = (1:length(tmp_dis))[-which(is.na(tmp_dis))]

        y_dis = tmp_dis[-which(is.na(tmp_dis))]}

      else {

        x_dis = 1:length(tmp_dis)

        y_dis = tmp_dis
      }

      tmp_silh = rep(NA, max_clusters)

      for (i in 2:max_clusters) { tmp_silh[i] = opt_cl[[i]]$avg_width_silhouette }

      SUM_sil = sum(stats::na.omit(tmp_silh))

      for (i in 2:max_clusters) { tmp_silh[i] = tmp_silh[i] / SUM_sil }             # the silhoutte widths are divided by the sum, so that they are in the range 0 to 1

      if (length(which(is.na(tmp_silh))) > 0) {

        x_sil = (1:length(tmp_silh))[-which(is.na(tmp_silh))]

        y_sil = tmp_silh[-which(is.na(tmp_silh))]}

      else {

        x_sil = 1:length(tmp_silh)

        y_sil = tmp_silh
      }

      tmp_VAL_ALL = as.vector(stats::na.omit(c(tmp_dis, tmp_silh)))

      y_MIN = min(tmp_VAL_ALL)

      y_MAX = max(tmp_VAL_ALL)

      graphics::plot(x = x_dis, y = y_dis, type = 'l', xlab = 'clusters', ylim = c(y_MIN, y_MAX), col = 'red', ylab = 'dissimilarity -- silhouette', axes = FALSE)

      graphics::axis(1, at = seq(1, length(tmp_dis) , by = 1))

      graphics::axis(2, at = round(seq(y_MIN, y_MAX + y_MAX/10 , by = 0.01), 2))

      graphics::abline(h = seq(0.0, y_MAX + y_MAX/10, 0.05), v = seq(1, length(tmp_dis) , by = 1), col = "gray", lty = 3)

      graphics::text(x = x_dis, y = y_dis, labels = round(y_dis,3), cex = 0.8)

      graphics::lines(x = x_sil[1:length(x_sil)], y = y_sil[1:length(y_sil)], type = 'l', col = 'blue')

      graphics::text(x = x_sil[1:length(x_sil)], y = y_sil[1:length(y_sil)], labels = round(y_sil[1:length(y_sil)],3), cex = 0.8)

      graphics::legend("topright", legend = c('avg. dissimilarity', 'avg. silhouette width'), col = c("red","blue"), lty = 1, text.font = 1)


      # consecutive plots [ the 'tryCatch' function is necessary otherwise it raises an error of the type : "Error in plot.new() : figure margins too large" ]

      try_c = tryCatch(function_interactive(opt_cl, max_clusters, inter_bool), error = function(e) e)

      if (inherits(try_c, "error")) {

        warning("The plot can not be created for the specified number of clusters. This means the output data do not fit in the figure (plot) margins.", call. = F)
      }

      else {

        try_c
      }
    }
  }

  return(opt_cl)
}



#' Plot of silhouette widths or dissimilarities
#'
#' @param evaluation_object the output of either a \emph{Cluster_Medoids} or \emph{Clara_Medoids} function
#' @param silhouette either TRUE or FALSE, indicating whether the silhouette widths or the dissimilarities should be plotted
#' @return TRUE if either the silhouette widths or the dissimilarities are plotted successfully, otherwise FALSE
#' @author Lampros Mouselimis
#' @details
#' This function takes the result-object of the \emph{Cluster_Medoids} or \emph{Clara_Medoids} function and depending on the argument \emph{silhouette} it plots either the dissimilarities or
#' the silhouette widths of the observations belonging to each cluster.
#'
#' @importFrom graphics par
#' @importFrom graphics barplot
#' @importFrom graphics title
#'
#' @export
#' @examples
#'
#' # data(soybean)
#'
#' # dat = soybean[, -ncol(soybean)]
#'
#' # cm = Cluster_Medoids(dat, clusters = 5, distance_metric = 'jaccard_coefficient')
#'
#' # plt_sd = Silhouette_Dissimilarity_Plot(cm, silhouette = TRUE)
#'


Silhouette_Dissimilarity_Plot = function(evaluation_object, silhouette = TRUE) {

  if (!'silhouette_plot' %in% names(evaluation_object)) {

    if (!inherits(evaluation_object, 'cluster medoids silhouette')) {

      stop("the evaluation_object parameter should be the output of a Cluster_Medoids or Clara_Medoids function")
    }
  }

  if (inherits(evaluation_object, 'cluster medoids silhouette')) {

    evaluation_object$silhouette_matrix = as.matrix(evaluation_object$silhouette_matrix)

    evaluation_object = split_rcpp_lst(evaluation_object)
  }

  if (!is.logical(silhouette)) stop('silhouette should be either TRUE or FALSE')

  # default graphics parameter setting

  if (grDevices::dev.cur() != 1) {

    grDevices::dev.off()                          # reset par()
  }

  success_plot_flag = rep(FALSE, length(evaluation_object$list_intra_dissm))

  len_object = length(evaluation_object$list_intra_dissm)

  op <- graphics::par(mfrow = c(len_object, 1),

                      oma = c(2,2,2.5,2) + 0.1,

                      mar = c(2,2,2,2) + 0.1,

                      mgp = c(2.0, 1.0, 0.0))


  # adjust y-axis using the max. number of the cluster obs.

  max_ylim = max(unlist(lapply(evaluation_object$list_intra_dissm, length)))


  # loop to plot the silhouette, dissimilarities

  for (i in 1:length(evaluation_object$list_intra_dissm)) {

    if (silhouette) {

      if (i == 1) {

        graphics::barplot(sort(as.vector(evaluation_object$list_silhouette[[i]]), decreasing = F), width = 2, horiz = T, xlim = c(-1.0, 1.0), ylim = c(0, max_ylim), axes = F)

        if (len_object < 8) {

          graphics::legend("topleft", legend = c(paste0("cluster ", i), paste0('silhouette : ', round(mean(evaluation_object$list_silhouette[[i]]), 3)),

                                       paste0('observations : ', length(evaluation_object$list_silhouette[[i]]))), text.font = 2, cex = 1.0)}

        else {

          graphics::legend("topleft", legend = paste0("cluster ", paste0(i,  paste0( ' , silhouette : ', paste0(round(mean(evaluation_object$list_silhouette[[i]]), 3)),

                                                                           paste0(' , observations : ', length(evaluation_object$list_silhouette[[i]]))))), text.font = 2, cex = 1.0)
        }

        graphics::title(main = paste0("Silhouette plot for the ", paste0(sum(unlist(lapply(evaluation_object$list_silhouette, length))), " data observations")),

              cex.main = 1.5, font.main = 4, col.main = "blue", outer = T)

        graphics::mtext(paste0("average silhouette width : ", round(evaluation_object$avg_width_silhouette, 3)), outer = T, side = 1,

              cex = 0.75, font = 2, col = "blue")}

      else if (i == length(evaluation_object$list_silhouette)) {

        graphics::barplot(sort(as.vector(evaluation_object$list_silhouette[[i]]), decreasing = F), width = 2, horiz = T, xlim = c(-1.0, 1.0), ylim = c(0, max_ylim), axes = T)

        if (len_object < 8) {

          graphics::legend("topleft", legend = c(paste0("cluster ", i), paste0('silhouette : ', round(mean(evaluation_object$list_silhouette[[i]]), 3)),

                                       paste0('observations : ', length(evaluation_object$list_silhouette[[i]]))), text.font = 2, cex = 1.0)}

        else {

          graphics::legend("topleft", legend = paste0("cluster ", paste0(i,  paste0( ' , silhouette : ', paste0(round(mean(evaluation_object$list_silhouette[[i]]), 3)),

                                                                           paste0(' , observations : ', length(evaluation_object$list_silhouette[[i]]))))), text.font = 2, cex = 1.0)
        }
      }

      else {

        graphics::barplot(sort(as.vector(evaluation_object$list_silhouette[[i]]), decreasing = F), width = 2, horiz = T, xlim = c(-1.0, 1.0), ylim = c(0, max_ylim), axes = F)

        if (len_object < 8) {

          graphics::legend("topleft", legend = c(paste0("cluster ", i), paste0('silhouette : ', round(mean(evaluation_object$list_silhouette[[i]]), 3)),

                                       paste0('observations : ', length(evaluation_object$list_silhouette[[i]]))), text.font = 2, cex = 1.0)}

        else {

          graphics::legend("topleft", legend = paste0("cluster ", paste0(i,  paste0( ' , silhouette : ', paste0(round(mean(evaluation_object$list_silhouette[[i]]), 3)),

                                                                           paste0(' , observations : ', length(evaluation_object$list_silhouette[[i]]))))), text.font = 2, cex = 1.0)
        }
      }
    }

    if (!silhouette) {

      # adjust the xlim using max. dissimilarity value

      max_dis = max(unlist(lapply(evaluation_object$list_intra_dissm, mean)))

      round_nearest_half = ceiling(max_dis / 0.5) * 0.5                                          # http://stackoverflow.com/questions/27518497/r-round-to-nearest-half [ modified using ceiling ]


      if (i == 1) {

        graphics::barplot(sort(as.vector(evaluation_object$list_intra_dissm[[i]]), decreasing = F), width = 2, horiz = T, xlim = c(-0.5, round_nearest_half + 0.1),

                ylim = c(0, max_ylim), axes = F)

        if (len_object < 8) {

          graphics::legend("topleft", legend = c(paste0("cluster ", i), paste0('dissimilarity : ', round(mean(evaluation_object$list_intra_dissm[[i]]), 3)),

                                       paste0('observations : ', length(evaluation_object$list_intra_dissm[[i]]))), text.font = 2, cex = 1.0)}

        else {

          graphics::legend("topleft", legend = paste0("cluster ", paste0(i,  paste0( ' , dissimilarity : ', paste0(round(mean(evaluation_object$list_intra_dissm[[i]]), 3)),

                                                                           paste0(' , observations : ', length(evaluation_object$list_intra_dissm[[i]]))))), text.font = 2, cex = 1.0)
        }

        graphics::title(main = paste0("Dissimilarity plot for the ", paste0(sum(unlist(lapply(evaluation_object$list_intra_dissm, length))), " data observations")),

              cex.main = 1.5, font.main = 4, col.main = "blue", outer = T)

        graphics::mtext(paste0("average dissimilarity : ", round(evaluation_object$avg_intra_clust_dissimilarity, 3)), outer = T, side = 1,

              cex = 0.75, font = 2, col = "blue")}

      else if (i == length(evaluation_object$list_intra_dissm)) {

        graphics::barplot(sort(as.vector(evaluation_object$list_intra_dissm[[i]]), decreasing = F), width = 2, horiz = T, xlim = c(-0.5, round_nearest_half + 0.1),

                ylim = c(0, max_ylim), axes = T)

        if (len_object < 8) {

          graphics::legend("topleft", legend = c(paste0("cluster ", i), paste0('dissimilarity : ', round(mean(evaluation_object$list_intra_dissm[[i]]), 3)),

                                       paste0('observations : ', length(evaluation_object$list_intra_dissm[[i]]))), text.font = 2, cex = 1.0)}

        else {

          graphics::legend("topleft", legend = paste0("cluster ", paste0(i,  paste0( ' , dissimilarity : ', paste0(round(mean(evaluation_object$list_intra_dissm[[i]]), 3)),

                                                                           paste0(' , observations : ', length(evaluation_object$list_intra_dissm[[i]]))))), text.font = 2, cex = 1.0)
        }
      }

      else {

        graphics::barplot(sort(as.vector(evaluation_object$list_intra_dissm[[i]]), decreasing = F), width = 2, horiz = T, xlim = c(-0.5, round_nearest_half + 0.1),

                ylim = c(0, max_ylim), axes = F)

        if (len_object < 8) {

          graphics::legend("topleft", legend = c(paste0("cluster ", i), paste0('dissimilarity : ', round(mean(evaluation_object$list_intra_dissm[[i]]), 3)),

                                       paste0('observations : ', length(evaluation_object$list_intra_dissm[[i]]))), text.font = 2, cex = 1.0)}

        else {

          graphics::legend("topleft", legend = paste0("cluster ", paste0(i,  paste0( ' , dissimilarity : ', paste0(round(mean(evaluation_object$list_intra_dissm[[i]]), 3)),

                                                                           paste0(' , observations : ', length(evaluation_object$list_intra_dissm[[i]]))))), text.font = 2, cex = 1.0)
        }
      }
    }

    success_plot_flag[i] = TRUE
  }

  if (sum(success_plot_flag) == length(evaluation_object$list_intra_dissm)) {

    return(T)}

  else {

    return(F)
  }
}



#' 2-dimensional plots
#'
#' @param data a 2-dimensional matrix or data frame
#' @param clusters numeric vector of length equal to the number of rows of the data, which is the result of a clustering method
#' @param centroids_medoids a matrix of centroids or medoids. The rows of the centroids_medoids should be equal to the length of the unique values of the clusters and
#' the columns should be equal to the columns of the data.
#' @return a plot
#' @author Lampros Mouselimis
#' @details
#' This function plots the clusters using 2-dimensional data and medoids or centroids.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_shape_manual scale_size_manual theme element_blank
#'
#' @export
#' @examples
#'
#' # data(dietary_survey_IBS)
#'
#' # dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' # dat = center_scale(dat)
#'
#' # pca_dat = stats::princomp(dat)$scores[, 1:2]
#'
#' # km = KMeans_rcpp(pca_dat, clusters = 2, num_init = 5, max_iters = 100)
#'
#' # plot_2d(pca_dat, km$clusters, km$centroids)
#'



plot_2d = function(data, clusters, centroids_medoids) {

  if (grDevices::dev.cur() != 1) {

    grDevices::dev.off()                          # reset par()
  }

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if ('data.frame' %in% class(centroids_medoids)) centroids_medoids = as.matrix(centroids_medoids)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (is.integer(clusters)) clusters = as.numeric(clusters)
  if (!is.vector(clusters) || !inherits(clusters, "numeric")) stop('The "clusters" parameter has to be a numeric vector!')
  if (!inherits(centroids_medoids, 'matrix') || nrow(centroids_medoids) != length(unique(clusters)) || ncol(centroids_medoids) != ncol(data))
    stop('centroids_medoids should be a matrix with number of rows equal to the unique labels of clusters and number of columns equal to the number of columns of the data')

  flag_non_finite = check_NaN_Inf(data)

  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  if (length(unique(as.vector(clusters))) > 26) stop("valid shape values are from 0 to 25, consider to reduce the number of class-levels")

  if (ncol(data) != 2) stop("the data should be 2-dimensional")


  # match cluster numbers to factor-levels

  levs = as.factor(paste0("cluster ", 1:length(unique(as.vector(clusters)))))[ match(as.vector(clusters), sort(unique(as.vector(clusters)))) ]


  # give to each cluster a different plot-symbol

  df_plot = data.frame(data, clusters = levs)

  colnames(df_plot) = c("x", "y", "clusters")

  add_points = data.frame(centroids_medoids)

  add_points$clusters = as.factor(paste0("centroid or medoid ", 1:nrow(add_points)))

  colnames(add_points) = c("x", "y", "clusters")


  # Change point shapes and colors

  ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = y, group = clusters)) + ggplot2::geom_point(ggplot2::aes(shape = clusters, color = clusters), size = 2.5) +

    ggplot2::geom_point(data = add_points, ggplot2::aes(x = x, y = y, shape = clusters, color = clusters), size = 5) +

    ggplot2::scale_shape_manual(values = 0:(length(unique(as.vector(clusters))) * 2)) +

    ggplot2::scale_size_manual(values = c(5,3,4)) +

    ggplot2::theme(legend.position = "right") + ggplot2::theme(legend.title = ggplot2::element_blank())
}




#' entropy formula (used in external_validation function)
#'
#' @keywords internal

entropy_formula = function(x_vec) {

  vec = rep(NA, length(x_vec))

  for (i in 1:length(x_vec)) {

    if (x_vec[i] == 0.0) {

      vec[i] = 0.0}

    else {

      vec[i] = ((x_vec[i]) * log2(x_vec[i]/sum(x_vec)))
    }
  }

  return(vec)
}




#' external clustering validation
#'
#' @param true_labels a numeric vector of length equal to the length of the clusters vector
#' @param clusters a numeric vector ( the result of a clustering method ) of length equal to the length of the true_labels
#' @param method one of \emph{rand_index},  \emph{adjusted_rand_index},  \emph{jaccard_index},  \emph{fowlkes_Mallows_index},  \emph{mirkin_metric},  \emph{purity},  \emph{entropy},  \emph{nmi} (normalized mutual information), \emph{var_info} (variation of information), and \emph{nvi} (normalized variation of information)
#' @param summary_stats besides the available methods the summary_stats parameter prints also the specificity, sensitivity, precision, recall and F-measure of the clusters
#' @return if summary_stats is FALSE the function returns a float number, otherwise it returns also a summary statistics table
#' @author Lampros Mouselimis
#' @details
#' This function uses external validation methods to evaluate the clustering results
#'
#' @importFrom gmp asNumeric chooseZ as.bigz
#'
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' X = center_scale(dat)
#'
#' km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'kmeans++')
#'
#' res = external_validation(dietary_survey_IBS$class, km$clusters, method = "adjusted_rand_index")
#'


external_validation = function(true_labels, clusters, method = "adjusted_rand_index", summary_stats = FALSE) {

  if (is.integer(true_labels)) true_labels = as.numeric(true_labels)
  if (is.integer(clusters)) clusters = as.numeric(clusters)
  if (!is.vector(true_labels) || !is.numeric(true_labels)) stop('true_labels should be a numeric vector')
  if (!is.vector(clusters) || !is.numeric(clusters)) stop('clusters should be a numeric vector')
  if (length(true_labels) != length(clusters)) stop('the length of the true_labels vector should equal the length of the clusters vector')
  if (!method %in% c('rand_index', 'adjusted_rand_index', 'jaccard_index', 'fowlkes_mallows_index', 'mirkin_metric', 'purity', 'entropy', 'nmi', 'var_info', 'nvi'))
    stop("supported methods are 'rand_index', 'adjusted_rand_index', 'jaccard_index', 'fowlkes_mallows_index', 'mirkin_metric', 'purity', 'entropy', 'nmi', 'var_info', 'nvi'")

  # if (method == 'nmi' || method == 'nvi') {
  #
  #   unq_true = unique(true_labels)
  #
  #   unq_clust = unique(clusters)
  #
  #   if (length(unq_true) == 1 && length(unq_clust) == 1) {                # account for the case where true-labels and clusters perfectly match, SEE comments of issue https://github.com/mlampros/ClusterR/issues/8, https://github.com/scikit-learn/scikit-learn/blob/a24c8b46/sklearn/metrics/cluster/supervised.py#L772
  #
  #     return(1.0)
  #   }
  # }

  tbl = table(clusters, true_labels)

  conv_df = as.data.frame.matrix(tbl)

  # Diagonal = rep(0, ncol(conv_df))
  #
  # for (i in 1:nrow(conv_df)) {
  #
  #   wh_idx = which.max(conv_df[i, ])
  #
  #   if (conv_df[i, wh_idx] > Diagonal[wh_idx]) {
  #
  #     Diagonal[wh_idx] = conv_df[i, wh_idx]
  #   }
  # }

  tp_plus_fp = sum(gmp::asNumeric(gmp::chooseZ(rowSums(conv_df), 2)))

  tp_plus_fn = sum(gmp::asNumeric(gmp::chooseZ(colSums(conv_df), 2)))

  tp = sum(gmp::asNumeric(gmp::chooseZ(as.vector(as.matrix(conv_df)), 2)))

  fp = tp_plus_fp - tp

  fn = tp_plus_fn - tp

  tn = gmp::asNumeric(gmp::chooseZ(sum(as.vector(as.matrix(conv_df))), 2)) - tp - fp - fn

  if (summary_stats || method == "adjusted_rand_index") {

    prod_comb = (tp_plus_fp * tp_plus_fn) / gmp::asNumeric(gmp::chooseZ(length(true_labels), 2))

    mean_comb = (tp_plus_fp + tp_plus_fn) / 2.0
  }

  if (summary_stats || method == 'purity') {

    tmp_pur = apply(conv_df, 1, max)

    res_purity = sum(tmp_pur)/length(true_labels)
  }

  if (summary_stats || method == 'entropy') {

    tmp_entropy = sum(apply(conv_df, 2, function(x) entropy_formula(x)))

    res_entropy = -(1/(sum(tbl) * log2(length(unique(true_labels))))) * tmp_entropy
  }

  if (summary_stats || method == 'nmi' || method == 'var_info' || method == 'nvi') {

    mutual_information = 0.0

    joint_entropy = 0.0

    for (i in 1:nrow(conv_df)) {

      for (j in 1:ncol(conv_df)) {

        if (conv_df[i,j] > 0.0) {

          joint_entropy = joint_entropy + (-((conv_df[i,j] / sum(tbl)) * log2(conv_df[i,j] / sum(tbl))))

          # mutual_information = mutual_information + ((conv_df[i,j] / sum(tbl)) * log2((sum(tbl) * conv_df[i,j]) / (sum(conv_df[i,]) * sum(conv_df[,j]))))       # SEE the comments of issue https://github.com/mlampros/ClusterR/issues/8

          mutual_information = mutual_information + ((conv_df[i,j] / sum(tbl)) * log2(as.numeric(gmp::as.bigz(as.numeric(sum(tbl)) * as.numeric(conv_df[i,j])) / gmp::as.bigz(as.numeric(sum(conv_df[i,])) * as.numeric(sum(conv_df[,j]))))))
        }
      }
    }

    entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))

    entr_class = sum(apply(conv_df, 2, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))

    if (summary_stats || method == 'nmi' || method == 'nvi') {

      unq_true = unique(true_labels)

      unq_clust = unique(clusters)

      if (length(unq_true) == 1 && length(unq_clust) == 1) {                # account for the case where true-labels and clusters perfectly match, SEE comments of issue https://github.com/mlampros/ClusterR/issues/8, https://github.com/scikit-learn/scikit-learn/blob/a24c8b46/sklearn/metrics/cluster/supervised.py#L772

        NMI = 1.0

        NVI = 1.0
      }

      else {

        NMI = (mutual_information / ((entr_cluster + entr_class) / 2.0))

        NVI = 1.0 - (mutual_information / joint_entropy)
      }
    }

    VAR_INFO = (entr_cluster + entr_class) - 2.0 * mutual_information
  }

  if (summary_stats) {

    prec = tp / (tp + fp)
    rec = tp / (tp + fn)

    cat('', '\n')
    cat('----------------------------------------', '\n')
    cat('purity                         :', round(res_purity, 4), '\n')
    cat('entropy                        :', round(res_entropy, 4), '\n')
    cat('normalized mutual information  :', round(NMI, 4), '\n')                      # between 0.0 and 1.0
    cat('variation of information       :', round(VAR_INFO, 4), '\n')                 # the lower the better [ non-negative ]
    cat('normalized var. of information :', round(NVI, 4), '\n')                      # between 0.0 and 1.0; the lower the better
    cat('----------------------------------------', '\n')
    cat('specificity                    :', round(tn / (tn + fp), 4), '\n')
    cat('sensitivity                    :', round(tp / (tp + fn), 4), '\n')
    cat('precision                      :', round(prec, 4), '\n')
    cat('recall                         :', round(rec, 4), '\n')
    cat('F-measure                      :', round(2.0 * ((prec * rec) / (prec + rec)), 4), '\n')
    cat('----------------------------------------', '\n')
    #cat('accuracy                      :', round(sum(Diagonal) / length(true_labels), 4), '\n')
    cat('accuracy OR rand-index         :', round((tp + tn) / (tp + fp + fn + tn), 4), '\n')
    cat('adjusted-rand-index            :', round((tp - prod_comb) / (mean_comb - prod_comb), 4), '\n')
    cat('jaccard-index                  :', round(tp / (tp + fp + fn), 4), '\n')
    cat('fowlkes-mallows-index          :', round(sqrt((tp / ((tp + fp))) * (tp / (tp + fn))), 4), '\n')
    cat('mirkin-metric                  :', round(2.0 * (fp + fn), 4), '\n')
    cat('----------------------------------------', '\n')
  }

  # if (method == "accuracy") {
  #
  #   return(sum(Diagonal) / length(true_labels))
  # }

  if (method == "rand_index") {                                # http://stats.stackexchange.com/questions/89030/rand-index-calculation
    return((tp + tn) / (tp + fp + fn + tn))
  }

  if (method == "adjusted_rand_index") {
    return((tp - prod_comb) / (mean_comb - prod_comb))         # https://github.com/scikit-learn/scikit-learn/blob/51a765a/sklearn/metrics/cluster/supervised.py#L90
  }

  if (method == "jaccard_index") {
    return(tp / (tp + fp + fn))                                # http://www.cs.ucsb.edu/~veronika/MAE/wagner07comparingclusterings.pdf
  }

  if (method == "fowlkes_mallows_index") {
    return(sqrt((tp / ((tp + fp))) * (tp / (tp + fn))))        # https://en.wikipedia.org/wiki/Fowlkes%E2%80%93Mallows_index
  }

  if (method == "mirkin_metric") {
    return(2.0 * (fp + fn))                                    # http://www.cs.ucsb.edu/~veronika/MAE/wagner07comparingclusterings.pdf
  }

  if (method == 'purity') {                                    # http://bioinformatics.oxfordjournals.org/content/23/12/1495.full.pdf+html [ page 1498 ]
    return(res_purity)
  }

  if (method == 'entropy') {                                   # http://bioinformatics.oxfordjournals.org/content/23/12/1495.full.pdf+html [ page 1498 ]
    return(res_entropy)
  }

  if (method == 'nmi') {                                       # http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html, http://stackoverflow.com/questions/35709562/how-to-calculate-clustering-entropy-a-working-example-or-software-code
    return(NMI)
  }

  if (method == 'var_info') {                                  # http://www.stat.washington.edu/mmp/Papers/compare-colt.pdf, http://www.cs.ucsb.edu/~veronika/MAE/wagner07comparingclusterings.pdf
    return(VAR_INFO)
  }

  if (method == 'nvi') {                                       # http://jmlr.csail.mit.edu/papers/volume11/vinh10a/vinh10a.pdf
    return(NVI)
  }
}




#' Function to scale and/or center the data
#'
#' @param data matrix or data frame
#' @param mean_center either TRUE or FALSE. If mean_center is TRUE then the mean of each column will be subtracted
#' @param sd_scale either TRUE or FALSE. See the details section for more information
#' @return a matrix
#' @details
#' If sd_scale is TRUE and mean_center is TRUE then each column will be divided by the standard deviation. If sd_scale is TRUE and mean_center is FALSE then each
#' column will be divided by sqrt( sum(x^2) / (n-1) ).
#' In case of missing values the function raises an error.
#' In case that the standard deviation equals zero then the standard deviation will be replaced with 1.0, so that NaN's can be avoided by division
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = center_scale(dat, mean_center = TRUE, sd_scale = TRUE)
#'


center_scale = function(data, mean_center = TRUE, sd_scale = TRUE) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!is.logical(mean_center)) stop('the mean_center parameter should be either TRUE or FALSE')
  if (!is.logical(sd_scale)) stop('the sd_scale parameter should be either TRUE or FALSE')

  flag_non_finite = check_NaN_Inf(data)

  if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  res = SCALE(data, mean_center, sd_scale)

  return(res)
}




#' Distance matrix calculation
#'
#' @param data matrix or data frame
#' @param method a string specifying the distance method. One of,  \emph{euclidean},  \emph{manhattan},  \emph{chebyshev},  \emph{canberra},  \emph{braycurtis},  \emph{pearson_correlation},  \emph{simple_matching_coefficient},  \emph{minkowski},  \emph{hamming},  \emph{jaccard_coefficient},  \emph{Rao_coefficient},  \emph{mahalanobis}, \emph{cosine}
#' @param upper either TRUE or FALSE specifying if the upper triangle of the distance matrix should be returned. If FALSE then the upper triangle will be filled with NA's
#' @param diagonal either TRUE or FALSE specifying if the diagonal of the distance matrix should be returned. If FALSE then the diagonal will be filled with NA's
#' @param minkowski_p a numeric value specifying the minkowski parameter in case that method = "minkowski"
#' @param threads the number of cores to run in parallel (if OpenMP is available)
#' @return a matrix
#' @export
#' @examples
#'
#' data(dietary_survey_IBS)
#'
#' dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
#'
#' dat = distance_matrix(dat, method = 'euclidean', upper = TRUE, diagonal = TRUE)
#'


distance_matrix = function(data, method = 'euclidean', upper = FALSE, diagonal = FALSE, minkowski_p = 1.0, threads = 1) {

  if ('data.frame' %in% class(data)) data = as.matrix(data)
  if (!inherits(data, 'matrix')) stop('data should be either a matrix or a data frame')
  if (!method %in% c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis", "pearson_correlation", "simple_matching_coefficient", "minkowski",
                     "hamming", "jaccard_coefficient", "Rao_coefficient", "mahalanobis", "cosine"))
    stop("the method should be one of 'euclidean', 'manhattan', 'chebyshev', 'canberra', 'braycurtis', 'pearson_correlation', 'simple_matching_coefficient',
         'minkowski', 'hamming', 'jaccard_coefficient', 'Rao_coefficient', 'mahalanobis', 'cosine'")
  if (!is.logical(upper)) stop('the upper parameter should be either TRUE or FALSE')
  if (!is.logical(diagonal)) stop('the diagonal parameter should be either TRUE or FALSE')
  if (method == 'minkowski' && minkowski_p == 0.0) stop('if distance metric is minkowski then the minkowski_p should be either a positive or a negative number but not 0.0')
  if (threads < 1) stop('the number of threads should be greater than 1')

  #flag_non_finite = check_NaN_Inf(data)                                     # from version 1.0.3 the "distance_matrix" function can accept data with missing values

  #if (!flag_non_finite) stop("the data includes NaN's or +/- Inf values")

  res = dissim_mat(data, method, minkowski_p, upper, diagonal, threads, 1.0e-6)

  return(res)
}





#' Affinity propagation clustering
#'
#' @param data a matrix. Either a similarity matrix (where number of rows equal to number of columns) or a 3-dimensional matrix where the 1st, 2nd and 3rd column correspond to (i-index, j-index, value) triplet of a similarity matrix.
#' @param p a numeric vector of size 1 or size equal to the number of rows of the input matrix. See the details section for more information.
#' @param maxits a numeric value specifying the maximum number of iterations (defaults to 1000)
#' @param convits a numeric value. If the estimated exemplars stay fixed for convits iterations, the affinity propagation algorithm terminates early (defaults to 100)
#' @param dampfact a float number specifying the update equation damping level in [0.5, 1). Higher values correspond to heavy damping, which may be needed if oscillations occur (defaults to 0.9)
#' @param details a boolean specifying if details should be printed in the console
#' @param nonoise a float number. The affinity propagation algorithm adds a small amount of noise to \emph{data} to prevent degenerate cases; this disables that.
#' @param time a boolean. If TRUE then the elapsed time will be printed in the console.
#' @export
#' @details
#'
#' The \emph{affinity propagation} algorithm automatically determines the number of clusters based on the input preference \emph{p}, a real-valued N-vector. p(i) indicates the preference that data point i be
#' chosen as an exemplar. Often a good choice is to set all preferences to median(data). The number of clusters identified can be adjusted by changing this value accordingly. If \emph{p} is a scalar, assumes all
#' preferences are that shared value.
#'
#' The number of clusters eventually emerges by iteratively passing messages between data points to update two matrices, A and R (Frey and Dueck 2007). The "responsibility" matrix R has values r(i, k)
#' that quantify how well suited point k is to serve as the exemplar for point i relative to other candidate exemplars for point i. The "availability" matrix A contains values a(i, k) representing how
#' "appropriate" point k would be as an exemplar for point i, taking into account other points' preferences for point k as an exemplar. Both matrices R and A are initialized with all zeros. The AP
#' algorithm then performs updates iteratively over the two matrices. First, "Responsibilities" r(i, k) are sent from data points to candidate exemplars to indicate how strongly each data point favors
#' the candidate exemplar over other candidate exemplars. "Availabilities" a(i, k) then are sent from candidate exemplars to data points to indicate the degree to which each candidate exemplar is
#' available to be a cluster center for the data point. In this case, the responsibilities and availabilities are messages that provide evidence about whether each data point should be an exemplar and,
#' if not, to what exemplar that data point should be assigned. For each iteration in the message-passing procedure, the sum of r(k; k) + a(k; k) can be used to identify exemplars. After the messages
#' have converged, two ways exist to identify exemplars. In the first approach, for data point i, if r(i, i) + a(i, i) > 0, then data point i is an exemplar. In the second approach, for data point i,
#' if r(i, i) + a(i, i) > r(i, j) + a(i, j) for all i not equal to j, then data point i is an exemplar. The entire procedure terminates after it reaches a predefined number of iterations or if the
#' determined clusters have remained constant for a certain number of iterations... ( https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650075/  -- See chapter 2 )
#'
#' Excluding the main diagonal of the similarity matrix when calculating the median as preference ('p') value can be considered as another option too.
#'
#' @references
#' https://www.psi.toronto.edu/index.php?q=affinity%20propagation
#'
#' https://www.psi.toronto.edu/affinitypropagation/faq.html
#'
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5650075/    ( SEE chapter 2 )
#'
#' @examples
#'
#' set.seed(1)
#' dat = matrix(sample(1:255, 2500, replace = TRUE), 100, 25)
#'
#' smt = 1.0 - distance_matrix(dat, method = 'euclidean', upper = TRUE, diagonal = TRUE)
#' diag(smt) = 0.0
#'
#' ap = AP_affinity_propagation(smt, p = median(as.vector(smt)))
#'
#' str(ap)
#'

AP_affinity_propagation = function(data, p, maxits = 1000, convits = 100, dampfact = 0.9, details = FALSE, nonoise = 0.0, time = FALSE) {

  if (!inherits(data, "matrix")) stop("The 'data' parameter should be a matrix!", call. = F)

  lst_res = affinity_propagation(data, p, maxits, convits, dampfact, details, nonoise, 2.2204e-16, time)

  vec_clust = rep(NA, nrow(data))
  ap_clust = lst_res$clusters
  nams_ap_clust = names(ap_clust)

  for (x in 1:length(ap_clust)) {
    vec_clust[ap_clust[[x]] + 1] = as.integer(nams_ap_clust[x])                  # add 1 to account for the difference in indexing between Rcpp and R
  }

  lst_res[['clusters_vectorized']] = vec_clust

  return(lst_res)
}



#' Affinity propagation preference range
#'
#' @param data a matrix. Either a similarity matrix (where number of rows equal to number of columns) or a 3-dimensional matrix where the 1st, 2nd and 3rd column correspond to (i-index, j-index, value) triplet of a similarity matrix.
#' @param method a character string specifying the preference range method to use. One of 'exact', 'bound'. See the details section for more information.
#' @param threads an integer specifying the number of cores to run in parallel ( applies only if \emph{method} is set to 'exact' which is more computationally intensive )
#' @export
#' @details
#'
#' Given a set of similarities, \emph{data}, this function computes a lower bound, pmin, on the value for the preference where the optimal number of clusters (exemplars) changes from 1 to 2,
#' and the exact value of the preference, pmax, where the optimal number of clusters changes from n-1 to n. For N data points, there may be as many as N^2-N pair-wise similarities (note that
#' the similarity of data point i to k need not be equal to the similarity of data point k to i). These may be passed in an NxN matrix of similarities, \emph{data}, where data(i,k) is the similarity of
#' point i to point k. In fact, only a smaller number of relevant similarities need to be provided, in which case the others are assumed to be -Inf. M similarity values are known, can be passed
#' in an Mx3 matrix \emph{data}, where each row of \emph{data} contains a pair of data point indices and a corresponding similarity value: data(j,3) is the similarity of data point data(j,1) to
#' data point data(j,2).
#'
#' A single-cluster solution may not exist, in which case pmin is set to NaN. The \emph{AP_preferenceRange} uses one of the methods below to compute pmin and pmax:
#'
#' \emph{exact} : Computes the exact values for pmin and pmax (Warning: This can be quite slow)
#' \emph{bound} : Computes the exact value for pmax, but estimates pmin using a bound (default)
#'
#' @references
#' https://www.psi.toronto.edu/affinitypropagation/preferenceRange.m
#'
#' @examples
#'
#' set.seed(1)
#' dat = matrix(sample(1:255, 2500, replace = TRUE), 100, 25)
#'
#' smt = 1.0 - distance_matrix(dat, method = 'euclidean', upper = TRUE, diagonal = TRUE)
#' diag(smt) = 0.0
#'
#' ap_range = AP_preferenceRange(smt, method = "bound")
#'


AP_preferenceRange = function(data, method = "bound", threads = 1) {

  if (!inherits(data, "matrix")) stop("The 'data' parameter should be a matrix!", call. = F)
  if (!method %in% c("bound", "exact")) stop("The 'method' parameter should be either 'bound' or 'exact'!", call. = F)

  return(preferenceRange(data, method, threads))
}

