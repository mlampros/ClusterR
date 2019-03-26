#============================================================

# data

data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

X = center_scale(dat)

# tbl = tibble::as.tibble(X)             # see line 356 [ it works, however I didn't add this test case because tibble has many dependencies -- I guess not sensible for a single test case ]


#=============================================================


context('k-means and mini-batch-k-means')


#############################
# error handling KMeans_arma
#############################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( KMeans_arma(tmp_x, clusters = 2, n_iter = 10, "random_subset", verbose = F) )
})


testthat::test_that("in case that the clusters parameter is not numeric, it returns an error", {
  
  tmp_m = data.frame(1)
  
  testthat::expect_error( KMeans_arma(X, clusters = tmp_m, n_iter = 10, "random_subset", verbose = F) )
})


testthat::test_that("in case that the length of the clusters parameter is not 1, it returns an error", {
  
  tmp_m = c(1,2)
  
  testthat::expect_error( KMeans_arma(X, clusters = tmp_m, n_iter = 10, "random_subset", verbose = F) )
})


testthat::test_that("in case that the clusters parameter is less than 1, it returns an error", {
  
  tmp_m = 0
  
  testthat::expect_error( KMeans_arma(X, clusters = tmp_m, n_iter = 10, "random_subset", verbose = F) )
})


testthat::test_that("in case that the n_iter parameter is less than 0, it returns an error", {
  
  testthat::expect_error( KMeans_arma(X, clusters = 2, n_iter = -1, "random_subset", verbose = F) )
})


testthat::test_that("in case that the seed_mode parameter is invalid, it returns an error", {
  
  testthat::expect_error( KMeans_arma(X, clusters = 2, n_iter = 5, "invalid", verbose = F) )
})


testthat::test_that("in case that the seed_mode parameter equals 'keep_existing' and the CENTROIDS has invalid columns, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X) - 1)), nrow = 2, ncol = ncol(X) - 1)
  
  testthat::expect_error( KMeans_arma(X, clusters = 2, n_iter = 5, "keep_existing", verbose = F, CENTROIDS = cntr) )
})

testthat::test_that("in case that the seed_mode parameter equals 'keep_existing' and the CENTROIDS has invalid rows, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 3, ncol = ncol(X))
  
  testthat::expect_error( KMeans_arma(X, clusters = 2, n_iter = 5, "keep_existing", verbose = F, CENTROIDS = cntr) )
})


testthat::test_that("in case that the seed_mode parameter equals 'keep_existing' and the CENTROIDS is NULL, it returns an error", {
  
  testthat::expect_error( KMeans_arma(X, clusters = 2, n_iter = 5, "keep_existing", verbose = F, CENTROIDS = NULL) )
})


testthat::test_that("in case that the seed_mode parameter does not equal 'keep_existing' and the CENTROIDS is not NULL, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  testthat::expect_error( KMeans_arma(X, clusters = 2, n_iter = 5, "static_subset", verbose = F, CENTROIDS = cntr) )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( KMeans_arma(X, clusters = 2, n_iter = 5, "static_subset", verbose = 'invalid') )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = NaN
  
  testthat::expect_error( KMeans_arma(tmp_dat, clusters = 2, n_iter = 10, "random_subset", verbose = F)  )
})



#######################
# KMeans_arma function
#######################


testthat::test_that("in case that the data is a matrix the result is a CENTROIDS-matrix and the class is 'k-means clustering' ", {
  
  km = KMeans_arma(X, clusters = 2, n_iter = 10, "random_subset", verbose = F)
  
  testthat::expect_true( is.matrix(km) && class(km) == "k-means clustering" )
})


testthat::test_that("in case that the data is a data frame the result is a CENTROIDS-matrix and the class is 'k-means clustering' ", {
  
  km = KMeans_arma(dat, clusters = 2, n_iter = 10, "random_subset", verbose = F)
  
  testthat::expect_true( is.matrix(km) && class(km) == "k-means clustering" )
})


testthat::test_that("it returns a CENTROID-matrix of class 'k-means clustering' when the CENTROIDS parameter is not NULL ", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  res = KMeans_arma(X, clusters = 2, n_iter = 5, "keep_existing", verbose = F, CENTROIDS = cntr)
  
  testthat::expect_true( is.matrix(res) && class(res) == 'k-means clustering' )
})


testthat::test_that("it returns a matrix of class 'k-means clustering' for different seed modes", {
  
  parms = c('static_subset','random_subset','static_spread','random_spread')
  
  res_vec = rep(NA, length(parms))
  
  for (i in 1:length(parms)) {
    
    tmp_km = KMeans_arma(X, clusters = 2, n_iter = 5, parms[i], verbose = F)
    
    res_vec[i] =  (is.matrix(tmp_km) && class(tmp_km) == 'k-means clustering')
  }
  
  testthat::expect_true( sum(res_vec) == length(parms) )
})



#############################
# error handling KMeans_rcpp
#############################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( KMeans_rcpp(tmp_x, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init') )
})


testthat::test_that("in case that the clusters parameter is not numeric, it returns an error", {
  
  tmp_m = data.frame(1)
  
  testthat::expect_error( KMeans_rcpp(X, clusters = tmp_m, num_init = 5, max_iters = 100, initializer = 'optimal_init') )
})


testthat::test_that("in case that the length of the clusters parameter is not 1, it returns an error", {
  
  tmp_m = c(1,2)
  
  testthat::expect_error( KMeans_rcpp(X, clusters = tmp_m, num_init = 5, max_iters = 100, initializer = 'optimal_init') )
})


testthat::test_that("in case that the clusters parameter is less than 1, it returns an error", {
  
  tmp_m = 0
  
  testthat::expect_error( KMeans_rcpp(X, clusters = tmp_m, num_init = 5, max_iters = 100, initializer = 'optimal_init') )
})


testthat::test_that("in case that the num_init parameter is less than 1, it returns an error", {
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 0, max_iters = 100, initializer = 'optimal_init') )
})


testthat::test_that("in case that the max_iters parameter is less than 1, it returns an error", {
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 1, max_iters = 0, initializer = 'optimal_init') )
})


testthat::test_that("in case that the initializer parameter is not one of c('kmeans++', 'random', 'optimal_init', 'quantile_init'), it returns an error", {
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 1, max_iters = 1, initializer = 'invalid') )
})


testthat::test_that("in case that the fuzzy parameter is not logical, it returns an error", {
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 1, max_iters = 1, initializer = 'optimal_init', fuzzy = 'invalid') )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 1, max_iters = 1, initializer = 'optimal_init', verbose = 'invalid') )
})


testthat::test_that("in case that CENTROIDS has invalid columns, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X) - 1)), nrow = 2, ncol = ncol(X) - 1)
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 1, max_iters = 1, initializer = 'optimal_init', CENTROIDS = cntr) )
})


testthat::test_that("in case that CENTROIDS has invalid rows, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 3, ncol = ncol(X))
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 1, max_iters = 1, initializer = 'optimal_init', CENTROIDS = cntr) )
})


testthat::test_that("in case that the tol parameter is less than or equal to 0.0, it returns an error", {
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 1, max_iters = 1, initializer = 'optimal_init', tol = 0.0) )
})


testthat::test_that("in case that the tol_optimal_init parameter is less than or equal to 0.0, it returns an error", {
  
  testthat::expect_error( KMeans_rcpp(X, clusters = 2, num_init = 1, max_iters = 1, initializer = 'optimal_init', tol_optimal_init = 0.0) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = Inf
  
  testthat::expect_error( KMeans_rcpp(tmp_dat, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')  )
})



#######################
# KMeans_rcpp function
#######################


testthat::test_that("in case that the data is a matrix the result is a list and the class is 'k-means clustering' ", {
  
  clust = 2
  
  km = KMeans_rcpp(X, clusters = clust, num_init = 5, max_iters = 100, initializer = 'optimal_init', fuzzy = TRUE)
  
  testthat::expect_true( all(names(km) %in% c("clusters", "fuzzy_clusters", "centroids", "total_SSE", "best_initialization", "WCSS_per_cluster", "obs_per_cluster", "between.SS_DIV_total.SS"))  && is.matrix(km$fuzzy_clusters) &&
                           
                           is.vector(km$clusters) && length(km$clusters) == nrow(X) && is.numeric(km$between.SS_DIV_total.SS) && length(km$between.SS_DIV_total.SS) == 1 && ncol(km$fuzzy_clusters) == clust &&
                           
                           nrow(km$centroids) == clust && ncol(km$centroids) == ncol(X) && length(km$total_SSE) == 1 && is.numeric(km$total_SSE) && length(km$best_initialization) == 1 && 
                           
                           is.numeric(km$best_initialization) && ncol(km$WCSS_per_cluster) == clust && ncol(km$obs_per_cluster) == clust && class(km) == "k-means clustering"  )
})


testthat::test_that("in case that the data is a data frame the result is a list and the class is 'k-means clustering' ", {
  
  clust = 2
  
  km = KMeans_rcpp(dat, clusters = clust, num_init = 5, max_iters = 100, initializer = 'optimal_init')
  
  testthat::expect_true( all(names(km) %in% c("clusters", "centroids", "total_SSE", "best_initialization", "WCSS_per_cluster", "obs_per_cluster", "between.SS_DIV_total.SS"))  
                         
                         && is.vector(km$clusters) && length(km$clusters) == nrow(X) && is.numeric(km$between.SS_DIV_total.SS) && length(km$between.SS_DIV_total.SS) == 1 &&
                           
                           nrow(km$centroids) == clust && ncol(km$centroids) == ncol(X) && length(km$total_SSE) == 1 && is.numeric(km$total_SSE) && length(km$best_initialization) == 1 && 
                           
                           is.numeric(km$best_initialization) && ncol(km$WCSS_per_cluster) == clust && ncol(km$obs_per_cluster) == clust && class(km) == "k-means clustering"  )
})


testthat::test_that("KMeans_rcpp returns the correct output for the initializers", {
  
  clust = 2
  
  res = rep(NA, 4)
  
  count = 1
  
  set.seed(1)
  
  for (i in c('kmeans++', 'random', 'optimal_init', 'quantile_init')) {
    
    km = KMeans_rcpp(X, clusters = clust, num_init = 5, max_iters = 10, initializer = i, tol_optimal_init = 0.2)
    
    res[count] = ( all(names(km) %in% c("clusters", "centroids", "total_SSE", "best_initialization", "WCSS_per_cluster", "obs_per_cluster", "between.SS_DIV_total.SS"))  
                  
                  && is.vector(km$clusters) && length(km$clusters) == nrow(X) && is.numeric(km$between.SS_DIV_total.SS) && length(km$between.SS_DIV_total.SS) == 1 &&
                    
                    nrow(km$centroids) == clust && ncol(km$centroids) == ncol(X) && length(km$total_SSE) == 1 && is.numeric(km$total_SSE) && length(km$best_initialization) == 1 && 
                    
                    is.numeric(km$best_initialization) && ncol(km$WCSS_per_cluster) == clust && ncol(km$obs_per_cluster) == clust && class(km) == "k-means clustering")
    
    count = count + 1
  }
  
  testthat::expect_true( all(res) )
})



testthat::test_that("KMeans_rcpp returns the correct output if CENTROIDS is user-defined ", {
  
  clust = 2
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  km = KMeans_rcpp(X, clusters = clust, num_init = 5, max_iters = 100, CENTROIDS = cntr)
  
  testthat::expect_true( all(names(km) %in% c("clusters", "centroids", "total_SSE", "best_initialization", "WCSS_per_cluster", "obs_per_cluster", "between.SS_DIV_total.SS"))
                         
                         && is.vector(km$clusters) && length(km$clusters) == nrow(X) && is.numeric(km$between.SS_DIV_total.SS) && length(km$between.SS_DIV_total.SS) == 1 &&
                           
                           nrow(km$centroids) == clust && ncol(km$centroids) == ncol(X) && length(km$total_SSE) == 1 && is.numeric(km$total_SSE) && length(km$best_initialization) == 1 && 
                           
                           is.numeric(km$best_initialization) && ncol(km$WCSS_per_cluster) == clust && ncol(km$obs_per_cluster) == clust && class(km) == "k-means clustering"  )
})



# testthat::test_that("in case that the data is of type 'tibble' the result is a list and the class is 'k-means clustering' ", {
# 
#   clust = 2
# 
#   km = KMeans_rcpp(tbl, clusters = clust, num_init = 5, max_iters = 100, initializer = 'optimal_init', fuzzy = TRUE)
# 
#   testthat::expect_true( names(km) %in% c("clusters", "centroids", "total_SSE", "best_initialization", "WCSS_per_cluster", "obs_per_cluster", "between.SS_DIV_total.SS")  && is.matrix(km$fuzzy_clusters) &&
# 
#                            is.vector(km$clusters) && length(km$clusters) == nrow(tbl) && is.numeric(km$between.SS_DIV_total.SS) && length(km$between.SS_DIV_total.SS) == 1 && ncol(km$fuzzy_clusters) == clust &&
# 
#                            nrow(km$centroids) == clust && ncol(km$centroids) == ncol(tbl) && length(km$total_SSE) == 1 && is.numeric(km$total_SSE) && length(km$best_initialization) == 1 &&
# 
#                            is.numeric(km$best_initialization) && ncol(km$WCSS_per_cluster) == clust && ncol(km$obs_per_cluster) == clust && class(km) == "k-means clustering"  )
# })




################################
# error handling predict_KMeans
################################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  tmp_x = list(X)
  
  testthat::expect_error( predict_KMeans(tmp_x, CENTROIDS = cntr) )
})


testthat::test_that("in case that the CENTROIDS is not a matrix, it returns an error", {
  
  cntr = as.data.frame(matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X)))
  
  testthat::expect_error( predict_KMeans(X, CENTROIDS = cntr) )
})


testthat::test_that("in case that the columns of the CENTROIDS is not equal to the columns of the data, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X) - 1)), nrow = 2, ncol = ncol(X) - 1)
  
  testthat::expect_error( predict_KMeans(X, CENTROIDS = cntr) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = Inf
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  testthat::expect_error( predict_KMeans(tmp_dat, CENTROIDS = cntr) )
})


##########################
# predict_KMeans function
##########################


testthat::test_that("predict_KMeans returns the correct output if the input is a data frame AND if the CENTROIDS is a matrix and has the correct dimensions ", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  km = predict_KMeans(dat, CENTROIDS = cntr)
  
  testthat::expect_true( length(km) == nrow(X) && class(km) == "k-means clustering"  )
})


testthat::test_that("predict_KMeans returns the correct output if the input is a matrix AND if the CENTROIDS is a matrix and has the correct dimensions ", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  km = predict_KMeans(X, CENTROIDS = cntr)
  
  testthat::expect_true( length(km) == nrow(X) && class(km) == "k-means clustering"  )
})


testthat::test_that("the predict_KMeans works using the CENTROIDS of the KMeans_rcpp function", {
  
  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')
  
  km_preds = predict_KMeans(X, CENTROIDS = km$centroids)
  
  testthat::expect_true( length(km_preds) == nrow(X) && class(km) == "k-means clustering"  )
})


testthat::test_that("the predict_KMeans works using the CENTROIDS of the KMeans_arma function", {
  
  km = KMeans_arma(X, clusters = 2, n_iter = 10, "random_subset", verbose = F)
  
  km_preds = predict_KMeans(X, CENTROIDS = km)
  
  testthat::expect_true( length(km_preds) == nrow(X) && class(km_preds) == "k-means clustering"  )
})



#########################################
# error handling Optimal_Clusters_KMeans
#########################################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( Optimal_Clusters_KMeans(tmp_x, max_clusters = 10, criterion = 'distortion_fK', plot_clusters = FALSE) )
})


testthat::test_that("in case that the max_clusters parameter is not numeric, it returns an error", {
  
  tmp_m = data.frame(1)
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = tmp_m, criterion = 'distortion_fK', plot_clusters = FALSE) )
})


testthat::test_that("in case that the max_clusters parameter is less than 1, it returns an error", {
  
  tmp_m = 0
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = tmp_m, criterion = 'distortion_fK', plot_clusters = FALSE) )
})


testthat::test_that("if the criterion is not one of c('variance_explained', 'WCSSE', 'dissimilarity', 'silhouette', 'distortion_fK', 'AIC', 'BIC', 'Adjusted_Rsquared'), it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, criterion = 'invalid', plot_clusters = FALSE) )
})


testthat::test_that("in case that the num_init parameter is less than 1, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, criterion = 'distortion_fK', num_init = 0, plot_clusters = FALSE) )
})


testthat::test_that("in case that the max_iters parameter is less than 1, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, criterion = 'distortion_fK', max_iters = 0, plot_clusters = FALSE) )
})


testthat::test_that("if the initializer is not one of c('kmeans++', 'random', 'optimal_init', 'quantile_init'), it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, initializer = 'invalid', plot_clusters = FALSE) )
})


testthat::test_that("in case that the threads parameter is less than 1, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, criterion = 'distortion_fK', threads = 0, plot_clusters = FALSE) )
})


testthat::test_that("in case that the tol parameter is less than or equal to 0.0, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, criterion = 'distortion_fK', tol = 0.0, plot_clusters = FALSE) )
})


testthat::test_that("in case that the plot_clusters parameter is not logical, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, criterion = 'distortion_fK', plot_clusters = 'FALSE') )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, criterion = 'distortion_fK', verbose = 'FALSE', plot_clusters = FALSE) )
})


testthat::test_that("in case that the tol_optimal_init parameter is less than or equal to 0.0, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_KMeans(X, max_clusters = 5, criterion = 'distortion_fK', tol_optimal_init = 0.0, plot_clusters = FALSE) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = -Inf
  
  testthat::expect_error( Optimal_Clusters_KMeans(tmp_dat, max_clusters = 5, criterion = 'distortion_fK', plot_clusters = FALSE) )
})


testthat::test_that("in case that the 'mini_batch_params' is not NULL and the named list is not valid it returns an error", {
  
  params_mbkm = list(invalid = 10, init_fraction = 0.3, early_stop_iter = 10)

  testthat::expect_error( Optimal_Clusters_KMeans(dat, max_clusters = 10, criterion = "distortion_fK",
                                                  
                                                  plot_clusters = FALSE, mini_batch_params = params_mbkm) )
})


testthat::test_that("in case that the 'mini_batch_params' is not NULL and the criterion is 'variance_explained' it returns an error", {
  
  params_mbkm = list(batch_size = 10, init_fraction = 0.3, early_stop_iter = 10)
  
  testthat::expect_error( Optimal_Clusters_KMeans(dat, max_clusters = 10, criterion = "variance_explained",
                                                  
                                                  plot_clusters = FALSE, mini_batch_params = params_mbkm) )
})



###################################
# Optimal_Clusters_KMeans function   [ in case that the 'max_clusters' parameter is of length 1 ]
###################################


testthat::test_that("Optimal_Clusters_KMeans returns the correct output if the input is a data frame ", {
  
  nr_clust = 10
  
  res =  Optimal_Clusters_KMeans(dat, max_clusters = nr_clust, criterion = 'distortion_fK', plot_clusters = FALSE, tol_optimal_init = 0.2)
  
  testthat::expect_true( length(res) == nr_clust && class(res) == "k-means clustering"  )
})



testthat::test_that("Optimal_Clusters_KMeans returns the correct output for different criteria", {
  
  vec = c('variance_explained', 'WCSSE', 'dissimilarity', 'silhouette', 'AIC', 'BIC', 'distortion_fK', 'Adjusted_Rsquared')
  
  out = rep(NA, length(vec))
  
  nr_clust = 5
  
  count = 1
  
  for (i in vec) {
    
    res =  Optimal_Clusters_KMeans(dat, max_clusters = nr_clust, criterion = i, plot_clusters = T, tol_optimal_init = 0.2)
    
    out[count] = (length(res) == nr_clust && class(res) == "k-means clustering")
    
    count = count + 1
  }
  
  testthat::expect_true( sum(out) == length(vec) )
})


testthat::test_that("Optimal_Clusters_KMeans returns the correct output if the 'mini_batch_params' is not NULL", {
  
  nr_clust = 10
  
  params_mbkm = list(batch_size = 10, init_fraction = 0.3, early_stop_iter = 10)
  
  res = Optimal_Clusters_KMeans(dat, max_clusters = nr_clust, criterion = "distortion_fK",
                                
                                plot_clusters = FALSE, mini_batch_params = params_mbkm)
  
  testthat::expect_true( length(res) == nr_clust && class(res) == "k-means clustering"  )
})


###################################
# Optimal_Clusters_KMeans function         [ in case that the 'max_clusters' parameter is a contiguous or non-contiguous vector ] [ here I tested only the 'KMeans_rcpp' function but the same applies to 'MiniBatchKmeans' ]
###################################


testthat::test_that("max_clusters-vector for 'variance_explained'", {

  subs = 2:3
  
  res1 =  Optimal_Clusters_KMeans(dat, max_clusters = 1:3, criterion = 'variance_explained')
  res2 =  Optimal_Clusters_KMeans(dat, max_clusters = subs, criterion = 'variance_explained')
  
  testthat::expect_true( all(res1[subs] == res2)  )
})


testthat::test_that("max_clusters-vector for 'WCSSE'", {
  
  subs = c(2,4)
  
  res1 =  Optimal_Clusters_KMeans(dat, max_clusters = 1:4, criterion = 'WCSSE')
  res2 =  Optimal_Clusters_KMeans(dat, max_clusters = subs, criterion = 'WCSSE')
  
  testthat::expect_true( all(res1[subs] == res2)  )
})


testthat::test_that("max_clusters-vector for 'dissimilarity'", {
  
  subs = c(1,3)
  
  res1 =  Optimal_Clusters_KMeans(dat, max_clusters = 1:4, criterion = 'dissimilarity')
  res2 =  Optimal_Clusters_KMeans(dat, max_clusters = subs, criterion = 'dissimilarity')
  
  testthat::expect_true( all(res1[subs] == res2)  )
})


testthat::test_that("max_clusters-vector for 'silhouette'", {
  
  subs = c(2,3)
  
  res1 =  Optimal_Clusters_KMeans(dat, max_clusters = 1:3, criterion = 'silhouette')
  res2 =  Optimal_Clusters_KMeans(dat, max_clusters = subs, criterion = 'silhouette')
  
  testthat::expect_true( all(res1[subs] == res2)  )
})


testthat::test_that("max_clusters-vector for 'AIC'", {
  
  subs = c(2,4)
  
  res1 =  Optimal_Clusters_KMeans(dat, max_clusters = 1:4, criterion = 'AIC')
  res2 =  Optimal_Clusters_KMeans(dat, max_clusters = subs, criterion = 'AIC')
  
  testthat::expect_true( all(res1[subs] == res2)  )
})


testthat::test_that("max_clusters-vector for 'BIC'", {
  
  subs = c(1,4)
  
  res1 =  Optimal_Clusters_KMeans(dat, max_clusters = 1:4, criterion = 'BIC')
  res2 =  Optimal_Clusters_KMeans(dat, max_clusters = subs, criterion = 'BIC')
  
  testthat::expect_true( all(res1[subs] == res2)  )
})


################################
# error handling MiniBatchKmeans
################################

testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( MiniBatchKmeans(tmp_x, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10) )
})


testthat::test_that("in case that the clusters parameter is not numeric, it returns an error", {
  
  tmp_m = data.frame(1)
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = tmp_m, batch_size = 20, num_init = 5, early_stop_iter = 10) )
})


testthat::test_that("in case that the length of the clusters parameter is not 1, it returns an error", {
  
  tmp_m = c(1,2)
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = tmp_m, batch_size = 20, num_init = 5, early_stop_iter = 10) )
})


testthat::test_that("in case that the clusters parameter is less than 1, it returns an error", {
  
  tmp_m = 0
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = tmp_m, batch_size = 20, num_init = 5, early_stop_iter = 10) )
})


testthat::test_that("in case that the batch_size parameter is less than 1, it returns an error", {
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 0, num_init = 5, early_stop_iter = 10) )
})



testthat::test_that("in case that the num_init parameter is less than 1, it returns an error", {
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 0, early_stop_iter = 10) )
})



testthat::test_that("in case that the max_iters parameter is less than 1, it returns an error", {
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, max_iters = 0, early_stop_iter = 10) )
})


testthat::test_that("in case that the init_fraction parameter is less than or equal to 0.0, it returns an error", {
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, init_fraction = 0.0, early_stop_iter = 10) )
})


testthat::test_that("in case that the initializer parameter is not one of c('kmeans++', 'random', 'optimal_init', 'quantile_init'), it returns an error", {
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, initializer = 'invalid', early_stop_iter = 10) )
})


testthat::test_that("in case that the early_stop_iter parameter is less than 1, it returns an error", {
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, early_stop_iter = 0) )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, early_stop_iter = 10, verbose = 'FALSE') )
})


testthat::test_that("in case that CENTROIDS has invalid columns, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X) - 1)), nrow = 2, ncol = ncol(X) - 1)
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, early_stop_iter = 10, CENTROIDS = cntr) )
})


testthat::test_that("in case that CENTROIDS has invalid rows, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 3, ncol = ncol(X))
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, early_stop_iter = 10, CENTROIDS = cntr) )
})


testthat::test_that("in case that the tol parameter is less than or equal to 0.0, it returns an error", {
  
  testthat::expect_error( MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, early_stop_iter = 10, tol = 0.0) )
})


testthat::test_that("in case that the tol_optimal_init parameter is less than or equal to 0.0, it returns an error", {
  
  testthat::expect_error(  MiniBatchKmeans(X, clusters = 2, batch_size = 10, num_init = 1, early_stop_iter = 10, tol_optimal_init = 0.0) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = Inf
  
  testthat::expect_error( MiniBatchKmeans(tmp_dat, clusters = 2, batch_size = 10, num_init = 1, early_stop_iter = 10) )
})


##########################
# MiniBatchKmeans function
##########################


testthat::test_that("in case that the data is a matrix the result is a list and the class is 'k-means clustering' ", {
  
  clust = 2
  
  numinit = 5
  
  km = MiniBatchKmeans(X, clusters = clust, batch_size = 20, num_init = numinit, early_stop_iter = 10, tol_optimal_init = 0.2)
  
  testthat::expect_true( all(names(km) %in% c("centroids", "WCSS_per_cluster", "best_initialization", "iters_per_initialization"))  && is.matrix(km$centroids) && nrow(km$centroids) == clust &&
                           
                           ncol(km$centroids) == ncol(X) && is.matrix(km$WCSS_per_cluster) && ncol(km$WCSS_per_cluster) == clust && is.numeric(km$best_initialization) && length(km$best_initialization) == 1 && 
                           
                           is.matrix(km$iters_per_initialization) && ncol(km$iters_per_initialization) == numinit && class(km) == "k-means clustering"  )
})


testthat::test_that("in case that the data is a matrix the result is a list and the class is 'k-means clustering' ", {
  
  clust = 2
  
  numinit = 5
  
  km = MiniBatchKmeans(dat, clusters = clust, batch_size = 20, num_init = numinit, early_stop_iter = 10, tol_optimal_init = 0.2)
  
  testthat::expect_true( all(names(km) %in% c("centroids", "WCSS_per_cluster", "best_initialization", "iters_per_initialization"))  && is.matrix(km$centroids) && nrow(km$centroids) == clust &&
                           
                           ncol(km$centroids) == ncol(X) && is.matrix(km$WCSS_per_cluster) && ncol(km$WCSS_per_cluster) == clust && is.numeric(km$best_initialization) && length(km$best_initialization) == 1 && 
                           
                           is.matrix(km$iters_per_initialization) && ncol(km$iters_per_initialization) == numinit && class(km) == "k-means clustering"  )
})



testthat::test_that("for different parameter settings it returns the correct output", {
  
  clust = 2
  
  numinit = 5
  
  inits = c('kmeans++', 'random', 'optimal_init', 'quantile_init')
  
  res = rep(NA, length(inits))
  
  for (i in 1:length(inits)) {
    
    km = MiniBatchKmeans(dat, clusters = clust, batch_size = 20, num_init = numinit, initializer = inits[i], early_stop_iter = 10, tol_optimal_init = 0.2)
    
    res[i] = ( all(names(km) %in% c("centroids", "WCSS_per_cluster", "best_initialization", "iters_per_initialization"))  && is.matrix(km$centroids) && nrow(km$centroids) == clust &&
                
                ncol(km$centroids) == ncol(X) && is.matrix(km$WCSS_per_cluster) && ncol(km$WCSS_per_cluster) == clust && is.numeric(km$best_initialization) && length(km$best_initialization) == 1 && 
                
                is.matrix(km$iters_per_initialization) && ncol(km$iters_per_initialization) == numinit && class(km) == "k-means clustering")
  }
  
  testthat::expect_true( sum(res) == length(inits) )
})



testthat::test_that("it returns the correct output if the CENTROIDS parameter is not NULL ", {
  
  clust = 2
  
  cntr = matrix(runif(clust * (ncol(X))), nrow = clust, ncol = ncol(dat))
  
  km = MiniBatchKmeans(dat, clusters = clust, batch_size = 20, early_stop_iter = 10, CENTROIDS = cntr, tol_optimal_init = 0.2)
  
  testthat::expect_true( all(names(km) %in% c("centroids", "WCSS_per_cluster", "best_initialization", "iters_per_initialization"))  && is.matrix(km$centroids) && nrow(km$centroids) == clust &&
                           
                           ncol(km$centroids) == ncol(X) && is.matrix(km$WCSS_per_cluster) && ncol(km$WCSS_per_cluster) == clust && is.numeric(km$best_initialization) && length(km$best_initialization) == 1 && 
                           
                           is.matrix(km$iters_per_initialization) && class(km) == "k-means clustering"  )
})



testthat::test_that("in case that the init_fraction is greater than 0.0 and the intializer equals to 'kmeans++' it returns the correct output ", {
  
  clust = 2
  
  numinit = 5
  
  km = MiniBatchKmeans(X, clusters = clust, batch_size = 20, num_init = numinit, early_stop_iter = 10, init_fraction = 0.4, initializer = 'kmeans++')
  
  testthat::expect_true( all(names(km) %in% c("centroids", "WCSS_per_cluster", "best_initialization", "iters_per_initialization"))  && is.matrix(km$centroids) && nrow(km$centroids) == clust &&
                           
                           ncol(km$centroids) == ncol(X) && is.matrix(km$WCSS_per_cluster) && ncol(km$WCSS_per_cluster) == clust && is.numeric(km$best_initialization) && length(km$best_initialization) == 1 && 
                           
                           is.matrix(km$iters_per_initialization) && ncol(km$iters_per_initialization) == numinit && class(km) == "k-means clustering"  )
})


testthat::test_that("in case that the init_fraction is greater than 0.0 and the intializer equals to 'quantile_init' it returns the correct output ", {
  
  clust = 2
  
  numinit = 5
  
  km = MiniBatchKmeans(X, clusters = clust, batch_size = 20, num_init = numinit, early_stop_iter = 10, init_fraction = 0.4, initializer = "quantile_init", tol_optimal_init = 0.2)
  
  testthat::expect_true( all(names(km) %in% c("centroids", "WCSS_per_cluster", "best_initialization", "iters_per_initialization"))  && is.matrix(km$centroids) && nrow(km$centroids) == clust &&
                           
                           ncol(km$centroids) == ncol(X) && is.matrix(km$WCSS_per_cluster) && ncol(km$WCSS_per_cluster) == clust && is.numeric(km$best_initialization) && length(km$best_initialization) == 1 && 
                           
                           is.matrix(km$iters_per_initialization) && ncol(km$iters_per_initialization) == numinit && class(km) == "k-means clustering"  )
})

#####################################
# error handling predict_MBatchKMeans
#####################################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {
  
  tmp_x = list(X)
  
  MbatchKm = MiniBatchKmeans(X, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10)
  
  testthat::expect_error( predict_MBatchKMeans(tmp_x, MbatchKm$centroids, fuzzy = FALSE) )
})


testthat::test_that("in case that the CENTROIDS is not a matrix, it returns an error", {
  
  cntr = as.data.frame(matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X)))
  
  testthat::expect_error( predict_MBatchKMeans(X, CENTROIDS = cntr, fuzzy = FALSE) )
})


testthat::test_that("in case that the columns of the CENTROIDS is not equal to the columns of the data, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X) - 1)), nrow = 2, ncol = ncol(X) - 1)
  
  testthat::expect_error( predict_MBatchKMeans(X, CENTROIDS = cntr) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = Inf
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  testthat::expect_error( predict_MBatchKMeans(tmp_dat, CENTROIDS = cntr) )
})


testthat::test_that("in case that the fuzzy parameter is not logical, it returns an error", {
  
  cntr = matrix(runif(2 * (ncol(X))), nrow = 2, ncol = ncol(X))
  
  testthat::expect_error( predict_MBatchKMeans(X, CENTROIDS = cntr, fuzzy = 'FALSE') )
})


################################
# predict_MBatchKMeans function
################################



testthat::test_that("in case that the data is a matrix (fuzzy = TRUE) the result is a list and the class is 'k-means clustering' ", {
  
  MbatchKm = MiniBatchKmeans(X, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10)
  
  km = predict_MBatchKMeans(X, MbatchKm$centroids, fuzzy = TRUE)
  
  testthat::expect_true( all(names(km) %in% c("clusters", "fuzzy_clusters"))  && is.matrix(km$fuzzy_clusters) && nrow(km$fuzzy_clusters) == nrow(X) && ncol(km$fuzzy_clusters) == 2 &&
                           
                           is.vector(km$clusters) && length(km$clusters) == nrow(X) && class(km) == "k-means clustering"  )
})



testthat::test_that("in case that the data is a data frame (fuzzy = TRUE) the result is a list and the class is 'k-means clustering' ", {
  
  MbatchKm = MiniBatchKmeans(dat, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10)
  
  km = predict_MBatchKMeans(dat, MbatchKm$centroids, fuzzy = TRUE)
  
  testthat::expect_true( all(names(km) %in% c("clusters", "fuzzy_clusters"))  && is.matrix(km$fuzzy_clusters) && nrow(km$fuzzy_clusters) == nrow(X) && ncol(km$fuzzy_clusters) == 2 &&
                           
                           is.vector(km$clusters) && length(km$clusters) == nrow(X) && class(km) == "k-means clustering"  )
})


testthat::test_that("in case that the data is a matrix (fuzzy = FALSE) the result is a vector and the class is 'k-means clustering' ", {
  
  MbatchKm = MiniBatchKmeans(X, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10)
  
  km = predict_MBatchKMeans(X, MbatchKm$centroids, fuzzy = FALSE)
  
  testthat::expect_true( is.numeric(km) && length(km) == nrow(X) && class(km) == "k-means clustering"  )
})

