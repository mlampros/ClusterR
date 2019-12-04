#============================================================

# data

data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

X = center_scale(dat)


#=============================================================


context('plot_2d - silhouette plot - external_validation - center_scale - distance_matrix')


########################
# error handling plot_2d
########################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {

  tmp_x = list(X)

  clust = sample(1:2, nrow(X), replace = T)

  centr = matrix(runif(ncol(X) * 2), nrow = 2, ncol = ncol(X))

  testthat::expect_error( plot_2d(tmp_x, clust, centr) )
})



testthat::test_that("in case that the data is not 2-dimensional, it returns an error", {

  clust = sample(1:2, nrow(X), replace = T)

  centr = data.frame(matrix(runif(ncol(X) * 2), nrow = 2, ncol = ncol(X)))

  testthat::expect_error( plot_2d(dat, clust, centr) )
})


testthat::test_that("in case that the unique levels of the clusters is greater than 26, it returns an error", {

  clust = sample(1:30, nrow(X), replace = T)

  centr = data.frame(matrix(runif(ncol(X) * length(unique(clust))), nrow = length(unique(clust)), ncol = ncol(X)))

  testthat::expect_error( plot_2d(dat, clust, centr) )
})


testthat::test_that("in case that the clusters parameter is not a numeric vector, it returns an error", {

  tmp_m = data.frame(1)

  centr = matrix(runif(ncol(X) * 2), nrow = 2, ncol = ncol(X))

  testthat::expect_error( plot_2d(X, tmp_m, centr) )
})



testthat::test_that("in case that the centroids/medoids is not a matrix or data frame, it returns an error", {

  clust = sample(1:2, nrow(X), replace = T)

  centr = list(matrix(runif(ncol(X) * 2), nrow = 2, ncol = ncol(X)))

  testthat::expect_error( plot_2d(X, clust, centr) )
})


testthat::test_that("in case that the centroids/medoids rows do not equal the length of the unique levels of the clusters vector, it returns an error", {

  clust = sample(1:2, nrow(X), replace = T)

  centr = matrix(runif(ncol(X) * 3), nrow = 3, ncol = ncol(X))

  testthat::expect_error( plot_2d(X, clust, centr) )
})


testthat::test_that("in case that the centroids/medoids columns do not equal the number of columns of the data, it returns an error", {

  clust = sample(1:2, nrow(X), replace = T)

  centr = matrix(runif((ncol(X) - 1) * 2), nrow = 2, ncol = ncol(X) - 1)

  testthat::expect_error( plot_2d(X, clust, centr) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {

  tmp_dat = X

  tmp_dat[1,1] = Inf

  clust = sample(1:2, nrow(X), replace = T)

  centr = matrix(runif(ncol(X) * 2), nrow = 2, ncol = ncol(X))

  testthat::expect_error( plot_2d(tmp_dat, clust, centr)  )
})



########################
# error handling plot_2d
########################


testthat::test_that("the function returns a plot of class 'gg', 'ggplot' ", {

  pca_dat = stats::princomp(X)$scores[, 1:2]

  km = KMeans_rcpp(pca_dat, clusters = 2, num_init = 5, max_iters = 100)

  testthat::expect_true( inherits(plot_2d(pca_dat, km$clusters, km$centroids), c("gg", "ggplot")) )
})



#####################################
# error handling external_validation
#####################################


testthat::test_that("in case that the true labels is not a numeric vector, it returns an error", {

  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')

  tmp_cl = list(dietary_survey_IBS$class)

  testthat::expect_error( external_validation(tmp_cl, km$clusters, method = "adjusted_rand_index") )
})


testthat::test_that("in case that the clusters is not a numeric vector, it returns an error", {

  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')

  tmp_cl = list(km$clusters)

  testthat::expect_error( external_validation(dietary_survey_IBS$class, tmp_cl, method = "adjusted_rand_index") )
})


testthat::test_that("in case that the clusters is not a numeric vector, it returns an error", {

  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')

  tmp_cl = sample(1:2, nrow(X) + 1, replace = T)

  testthat::expect_error( external_validation(dietary_survey_IBS$class, tmp_cl, method = "adjusted_rand_index") )
})


testthat::test_that("in case that the method is not valid, it returns an error", {

  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')

  testthat::expect_error( external_validation(dietary_survey_IBS$class, km$clusters, method = "invalid") )
})



###############################
# external_validation function
###############################


testthat::test_that("in case that the true labels is an integer vector they will be converted to numeric", {

  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')

  tmp_cl = as.integer(dietary_survey_IBS$class)

  res = external_validation(tmp_cl, km$clusters, method = "adjusted_rand_index")

  testthat::expect_true( is.numeric(res) && length(res) == 1 )
})


testthat::test_that("in case that the clusters is an integer vector they will be converted to numeric", {

  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')

  tmp_cl = as.integer(km$clusters)

  res = external_validation(dietary_survey_IBS$class, tmp_cl, method = "adjusted_rand_index")

  testthat::expect_true( is.numeric(res) && length(res) == 1 )
})


testthat::test_that("the function for the different methods returns the correct output", {

  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')

  mth = c('rand_index', 'adjusted_rand_index', 'jaccard_index', 'fowlkes_mallows_index', 'mirkin_metric', 'purity', 'entropy', 'nmi', 'var_info', 'nvi')

  res = rep(NA, length(mth))

  for (i in 1:length(mth)) {

    tmp_res = external_validation(dietary_survey_IBS$class, km$clusters, method = mth[i])

    res[i] = (is.numeric(tmp_res) && length(tmp_res) == 1)
  }

  testthat::expect_true( sum(res) == length(mth) )
})



testthat::test_that("if summary_stats is TRUE the function returns output", {

  km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'optimal_init')

  res = external_validation(dietary_survey_IBS$class, km$clusters, method = "adjusted_rand_index", summary_stats = T)

  testthat::expect_output( str(res), 'num' )
})



#############################
# error handling center_scale
#############################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {

  tmp_x = list(X)

  testthat::expect_error( center_scale(tmp_x, mean_center = TRUE, sd_scale = TRUE) )
})


testthat::test_that("in case that the mean_center parameter is not logical, it returns an error", {

  testthat::expect_error( center_scale(X, mean_center = 'TRUE', sd_scale = TRUE) )
})


testthat::test_that("in case that the sd_scale parameter is not logical, it returns an error", {

  testthat::expect_error( center_scale(X, mean_center = TRUE, sd_scale = 'TRUE') )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {

  tmp_dat = X

  tmp_dat[1,1] = Inf

  testthat::expect_error( center_scale(tmp_dat, mean_center = TRUE, sd_scale = TRUE)  )
})



#######################
# center_scale function
#######################


testthat::test_that("in case that the data is a matrix the function returns a matrix", {

  res = center_scale(X, mean_center = TRUE, sd_scale = TRUE)

  testthat::expect_true( is.matrix(res) && nrow(X) == nrow(res) && ncol(X) == ncol(res) )
})


testthat::test_that("in case that the data is a data frame the function returns a matrix", {

  res = center_scale(dat, mean_center = TRUE, sd_scale = TRUE)

  testthat::expect_true( is.matrix(res) && nrow(X) == nrow(res) && ncol(X) == ncol(res) )
})


testthat::test_that("in case that the data is a matrix and mean_center = TRUE and sd_scale = FALSE, the function returns a matrix", {

  res = center_scale(X, mean_center = TRUE, sd_scale = FALSE)

  testthat::expect_true( is.matrix(res) && nrow(X) == nrow(res) && ncol(X) == ncol(res) )
})


testthat::test_that("in case that the data is a matrix and mean_center = FALSE and sd_scale = TRUE, the function returns a matrix", {

  res = center_scale(X, mean_center = FALSE, sd_scale = TRUE)

  testthat::expect_true( is.matrix(res) && nrow(X) == nrow(res) && ncol(X) == ncol(res) )
})


testthat::test_that("in case that the data is a matrix and mean_center = FALSE and sd_scale = FALSE, the function returns a matrix", {

  res = center_scale(X, mean_center = FALSE, sd_scale = FALSE)

  testthat::expect_true( is.matrix(res) && nrow(X) == nrow(res) && ncol(X) == ncol(res) )
})


#################################
# error handling distance_matrix
#################################


testthat::test_that("in case that the data is not a matrix or data frame, it returns an error", {

  tmp_x = list(X)

  testthat::expect_error( distance_matrix(tmp_x, method = 'euclidean', upper = TRUE, diagonal = TRUE) )
})


testthat::test_that("in case that the method parameter is invalid, it returns an error", {

  testthat::expect_error( distance_matrix(X, method = 'invalid', upper = TRUE, diagonal = TRUE) )
})


testthat::test_that("in case that the upper parameter is not logical, it returns an error", {

  testthat::expect_error( distance_matrix(X, method = 'euclidean', upper = 'TRUE', diagonal = TRUE) )
})



testthat::test_that("in case that the diagonal parameter is not logical, it returns an error", {

  testthat::expect_error( distance_matrix(X, method = 'euclidean', upper = TRUE, diagonal = 'TRUE') )
})


testthat::test_that("in case that the method parameter is minkowski and the minkowski_p parameter is 0.0, it returns an error", {

  testthat::expect_error( distance_matrix(X, method = 'minkowski', minkowski_p = 0.0) )
})


testthat::test_that("in case that the threads parameter is less than 1, it returns an error", {

  testthat::expect_error( distance_matrix(X, method = 'euclidean', threads = 0) )
})


# from version 1.0.3 the "distance_matrix" function can accept data with missing values
#
# testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
#
#   tmp_dat = X
#
#   tmp_dat[1,1] = Inf
#
#   testthat::expect_error( distance_matrix(tmp_dat, method = 'euclidean', upper = TRUE, diagonal = TRUE)  )
# })


##########################
# distance_matrix function
##########################



testthat::test_that("in case that the data is a matrix the function returns a matrix", {

  res = distance_matrix(X, method = 'euclidean', upper = TRUE, diagonal = TRUE)

  testthat::expect_true( is.matrix(res) )
})


testthat::test_that("in case that the data is a data frame the function returns a matrix", {

  res = distance_matrix(dat, method = 'euclidean', upper = TRUE, diagonal = TRUE)

  testthat::expect_true( is.matrix(res) )
})



###############################################
# error handling Silhouette_Dissimilarity_Plot
###############################################


testthat::test_that("in case that the evaluation_object is invalid, it returns an error", {

  km = KMeans_arma(X, clusters = 2, n_iter = 10, "random_subset", verbose = F)

  testthat::expect_error( Silhouette_Dissimilarity_Plot(km, silhouette = TRUE) )
})


testthat::test_that("in case that the silhouette parameter is not logical, it returns an error", {

  cm = Cluster_Medoids(X, clusters = 2, distance_metric = 'euclidean')

  testthat::expect_error( Silhouette_Dissimilarity_Plot(cm, silhouette = 'TRUE') )
})


#########################################
# Silhouette_Dissimilarity_Plot function
#########################################


testthat::test_that("in case of 'Cluster_Medoids' function AND silhouette = TRUE, it returns TRUE", {

  cm = Cluster_Medoids(X, clusters = 4, distance_metric = 'euclidean')

  plt_sd = Silhouette_Dissimilarity_Plot(cm, silhouette = TRUE)

  testthat::expect_true( plt_sd == T )
})


testthat::test_that("in case of 'Cluster_Medoids' function AND silhouette = FALSE, it returns TRUE", {

  cm = Cluster_Medoids(X, clusters = 8, distance_metric = 'euclidean')

  plt_sd = Silhouette_Dissimilarity_Plot(cm, silhouette = FALSE)

  testthat::expect_true( plt_sd == T )
})


testthat::test_that("in case of 'Clara_Medoids' function AND silhouette = TRUE, it returns TRUE", {

  cm = Clara_Medoids(X, clusters = 8, samples = 5, sample_size = 0.2, swap_phase = TRUE, fuzzy = T)

  plt_sd = Silhouette_Dissimilarity_Plot(cm, silhouette = TRUE)

  testthat::expect_true( plt_sd == T )
})


testthat::test_that("in case of 'Clara_Medoids' function AND silhouette = FALSE, it returns TRUE", {

  cm = Clara_Medoids(X, clusters = 4, samples = 5, sample_size = 0.2, swap_phase = TRUE, fuzzy = T)

  plt_sd = Silhouette_Dissimilarity_Plot(cm, silhouette = FALSE)

  testthat::expect_true( plt_sd == T )
})
