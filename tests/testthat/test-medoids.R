#============================================================

# data

data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

X = center_scale(dat)


#=============================================================


context('clustering using medoids')


################################
# error handling Cluster_Medoids
################################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {

  tmp_x = list(X)

  testthat::expect_error( Cluster_Medoids(tmp_x, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the clusters parameter is not numeric, it returns an error", {

  tmp_m = data.frame(1)

  testthat::expect_error( Cluster_Medoids(X, clusters = tmp_m, distance_metric = 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the length of the clusters parameter is not 1, it returns an error", {

  tmp_m = c(1,2)

  testthat::expect_error( Cluster_Medoids(X, clusters = tmp_m, distance_metric = 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the clusters parameter is less than 1, it returns an error", {

  tmp_m = 0

  testthat::expect_error( Cluster_Medoids(X, clusters = tmp_m, distance_metric = 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the distance_metric parameter is invalid, it returns an error", {

  testthat::expect_error( Cluster_Medoids(X, clusters = 2, distance_metric = 'invalid', swap_phase = TRUE) )
})


testthat::test_that("in case that the distance_metric parameter is minkowski and the minkowski_p parameter is 0.0, it returns an error", {

  testthat::expect_error( Cluster_Medoids(X, clusters = 2, distance_metric = 'minkowski', minkowski_p = 0.0, swap_phase = TRUE) )
})


testthat::test_that("in case that the threads parameter is less than 1, it returns an error", {

  testthat::expect_error( Cluster_Medoids(X, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, threads = 0) )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {

  testthat::expect_error( Cluster_Medoids(X, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, verbose = 'invalid') )
})


testthat::test_that("in case that the swap_phase parameter is not logical, it returns an error", {

  testthat::expect_error( Cluster_Medoids(X, clusters = 2, distance_metric = 'euclidean', swap_phase = 'invalid') )
})


testthat::test_that("in case that the fuzzy parameter is not logical, it returns an error", {

  testthat::expect_error( Cluster_Medoids(X, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, fuzzy = 'invalid') )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {

  tmp_dat = X

  tmp_dat[1,1] = NaN

  testthat::expect_error( Cluster_Medoids(tmp_dat, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE)  )
})


##########################
# Cluster_Medoids function
##########################


testthat::test_that("in case that the data is a matrix, it returns the correct output", {

  cm = Cluster_Medoids(X, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, fuzzy = T)

  testthat::expect_true( sum(names(cm) %in% c("medoids", "medoid_indices", "best_dissimilarity", "dissimilarity_matrix",
                                          "clusters", "silhouette_matrix", "fuzzy_probs", "clustering_stats")) == 8 && class(cm) == "cluster medoids silhouette" && is.matrix(cm$medoids) && is.vector(cm$medoid_indices) &&
                           is.numeric(cm$best_dissimilarity) && is.vector(cm$clusters) && is.data.frame(cm$silhouette_matrix) && is.matrix(cm$fuzzy_probs) && is.data.frame(cm$clustering_stats) )
})



testthat::test_that("in case that the data is a dissimilarity matrix, it returns the correct output", {

  dism_mat = distance_matrix(X, method = 'euclidean', upper = TRUE, diagonal = TRUE)

  cm = Cluster_Medoids(dism_mat, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, fuzzy = T)

  testthat::expect_true( sum(names(cm) %in% c("medoids", "medoid_indices", "best_dissimilarity", "dissimilarity_matrix",
                                              "clusters", "silhouette_matrix", "fuzzy_probs", "clustering_stats")) == 8 && class(cm) == "cluster medoids silhouette" && is.vector(cm$medoids) && is.vector(cm$medoid_indices) &&
                           is.numeric(cm$best_dissimilarity) && is.vector(cm$clusters) && is.data.frame(cm$silhouette_matrix) && is.matrix(cm$fuzzy_probs) && is.data.frame(cm$clustering_stats) )
})



testthat::test_that("in case that the data is a matrix, it returns the correct output", {

  cm = Cluster_Medoids(dat, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, fuzzy = T)

  testthat::expect_true( sum(names(cm) %in% c("medoids", "medoid_indices", "best_dissimilarity", "dissimilarity_matrix",
                                              "clusters", "silhouette_matrix", "fuzzy_probs", "clustering_stats")) == 8 && class(cm) == "cluster medoids silhouette" && is.matrix(cm$medoids) && is.vector(cm$medoid_indices) &&
                           is.numeric(cm$best_dissimilarity) && is.vector(cm$clusters) && is.data.frame(cm$silhouette_matrix) && is.matrix(cm$fuzzy_probs) && is.data.frame(cm$clustering_stats) )
})


testthat::test_that("in case that clusters = 1, it returns the correct output", {

  cm = Cluster_Medoids(dat, clusters = 1, distance_metric = 'euclidean', swap_phase = TRUE, fuzzy = T)

  testthat::expect_true( sum(names(cm) %in% c("medoids", "medoid_indices", "best_dissimilarity", "dissimilarity_matrix",
                                              "clusters", "silhouette_matrix", "fuzzy_probs", "clustering_stats")) == 8 && class(cm) == "cluster medoids silhouette" && is.vector(cm$medoids) && is.vector(cm$medoid_indices) &&
                           is.numeric(cm$best_dissimilarity) && is.vector(cm$clusters) && is.null(cm$silhouette_matrix) && is.matrix(cm$fuzzy_probs) && is.null(cm$clustering_stats) )
})


##############################
# error handling Clara_Medoids
##############################



testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {

  tmp_x = list(X)

  testthat::expect_error( Clara_Medoids(tmp_x, clusters = 2, samples = 5, sample_size = 0.2, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the clusters parameter is not numeric, it returns an error", {

  tmp_m = data.frame(1)

  testthat::expect_error( Clara_Medoids(X, clusters = tmp_m, samples = 5, sample_size = 0.2, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the length of the clusters parameter is not 1, it returns an error", {

  tmp_m = c(1,2)

  testthat::expect_error( Clara_Medoids(X, clusters = tmp_m, samples = 5, sample_size = 0.2, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the clusters parameter is less than 1, it returns an error", {

  tmp_m = 0

  testthat::expect_error( Clara_Medoids(X, clusters = tmp_m, samples = 5, sample_size = 0.2, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the samples parameter is not numeric, it returns an error", {

  tmp_s = '0'

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = tmp_s, sample_size = 0.2, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the length of the samples parameter is not 1, it returns an error", {

  tmp_s = c(0,1)

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = tmp_s, sample_size = 0.2, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the samples parameter is not less than 1, it returns an error", {

  tmp_s = 0

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = tmp_s, sample_size = 0.2, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the sample_size parameter is not numeric, it returns an error", {

  tmp_s = '0'

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = 5, sample_size = tmp_s, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the sample_size parameter is less than 0.0, it returns an error", {

  tmp_s = -1.0

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = 5, sample_size = tmp_s, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the samples parameter is not less than 1, it returns an error", {

  tmp_s = 2.0

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = 5, sample_size = tmp_s, 'euclidean', swap_phase = TRUE) )
})


testthat::test_that("in case that the distance_metric parameter is invalid, it returns an error", {

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = 5, sample_size = 0.2, 'invalid', swap_phase = TRUE) )
})


testthat::test_that("in case that the distance_metric parameter is minkowski and the minkowski_p parameter is 0.0, it returns an error", {

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = 5, sample_size = 0.2, 'minkowski', minkowski_p = 0.0, swap_phase = TRUE) )
})


testthat::test_that("in case that the threads parameter is less than 1, it returns an error", {

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = 5, sample_size = 0.2, threads = 0, swap_phase = TRUE) )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = 5, sample_size = 0.2, verbose = 'invalid', swap_phase = TRUE) )
})


testthat::test_that("in case that the swap_phase parameter is not logical, it returns an error", {

  testthat::expect_error(  Clara_Medoids(X, clusters = 2, samples = 5, sample_size = 0.2, swap_phase = 'TRUE') )
})


testthat::test_that("in case that the fuzzy parameter is not logical, it returns an error", {

  testthat::expect_error( Clara_Medoids(X, clusters = 2, samples = 5, sample_size = 0.2, fuzzy = 'TRUE') )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {

  tmp_dat = X

  tmp_dat[1,1] = NaN

  testthat::expect_error( Clara_Medoids(tmp_dat, clusters = 2, samples = 5, sample_size = 0.2, swap_phase = TRUE) )
})



#########################
# Clara_Medoids function
#########################



testthat::test_that("in case that the data is a matrix, it returns the correct output", {

  cm = Clara_Medoids(X, clusters = 2, samples = 5, sample_size = 0.2, swap_phase = TRUE, fuzzy = T)

  testthat::expect_true( sum(names(cm) %in% c("medoids", "medoid_indices", "sample_indices", "best_dissimilarity",
                                              "clusters", "silhouette_matrix", "fuzzy_probs", "clustering_stats",
                                              "dissimilarity_matrix")) == 9 && class(cm) == "cluster medoids silhouette" && is.matrix(cm$medoids) && is.vector(cm$medoid_indices) && is.vector(cm$sample_indices) &&
                           is.numeric(cm$best_dissimilarity) && is.vector(cm$clusters) && is.data.frame(cm$silhouette_matrix) && is.matrix(cm$fuzzy_probs) && is.data.frame(cm$clustering_stats) )
})


testthat::test_that("in case that the data is a data frame, it returns the correct output", {

  cm = Clara_Medoids(dat, clusters = 2, samples = 5, sample_size = 0.2, swap_phase = TRUE, fuzzy = T)

  testthat::expect_true( sum(names(cm) %in% c("medoids", "medoid_indices", "sample_indices", "best_dissimilarity",
                                              "clusters", "silhouette_matrix", "fuzzy_probs", "clustering_stats",
                                              "dissimilarity_matrix")) == 9 && class(cm) == "cluster medoids silhouette" && is.matrix(cm$medoids) && is.vector(cm$medoid_indices) && is.vector(cm$sample_indices) &&
                           is.numeric(cm$best_dissimilarity) && is.vector(cm$clusters) && is.data.frame(cm$silhouette_matrix) && is.matrix(cm$fuzzy_probs) && is.data.frame(cm$clustering_stats) )
})



testthat::test_that("in case that the clusters parameter is 1, it returns the correct output", {

  cm = Clara_Medoids(dat, clusters = 1, samples = 5, sample_size = 0.2, swap_phase = TRUE, fuzzy = T)

  testthat::expect_true( sum(names(cm) %in% c("medoids", "medoid_indices", "sample_indices", "best_dissimilarity",
                                              "clusters", "silhouette_matrix", "fuzzy_probs", "clustering_stats",
                                              "dissimilarity_matrix")) == 9 && class(cm) == "cluster medoids silhouette" && is.matrix(cm$medoids) && is.vector(cm$medoid_indices) && is.vector(cm$sample_indices) &&
                           is.numeric(cm$best_dissimilarity) && is.vector(cm$clusters) && is.null(cm$silhouette_matrix) && is.matrix(cm$fuzzy_probs) && is.data.frame(cm$clustering_stats) )
})



#################################
# error handling predict_Medoids
#################################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {

  cm = Cluster_Medoids(X, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, fuzzy = T)

  tmp_x = list(X)

  testthat::expect_error( predict_Medoids(tmp_x, MEDOIDS = cm$medoids, 'euclidean', fuzzy = TRUE) )
})


testthat::test_that("in case that the MEDOIDS is NULL, it returns an error", {

  testthat::expect_error( predict_Medoids(X, MEDOIDS = NULL, 'euclidean', fuzzy = TRUE) )
})


testthat::test_that("in case that the MEDOIDS is NULL, it returns an error", {

  tmp_cm = data.frame(matrix(runif(82), ncol = 41, nrow = 2))

  testthat::expect_error( predict_Medoids(X, MEDOIDS = tmp_cm, 'euclidean', fuzzy = TRUE) )
})



testthat::test_that("in case that the distance_metric parameter is invalid, it returns an error", {

  tmp_cm = data.frame(matrix(runif(84), ncol = 42, nrow = 2))

  testthat::expect_error( predict_Medoids(X, MEDOIDS = tmp_cm, 'invalid', fuzzy = TRUE) )
})



testthat::test_that("in case that the distance_metric parameter is minkowski and the minkowski_p parameter is 0.0, it returns an error", {

  tmp_cm = data.frame(matrix(runif(84), ncol = 42, nrow = 2))

  testthat::expect_error( predict_Medoids(X, MEDOIDS = tmp_cm, 'minkowski', minkowski_p = 0.0) )
})


testthat::test_that("in case that the threads parameter is less than 1, it returns an error", {

  tmp_cm = data.frame(matrix(runif(84), ncol = 42, nrow = 2))

  testthat::expect_error( predict_Medoids(X, MEDOIDS = tmp_cm, threads = 0) )
})


testthat::test_that("in case that the fuzzy parameter is not logical, it returns an error", {

  tmp_cm = data.frame(matrix(runif(84), ncol = 42, nrow = 2))

  testthat::expect_error( predict_Medoids(X, MEDOIDS = tmp_cm, fuzzy = 'TRUE') )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {

  tmp_cm = data.frame(matrix(runif(84), ncol = 42, nrow = 2))

  tmp_dat = X

  tmp_dat[1,1] = NaN

  testthat::expect_error( predict_Medoids(tmp_dat, MEDOIDS = tmp_cm) )
})



###########################
# predict_Medoids function
###########################


testthat::test_that("in case that the data is a matrix, it returns the correct output", {

  cm = Cluster_Medoids(X, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, fuzzy = T)

  pr = predict_Medoids(X, MEDOIDS = cm$medoids, 'euclidean', fuzzy = TRUE)

  testthat::expect_true( sum(names(pr) %in% c("clusters", "fuzzy_clusters", "dissimilarity")) == 3 && class(pr) == "cluster medoids silhouette" && is.vector(pr$clusters) &&
                           is.matrix(pr$fuzzy_clusters) && is.numeric(pr$dissimilarity) )
})


testthat::test_that("in case that the data is a data frame, it returns the correct output", {

  cm = Cluster_Medoids(dat, clusters = 2, distance_metric = 'euclidean', swap_phase = TRUE, fuzzy = T)

  pr = predict_Medoids(dat, MEDOIDS = cm$medoids, 'euclidean', fuzzy = TRUE)

  testthat::expect_true( sum(names(pr) %in% c("clusters", "fuzzy_clusters", "dissimilarity")) == 3 && class(pr) == "cluster medoids silhouette" && is.vector(pr$clusters) &&
                           is.matrix(pr$fuzzy_clusters) && is.numeric(pr$dissimilarity) )
})


testthat::test_that("in case that MEDOIDS is a data frame, it returns the correct output", {

  tmp_cm = data.frame(matrix(runif(84), ncol = 42, nrow = 2))

  pr = predict_Medoids(X, MEDOIDS = tmp_cm, 'euclidean', fuzzy = TRUE)

  testthat::expect_true( sum(names(pr) %in% c("clusters", "fuzzy_clusters", "dissimilarity")) == 3 && class(pr) == "cluster medoids silhouette" && is.vector(pr$clusters) &&
                           is.matrix(pr$fuzzy_clusters) && is.numeric(pr$dissimilarity) )
})


#########################################
# error handling Optimal_Clusters_Medoids
#########################################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {

  tmp_x = list(X)

  testthat::expect_error( Optimal_Clusters_Medoids(tmp_x, max_clusters = 10, 'euclidean', 'dissimilarity', plot_clusters = FALSE) )
})


testthat::test_that("in case that the max_clusters parameter is not numeric, it returns an error", {

  tmp_m = data.frame(1)

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = tmp_m, 'euclidean', 'dissimilarity', plot_clusters = FALSE) )
})


testthat::test_that("in case that the length of the max_clusters parameter is not 1, it returns an error", {

  tmp_m = list(1,2)

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = tmp_m, 'euclidean', 'dissimilarity', plot_clusters = FALSE) )
})


testthat::test_that("in case that the max_clusters parameter is less than 1, it returns an error", {

  tmp_m = 0

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = tmp_m, 'euclidean', 'dissimilarity', plot_clusters = FALSE) )
})


testthat::test_that("in case that the distance_metric parameter is invalid, it returns an error", {

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, 'invalid', 'dissimilarity', plot_clusters = FALSE) )
})


testthat::test_that("in case that the distance_metric parameter is minkowski and the minkowski_p parameter is 0.0, it returns an error", {

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, 'minkowski', minkowski_p = 0.0, 'dissimilarity', plot_clusters = FALSE) )
})


testthat::test_that("in case that the threads parameter is less than 1, it returns an error", {

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, 'euclidean', 'dissimilarity', plot_clusters = FALSE, threads = 0) )
})


testthat::test_that("in case that the criterion parameter is not one of 'silhouette', 'dissimilarity', it returns an error", {

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, 'euclidean', 'invalid', plot_clusters = FALSE) )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, 'euclidean', 'dissimilarity', plot_clusters = FALSE, verbose = 0)  )
})


testthat::test_that("in case that the swap_phase parameter is not logical, it returns an error", {

  testthat::expect_error(  Optimal_Clusters_Medoids(X, max_clusters = 5, 'euclidean', 'dissimilarity', plot_clusters = FALSE, swap_phase = 0) )
})


testthat::test_that("in case that the plot_clusters parameter is not logical, it returns an error", {

  testthat::expect_error(  Optimal_Clusters_Medoids(X, max_clusters = 5, 'euclidean', 'dissimilarity', plot_clusters = 0) )
})


testthat::test_that("in case that the one of clara_samples, clara_sample_size is 0, it returns an error", {

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, clara_samples = 5, clara_sample_size = 0.0, 'euclidean', 'dissimilarity', plot_clusters = F) )
})


testthat::test_that("in case that the one of clara_samples, clara_sample_size is 0, it returns an error", {

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, clara_samples = 0, clara_sample_size = 0.5, 'euclidean', 'dissimilarity', plot_clusters = F) )
})


testthat::test_that("in case that the clara_samples parameter is not numeric, it returns an error", {

  tmp_s = '0'

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, clara_samples = tmp_s, clara_sample_size = 0.2, 'euclidean', 'dissimilarity', plot_clusters = F) )
})


testthat::test_that("in case that the clara_samples parameter is not less than 1, it returns an error", {

  tmp_s = 0

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, clara_samples = tmp_s, clara_sample_size = 0.2, 'euclidean', 'dissimilarity', plot_clusters = F) )
})


testthat::test_that("in case that the clara_sample_size parameter is not numeric, it returns an error", {

  tmp_s = '0'

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, clara_samples = 5, clara_sample_size = tmp_s, 'euclidean', 'dissimilarity', plot_clusters = F) )
})


testthat::test_that("in case that the clara_sample_size parameter is less than 0.0, it returns an error", {

  tmp_s = -1.0

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, clara_samples = tmp_s, clara_sample_size = 0.2, 'euclidean', 'dissimilarity', plot_clusters = F) )
})


testthat::test_that("in case that the clara_samples parameter is not less than 1, it returns an error", {

  tmp_s = 2.0

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = 5, clara_samples = 5, clara_sample_size = tmp_s, 'euclidean', 'dissimilarity', plot_clusters = F) )
})


testthat::test_that("in case that the clara_samples and clara_sample_size parameters is greater than 0.0 and the data is a dissimilarity matrix, it returns an error", {

  vec = rep(0, 10)

  mat_diag = diag(vec)

  testthat::expect_error( Optimal_Clusters_Medoids(mat_diag, max_clusters = 5, clara_samples = 5, clara_sample_size = 0.3, 'euclidean', 'dissimilarity', plot_clusters = F) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {

  tmp_cm = data.frame(matrix(runif(84), ncol = 42, nrow = 2))

  tmp_dat = X

  tmp_dat[1,1] = NaN

  testthat::expect_error( Optimal_Clusters_Medoids(tmp_dat, max_clusters = 5, 'euclidean', 'dissimilarity', plot_clusters = F)  )
})


testthat::test_that("in case that the max_clusters parameter includes a 0, it returns an error", {

  testthat::expect_error( Optimal_Clusters_Medoids(X, max_clusters = c(0,3,5), clara_samples = 5, clara_sample_size = 0.5, 'euclidean', 'dissimilarity', plot_clusters = F) )
})



####################################
# Optimal_Clusters_Medoids function
####################################


testthat::test_that("in case that the data is a matrix, it returns the correct output", {

  opt_md = Optimal_Clusters_Medoids(X, max_clusters = 10, 'euclidean', 'dissimilarity', plot_clusters = F)

  testthat::expect_true( class(opt_md) == "cluster medoids silhouette" && mean(unlist(lapply(opt_md, length))[-1]) == 6 )
})


testthat::test_that("in case that the data is a matrix, it returns the correct output [ non-contiguous vector ]", {

  opt_md = Optimal_Clusters_Medoids(X, max_clusters = 6, 'euclidean', 'dissimilarity', plot_clusters = F)
  opt_md1 = Optimal_Clusters_Medoids(X, max_clusters = c(2,4,6), 'euclidean', 'dissimilarity', plot_clusters = F)

  clust_two = opt_md[[2]]$avg_intra_clust_dissimilarity == opt_md1[[1]]$avg_intra_clust_dissimilarity
  clust_four = opt_md[[4]]$sum_intra_dissim == opt_md1[[2]]$sum_intra_dissim
  clust_six = opt_md[[6]]$avg_width_silhouette == opt_md1[[3]]$avg_width_silhouette

  testthat::expect_true( class(opt_md) == "cluster medoids silhouette" && mean(unlist(lapply(opt_md, length))[-1]) == 6 && all(clust_two, clust_four, clust_six) )
})


testthat::test_that("in case that the data is a data frame, it returns the correct output", {

  opt_md = Optimal_Clusters_Medoids(dat, max_clusters = 10, 'euclidean', 'dissimilarity', plot_clusters = FALSE)

  testthat::expect_true( class(opt_md) == "cluster medoids silhouette" && mean(unlist(lapply(opt_md, length))[-1]) == 6 )
})


testthat::test_that("in case of Cluster_Medoids for different parameter settings, it returns the correct output", {

  tmp = c("silhouette","dissimilarity")

  res = rep(NA, length(tmp))

  for (i in 1:length(tmp)) {

    opt_md = Optimal_Clusters_Medoids(dat, max_clusters = 10, 'euclidean', tmp[i], plot_clusters = F)

    res[i] = (inherits(opt_md, 'cluster medoids silhouette') && mean(unlist(lapply(opt_md, length))[-1]) == 6)
  }

  testthat::expect_true( sum(res) == length(tmp) )
})


testthat::test_that("in case of Clara_Medoids for different parameter settings, it returns the correct output", {

  tmp = c("silhouette","dissimilarity")

  res = rep(NA, length(tmp))

  for (i in 1:length(tmp)) {

    opt_md = Optimal_Clusters_Medoids(dat, max_clusters = 10, 'euclidean', tmp[i], clara_samples = 5, clara_sample_size = 0.2, plot_clusters = F)

    res[i] = (inherits(opt_md, 'cluster medoids silhouette') && mean(unlist(lapply(opt_md, length))[-1]) == 6)
  }

  testthat::expect_true( sum(res) == length(tmp) )
})

