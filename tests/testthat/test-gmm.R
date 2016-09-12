#============================================================

# data

data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

X = center_scale(dat)


# parameters for the predict_GMM function

CENTROIDS = matrix(runif(84), nrow = 2, ncol = 42)

COVARIANCE = matrix(runif(84), nrow = 2, ncol = 42)

WEIGHTS = c(0.5, 0.5)

#=============================================================


context('gaussian mixture models')


##############################
# error handling GMM function
##############################

testthat::test_that("in case that the data is not a matrix or data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( GMM(tmp_x, 2, "maha_dist", "random_subset", 10, 10)  )
})


testthat::test_that("in case that the number of gaussian mixture components is less than or equal to 0, it returns an error", {
  
  testthat::expect_error( GMM(X, 0, "maha_dist", "random_subset", 10, 10)  )
})


testthat::test_that("in case that the dist_mode is not one of 'eucl_dist', 'maha_dist', it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "unknown_dist", "random_subset", 10, 10)  )
})


testthat::test_that("in case that the seed_mode is not one of 'static_subset','random_subset','static_spread','random_spread', it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "unknown_subset", 10, 10)  )
})


testthat::test_that("in case that the km_iter is negative, it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "random_subset", -1, 10)  )
})


testthat::test_that("in case that the em_iter is negative, it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "random_subset", 10, -1)  )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "random_subset", 10, 10, verbose = 'NA')  )
})


testthat::test_that("in case that the var_floor parameter is negative, it returns an error", {
  
  testthat::expect_error( GMM(X, 2, "maha_dist", "random_subset", 10, 10, var_floor = -1)  )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = NaN

  testthat::expect_error( GMM(tmp_dat, 2, "maha_dist", "random_subset", 10, 10)  )
})


#################
# GMM function
#################


testthat::test_that("in case that the data is a matrix the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(X, 2, "maha_dist", "random_subset", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})


testthat::test_that("in case that the data is a data frame the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(dat, 2, "eucl_dist", "static_subset", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})


testthat::test_that("in case that the data is a matrix the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(X, 2, "maha_dist", "random_subset", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})


testthat::test_that("in case that the data is a data frame the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(dat, 2, "eucl_dist", "static_subset", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})



testthat::test_that("in case that the data is a matrix the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(X, 2, "maha_dist", "random_spread", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})


testthat::test_that("in case that the data is a data frame the result is a list of length 4 and the class is 'Gaussian Mixture Models' ", {
  
  res = GMM(dat, 2, "eucl_dist", "static_spread", 10, 10)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == 4 && sum(names(res) %in% c("centroids", "covariance_matrices", "weights", "Log_likelihood")) == 4 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 4 && class(res) == "Gaussian Mixture Models" )
  }
})



######################################
# error handling predict_GMM function
######################################


testthat::test_that("in case that the data is not a matrix or data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( predict_GMM(tmp_x, CENTROIDS, COVARIANCE, WEIGHTS)  )
})


testthat::test_that("in case that the data is a matrix AND the CENTROIDS is NOT a matrix or data frame, it returns an error", {
  
  tmp_c = list(CENTROIDS)

  testthat::expect_error( predict_GMM(X, tmp_c, COVARIANCE, WEIGHTS) )
})


testthat::test_that("in case that the data is a matrix AND the COVARIANCE is NOT a matrix or data frame, it returns an error", {
  
  tmp_c = list(COVARIANCE)
  
  testthat::expect_error( predict_GMM(X, CENTROIDS, tmp_c, WEIGHTS) )
})


testthat::test_that("in case that the number of columns of the data do not equal the number of columns of the COVARIANCE or CENTROIDS matrices, it returns an error", {

  testthat::expect_error( predict_GMM(X[, -1], CENTROIDS, COVARIANCE, WEIGHTS) )
})


testthat::test_that("in case that the length of WEIGHTS vector does not equal the number of rows of the COVARIANCE or CENTROIDS matrices, it returns an error", {
  
  testthat::expect_error( predict_GMM(X, CENTROIDS, COVARIANCE, WEIGHTS[-1]) )
})


testthat::test_that("in case that the class of the WEIGHTS vector is not numeric, it returns an error", {
  
  tmp_w = matrix(WEIGHTS, nrow = 1)
  
  testthat::expect_error( predict_GMM(X, CENTROIDS, COVARIANCE, tmp_w) )
})


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = Inf
  
  testthat::expect_error( predict_GMM(tmp_dat, CENTROIDS, COVARIANCE, WEIGHTS) )
})


#######################
# predict_GMM function
#######################


testthat::test_that("in case that the data is a matrix the result is a list of length 3 and the class is 'Gaussian Mixture Models' ", {
  
  res = predict_GMM(X, CENTROIDS, COVARIANCE, WEIGHTS)
  
  testthat::expect_true( length(res) == 3 && sum(names(res) %in% c("log_likelihood", "cluster_proba", "cluster_labels")) == 3 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 3 && class(res) == "Gaussian Mixture Models")
})


testthat::test_that("in case that the data is a data frame the result is a list of length 3 and the class is 'Gaussian Mixture Models' ", {
  
  res = predict_GMM(dat, CENTROIDS, COVARIANCE, WEIGHTS)
  
  testthat::expect_true( length(res) == 3 && sum(names(res) %in% c("log_likelihood", "cluster_proba", "cluster_labels")) == 3 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 3 && class(res) == "Gaussian Mixture Models" )
})


testthat::test_that("in case that the data is a matrix AND the CENTROIDS is a data frame the result is a list of length 3 and the class is 'Gaussian Mixture Models' ", {
  
  tmp_c = data.frame(CENTROIDS)
  
  res = predict_GMM(X, tmp_c, COVARIANCE, WEIGHTS)
  
  testthat::expect_true( length(res) == 3 && sum(names(res) %in% c("log_likelihood", "cluster_proba", "cluster_labels")) == 3 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 3 && class(res) == "Gaussian Mixture Models")
})


testthat::test_that("in case that the data is a matrix AND the COVARIANCE is a data frame the result is a list of length 3 and the class is 'Gaussian Mixture Models' ", {
  
  tmp_c = data.frame(COVARIANCE)
  
  res = predict_GMM(X, CENTROIDS, tmp_c, WEIGHTS)
  
  testthat::expect_true( length(res) == 3 && sum(names(res) %in% c("log_likelihood", "cluster_proba", "cluster_labels")) == 3 &&
                           
                           sum(unlist(lapply(res, function(x) !is.null(x)))) == 3 && class(res) == "Gaussian Mixture Models")
})



###############################################
# error handling Optimal_Clusters_GMM function
###############################################


testthat::test_that("in case that the data is not a matrix or a data frame, it returns an error", {
  
  tmp_x = list(X)
  
  testthat::expect_error( Optimal_Clusters_GMM(tmp_x, Nr_clusters, criterion = "BIC", plot_data = FALSE) )
})


testthat::test_that("in case that the max_clusters parameter is less than 2 and plot_data is TRUE, it returns an error", {

  testthat::expect_error( Optimal_Clusters_GMM(X, 1, criterion = "BIC", plot_data = T) )
})


testthat::test_that("in case that the max_clusters parameter is not numeric, it returns an error", {
  
  tmp_m = data.frame(1)
  
  testthat::expect_error( Optimal_Clusters_GMM(X, tmp_m, criterion = "BIC", plot_data = F) )
})


testthat::test_that("in case that the length of the max_clusters parameter is not 1, it returns an error", {
  
  tmp_m = c(1,2)
  
  testthat::expect_error( Optimal_Clusters_GMM(X, tmp_m, criterion = "BIC", plot_data = F) )
})


testthat::test_that("in case that the max_clusters parameter is less than 1, it returns an error", {
  
  tmp_m = 0
  
  testthat::expect_error( Optimal_Clusters_GMM(X, tmp_m, criterion = "BIC", plot_data = F) )
})


testthat::test_that("in case that the criterion parameter is not 'AIC', 'BIC', it returns an error", {

  testthat::expect_error( Optimal_Clusters_GMM(X, 5, criterion = "invalid", plot_data = F) )
})



testthat::test_that("in case that the dist_mode parameter is not 'eucl_dist', 'maha_dist', it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'invalid', criterion = "BIC", plot_data = F) )
})



testthat::test_that("in case that the seed_mode parameter is not 'static_subset','random_subset','static_spread','random_spread', it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'invalid', plot_data = F) )
})



testthat::test_that("in case that the km_iter is less than 0, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = -1, plot_data = F) )
})



testthat::test_that("in case that the em_iter is less than 0, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = 5, 
                                               
                                               em_iter = -1, plot_data = F) )
})



testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = 5, 
                                               
                                               em_iter = 5, plot_data = F, verbose = 'invalid') )
})


testthat::test_that("in case that the verbose parameter is not logical, it returns an error", {
  
  testthat::expect_error( Optimal_Clusters_GMM(X, 5, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = 5, 
                                               
                                               em_iter = 5, plot_data = F, var_floor = -1) )
})


# testthat::test_that("in case that max_clusters is greater than the number of the columns and verbose is TRUE, it returns a warning", {
#   
#   testthat::expect_warning( Optimal_Clusters_GMM(X, ncol(X) + 1, dist_mode = 'eucl_dist', criterion = "BIC", seed_mode = 'random_subset', km_iter = 5, 
#                                                  
#                                                  em_iter = 5, plot_data = F, verbose = T) )
# })


testthat::test_that("in case that the data includes NaN or Inf values, it returns an error", {
  
  tmp_dat = X
  
  tmp_dat[1,1] = -Inf
  
  Nr_clusters = 5
  
  testthat::expect_error( Optimal_Clusters_GMM(tmp_dat, Nr_clusters, criterion = "BIC", plot_data = FALSE) )
})


################################
# Optimal_Clusters_GMM function
################################


testthat::test_that("in case that the data is a matrix the result is a vector and the class is 'Gaussian Mixture Models' ", {

  Nr_clusters = 5

  res = Optimal_Clusters_GMM(X, Nr_clusters, criterion = "BIC", plot_data = FALSE)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == Nr_clusters && class(res) == "Gaussian Mixture Models" )
  }
})



testthat::test_that("in case that the data is a data frame the result is a vector and the class is 'Gaussian Mixture Models' ", {

  Nr_clusters = 5

  res = Optimal_Clusters_GMM(dat, Nr_clusters, criterion = "BIC", plot_data = T)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == Nr_clusters && class(res) == "Gaussian Mixture Models" )
  }
})



testthat::test_that("in case of different parameters the result is a vector and the class is 'Gaussian Mixture Models' ", {

  Nr_clusters = 5

  res = Optimal_Clusters_GMM(dat, Nr_clusters, criterion = "AIC", dist_mode = 'maha_dist', seed_mode = 'static_spread', plot_data = T)
  
  if ('Error' %in% names(res)) {
    
    testthat::expect_true( length(res) == 2)}
  
  else {
    
    testthat::expect_true( length(res) == Nr_clusters && class(res) == "Gaussian Mixture Models" )
  }
})

