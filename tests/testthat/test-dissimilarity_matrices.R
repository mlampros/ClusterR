#============================================================

# data

set.seed(0)
binary_data = do.call(cbind, lapply(1:10, function(x) sample(0:1, 100, replace = T)))

set.seed(1)
continuous_data = matrix(runif(1000), ncol = 10, nrow = 100)

set.seed(2)
MEDOIDS_continuous = matrix(runif(20), ncol = ncol(binary_data), nrow = 2)

set.seed(3)
MEDOIDS_binary = matrix(sample(0:1, 20, replace = T), ncol = ncol(binary_data), nrow = 2)

#=============================================================


context('dissimilarity - matrices')



######################
# dissim_mat function
######################

testthat::test_that("in case that the data is binary it returns the correct output for the binary methods", {

  binary_methods = c("simple_matching_coefficient", "hamming", "jaccard_coefficient", "Rao_coefficient")
  
  res = rep(NA, length(binary_methods))
  
  for (i in 1:length(binary_methods)) {
    
    out = dissim_mat(binary_data, binary_methods[i], upper = T, diagonal = T)
    
    res[i] = is.matrix(out) && nrow(out) == ncol(out) && nrow(out) == nrow(binary_data) && ncol(out) == nrow(binary_data)
  }
  
  testthat::expect_true( length(binary_methods) == sum(res) )
})


testthat::test_that("in case that the data is numeric it returns the correct output for the numeric methods", {
  
  continuous_methods = c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis", "pearson_correlation", "mahalanobis")
  
  res = rep(NA, length(continuous_methods))
  
  for (i in 1:length(continuous_methods)) {
    
    out = dissim_mat(continuous_data, continuous_methods[i], upper = T, diagonal = T)
    
    res[i] = is.matrix(out) && nrow(out) == ncol(out) && nrow(out) == nrow(continuous_data) && ncol(out) == nrow(continuous_data)
  }
  
  testthat::expect_true( length(continuous_methods) == sum(res) )
})


testthat::test_that("in case that the data is numeric it returns the correct output for the minkowski method", {
  
  out = dissim_mat(continuous_data, "minkowski", upper = T, diagonal = T, minkowski_p = 1.0)
  
  testthat::expect_true( is.matrix(out) && nrow(out) == ncol(out) && nrow(out) == nrow(continuous_data) && ncol(out) == nrow(continuous_data) )
})



##########################
# dissim_MEDOIDS function
##########################


testthat::test_that("in case that the data is binary it returns the correct output for the binary methods", {
  
  binary_methods = c("simple_matching_coefficient", "hamming", "jaccard_coefficient", "Rao_coefficient")
  
  res = rep(NA, length(binary_methods))
  
  for (i in 1:length(binary_methods)) {
    
    out = dissim_MEDOIDS(binary_data, binary_methods[i], MEDOIDS_binary)
    
    res[i] = is.matrix(out) && nrow(out) == nrow(binary_data) && ncol(out) == nrow(MEDOIDS_binary)
  }
  
  testthat::expect_true( length(binary_methods) == sum(res) )
})


testthat::test_that("in case that the data is numeric it returns the correct output for the numeric methods", {
  
  continuous_methods = c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis", "pearson_correlation", "mahalanobis")
  
  res = rep(NA, length(continuous_methods))
  
  for (i in 1:length(continuous_methods)) {
    
    out = dissim_MEDOIDS(continuous_data, continuous_methods[i], MEDOIDS_continuous)
    
    res[i] = is.matrix(out) && nrow(out) == nrow(continuous_data) && ncol(out) == nrow(MEDOIDS_continuous)
  }
  
  testthat::expect_true( length(continuous_methods) == sum(res) )
})


testthat::test_that("in case that the data is numeric it returns the correct output for the minkowski method", {
  
  out = dissim_MEDOIDS(continuous_data, "minkowski", MEDOIDS_continuous, minkowski_p = 1.0)
  
  testthat::expect_true( is.matrix(out) && nrow(out) == nrow(continuous_data) && ncol(out) == nrow(MEDOIDS_continuous) )
})

