
#------
# data
#------

set.seed(1)
sm_data = matrix(runif(1000), ncol = 10, nrow = 100)
sm_mt = 1 - distance_matrix(sm_data, upper = T, diagonal = T)
diag(sm_mt) = 0.0

context('affinity propagation functions')


#--------------------------
# affinity propagation (AP)
#--------------------------

testthat::test_that("the AP function returns a error if the input data is not of type matrix", {
  
  testthat::expect_error( AP_affinity_propagation(data = list(sm_mt), p = median(as.vector(sm_mt))) )
})


testthat::test_that("the AP function returns error-free", {
  
  ap = AP_affinity_propagation(data = sm_mt, p = median(as.vector(sm_mt)))
  
  testthat::expect_true( all(names(ap) %in% c("K", "N", "netsim", "dpsim", "expref", "iterations", "exemplars", "idx", "clusters", "clusters_vectorized")) )
})



#-----------------------------------------------------
# preference range value for affinity propagation (AP)
#-----------------------------------------------------

testthat::test_that("the AP_preferenceRange function returns a error if the input data is not of type matrix", {
  
  testthat::expect_error( AP_preferenceRange(data = list(sm_mt), method = "bound") )
})


testthat::test_that("the AP_preferenceRange function returns a error if method is not one of 'bound', 'exact'", {
  
  testthat::expect_error( AP_preferenceRange(data = sm_mt, method = "invalid") )
})


testthat::test_that("the AP_preferenceRange function returns the correct output if method is 'bound'", {
  
  vec_ap = AP_preferenceRange(data = sm_mt, method = "bound")
  
  testthat::expect_true( inherits(vec_ap, 'numeric') && length(vec_ap) == 2 )
})


testthat::test_that("the AP_preferenceRange function returns the correct output if method is 'exact' (both for single and multi-threaded version)", {
  
  vec_ap1 = AP_preferenceRange(data = sm_mt, method = "exact", threads = 1)
  vec_ap2 = AP_preferenceRange(data = sm_mt, method = "exact", threads = 2)
  
  testthat::expect_true( inherits(vec_ap1, 'numeric') && length(vec_ap1) == 2 && inherits(vec_ap2, 'numeric') && length(vec_ap2) == 2 && all(vec_ap1 == vec_ap2) )
})
