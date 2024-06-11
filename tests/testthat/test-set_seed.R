data(dietary_survey_IBS)
dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
X = center_scale(dat)

test_that("Default seed is 1", {
  km_seed_default = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F)
  runif_seed_default = runif(3L)

  km_seed_1 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = 1)
  runif_seed_1 = runif(3L)

  expect_identical(km_seed_default, km_seed_1)
  expect_identical(runif_seed_default, runif_seed_1)
})

test_that("set.seed() is skipped if seed=NA is given", {
  km_seed_1 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F)
  runif_seed_1 = runif(3L)

  km_seed_2 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = 2)
  runif_seed_2 = runif(3L)

  set.seed(2L)
  km_seed_na = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = NA)
  runif_seed_na = runif(3L)

  expect_false(identical(km_seed_1, km_seed_2))
  expect_false(identical(runif_seed_1, runif_seed_2))
  expect_identical(km_seed_2, km_seed_na)
  expect_identical(runif_seed_2, runif_seed_na)
})

test_that("reproducible results are obtained with seed=NA", {
  set.seed(2L)
  km_a1 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = NA)
  km_a2 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = NA)
  km_a3 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = NA)
  runif_a = runif(3L)

  set.seed(2L)
  km_b1 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = NA)
  km_b2 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = NA)
  km_b3 = KMeans_arma(X, clusters = 2, n_iter = 2, "random_subset", verbose = F, seed = NA)
  runif_b = runif(3L)

  expect_false(identical(km_a1, km_a2))
  expect_false(identical(km_a2, km_a3))
  expect_false(identical(km_a3, km_a1))
  expect_identical(km_a1, km_b1)
  expect_identical(km_a2, km_b2)
  expect_identical(km_a3, km_b3)
  expect_identical(runif_a, runif_a)
})
