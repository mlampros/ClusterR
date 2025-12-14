# Gaussian Mixture Model clustering

Gaussian Mixture Model clustering

## Usage

``` r
GMM(
  data,
  gaussian_comps = 1,
  dist_mode = "eucl_dist",
  seed_mode = "random_subset",
  km_iter = 10,
  em_iter = 5,
  verbose = FALSE,
  var_floor = 1e-10,
  seed = 1,
  full_covariance_matrices = FALSE
)
```

## Arguments

- data:

  matrix or data frame

- gaussian_comps:

  the number of gaussian mixture components

- dist_mode:

  the distance used during the seeding of initial means and k-means
  clustering. One of, *eucl_dist*, *maha_dist*.

- seed_mode:

  how the initial means are seeded prior to running k-means and/or EM
  algorithms. One of, *static_subset*, *random_subset*, *static_spread*,
  *random_spread*.

- km_iter:

  the number of iterations of the k-means algorithm

- em_iter:

  the number of iterations of the EM algorithm

- verbose:

  either TRUE or FALSE; enable or disable printing of progress during
  the k-means and EM algorithms

- var_floor:

  the variance floor (smallest allowed value) for the diagonal
  covariances

- seed:

  integer value for random number generator (RNG)

- full_covariance_matrices:

  a boolean. If FALSE "diagonal" covariance matrices (i.e. in each
  covariance matrix, all entries outside the main diagonal are assumed
  to be zero) otherwise "full" covariance matrices will be returned. Be
  aware in case of "full" covariance matrices a cube (3-dimensional)
  rather than a matrix for the output "covariance_matrices" value will
  be returned.

## Value

a list consisting of the centroids, covariance matrix ( where each row
of the matrix represents a diagonal covariance matrix), weights and the
log-likelihoods for each gaussian component. In case of Error it returns
the error message and the possible causes.

## Details

This function is an R implementation of the 'gmm_diag' class of the
Armadillo library. The only exception is that user defined parameter
settings are not supported, such as seed_mode = 'keep_existing'. For
probabilistic applications, better model parameters are typically
learned with dist_mode set to maha_dist. For vector quantisation
applications, model parameters should be learned with dist_mode set to
eucl_dist, and the number of EM iterations set to zero. In general, a
sufficient number of k-means and EM iterations is typically about 10.
The number of training samples should be much larger than the number of
Gaussians. Seeding the initial means with static_spread and
random_spread can be much more time consuming than with static_subset
and random_subset. The k-means and EM algorithms will run faster on
multi-core machines when OpenMP is enabled in your compiler (eg.
-fopenmp in GCC)

## References

http://arma.sourceforge.net/docs.html

## Examples

``` r
data(dietary_survey_IBS)

dat = as.matrix(dietary_survey_IBS[, -ncol(dietary_survey_IBS)])

dat = center_scale(dat)

gmm = GMM(dat, 2, "maha_dist", "random_subset", 10, 10)
```
