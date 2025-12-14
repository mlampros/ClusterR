# Optimal number of Clusters for the gaussian mixture models

Optimal number of Clusters for the gaussian mixture models

## Usage

``` r
Optimal_Clusters_GMM(
  data,
  max_clusters,
  criterion = "AIC",
  dist_mode = "eucl_dist",
  seed_mode = "random_subset",
  km_iter = 10,
  em_iter = 5,
  verbose = FALSE,
  var_floor = 1e-10,
  plot_data = TRUE,
  seed = 1,
  full_covariance_matrices = FALSE
)
```

## Arguments

- data:

  matrix or data frame

- max_clusters:

  either a numeric value, a contiguous or non-continguous numeric vector
  specifying the cluster search space

- criterion:

  one of 'AIC' or 'BIC'

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

- plot_data:

  either TRUE or FALSE indicating whether the results of the function
  should be plotted

- seed:

  integer value for random number generator (RNG)

- full_covariance_matrices:

  a boolean. If FALSE "diagonal" covariance matrices (i.e. in each
  covariance matrix, all entries outside the main diagonal are assumed
  to be zero) otherwise "full" covariance matrices will be used. Note:
  when using full covariance matrices, the AIC/BIC calculation accounts
  for the increased number of parameters.

## Value

a vector with either the AIC or BIC for each iteration. In case of Error
it returns the error message and the possible causes.

## Details

**AIC** : the Akaike information criterion

**BIC** : the Bayesian information criterion

In case that the *max_clusters* parameter is a contiguous or
non-contiguous vector then plotting is disabled. Therefore, plotting is
enabled only if the *max_clusters* parameter is of length 1.

When *full_covariance_matrices* is TRUE, the AIC/BIC values will be
different from when it is FALSE because full covariance matrices have
more free parameters (k\*(d + d\*(d+1)/2)) compared to diagonal
covariance matrices (k\*2\*d), where k is the number of clusters and d
is the number of dimensions.

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)

opt_gmm = Optimal_Clusters_GMM(dat, 10, criterion = "AIC", plot_data = FALSE)


#----------------------------
# non-contiguous search space
#----------------------------

search_space = c(2,5)

opt_gmm = Optimal_Clusters_GMM(dat, search_space, criterion = "AIC", plot_data = FALSE)
```
