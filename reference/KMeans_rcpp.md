# k-means using RcppArmadillo

k-means using RcppArmadillo

## Usage

``` r
KMeans_rcpp(
  data,
  clusters,
  num_init = 1,
  max_iters = 100,
  initializer = "kmeans++",
  fuzzy = FALSE,
  verbose = FALSE,
  CENTROIDS = NULL,
  tol = 1e-04,
  tol_optimal_init = 0.3,
  seed = 1
)
```

## Arguments

- data:

  matrix or data frame

- clusters:

  the number of clusters

- num_init:

  number of times the algorithm will be run with different centroid
  seeds

- max_iters:

  the maximum number of clustering iterations

- initializer:

  the method of initialization. One of, *optimal_init*, *quantile_init*,
  *kmeans++* and *random*. See details for more information

- fuzzy:

  either TRUE or FALSE. If TRUE, then prediction probabilities will be
  calculated using the distance between observations and centroids

- verbose:

  either TRUE or FALSE, indicating whether progress is printed during
  clustering.

- CENTROIDS:

  a matrix of initial cluster centroids. The rows of the CENTROIDS
  matrix should be equal to the number of clusters and the columns
  should be equal to the columns of the data.

- tol:

  a float number. If, in case of an iteration (iteration \> 1 and
  iteration \< max_iters) 'tol' is greater than the squared norm of the
  centroids, then kmeans has converged

- tol_optimal_init:

  tolerance value for the 'optimal_init' initializer. The higher this
  value is, the far appart from each other the centroids are.

- seed:

  integer value for random number generator (RNG)

## Value

a list with the following attributes: clusters, fuzzy_clusters (if fuzzy
= TRUE), centroids, total_SSE, best_initialization, WCSS_per_cluster,
obs_per_cluster, between.SS_DIV_total.SS

## Details

This function has the following features in comparison to the
KMeans_arma function:

Besides optimal_init, quantile_init, random and kmeans++ initilizations
one can specify the centroids using the CENTROIDS parameter.

The running time and convergence of the algorithm can be adjusted using
the num_init, max_iters and tol parameters.

If num_init \> 1 then KMeans_rcpp returns the attributes of the best
initialization using as criterion the
within-cluster-sum-of-squared-error.

—————initializers———————-

**optimal_init** : this initializer adds rows of the data incrementally,
while checking that they do not already exist in the centroid-matrix \[
experimental \]

**quantile_init** : initialization of centroids by using the cummulative
distance between observations and by removing potential duplicates \[
experimental \]

**kmeans++** : kmeans++ initialization. Reference :
http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf AND
http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work

**random** : random selection of data rows as initial centroids

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)

km = KMeans_rcpp(dat, clusters = 2, num_init = 5, max_iters = 100, initializer = 'kmeans++')
```
