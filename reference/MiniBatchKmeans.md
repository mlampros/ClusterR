# Mini-batch-k-means using RcppArmadillo

Mini-batch-k-means using RcppArmadillo

## Usage

``` r
MiniBatchKmeans(
  data,
  clusters,
  batch_size = 10,
  num_init = 1,
  max_iters = 100,
  init_fraction = 1,
  initializer = "kmeans++",
  early_stop_iter = 10,
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

- batch_size:

  the size of the mini batches

- num_init:

  number of times the algorithm will be run with different centroid
  seeds

- max_iters:

  the maximum number of clustering iterations

- init_fraction:

  percentage of data to use for the initialization centroids (applies if
  initializer is *kmeans++* or *optimal_init*). Should be a float number
  between 0.0 and 1.0.

- initializer:

  the method of initialization. One of, *optimal_init*, *quantile_init*,
  *kmeans++* and *random*. See details for more information

- early_stop_iter:

  continue that many iterations after calculation of the best
  within-cluster-sum-of-squared-error

- verbose:

  either TRUE or FALSE, indicating whether progress is printed during
  clustering

- CENTROIDS:

  a matrix of initial cluster centroids. The rows of the CENTROIDS
  matrix should be equal to the number of clusters and the columns
  should be equal to the columns of the data

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

a list with the following attributes: centroids, WCSS_per_cluster,
best_initialization, iters_per_initialization

## Details

This function performs k-means clustering using mini batches.

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

## References

http://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf,
https://github.com/siddharth-agrawal/Mini-Batch-K-Means

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)

MbatchKm = MiniBatchKmeans(dat, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10)
```
