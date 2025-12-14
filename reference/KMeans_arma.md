# k-means using the Armadillo library

k-means using the Armadillo library

## Usage

``` r
KMeans_arma(
  data,
  clusters,
  n_iter = 10,
  seed_mode = "random_subset",
  verbose = FALSE,
  CENTROIDS = NULL,
  seed = 1
)
```

## Arguments

- data:

  matrix or data frame

- clusters:

  the number of clusters

- n_iter:

  the number of clustering iterations (about 10 is typically sufficient)

- seed_mode:

  how the initial centroids are seeded. One of, *keep_existing*,
  *static_subset*, *random_subset*, *static_spread*, *random_spread*.

- verbose:

  either TRUE or FALSE, indicating whether progress is printed during
  clustering

- CENTROIDS:

  a matrix of initial cluster centroids. The rows of the CENTROIDS
  matrix should be equal to the number of clusters and the columns
  should be equal to the columns of the data. CENTROIDS should be used
  in combination with seed_mode 'keep_existing'.

- seed:

  integer value for random number generator (RNG)

## Value

the centroids as a matrix. In case of Error it returns the error
message, whereas in case of an empty centroids-matrix it returns a
warning-message.

## Details

This function is an R implementation of the 'kmeans' class of the
Armadillo library. It is faster than the KMeans_rcpp function but it
lacks some features. For more info see the details section of the
KMeans_rcpp function. The number of columns should be larger than the
number of clusters or CENTROIDS. If the clustering fails, the means
matrix is reset and a bool set to false is returned. The clustering will
run faster on multi-core machines when OpenMP is enabled in your
compiler (eg. -fopenmp in GCC)

## References

http://arma.sourceforge.net/docs.html

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)

km = KMeans_arma(dat, clusters = 2, n_iter = 10, "random_subset")
```
