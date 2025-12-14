# Prediction function for the k-means

Prediction function for the k-means

## Usage

``` r
predict_KMeans(data, CENTROIDS, threads = 1, fuzzy = FALSE)

# S3 method for class 'KMeansCluster'
predict(object, newdata, fuzzy = FALSE, threads = 1, ...)
```

## Arguments

- data:

  matrix or data frame

- CENTROIDS:

  a matrix of initial cluster centroids. The rows of the CENTROIDS
  matrix should be equal to the number of clusters and the columns
  should be equal to the columns of the data.

- threads:

  an integer specifying the number of cores to run in parallel

- fuzzy:

  either TRUE or FALSE. If TRUE, then probabilities for each cluster
  will be returned based on the distance between observations and
  centroids.

- object, newdata, ...:

  arguments for the \`predict\` generic

## Value

a vector (clusters)

## Details

This function takes the data and the output centroids and returns the
clusters.

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)

km = KMeans_rcpp(dat, clusters = 2, num_init = 5, max_iters = 100, initializer = 'kmeans++')

pr = predict_KMeans(dat, km$centroids, threads = 1)
```
