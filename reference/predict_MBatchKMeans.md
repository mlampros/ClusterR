# Prediction function for Mini-Batch-k-means

Prediction function for Mini-Batch-k-means

## Usage

``` r
predict_MBatchKMeans(data, CENTROIDS, fuzzy = FALSE, updated_output = FALSE)

# S3 method for class 'MBatchKMeans'
predict(object, newdata, fuzzy = FALSE, ...)
```

## Arguments

- data:

  matrix or data frame

- CENTROIDS:

  a matrix of initial cluster centroids. The rows of the CENTROIDS
  matrix should be equal to the number of clusters and the columns
  should equal the columns of the data.

- fuzzy:

  either TRUE or FALSE. If TRUE then prediction probabilities will be
  calculated using the distance between observations and centroids.

- updated_output:

  either TRUE or FALSE. If TRUE then the 'predict_MBatchKMeans' function
  will follow the same output object behaviour as the 'predict_KMeans'
  function (if fuzzy is TRUE it will return probabilities otherwise it
  will return the hard clusters). This parameter will be removed in
  version 1.4.0 because this will become the default output format.

- object, newdata, ...:

  arguments for the \`predict\` generic

## Value

if fuzzy = TRUE the function returns a list with two attributes: a
vector with the clusters and a matrix with cluster probabilities.
Otherwise, it returns a vector with the clusters.

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

MbatchKm = MiniBatchKmeans(dat, clusters = 2, batch_size = 20, num_init = 5, early_stop_iter = 10)

pr = predict_MBatchKMeans(dat, MbatchKm$centroids, fuzzy = FALSE)
#> Warning: `predict_MBatchKMeans()` was deprecated in ClusterR 1.3.0.
#> ℹ Beginning from version 1.4.0, if the fuzzy parameter is TRUE the function
#>   'predict_MBatchKMeans' will return only the probabilities, whereas currently
#>   it also returns the hard clusters
```
