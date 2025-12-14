# Predictions for the Medoid functions

Predictions for the Medoid functions

## Usage

``` r
predict_Medoids(
  data,
  MEDOIDS = NULL,
  distance_metric = "euclidean",
  fuzzy = FALSE,
  minkowski_p = 1,
  threads = 1
)

# S3 method for class 'MedoidsCluster'
predict(object, newdata, fuzzy = FALSE, threads = 1, ...)
```

## Arguments

- data:

  matrix or data frame

- MEDOIDS:

  a matrix of initial cluster medoids (data observations). The rows of
  the MEDOIDS matrix should be equal to the number of clusters and the
  columns of the MEDOIDS matrix should be equal to the columns of the
  data.

- distance_metric:

  a string specifying the distance method. One of, *euclidean*,
  *manhattan*, *chebyshev*, *canberra*, *braycurtis*,
  *pearson_correlation*, *simple_matching_coefficient*, *minkowski*,
  *hamming*, *jaccard_coefficient*, *Rao_coefficient*, *mahalanobis*,
  *cosine*

- fuzzy:

  either TRUE or FALSE. If TRUE, then probabilities for each cluster
  will be returned based on the distance between observations and
  medoids.

- minkowski_p:

  a numeric value specifying the minkowski parameter in case that
  distance_metric = "minkowski"

- threads:

  an integer specifying the number of cores to run in parallel. Openmp
  will be utilized to parallelize the number of initializations
  (num_init)

- object, newdata, ...:

  arguments for the \`predict\` generic

## Value

a list with the following attributes will be returned : clusters,
fuzzy_clusters (if fuzzy = TRUE), dissimilarity.

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)

cm = Cluster_Medoids(dat, clusters = 3, distance_metric = 'euclidean', swap_phase = TRUE)

pm = predict_Medoids(dat, MEDOIDS = cm$medoids, 'euclidean', fuzzy = TRUE)
```
