# Silhouette width based on pre-computed clusters

Silhouette width based on pre-computed clusters

## Usage

``` r
silhouette_of_clusters(data, clusters)
```

## Arguments

- data:

  a matrix or a data frame

- clusters:

  a numeric vector which corresponds to the pre-computed clusters (see
  the example section for more details). The size of the 'clusters'
  vector must be equal to the number of rows of the input data

## Value

a list object where the first sublist is the 'silhouette summary', the
second sublist is the 'silhouette matrix' and the third sublist is the
'global average silhouette' (based on the silhouette values of all
observations)

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)
dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
dat = center_scale(dat)

clusters = 2

# compute k-means
km = KMeans_rcpp(dat, clusters = clusters, num_init = 5, max_iters = 100, initializer = 'kmeans++')

# compute the silhouette width
silh_km = silhouette_of_clusters(data = dat, clusters = km$clusters)

# silhouette summary
silh_summary = silh_km$silhouette_summary

# silhouette matrix (including cluster & dissimilarity)
silh_mtrx = silh_km$silhouette_matrix

# global average silhouette
glob_avg = silh_km$silhouette_global_average
```
