# Plot of silhouette widths or dissimilarities

Plot of silhouette widths or dissimilarities

## Usage

``` r
Silhouette_Dissimilarity_Plot(evaluation_object, silhouette = TRUE)
```

## Arguments

- evaluation_object:

  the output of either a *Cluster_Medoids* or *Clara_Medoids* function

- silhouette:

  either TRUE or FALSE, indicating whether the silhouette widths or the
  dissimilarities should be plotted

## Value

TRUE if either the silhouette widths or the dissimilarities are plotted
successfully, otherwise FALSE

## Details

This function takes the result-object of the *Cluster_Medoids* or
*Clara_Medoids* function and depending on the argument *silhouette* it
plots either the dissimilarities or the silhouette widths of the
observations belonging to each cluster.

## Author

Lampros Mouselimis

## Examples

``` r
# data(soybean)

# dat = soybean[, -ncol(soybean)]

# cm = Cluster_Medoids(dat, clusters = 5, distance_metric = 'jaccard_coefficient')

# plt_sd = Silhouette_Dissimilarity_Plot(cm, silhouette = TRUE)
```
