# Compute the cost and clusters based on an input dissimilarity matrix and medoids

Compute the cost and clusters based on an input dissimilarity matrix and
medoids

## Usage

``` r
cost_clusters_from_dissim_medoids(data, medoids)
```

## Arguments

- data:

  a dissimilarity matrix, where the main diagonal equals 0.0 and the
  number of rows equals the number of columns

- medoids:

  a vector of output medoids of the 'Cluster_Medoids', 'Clara_Medoids'
  or any other 'partition around medoids' function

## Value

a list object that includes the cost and the clusters

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)
dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]
dat = center_scale(dat)

cm = Cluster_Medoids(dat, clusters = 3, distance_metric = 'euclidean', swap_phase = TRUE)
res = cost_clusters_from_dissim_medoids(data = cm$dissimilarity_matrix, medoids = cm$medoid_indices)

# cm$best_dissimilarity == res$cost
# table(cm$clusters, res$clusters)
```
