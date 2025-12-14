# Partitioning around medoids

Partitioning around medoids

## Usage

``` r
Cluster_Medoids(
  data,
  clusters,
  distance_metric = "euclidean",
  minkowski_p = 1,
  threads = 1,
  swap_phase = TRUE,
  fuzzy = FALSE,
  verbose = FALSE,
  seed = 1
)
```

## Arguments

- data:

  matrix or data frame. The data parameter can be also a dissimilarity
  matrix, where the main diagonal equals 0.0 and the number of rows
  equals the number of columns

- clusters:

  the number of clusters

- distance_metric:

  a string specifying the distance method. One of, *euclidean*,
  *manhattan*, *chebyshev*, *canberra*, *braycurtis*,
  *pearson_correlation*, *simple_matching_coefficient*, *minkowski*,
  *hamming*, *jaccard_coefficient*, *Rao_coefficient*, *mahalanobis*,
  *cosine*

- minkowski_p:

  a numeric value specifying the minkowski parameter in case that
  distance_metric = "minkowski"

- threads:

  an integer specifying the number of cores to run in parallel

- swap_phase:

  either TRUE or FALSE. If TRUE then both phases ('build' and 'swap')
  will take place. The 'swap_phase' is considered more computationally
  intensive.

- fuzzy:

  either TRUE or FALSE. If TRUE, then probabilities for each cluster
  will be returned based on the distance between observations and
  medoids

- verbose:

  either TRUE or FALSE, indicating whether progress is printed during
  clustering

- seed:

  \`r lifecycle::badge("deprecated")\` \`seed\` (integer value for
  random number generator (RNG)) is no longer supported and will be
  removed in version 1.4.0

## Value

a list with the following attributes: medoids, medoid_indices,
best_dissimilarity, dissimilarity_matrix, clusters, fuzzy_probs (if
fuzzy = TRUE), silhouette_matrix, clustering_stats

## Details

Due to the fact that I didn't have access to the book 'Finding Groups in
Data, Kaufman and Rousseeuw, 1990' (which includes the exact algorithm)
I implemented the 'Cluster_Medoids' function based on the paper
'Clustering in an Object-Oriented Environment' (see 'References').
Therefore, the 'Cluster_Medoids' function is an approximate
implementation and not an exact one. Furthermore, in comparison to
k-means clustering, the function 'Cluster_Medoids' is more robust,
because it minimizes the sum of unsquared dissimilarities. Moreover, it
doesn't need initial guesses for the cluster centers.

## References

Anja Struyf, Mia Hubert, Peter J. Rousseeuw, (Feb. 1997), Clustering in
an Object-Oriented Environment, Journal of Statistical Software, Vol 1,
Issue 4

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)

cm = Cluster_Medoids(dat, clusters = 3, distance_metric = 'euclidean', swap_phase = TRUE)
#> Warning: The `seed` argument of `Cluster_Medoids()` is deprecated as of ClusterR 1.2.6.
#> ℹ The 'seed' parameter will be removed in version 1.4.0
```
