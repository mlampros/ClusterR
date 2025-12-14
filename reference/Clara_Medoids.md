# Clustering large applications

Clustering large applications

## Usage

``` r
Clara_Medoids(
  data,
  clusters,
  samples,
  sample_size,
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

  matrix or data frame

- clusters:

  the number of clusters

- samples:

  number of samples to draw from the data set

- sample_size:

  fraction of data to draw in each sample iteration. It should be a
  float number greater than 0.0 and less or equal to 1.0

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

  an integer specifying the number of cores to run in parallel. Openmp
  will be utilized to parallelize the number of the different sample
  draws

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

  integer value for random number generator (RNG)

## Value

a list with the following attributes : medoids, medoid_indices,
sample_indices, best_dissimilarity, clusters, fuzzy_probs (if fuzzy =
TRUE), clustering_stats, dissimilarity_matrix, silhouette_matrix

## Details

The Clara_Medoids function is implemented in the same way as the 'clara'
(clustering large applications) algorithm (Kaufman and Rousseeuw(1990)).
In the 'Clara_Medoids' the 'Cluster_Medoids' function will be applied to
each sample draw.

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

clm = Clara_Medoids(dat, clusters = 3, samples = 5, sample_size = 0.2, swap_phase = TRUE)
```
