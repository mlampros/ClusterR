# Optimal number of Clusters for the partitioning around Medoids functions

Optimal number of Clusters for the partitioning around Medoids functions

## Usage

``` r
Optimal_Clusters_Medoids(
  data,
  max_clusters,
  distance_metric,
  criterion = "dissimilarity",
  clara_samples = 0,
  clara_sample_size = 0,
  minkowski_p = 1,
  swap_phase = TRUE,
  threads = 1,
  verbose = FALSE,
  plot_clusters = TRUE,
  seed = 1
)
```

## Arguments

- data:

  matrix or data.frame. If both clara_samples and clara_sample_size
  equal 0, then the data parameter can be also a dissimilarity matrix,
  where the main diagonal equals 0.0 and the number of rows equals the
  number of columns

- max_clusters:

  either a numeric value, a contiguous or non-continguous numeric vector
  specifying the cluster search space

- distance_metric:

  a string specifying the distance method. One of, *euclidean*,
  *manhattan*, *chebyshev*, *canberra*, *braycurtis*,
  *pearson_correlation*, *simple_matching_coefficient*, *minkowski*,
  *hamming*, *jaccard_coefficient*, *Rao_coefficient*, *mahalanobis*,
  *cosine*

- criterion:

  one of 'dissimilarity' or 'silhouette'

- clara_samples:

  number of samples to draw from the data set in case of clustering
  large applications (clara)

- clara_sample_size:

  fraction of data to draw in each sample iteration in case of
  clustering large applications (clara). It should be a float number
  greater than 0.0 and less or equal to 1.0

- minkowski_p:

  a numeric value specifying the minkowski parameter in case that
  distance_metric = "minkowski"

- swap_phase:

  either TRUE or FALSE. If TRUE then both phases ('build' and 'swap')
  will take place. The 'swap_phase' is considered more computationally
  intensive.

- threads:

  an integer specifying the number of cores to run in parallel. Openmp
  will be utilized to parallelize the number of sample draws

- verbose:

  either TRUE or FALSE, indicating whether progress is printed during
  clustering

- plot_clusters:

  TRUE or FALSE, indicating whether the iterative results should be
  plotted. See the details section for more information

- seed:

  integer value for random number generator (RNG)

## Value

a list of length equal to the max_clusters parameter (the first sublist
equals NULL, as dissimilarities and silhouette widths can be calculated
if the number of clusters \> 1). If plot_clusters is TRUE then the
function plots also the results.

## Details

In case of plot_clusters = TRUE, the first plot will be either a plot of
dissimilarities or both dissimilarities and silhouette widths giving an
indication of the optimal number of the clusters. Then, the user will be
asked to give an optimal value for the number of the clusters and after
that the second plot will appear with either the dissimilarities or the
silhouette widths belonging to each cluster.

In case that the *max_clusters* parameter is a contiguous or
non-contiguous vector then plotting is disabled. Therefore, plotting is
enabled only if the *max_clusters* parameter is of length 1.

## Author

Lampros Mouselimis

## Examples

``` r
if (FALSE) { # \dontrun{
data(soybean)

dat = soybean[, -ncol(soybean)]

opt_md = Optimal_Clusters_Medoids(dat, 10, 'jaccard_coefficient', plot_clusters = FALSE)


#----------------------------
# non-contiguous search space
#----------------------------

search_space = c(2,5)

opt_md = Optimal_Clusters_Medoids(dat, search_space, 'jaccard_coefficient', plot_clusters = FALSE)

} # }
```
