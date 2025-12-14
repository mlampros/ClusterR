# Optimal number of Clusters for Kmeans or Mini-Batch-Kmeans

Optimal number of Clusters for Kmeans or Mini-Batch-Kmeans

## Usage

``` r
Optimal_Clusters_KMeans(
  data,
  max_clusters,
  criterion = "variance_explained",
  fK_threshold = 0.85,
  num_init = 1,
  max_iters = 200,
  initializer = "kmeans++",
  tol = 1e-04,
  plot_clusters = TRUE,
  verbose = FALSE,
  tol_optimal_init = 0.3,
  seed = 1,
  mini_batch_params = NULL
)
```

## Arguments

- data:

  matrix or data frame

- max_clusters:

  either a numeric value, a contiguous or non-continguous numeric vector
  specifying the cluster search space

- criterion:

  one of *variance_explained*, *WCSSE*, *dissimilarity*, *silhouette*,
  *distortion_fK*, *AIC*, *BIC* and *Adjusted_Rsquared*. See details for
  more information.

- fK_threshold:

  a float number used in the 'distortion_fK' criterion

- num_init:

  number of times the algorithm will be run with different centroid
  seeds

- max_iters:

  the maximum number of clustering iterations

- initializer:

  the method of initialization. One of, *optimal_init*, *quantile_init*,
  *kmeans++* and *random*. See details for more information

- tol:

  a float number. If, in case of an iteration (iteration \> 1 and
  iteration \< max_iters) 'tol' is greater than the squared norm of the
  centroids, then kmeans has converged

- plot_clusters:

  either TRUE or FALSE, indicating whether the results of the
  *Optimal_Clusters_KMeans* function should be plotted

- verbose:

  either TRUE or FALSE, indicating whether progress is printed during
  clustering

- tol_optimal_init:

  tolerance value for the 'optimal_init' initializer. The higher this
  value is, the far appart from each other the centroids are.

- seed:

  integer value for random number generator (RNG)

- mini_batch_params:

  either NULL or a list of the following parameters : *batch_size*,
  *init_fraction*, *early_stop_iter*. If not NULL then the optimal
  number of clusters will be found based on the Mini-Batch-Kmeans. See
  the details and examples sections for more information.

## Value

a vector with the results for the specified criterion. If plot_clusters
is TRUE then it plots also the results.

## Details

—————criteria————————–

**variance_explained** : the sum of the
within-cluster-sum-of-squares-of-all-clusters divided by the total sum
of squares

**WCSSE** : the sum of the within-cluster-sum-of-squares-of-all-clusters

**dissimilarity** : the average intra-cluster-dissimilarity of all
clusters (the distance metric defaults to euclidean)

**silhouette** : the average silhouette width where first the average
per cluster silhouette is computed and then the global average (the
distance metric defaults to euclidean). To compute the silhouette width
for each cluster separately see the 'silhouette_of_clusters()' function

**distortion_fK** : this criterion is based on the following paper,
'Selection of K in K-means clustering'
(https://www.ee.columbia.edu/~dpwe/papers/PhamDN05-kmeans.pdf)

**AIC** : the Akaike information criterion

**BIC** : the Bayesian information criterion

**Adjusted_Rsquared** : the adjusted R^2 statistic

—————initializers———————-

**optimal_init** : this initializer adds rows of the data incrementally,
while checking that they do not already exist in the centroid-matrix \[
experimental \]

**quantile_init** : initialization of centroids by using the cummulative
distance between observations and by removing potential duplicates \[
experimental \]

**kmeans++** : kmeans++ initialization. Reference :
http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf AND
http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work

**random** : random selection of data rows as initial centroids

If the *mini_batch_params* parameter is not NULL then the optimal number
of clusters will be found based on the Mini-batch-Kmeans algorithm,
otherwise based on the Kmeans. The higher the *init_fraction* parameter
is the more close the results between Mini-Batch-Kmeans and Kmeans will
be.

In case that the *max_clusters* parameter is a contiguous or
non-contiguous vector then plotting is disabled. Therefore, plotting is
enabled only if the *max_clusters* parameter is of length 1. Moreover,
the *distortion_fK* criterion can't be computed if the *max_clusters*
parameter is a contiguous or non-continguous vector ( the
*distortion_fK* criterion requires consecutive clusters ). The same
applies also to the *Adjusted_Rsquared* criterion which returns
incorrect output.

## References

https://www.ee.columbia.edu/~dpwe/papers/PhamDN05-kmeans.pdf

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)


#-------
# kmeans
#-------

opt_km = Optimal_Clusters_KMeans(dat, max_clusters = 10, criterion = "distortion_fK",

                                 plot_clusters = FALSE)

#------------------
# mini-batch-kmeans
#------------------


params_mbkm = list(batch_size = 10, init_fraction = 0.3, early_stop_iter = 10)

opt_mbkm = Optimal_Clusters_KMeans(dat, max_clusters = 10, criterion = "distortion_fK",

                                   plot_clusters = FALSE, mini_batch_params = params_mbkm)


#----------------------------
# non-contiguous search space
#----------------------------

search_space = c(2,5)

opt_km = Optimal_Clusters_KMeans(dat, max_clusters = search_space,

                                 criterion = "variance_explained",

                                 plot_clusters = FALSE)
```
