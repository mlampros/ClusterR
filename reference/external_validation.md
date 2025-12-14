# external clustering validation

external clustering validation

## Usage

``` r
external_validation(
  true_labels,
  clusters,
  method = "adjusted_rand_index",
  summary_stats = FALSE
)
```

## Arguments

- true_labels:

  a numeric vector of length equal to the length of the clusters vector

- clusters:

  a numeric vector ( the result of a clustering method ) of length equal
  to the length of the true_labels

- method:

  one of *rand_index*, *adjusted_rand_index*, *jaccard_index*,
  *fowlkes_Mallows_index*, *mirkin_metric*, *purity*, *entropy*, *nmi*
  (normalized mutual information), *var_info* (variation of
  information), and *nvi* (normalized variation of information)

- summary_stats:

  besides the available methods the summary_stats parameter prints also
  the specificity, sensitivity, precision, recall and F-measure of the
  clusters

## Value

if summary_stats is FALSE the function returns a float number, otherwise
it returns also a summary statistics table

## Details

This function uses external validation methods to evaluate the
clustering results

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

X = center_scale(dat)

km = KMeans_rcpp(X, clusters = 2, num_init = 5, max_iters = 100, initializer = 'kmeans++')

res = external_validation(dietary_survey_IBS$class, km$clusters, method = "adjusted_rand_index")
```
