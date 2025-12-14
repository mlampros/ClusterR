# Distance matrix calculation

Distance matrix calculation

## Usage

``` r
distance_matrix(
  data,
  method = "euclidean",
  upper = FALSE,
  diagonal = FALSE,
  minkowski_p = 1,
  threads = 1
)
```

## Arguments

- data:

  matrix or data frame

- method:

  a string specifying the distance method. One of, *euclidean*,
  *manhattan*, *chebyshev*, *canberra*, *braycurtis*,
  *pearson_correlation*, *simple_matching_coefficient*, *minkowski*,
  *hamming*, *jaccard_coefficient*, *Rao_coefficient*, *mahalanobis*,
  *cosine*

- upper:

  either TRUE or FALSE specifying if the upper triangle of the distance
  matrix should be returned. If FALSE then the upper triangle will be
  filled with NA's

- diagonal:

  either TRUE or FALSE specifying if the diagonal of the distance matrix
  should be returned. If FALSE then the diagonal will be filled with
  NA's

- minkowski_p:

  a numeric value specifying the minkowski parameter in case that method
  = "minkowski"

- threads:

  the number of cores to run in parallel (if OpenMP is available)

## Value

a matrix

## Examples

``` r
data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = distance_matrix(dat, method = 'euclidean', upper = TRUE, diagonal = TRUE)
```
