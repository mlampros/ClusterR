# Prediction function for a Gaussian Mixture Model object

Prediction function for a Gaussian Mixture Model object

## Usage

``` r
predict_GMM(data, CENTROIDS, COVARIANCE, WEIGHTS)

# S3 method for class 'GMMCluster'
predict(object, newdata, ...)
```

## Arguments

- data:

  matrix or data frame

- CENTROIDS:

  matrix or data frame containing the centroids (means), stored as row
  vectors

- COVARIANCE:

  matrix or data frame (for diagonal covariance) or 3D array (for full
  covariance matrices)

- WEIGHTS:

  vector containing the weights

- object, newdata, ...:

  arguments for the \`predict\` generic

## Value

a list consisting of the log-likelihoods, cluster probabilities and
cluster labels.

## Details

This function takes the centroids, covariance matrix and weights from a
trained model and returns the log-likelihoods, cluster probabilities and
cluster labels for new data. The function handles both diagonal
covariance matrices (2D matrix) and full covariance matrices (3D
array/cube).

## Author

Lampros Mouselimis

## Examples

``` r
data(dietary_survey_IBS)

dat = as.matrix(dietary_survey_IBS[, -ncol(dietary_survey_IBS)])

dat = center_scale(dat)

gmm = GMM(dat, 2, "maha_dist", "random_subset", 10, 10)

# pr = predict_GMM(dat, gmm$centroids, gmm$covariance_matrices, gmm$weights)
```
