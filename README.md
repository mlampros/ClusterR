[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ClusterR)](http://cran.r-project.org/package=ClusterR)
[![Travis-CI Build Status](https://travis-ci.org/mlampros/ClusterR.svg?branch=master)](https://travis-ci.org/mlampros/ClusterR)
[![codecov.io](https://codecov.io/github/mlampros/ClusterR/coverage.svg?branch=master)](https://codecov.io/github/mlampros/ClusterR?branch=master)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ClusterR?color=blue)](http://www.r-pkg.org/pkg/ClusterR)


## ClusterR
<br>

The ClusterR package consists of Gaussian mixture models, k-means, mini-batch-kmeans and k-medoids clustering algorithms with the option to plot, validate, predict (new data) and find the optimal number of clusters. The package takes advantage of 'RcppArmadillo' to speed up the computationally intensive parts of the functions. More details on the functionality of ClusterR can be found in the [blog-post](http://mlampros.github.io/2016/09/12/clusterR_package/) and in the package Vignette. ClusterR can be installed, currently, in the following OS's: Linux, Mac and Windows.
<br><br>

To install the package from CRAN use, 

```R

install.packages("ClusterR")


```
<br>

and to download the latest version from Github use the *install_github* function of the devtools package,
<br><br>

```R

devtools::install_github('mlampros/ClusterR')


```
<br>

Use the following link to report bugs/issues,
<br><br>

[https://github.com/mlampros/ClusterR/issues](https://github.com/mlampros/ClusterR/issues)
