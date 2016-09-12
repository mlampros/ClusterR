
## ClusterR
<br>

The ClusterR package consists of Gaussian mixture models, k-means, mini-batch-kmeans and k-medoids clustering algorithms with the option to plot, validate, predict (new data) and find the optimal number of clusters. The package takes advantage of 'RcppArmadillo' to speed up the computationally intensive parts of the functions. More details on the functionality of ClusterR can be found in the [blog-post](http://mlampros.github.io/2016/09/12/clusterR_package/) and in the package Vignette.
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
