[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ClusterR)](http://cran.r-project.org/package=ClusterR)
[![Travis-CI Build Status](https://travis-ci.org/mlampros/ClusterR.svg?branch=master)](https://travis-ci.org/mlampros/ClusterR)
[![codecov.io](https://codecov.io/github/mlampros/ClusterR/coverage.svg?branch=master)](https://codecov.io/github/mlampros/ClusterR?branch=master)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ClusterR?color=blue)](http://www.r-pkg.org/pkg/ClusterR)
<a href="https://www.buymeacoffee.com/VY0x8snyh" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" height="21px" ></a>
[![](https://img.shields.io/docker/automated/mlampros/clusterr.svg)](https://hub.docker.com/r/mlampros/clusterr)


## ClusterR
<br>

The ClusterR package consists of Gaussian mixture models, k-means, mini-batch-kmeans, k-medoids and affinity propagation clustering algorithms with the option to plot, validate, predict (new data) and find the optimal number of clusters. The package takes advantage of 'RcppArmadillo' to speed up the computationally intensive parts of the functions. More details on the functionality of ClusterR can be found in the [blog-post](http://mlampros.github.io/2016/09/12/clusterR_package/), Vignette and in the package Documentation.
<br><br>

**UPDATE 16-08-2018**


As of version 1.1.4 the *ClusterR* package allows R package maintainers to perform **linking between packages at a C++ code (Rcpp) level**. This means that the Rcpp functions of the *ClusterR* package can be called in the C++ files of another package. In the next lines I'll give detailed explanations on how this can be done:

<br>

Assumming that an R package ('PackageA') calls one of the *ClusterR* Rcpp functions. Then the maintainer of 'PackageA' has to :

<br>

* **1st.** install the *ClusterR* package to take advantage of the new functionality either from CRAN using,

<br>


```R

install.packages("ClusterR")
 

```

<br>

or download the latest version from Github using the *remotes* package,

<br>

```R

remotes::install_github('mlampros/ClusterR', upgrade = 'always', dependencies = TRUE, repos = 'https://cloud.r-project.org/')
 

```

<br>

* **2nd.** update the **DESCRIPTION** file of 'PackageA' and especially the *LinkingTo* field by adding the *ClusterR* package (besides any other packages),

<br>

```R

LinkingTo: ClusterR

```

<br>

* **3rd.** open a **new C++ file** (for instance in Rstudio) and at the top of the file add the following 'headers', 'depends' and 'plugins',

<br>

```R

# include <RcppArmadillo.h>
# include <ClusterRHeader.h>
# include <affinity_propagation.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(ClusterR)]]
// [[Rcpp::plugins(cpp11)]]


```
<br>

The available functions can be found in the following files: **inst/include/ClusterRHeader.h** and **inst/include/affinity_propagation.h**

<br>

A *complete minimal example* would be :

<br>

```R
# include <RcppArmadillo.h>
# include <ClusterRHeader.h>
# include <affinity_propagation.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(ClusterR)]]
// [[Rcpp::plugins(cpp11)]]


using namespace clustR;


// [[Rcpp::export]]
Rcpp::List mini_batch_kmeans(arma::mat& data, int clusters, int batch_size, int max_iters, int num_init = 1, 

                            double init_fraction = 1.0, std::string initializer = "kmeans++",
                            
                            int early_stop_iter = 10, bool verbose = false, 
                            
                            Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, 
                            
                            double tol = 1e-4, double tol_optimal_init = 0.5, int seed = 1) {

  ClustHeader clust_header;

  return clust_header.mini_batch_kmeans(data, clusters, batch_size, max_iters, num_init, init_fraction, 
  
                                        initializer, early_stop_iter, verbose, CENTROIDS, tol, 
                                        
                                        tol_optimal_init, seed);
}


```

<br>

Then, by opening an R file a user can call the *mini_batch_kmeans* function using,

<br>

```R

Rcpp::sourceCpp('example.cpp')              # assuming that the previous Rcpp code is included in 'example.cpp' 
             
set.seed(1)
dat = matrix(runif(100000), nrow = 1000, ncol = 100)

mbkm = mini_batch_kmeans(dat, clusters = 3, batch_size = 50, max_iters = 100, num_init = 2, 

                         init_fraction = 1.0, initializer = "kmeans++", early_stop_iter = 10, 
                         
                         verbose = T, CENTROIDS = NULL, tol = 1e-4, tol_optimal_init = 0.5, seed = 1)
                         
str(mbkm)


```

<br>


Use the following link to report bugs/issues,
<br><br>

[https://github.com/mlampros/ClusterR/issues](https://github.com/mlampros/ClusterR/issues)

<br>


**UPDATE 15-08-2019**

<br>

I've added a **Dockerfile** to the respository so that users can build and run a docker image either on their *personal computer*, a *cloud instance* or a *cluster of computers*. It requires that [docker](https://www.docker.com/) is installed. I've included a [web-link](https://docs.docker.com/install/linux/docker-ce/ubuntu/) in the references on how to install docker on ubuntu but there are similar instructions for windows and Macintosh (either in the official docker website or on the web generally). The whole process was tested on Ubuntu 18.04

<br>

To build the docker image on Ubuntu use,

<br>

```R

sudo docker build -f Dockerfile --rm --force-rm -t rstudio/clusterr .


```

<br>

the name of the image (both in *build* and *run* steps) can be changed depending on your preference (this is just an example usage). Then to run the image use,

<br>

```R

docker run --rm -d --name dckr_rstd_clust -e USER=rstudio -e PASSWORD=give_here_your_password -p 8787:8787 rstudio/clusterr


```

<br>

the **USER** defaults to *rstudio* but you have to give your **PASSWORD** of preference (see also the [second reference web-link](https://www.rocker-project.org/) for more information).

<br>

Open your web-browser and depending where the docker image was *build / run* give, 

<br>

**1st. Option** on your personal computer,

<br>

```R
http://0.0.0.0:8787 

```

<br>

**2nd. Option** on a cloud instance, 

<br>

```R
http://Public DNS:8787

```

<br>

to access the Rstudio console in order to give your username and password.

<br>

*References:*

* [https://docs.docker.com/install/linux/docker-ce/ubuntu/](https://docs.docker.com/install/linux/docker-ce/ubuntu/)
* [https://www.rocker-project.org/](https://www.rocker-project.org/)

<br>
