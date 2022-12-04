
[![tic](https://github.com/mlampros/ClusterR/workflows/tic/badge.svg?branch=master)](https://github.com/mlampros/ClusterR/actions)
[![codecov.io](https://codecov.io/github/mlampros/ClusterR/coverage.svg?branch=master)](https://codecov.io/github/mlampros/ClusterR?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ClusterR)](http://cran.r-project.org/package=ClusterR)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ClusterR?color=blue)](http://www.r-pkg.org/pkg/ClusterR)
<a href="https://www.buymeacoffee.com/VY0x8snyh" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" height="21px" ></a>
[![](https://img.shields.io/docker/automated/mlampros/clusterr.svg)](https://hub.docker.com/r/mlampros/clusterr)
[![Dependencies](https://tinyverse.netlify.com/badge/ClusterR)](https://cran.r-project.org/package=ClusterR)


## ClusterR
<br>

The ClusterR package consists of Gaussian mixture models, k-means, mini-batch-kmeans, k-medoids and affinity propagation clustering algorithms with the option to plot, validate, predict (new data) and find the optimal number of clusters. The package takes advantage of 'RcppArmadillo' to speed up the computationally intensive parts of the functions. More details on the functionality of ClusterR can be found in the [blog-post](http://mlampros.github.io/2016/09/12/clusterR_package/), Vignette and in the package Documentation ( *scroll down for information on how to use the* **docker image** )
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


**UPDATE 28-11-2019**

<br>

**Docker images** of the *ClusterR* package are available to download from my [dockerhub](https://hub.docker.com/r/mlampros/clusterr) account. The images come with *Rstudio* and the *R-development* version (latest) installed. The whole process was tested on Ubuntu 18.04. To **pull** & **run** the image do the following,

<br>

```R

docker pull mlampros/clusterr:rstudiodev

docker run -d --name rstudio_dev -e USER=rstudio -e PASSWORD=give_here_your_password --rm -p 8787:8787 mlampros/clusterr:rstudiodev

```

<br>

The user can also **bind** a home directory / folder to the image to use its files by specifying the **-v** command,

<br>

```R

docker run -d --name rstudio_dev -e USER=rstudio -e PASSWORD=give_here_your_password --rm -p 8787:8787 -v /home/YOUR_DIR:/home/rstudio/YOUR_DIR mlampros/clusterr:rstudiodev


```

<br>

In the latter case you might have first give permission privileges for write access to **YOUR_DIR** directory (not necessarily) using,

<br>

```R

chmod -R 777 /home/YOUR_DIR


```

<br>

The **USER** defaults to *rstudio* but you have to give your **PASSWORD** of preference (see [https://rocker-project.org/](https://rocker-project.org/) for more information).

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

### **Citation:**

If you use the code of this repository in your paper or research please cite both **ClusterR** and the **original articles / software** `https://CRAN.R-project.org/package=ClusterR`:

<br>

```R
@Manual{,
  title = {{ClusterR}: Gaussian Mixture Models, K-Means, Mini-Batch-Kmeans, K-Medoids and Affinity Propagation Clustering},
  author = {Lampros Mouselimis},
  year = {2022},
  note = {R package version 1.2.8},
  url = {https://CRAN.R-project.org/package=ClusterR},
}
```

<br>

