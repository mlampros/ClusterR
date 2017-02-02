
## ClusterR 1.0.4

I removed the warnings, which occured during compilation.
I corrected the UBSAN memory errors which occured due to a mistake in the *check_medoids()* function of the *utils_rcpp.cpp* file.
I also modified the *quantile_init_rcpp()* function of the *utils_rcpp.cpp* file to print a warning if duplicates are present in the initial centroid matrix.


## ClusterR 1.0.3

* I updated the dissimilarity functions to accept data with missing values.
* I added an error exception in the predict_GMM() function in case that the determinant is equal to zero. The latter is possible if the data includes highly correlated variables or variables with low variance.
* I replaced all unsigned int's in the rcpp files with int data types


## ClusterR 1.0.2

I modified the RcppArmadillo functions so that ClusterR passes the Windows and OSX OS package check results


## ClusterR 1.0.1

I modified the RcppArmadillo functions so that ClusterR passes the Windows and OSX OS package check results


## ClusterR 1.0.0




