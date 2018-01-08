
## ClusterR 1.1.0

* I added the *DARMA_64BIT_WORD* flag in the Makevars file to allow the package processing big datasets


## ClusterR 1.0.9

* I modified the *Optimal_Clusters_KMeans* function to return a vector with the *distortion_fK* values if criterion is *distortion_fK* (instead of the *WCSSE* values).
* I added the 'Moore-Penrose pseudo-inverse' for the case of the 'mahalanobis' distance calculation.


## ClusterR 1.0.8

* I modified the *OpenMP* clauses of the .cpp files to address the ASAN errors.
* I removed the *threads* parameter from the *KMeans_rcpp* function, to address the ASAN errors ( negligible performance difference between threaded and non-threaded version especially if the *num_init* parameter is less than 10 ). The *threads* parameter was removed also from the *Optimal_Clusters_KMeans* function as it utilizes the *KMeans_rcpp* function to find the optimal clusters for the various methods.


## ClusterR 1.0.7

I modified the *kmeans_miniBatchKmeans_GMM_Medoids.cpp* file in the following lines in order to fix the clang-ASAN errors (without loss in performance):

* lines 1156-1160 : I commented the second OpenMp parallel-loop and I replaced the *k* variable with the *i* variable in the second for-loop [in the *dissim_mat()* function]
* lines 1739-1741 : I commented the second OpenMp parallel-loop [in the *silhouette_matrix()* function]
* I replaced (all) the *silhouette_matrix* (arma::mat) variable names with *Silhouette_matrix*, because the name overlapped with the name of the Rcpp function [in the *silhouette_matrix* function]
* I replaced all *sorted_medoids.n_elem* with the variable *unsigned int sorted_medoids_elem* [in the *silhouette_matrix* function]

I modified the following *functions* in the *clustering_functions.R* file:

* *KMeans_rcpp()* : I added an *experimental* note in the details for the *optimal_init* and *quantile_init* initializers.
* *Optimal_Clusters_KMeans()* : I added an *experimental* note in the details for the *optimal_init* and *quantile_init* initializers.
* *MiniBatchKmeans()* : I added an *experimental* note in the details for the *optimal_init* and *quantile_init* initializers.


## ClusterR 1.0.6

The *normalized variation of information* was added in the *external_validation* function (https://github.com/mlampros/ClusterR/pull/1)


## ClusterR 1.0.5

I fixed the valgrind memory errors


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




