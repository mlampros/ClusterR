#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]


# include "ClusterRHeader.h"
# include "affinity_propagation.h"

using namespace clustR;


//============================================ utility functions ==========================================================================


// it returns TRUE if the matrix does not include NaN's or +/- Inf
// it returns FALSE if at least one value is NaN or +/- Inf

// [[Rcpp::export]]
bool check_NaN_Inf(arma::mat x) {

  ClustHeader CRH;

  return CRH.check_NaN_Inf(x);
}


// this function takes the resulted centers from the 'KMEANS_rcpp' OR 'KMEANS_arma' function and
// returns the clusters. The returned clusters should match the clusters in case of the KMEANS_rcpp function.
//

// [[Rcpp::export]]
arma::rowvec validate_centroids(arma::mat& data, arma::mat init_centroids, int threads) {

  ClustHeader CRH;

  return CRH.validate_centroids(data, init_centroids, threads);
}


// this function scales and/or centers the data
//

// [[Rcpp::export]]
arma::mat SCALE(arma::mat data, bool mean_center = true, bool sd_scale = true) {

  ClustHeader CRH;

  return CRH.SCALE(data, mean_center, sd_scale);
}



//============================================ k-means ====================================================================================


// simple k-means
// It uses 4 different initializers and the CENTROIDS parameter
//
// fuzzy returns prediction probabilities in a fast way based on the distance between observations and centroids and it isn't similar to fuzzy-c-means
//
// https://github.com/scikit-learn/scikit-learn/blob/51a765a/sklearn/cluster/k_means_.py#L660
//

// [[Rcpp::export]]
Rcpp::List KMEANS_rcpp(arma::mat& data, int clusters, int num_init = 1, int max_iters = 200, std::string initializer = "kmeans++", bool fuzzy = false, bool verbose = false,

                       Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, double tol = 1e-4, double eps = 1.0e-6, double tol_optimal_init = 0.5, int seed = 1) {

  ClustHeader CRH;

  return CRH.KMEANS_rcpp(data, clusters, num_init, max_iters, initializer, fuzzy, verbose, R_NilValue, tol, eps, tol_optimal_init, seed);
}



// the KMEANS_arma returns only the centroids and if openmp exists then it uses all available threads
// the number of columns of the data should be much larger than the number of clusters  (however it works for number_columns == clusters)
// seed_mode is one of : "keep_existing" (I've to give the CENTROIDS), "static_subset", "random_subset", "static_spread", "random_spread"
// in comparison to my 'KMEANS_rcpp' the 'kmeans_arma' does a single initialization and it suggests maximum 10 iterations
//

// [[Rcpp::export]]
arma::mat KMEANS_arma(arma::mat& data, int clusters, int n_iter, bool verbose, std::string seed_mode = "random_subset",

                      Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, int seed = 1) {

  ClustHeader CRH;

  return CRH.KMEANS_arma(data, clusters, n_iter, verbose, seed_mode, CENTROIDS, seed);
}



// optimal number of clusters using distortion in kmeans
// https://www.ee.columbia.edu/~dpwe/papers/PhamDN05-kmeans.pdf
//

// [[Rcpp::export]]
Rcpp::List opt_clust_fK(arma::vec sum_distortion, int data_num_cols, double threshold = 0.85) {

  ClustHeader CRH;

  return CRH.opt_clust_fK(sum_distortion, data_num_cols, threshold);
}



//------------------------------------------------------------------------------------
// dissimilarity and silhouette estimation for k-means based on the euclidean-distance
//------------------------------------------------------------------------------------

// the dissimilarities and silhouette widths of k-means differ from those of the k-medoids, because the latter are based on
// the dissimilarity matrix [ observations of the data are used as medoids ] whereas the centroids of the k-means algorithm
// depend on random initialization and the output-centroids are averages of multiple runs


// dissimilarity or silhouette evaluation of clustering [ dissimilarity returns faster as it only requires intra-distances to be calculated ]
// http://blog.data-miners.com/2011/03/cluster-silhouettes.html
//

// [[Rcpp::export]]
Rcpp::List evaluation_rcpp(arma::mat& data, arma::vec CLUSTER, bool silhouette = false) {

  ClustHeader CRH;

  return CRH.evaluation_rcpp(data, CLUSTER, silhouette);
}



//============================================ mini-batch-k-means ===================================================================================


// mini-batch-kmeans
//
// in each iteration only a batch will be used to update the centroids
//
// PARAMETERS:
// max_iters : maximum number of iterations in each initialization
// num_init  : number of initializations
// init_fraction : percentage of data to use in the initialization phase, i.e. use only a fraction of the data , 'init_fraction', to build the initialization centroids [ applies only if initializer = 'kmeans++']
// early_stop_iter : continue that many iterations after calculation of the best WCSS
//
// http://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf
// https://github.com/siddharth-agrawal/Mini-Batch-K-Means
// https://github.com/scikit-learn/scikit-learn/blob/51a765a/sklearn/cluster/k_means_.py#L1113
//

// [[Rcpp::export]]
Rcpp::List mini_batch_kmeans(arma::mat& data, int clusters, int batch_size, int max_iters, int num_init = 1, double init_fraction = 1.0, std::string initializer = "kmeans++",

                             int early_stop_iter = 10, bool verbose = false, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, double tol = 1e-4, double tol_optimal_init = 0.5, int seed = 1) {

  ClustHeader CRH;

  return CRH.mini_batch_kmeans(data, clusters, batch_size, max_iters, num_init, init_fraction, initializer, early_stop_iter, verbose, CENTROIDS, tol, tol_optimal_init, seed);
}




// predict function for mini-batch-kmeans, which takes a matrix of centroids
//

// [[Rcpp::export]]
Rcpp::List Predict_mini_batch_kmeans(arma::mat& data, arma::mat& CENTROIDS, bool fuzzy = false, double eps = 1.0e-6) {

  ClustHeader CRH;

  return CRH.Predict_mini_batch_kmeans(data, CENTROIDS, fuzzy, eps);
}




//============================================ Gaussian Mixture Models (GMM) =========================================================================


// the GMM_arma returns the centroids, covariance matrix and the weights. If openmp exists then it uses all available threads
// distance is one of : "eucl_dist", "maha_dist"
// seed_mode is one of : "static_subset", "random_subset", "static_spread", "random_spread"  [ I excluded the "keep_existing" seed_mode ] -- user-defined parameter setting is not enabled
//

// [[Rcpp::export]]
Rcpp::List GMM_arma(arma::mat& data, int gaussian_comps, std::string dist_mode, std::string seed_mode, int km_iter, int em_iter,

                    bool verbose, double var_floor = 1e-10, int seed = 1) {

  ClustHeader CRH;

  return CRH.GMM_arma(data, gaussian_comps, dist_mode, seed_mode, km_iter, em_iter, verbose, var_floor, seed);
}




// multivariate normal (gaussian) distribution probability density function
//
// use this function to predict new-data. It takes the data, centroids, covariance matrix and weights (from a trained model) and it returns log-likelihood,
// probabilities and labels for unknown data.
//
// formula for multivariate normal distribution probability density function : SEE book "Fundamentals of Stream Processing: Application Design, Systems, and Analytics", page 406
//
// formula for log-likelihood function : SEE "http://jonathantemplin.com/files/multivariate/mv12uga/mv12uga_section05.pdf", page 15 and 26
//
// The predictions are based on centroids, covariance matrix and weights using the formulas and not on the log-likelihoods.
//

// [[Rcpp::export]]
Rcpp::List predict_MGausDPDF(arma::mat data, arma::mat CENTROIDS, arma::mat COVARIANCE, arma::vec WEIGHTS, double eps = 1.0e-8) {

  ClustHeader CRH;

  return CRH.predict_MGausDPDF(data, CENTROIDS, COVARIANCE, WEIGHTS, eps);
}




// function to calculate bic-aic
//

// [[Rcpp::export]]
arma::rowvec GMM_arma_AIC_BIC(arma::mat& data, arma::rowvec max_clusters, std::string dist_mode, std::string seed_mode,

                              int km_iter, int em_iter, bool verbose, double var_floor = 1e-10, std::string criterion = "AIC", int seed = 1) {

  ClustHeader CRH;

  return CRH.GMM_arma_AIC_BIC(data, max_clusters, dist_mode, seed_mode, km_iter, em_iter, verbose, var_floor, criterion, seed);
}




//============================================ cluster Medoids ===================================================================================


// dissimilarity matrix using various distance metric methods [ the function can handle missing values by using pair-wise deletion ]
//

// [[Rcpp::export]]
arma::mat dissim_mat(arma::mat& data, std::string& method, double minkowski_p = 1.0, bool upper = true, bool diagonal = true, int threads = 1, double eps = 1.0e-6) {

  ClustHeader CRH;

  return CRH.dissim_mat(data, method, minkowski_p, upper, diagonal, threads, eps);
}


// the 'ClusterMedoids' function should accept either a matrix or a dissimilarity matrix
// the dissimilarity matrix in comparison to a single matrix will have nrows == ncols AND diagonal == 0.0 [ it should be also a matrix ]
//

// [[Rcpp::export]]
Rcpp::List ClusterMedoids(arma::mat& data, int clusters, std::string method, double minkowski_p = 1.0, int threads = 1, bool verbose = false, bool swap_phase = false,

                           bool fuzzy = false, int seed = 1) {

  ClustHeader CRH;

  return CRH.ClusterMedoids(data, clusters, method, minkowski_p, threads, verbose, swap_phase, fuzzy, seed);
}



// calculate global dissimilarities for claraMedoids
// the function can handle missing values by using pair-wise deletion
//

// [[Rcpp::export]]
arma::mat dissim_MEDOIDS(arma::mat& data, std::string& method, arma::mat MEDOIDS, double minkowski_p = 1.0, int threads = 1, double eps = 1.0e-6) {

  ClustHeader CRH;

  return CRH.dissim_MEDOIDS(data, method, MEDOIDS, minkowski_p, threads, eps);
}



// http://www.sthda.com/english/wiki/partitioning-cluster-analysis-quick-start-guide-unsupervised-machine-learning#pam-partitioning-around-medoids

// 1. Split randomly the data sets in multiple subsets with fixed size
// 2. Compute PAM algorithm on each subset and choose the corresponding k representative objects (medoids). Assign each observation of the entire dataset to the nearest medoid.
// 3. Calculate the mean (or the sum) of the dissimilarities (for instance 'euclidean' distance) of the observations to their closest medoid. This is used as a measure of the goodness of the clustering.
// 4. Retain (keep) the sub-dataset for which the mean (or sum) is minimal. A further analysis is carried out on the final partition.

//--------------------------------------------------------------------------------------------------------------------------------------------------
// Important: if a medoid is significant --> add this medoid to the next sample by subsetting the data and appending the corresponding observations
//--------------------------------------------------------------------------------------------------------------------------------------------------

// Instead of finding representative objects for the entire data set, 'ClaraMedoids' draws a sample of the data set, applies the 'ClusterMedoids' function
// on the sample and finds the medoids of the sample. The point is that if the sample is drawn in a sufficiently random way the medoids of the sample would
// approximate the medoids of the entire data set. To come up with better approximation, 'ClaraMedoids' draws mulitple samples and GIVES THE BEST CLUSTERING
// as the OUTPUT. Here, for ACCURACY the QUALITY of a clustering is measured based on the AVERAGE DISSIMILARITY of all objects in the entire data set and NOT
// ONLY of those objects IN THE SAMPLES. Experiments, indicate that 5 samples of size 40 + 2 * k (i.e. 40 + 2 * number_of_clusters) give satisfactory results.
//


// [[Rcpp::export]]
Rcpp::List ClaraMedoids(arma::mat& data, int clusters, std::string method, int samples, double sample_size, double minkowski_p = 1.0,

                        int threads = 1, bool verbose = false, bool swap_phase = false, bool fuzzy = false, int seed = 1) {

  ClustHeader CRH;

  return CRH.ClaraMedoids(data, clusters, method, samples, sample_size, minkowski_p, threads, verbose, swap_phase, fuzzy, seed);
}



// prediction function for the k-medoids
//

// [[Rcpp::export]]
Rcpp::List predict_medoids(arma::mat& data, std::string method, arma::mat MEDOIDS, double minkowski_p = 1.0, int threads = 1, bool fuzzy = false, double eps = 1.0e-6) {

  ClustHeader CRH;

  return CRH.predict_medoids(data, method, MEDOIDS, minkowski_p, threads, fuzzy, eps);
}



// function which takes an Rcpp List object of matrices from Cluster_Medoids OR Clara_Medoids and returns
// those matrices split on the specific clusters


// [[Rcpp::export]]
Rcpp::List split_rcpp_lst(Rcpp::List lst) {

  ClustHeader CRH;

  return CRH.split_rcpp_lst(lst);
}




// optimal number of clusters for the k-medoids
//

// [[Rcpp::export]]
Rcpp::List OptClust(arma::mat& data, arma::rowvec iter_clust, std::string method, bool clara = false, int samples = 5, double sample_size = 0.001, double minkowski_p = 1.0,

                    std::string criterion = "dissimilarity", int threads = 1, bool swap_phase = false, bool verbose = false, int seed = 1) {

  ClustHeader CRH;

  return CRH.OptClust(data, iter_clust, method, clara, samples, sample_size, minkowski_p, criterion, threads, swap_phase, verbose, seed);
}


//============================================ affinity propagation ===================================================================================

// affinity propagation algorithm
//

// [[Rcpp::export]]
Rcpp::List affinity_propagation(arma::mat &s, std::vector<double> p, int maxits = 1000, int convits = 100, double dampfact = 0.9,
                                bool details = false, double nonoise = 0.0, double eps = 2.2204e-16, bool time = false) {
  
  Affinity_Propagation AFN;
  return AFN.affinity_propagation(s, p, maxits, convits, dampfact, details, nonoise, eps, time);
}


// preferenceRange function
//

// [[Rcpp::export]]
std::vector<double> preferenceRange(arma::mat &s, std::string method = "bound", int threads = 1) {

  Affinity_Propagation AFN;
  return AFN.preferenceRange(s, method, threads);
}

//=====================================================================================================================================================
