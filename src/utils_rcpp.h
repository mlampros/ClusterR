#ifndef __utils_rcpp__
#define __utils_rcpp__

arma::rowvec sample_vec(int num_elem, int start, int end, bool replace);
double squared_norm(arma::mat x);
int MinMat(arma::vec x);
arma::vec WCSS(arma::rowvec vec, arma::mat centroids);
arma::rowvec validate_centroids(arma::mat& data, arma::mat init_medoids);
double kmeans_pp_dist(arma::rowvec vec, arma::rowvec centroid);
arma::mat kmeans_pp_init(arma::mat& data, int clusters, bool medoids);
arma::rowvec norm_fuzzy(arma::rowvec vec, double eps);
Rcpp::NumericVector quantile_value(arma::rowvec x, int clusters);
arma::mat quantile_init_rcpp(arma::mat data, int sample_rows, int clusters, bool medoids);
arma::mat check_medoids(arma::mat data, int clust, double tol, bool medoids);
arma::mat SCALE(arma::mat data, bool mean_center, bool sd_scale);
double calc_silhouette(double inter, double outer);
bool check_NaN_Inf(arma::mat x);
Rcpp::List cluster_indices(arma::vec CLUSTER);
void set_seed(int seed);

#endif
