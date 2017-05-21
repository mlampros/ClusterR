
# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils_rcpp.h"


// use base R's set.seed() in Rcpp for RNG
// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case
//

// [[Rcpp::export]]
void set_seed(int seed) {

  Rcpp::Environment base_env("package:base");

  Rcpp::Function set_seed_r = base_env["set.seed"];

  set_seed_r(seed);
}


// secondary function which takes a vector with the clusters and returns a list with the cluster indices
//

// [[Rcpp::export]]
Rcpp::List cluster_indices(arma::vec CLUSTER) {

  arma::vec unq_values = arma::unique(CLUSTER);

  arma::vec count_clust(unq_values.n_elem, arma::fill::zeros);

  for (unsigned int j = 0; j < unq_values.n_elem; j++) {

    int midl_val = 0;

    for (unsigned int i = 0; i < CLUSTER.n_elem; i++) {

      if (unq_values(j) == CLUSTER(i)) {

        midl_val += 1;
      }
    }

    count_clust(j) = midl_val;
  }

  Rcpp::List tmp_clust(unq_values.n_elem);                            // initializes the list to store the indices

  int count_init = 0;

  for (unsigned int k = 0; k < unq_values.n_elem; k++) {

    Rcpp::NumericVector tmp_VEC(count_clust(k));

    for (unsigned int j = 0; j < CLUSTER.n_elem; j++) {

      if (unq_values(k) == CLUSTER(j)) {

        tmp_VEC[count_init] = j;

        count_init += 1;
      }
    }

    count_init = 0;

    tmp_clust[k] = tmp_VEC;
  }

  return tmp_clust;
}


// it returns TRUE if the matrix does not include NaN's or +/- Inf
// it returns FALSE if at least one value is NaN or +/- Inf

// [[Rcpp::export]]
bool check_NaN_Inf(arma::mat x) {

  return x.is_finite();
}


// silhouette formula for the intra-, outer-cluster-dissimilarities
// https://en.wikipedia.org/wiki/Silhouette_(clustering)
//

// [[Rcpp::export]]
double calc_silhouette(double intra, double outer) {

  if (outer > intra) {

    return 1.0 - (intra / outer);}

  else if (outer < intra) {

    return (outer / intra) - 1.0;}

  else {                             // outer == intra

    return 0.0;
  }
}



// sample with or without replacement [ 'num_elem' is an integer not a vector -- create a similar function for a vector ]
//

// [[Rcpp::export]]
arma::rowvec sample_vec(int num_elem, int start, int end, bool replace) {

  if (replace) {

    return(arma::randi<arma::rowvec>( num_elem, arma::distr_param(start,end) ));}

  else {

    arma::rowvec tmp = arma::conv_to< arma::rowvec >::from(arma::shuffle(arma::regspace(start,end)));

    return(tmp.subvec(0,num_elem-1));
  }
}



// early-stopping criterium

// [[Rcpp::export]]
double squared_norm(arma::mat x) {

  return std::sqrt(arma::accu(arma::pow(x, 2)));
}



// this function gives the index of the minimum value in a vector
//

// [[Rcpp::export]]
int MinMat(arma::vec x) {

  double out = arma::datum::inf;

  int idx = 0;

  for (unsigned int i = 0; i < x.n_elem; i++) {

    if (x(i) < out) {

      out = x(i);

      idx = i;
    }
  }

  return(idx);
}




// within-cluster-sum-of-squares
//

// [[Rcpp::export]]
arma::vec WCSS(arma::rowvec vec, arma::mat centroids) {

  arma::vec tmp_c(centroids.n_rows);

  for (unsigned int i = 0; i < centroids.n_rows; i++) {

    tmp_c(i) = arma::as_scalar(arma::accu(arma::pow(vec - centroids.row(i), 2)));
  }

  return tmp_c;
}




// this function takes the resulted centers from the 'KMEANS_rcpp' OR 'KMEANS_arma' function and
// returns the clusters. The returned clusters should match the clusters in case of the KMEANS_rcpp function.
//

// [[Rcpp::export]]
arma::rowvec validate_centroids(arma::mat& data, arma::mat init_centroids) {

  arma::rowvec tmp_idx(data.n_rows);

  for (unsigned int k = 0; k < data.n_rows; k++) {

    arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data.row(k)), init_centroids);

    tmp_idx(k) = MinMat(tmp_vec);
  }

  return tmp_idx;
}



// kmeans++ distance
//

// [[Rcpp::export]]
double kmeans_pp_dist(arma::rowvec vec, arma::rowvec centroid) {

  return arma::as_scalar(arma::accu(arma::pow(vec - centroid, 2)));
}



// kmeans++ algorithm to build initialization-centroids [ it works for 1:data.n_row-2 ]
//
// http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work
// https://datasciencelab.wordpress.com/2014/01/15/improved-seeding-for-clustering-with-k-means/
// https://en.wikipedia.org/wiki/K-means%2B%2B
//

// [[Rcpp::export]]
arma::mat kmeans_pp_init(arma::mat& data, int clusters, bool medoids = false) {

  arma::rowvec centroid_data(data.n_rows);                                                 // begin with inf values as I'll pick in every iteration the minimum

  centroid_data.fill(arma::datum::inf);

  arma::mat centroids(clusters, data.n_cols);                                              // save the end-centroids

  arma::mat medoids_vec(1, clusters);                                                      // in case of medoids

  arma::rowvec first = sample_vec(1, 0, data.n_rows - 1, false);                           // choose the first observation at random

  arma::rowvec indices(clusters);

  int idx = first(0);                                                                      // update idx in every iteration

  for (int clust = 0; clust < clusters; clust++) {

    indices(clust) = idx;

    arma::rowvec choose_row = data.row(idx);

    arma::rowvec tmp_dist(data.n_rows);

    for (unsigned int i = 0; i < data.n_rows; i++) {                                                // iterate over the data-rows and calculate the distance between centroid and data

      double tmp_val = kmeans_pp_dist(arma::conv_to< arma::rowvec >::from(data.row(i)), choose_row);

      arma::vec tmp_cent = { centroid_data(i), tmp_val };                                  // then pick the minimum between the current and the previous iters values

      tmp_dist(i) = arma::min(tmp_cent);
    }

    centroid_data = tmp_dist;

    arma::rowvec prob = tmp_dist / arma::sum(tmp_dist);

    arma::rowvec cum_sum = arma::cumsum(prob);                                             // calculate the cumulative probability

    cum_sum = arma::unique(cum_sum);

    arma::rowvec rr = arma::conv_to< arma::rowvec >::from(arma::randu(cum_sum.n_elem));

    double r = rr(0);                                                                      // pick a random value using the randu() function of armadillo

    for (unsigned int j = 0; j < cum_sum.n_elem; j++) {

      if (cum_sum(j) > r) {

        idx = j;                                                                           // update idx

        centroids.row(clust) = arma::conv_to< arma::rowvec >::from(data.row(j));           // add the centroids

        if (medoids) { medoids_vec.col(clust) = j; }

        break;
      }
    }
  }

  if (medoids) {

    return medoids_vec;}

  else {

    return centroids;
  }
}



// calculate fuzzy (soft) clusters
//

// [[Rcpp::export]]
arma::rowvec norm_fuzzy(arma::rowvec vec, double eps) {

  vec = arma::abs(vec);

  vec += eps;

  arma::rowvec d = arma::accu(vec) / vec;

  return d / arma::accu(d);
}




// secondary function for the 'quantile_init_rcpp' function
//

// [[Rcpp::export]]
Rcpp::NumericVector quantile_value(arma::rowvec x, int clusters) {

  arma::vec IDX = arma::regspace<arma::vec>(0.0, arma::max(x) / (clusters - 1), 1.0);

  Rcpp::Environment stats("package:stats");

  Rcpp::Function quantile = stats["quantile"];

  Rcpp::NumericVector ans(IDX.n_elem);

  for(unsigned int i = 0; i < IDX.n_elem; i++){

    ans[i] = Rcpp::as<double>(quantile(x, IDX(i)));
  }

  return ans;
}


// retrurns 'true' if duplicated indices are present in the centroids matrix
//

// [[Rcpp::export]]
bool duplicated_flag(arma::uvec x) {
  
  std::map<double, int> counts;
  
  for (unsigned int i = 0; i < x.n_elem; i++) {
    
    counts[x(i)]++;
  }
  
  arma::rowvec out_vec(counts.size());
  
  int iter = 0;
  
  for (auto& kv : counts) {
    
    out_vec(iter) = kv.second;
    
    iter++;
  }
  
  return arma::accu(out_vec - 1) > 0;
}




// initialization of centroids using the cumulative distance between observations and by removing duplicates
//

// [[Rcpp::export]]
arma::uvec quantile_init_rcpp(arma::mat data, int sample_rows, int clusters) {

  arma::rowvec tmp_idx = sample_vec(sample_rows, 0, data.n_rows - 1, false);

  arma::rowvec vec(sample_rows, arma::fill::zeros);

  for (int i = 0; i < sample_rows; i++) {

    int tmp_val = floorf(std::abs(arma::accu(data.row(tmp_idx(0)) - data.row(tmp_idx(i)))) * 10000000.0) / 10000000;      // use floorf with a float number

    vec(i) = tmp_val;
  }

  Rcpp::NumericVector tmp_rcpp = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(vec));

  arma::rowvec out = Rcpp::as<arma::rowvec>(Rcpp::duplicated(tmp_rcpp));

  arma::uvec tmp_out = arma::find(out == 0);

  arma::rowvec vec1 = arma::conv_to< arma::rowvec >::from(vec(tmp_out));

  arma::uvec vec_sort_idx = arma::sort_index(vec1);

  arma::rowvec vec_sorted = arma::conv_to< arma::rowvec >::from(vec1(vec_sort_idx));

  arma::rowvec idx1 = arma::conv_to< arma::rowvec >::from(tmp_idx(tmp_out));

  arma::rowvec idx1_sorted = arma::conv_to< arma::rowvec >::from(idx1(vec_sort_idx));

  arma::rowvec vec_cumsum = arma::cumsum(vec_sorted) / arma::accu(vec_sorted);

  arma::rowvec quant = Rcpp::as<arma::rowvec> (quantile_value(vec_cumsum, clusters));

  arma::uvec idx_out(clusters, arma::fill::zeros);              // medoids

  for (unsigned int j = 0; j < quant.n_elem; j++) {

    if (j == 0) {

      idx_out(j) = idx1_sorted(0);}

    else {

      arma::uvec tmp_find = arma::find(vec_cumsum <= quant[j]);

      idx_out(j) = idx1_sorted(tmp_find(tmp_find.n_elem - 1));
    }
  }
  
  return idx_out;
}




// it adds medoids/centroids incrementally while checking that they do not already exist in the matrix
// [ only in case that the data set includes duplicated rows is possible that medoids returned will be NaN's ]
//
// the 'tol' parameter is important especially in case of duplicated rows
//

// [[Rcpp::export]]
arma::vec check_medoids(arma::mat data, int clust, double tol = 0.5) {
  
  arma::vec medoids_out(clust);                         // medoids

  medoids_out.fill(arma::datum::nan);

  arma::mat medoids_dat(clust, data.n_cols);

  arma::uvec idx = arma::regspace<arma::uvec>(0, 1, data.n_rows - 1);

  arma::uvec shufl_idx = arma::shuffle(idx);

  arma::mat shufl_dat = data.rows(shufl_idx);

  int count_m = 0;

  for (unsigned int i = 0; i < shufl_dat.n_rows; i++) {

    if (i == 0) {

      medoids_dat.row(count_m) = arma::conv_to< arma::rowvec >::from(shufl_dat.row(i));

      medoids_out(count_m) = shufl_idx(i);

      count_m += 1;

      if (count_m == clust) {

        break;
      }
    }

    else {

      if (count_m == clust) {

        break;
      }

      unsigned int count_diff = 0;

      arma::mat sub_mat = medoids_dat.submat(0, 0, count_m - 1, medoids_dat.n_cols - 1);

      arma::rowvec pot_row = arma::conv_to< arma::rowvec >::from(shufl_dat.row(i));

      for (unsigned int j = 0; j < sub_mat.n_rows; j++) {

        double tmp_diff = arma::accu(arma::abs(arma::conv_to< arma::rowvec >::from(sub_mat.row(j)) - pot_row));

        if (tmp_diff > tol) {

          count_diff += 1;
        }
      }

      if (count_diff == sub_mat.n_rows) {

        medoids_dat.row(count_m) = pot_row;

        medoids_out(count_m) = shufl_idx(i);

        count_m += 1;
      }
    }
  }
  
  return medoids_out;
}



// this function scales and/or centers the data
//

// [[Rcpp::export]]
arma::mat SCALE(arma::mat data, bool mean_center = true, bool sd_scale = true) {

  arma::mat mat_out(data.n_rows, data.n_cols, arma::fill::zeros);

  for (unsigned int i = 0; i < data.n_cols; i++) {

    arma::vec tmp_vec = arma::conv_to< arma::vec >::from(data.col(i));

    if (mean_center) {

      double tmp_mean = arma::as_scalar(arma::mean(tmp_vec));

      tmp_vec -= tmp_mean;
    }

    if (mean_center && sd_scale) {

      double sd = arma::as_scalar(arma::stddev(tmp_vec));

      if (sd == 0.0) { sd = 1.0; }

      tmp_vec /= sd;
    }

    if (!mean_center && sd_scale) {

      double sd = std::sqrt(arma::as_scalar(arma::accu(arma::pow(tmp_vec, 2)) / (tmp_vec.n_elem - 1)));

      if (sd == 0.0) { sd = 1.0; }

      tmp_vec /= sd;
    }

    mat_out.col(i) = tmp_vec;
  }

  if (!mean_center && !sd_scale) {

    return data;}

  else {

    return mat_out;
  }
}

