#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils_rcpp.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>


//using Rcpp::Environment;

//============================================ k-means ====================================================================================



// center, square and take the sum of the data to calculate the total sum of squares
// http://stackoverflow.com/questions/8637460/k-means-return-value-in-r
//

// [[Rcpp::export]]
double tot_ss_data(arma::mat x) {
  
  double tot_ss = 0.0;
  
  for (unsigned int i = 0; i < x.n_cols; i++) {
    
    arma::vec tmp_vec = arma::conv_to< arma::vec >::from(x.col(i));
    
    double tmp_mean = arma::as_scalar(arma::mean(tmp_vec));
    
    tot_ss += arma::accu(arma::pow(tmp_vec - tmp_mean, 2));
  }
  
  return tot_ss;
}



// simple k-means 
// It uses 4 different initializers and the CENTROIDS parameter
//
// fuzzy returns prediction probabilities in a fast way based on the distance between observations and centroids and it isn't similar to fuzzy-c-means
//
// https://github.com/scikit-learn/scikit-learn/blob/51a765a/sklearn/cluster/k_means_.py#L660
//

// [[Rcpp::export]]
Rcpp::List KMEANS_rcpp(arma::mat& data, unsigned int clusters, int num_init = 1, int max_iters = 200, std::string initializer = "kmeans++", bool fuzzy = false, int threads = 1, 
                       
                       bool verbose = false, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, double tol = 1e-4, double eps = 1.0e-6, double tol_optimal_init = 0.5, int seed = 1) {
  
  if (clusters > data.n_rows - 2 || clusters < 1) { Rcpp::stop("the number of clusters should be at most equal to nrow(data) - 2 and never less than 1"); }

  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif
  
  set_seed(seed);             // R's RNG
  
  bool flag = false;
  
  arma::mat CENTROIDS1;
  
  if (CENTROIDS.isNotNull()) {
    
    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);
    
    num_init = 1;
    
    flag = true;
  }
  
  arma::rowvec lst_out;
  
  arma::mat centers_out, lst_fuzzy_out;
  
  arma::rowvec bst_obs(clusters, arma::fill::zeros);
  
  arma::rowvec bst_WCSS(clusters);
  
  bst_WCSS.fill(arma::datum::inf);           // initialize WCSS to Inf, so that in first iteration it can be compared with the minimum of the 'WSSE'
  
  int end_init = 0;
  
  if (verbose && threads == 1) { Rcpp::Rcout << " " << std::endl; }
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (int init_count = 0; init_count < num_init; init_count++) {
    
    arma::mat conv;
    
    if (!flag) {
      
      if (initializer == "kmeans++") {
        
        conv = kmeans_pp_init(data, clusters, false);}
      
      if (initializer == "random") {
        
        arma::uvec samp = arma::conv_to< arma::uvec >::from(sample_vec(clusters, 0, data.n_rows - 1, false));
        
        conv = data.rows(samp);
      }
      
      if (initializer == "quantile_init") {
        
        int samp_ROWS = data.n_rows / 10;
        
        conv = quantile_init_rcpp(data, samp_ROWS, clusters, false);
      }
      
      if (initializer == "optimal_init") {
       
        conv = check_medoids(data, clusters, tol_optimal_init, false);       // tolerance parameter 'tol' here by default equals to 0.5 [ this parameter is important in case of duplicated rows ]
      }
    }
    
    else {
      
      conv = CENTROIDS1;
    }
    
    arma::mat bst_centers, fuzzy_OUT;
    
    arma::rowvec CLUSTERS_OUT, OBS(clusters, arma::fill::zeros), SSE(clusters, arma::fill::zeros);
    
    int iter = 0; 
    
    while(true) {
      
      arma::mat new_centroids(clusters, data.n_cols, arma::fill::zeros), soft_CLUSTERS(data.n_rows, clusters);
      
      arma::rowvec total_WSSE(clusters, arma::fill::zeros), num_obs(clusters, arma::fill::zeros), CLUSTERS(data.n_rows);
      
      for (unsigned int i = 0; i < data.n_rows; i++) {
        
        arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data.row(i)), conv);              // returns a rowvec with the WSSE for each cluster
        
        if (fuzzy) {
          
          soft_CLUSTERS.row(i) = arma::conv_to< arma::rowvec >::from(tmp_vec);
        }
        
        int tmp_idx = MinMat(tmp_vec);                         // returns the index of the tmp_vec with the lowest WSSE
        
        arma::rowvec tmp_row = data.row(i);
        
        total_WSSE(tmp_idx) += tmp_vec(tmp_idx);                // assigns to total_WSSE the minimum cost
        
        num_obs(tmp_idx) += 1;                                 // number of observations in each cluster
        
        new_centroids.row(tmp_idx) +=  tmp_row;                // adds to the corresponding row of the tmp_clusters matrix the row of the data, so that the new centroids can be calculated
        
        CLUSTERS(i) = tmp_idx;
      }
      
      for (unsigned int j = 0; j < clusters; j++) {
        
        new_centroids.row(j) /= arma::as_scalar(num_obs(j));
      }
      
      double tmp_norm = squared_norm(conv - new_centroids);
      
      conv = new_centroids;
      
      if (verbose && threads == 1) { Rcpp::Rcout << "iteration: " << iter + 1 << " --> total WCSS: " << arma::accu(total_WSSE) << "  -->  squared norm: " << tmp_norm << std::endl; }
      
      if (tmp_norm < tol || iter == max_iters - 1) {            // break, if the squared_norm is less than tol or the iter equals to max_iters
        
        CLUSTERS_OUT = CLUSTERS;                                // update clusters
        
        if (fuzzy) {
          
          fuzzy_OUT = soft_CLUSTERS;                            // update soft clustering
        }
        
        SSE = total_WSSE;                                       // update WSSE
        
        OBS = num_obs;                                          // update number of observations
        
        bst_centers = conv;                                     // update centers
        
        break;
      }
      
      iter += 1;
    }
    
    if (arma::accu(SSE) < arma::accu(bst_WCSS)) {
      
      end_init = init_count + 1;
      
      bst_WCSS = SSE;
      
      bst_obs = OBS;
      
      centers_out = bst_centers;
      
      lst_out = CLUSTERS_OUT;
      
      if (fuzzy) {
        
        lst_fuzzy_out = fuzzy_OUT;
      }
    }
    
    if (verbose && threads == 1) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << "===================== end of initialization " << init_count + 1 << " =====================" << std::endl; Rcpp::Rcout << " " << std::endl; }
  }
  
  double tmp_sse = tot_ss_data(data);
  
  if (fuzzy) {
    
    arma::mat fuzzy_mat(lst_fuzzy_out.n_rows, lst_fuzzy_out.n_cols);
    
    for (unsigned int i = 0; i < lst_fuzzy_out.n_rows; i++) {
      
      fuzzy_mat.row(i) = norm_fuzzy(arma::conv_to< arma::rowvec >::from(lst_fuzzy_out.row(i)), eps);
    }
    
    return Rcpp::List::create(Rcpp::Named("clusters") = lst_out, Rcpp::Named("fuzzy_clusters") = fuzzy_mat, Rcpp::Named("centers") = centers_out, Rcpp::Named("total_SSE") = tmp_sse,
                                          
                                          Rcpp::Named("best_initialization") = end_init, Rcpp::Named("WCSS_per_cluster") = bst_WCSS, Rcpp::Named("obs_per_cluster") = bst_obs);
  }
  
  else {
    
    return Rcpp::List::create(Rcpp::Named("clusters") = lst_out, Rcpp::Named("centers") = centers_out, Rcpp::Named("total_SSE") = tmp_sse, 
                                          
                                          Rcpp::Named("best_initialization") = end_init, Rcpp::Named("WCSS_per_cluster") = bst_WCSS, Rcpp::Named("obs_per_cluster") = bst_obs);
  }
}



// the KMEANS_arma returns only the centroids and if openmp exists then it uses all available threads
// the number of columns of the data should be much larger than the number of clusters  (however it works for number_columns == clusters)
// seed_mode is one of : "keep_existing" (I've to give the CENTROIDS), "static_subset", "random_subset", "static_spread", "random_spread"
// in comparison to my 'KMEANS_rcpp' the 'kmeans_arma' does a single initialization and it suggests maximum 10 iterations
//

// [[Rcpp::export]]
arma::mat KMEANS_arma(arma::mat& data, unsigned int clusters, int n_iter, bool verbose, std::string seed_mode = "random_subset",
                      
                      Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, int seed = 1) {
  
  set_seed(seed);             // R's RNG
  
  arma::mat means;
  
  bool status;
  
  if (CENTROIDS.isNotNull() && seed_mode == "keep_existing") {
    
    means = Rcpp::as<arma::mat>(CENTROIDS);
    
    means = means.t();
    
    status = kmeans(means, data.t(), clusters, arma::keep_existing, n_iter, verbose);}
  
  else if (seed_mode == "static_subset") {
    
    status = kmeans(means, data.t(), clusters, arma::static_subset, n_iter, verbose);}
  
  else if (seed_mode == "random_subset") {
    
    status = kmeans(means, data.t(), clusters, arma::random_subset, n_iter, verbose);}
  
  else if (seed_mode == "static_spread") {
    
    status = kmeans(means, data.t(), clusters, arma::static_spread, n_iter, verbose);}
  
  else if (seed_mode == "random_spread") {
    
    status = kmeans(means, data.t(), clusters, arma::random_spread, n_iter, verbose);}
  
  else {
    
    Rcpp::stop("invalid seed_mode");
  }
  
  // if(status == false) {
  //   
  //   Rcpp::Rcout << "clustering failed" << std::endl;
  // }
  
  return means.t();
}



// optimal number of clusters using distortion in kmeans
// https://www.ee.columbia.edu/~dpwe/papers/PhamDN05-kmeans.pdf
//


// [[Rcpp::export]]
Rcpp::List opt_clust_fK(arma::vec sum_distortion, int data_num_cols, double threshold = 0.85) {
  
  arma::rowvec f_K(sum_distortion.n_elem, arma::fill::zeros), a_k_vec(sum_distortion.n_elem - 1, arma::fill::zeros);
  
  for (unsigned int i = 1; i < sum_distortion.n_elem + 1; i++) {
    
    if (i == 1) {
      
      f_K(i - 1) = 1.0;}
    
    else if (i > 1 && sum_distortion(i - 2) != 0.0) {
      
      double a_k = 0.0;
      
      if (i == 2 && data_num_cols > 1) {
        
        a_k = 1.0 - (3.0 / (4.0 * data_num_cols));
      }
      
      if (i > 2 && data_num_cols > 1) {
        
        a_k = a_k_vec(i - 3) + ((1.0 - a_k_vec(i - 3)) / 6.0);
      }
      
      a_k_vec(i - 2) = a_k;
      
      f_K(i - 1) = sum_distortion(i - 1) / (a_k * sum_distortion(i - 2));
    }
    
    else if (i > 1 && sum_distortion(i - 2) == 0.0) {
      
      f_K(i - 1) = 1.0;
    }
    
    else {
      
      Rcpp::stop("iteration " + std::to_string(i) + " out of range");
    }
  }
  
  arma::rowvec opt_cl = arma::conv_to< arma::rowvec >::from(arma::find(f_K <= threshold) + 1);
  
  return Rcpp::List::create(Rcpp::Named("fK_evaluation") = f_K, Rcpp::Named("ak_weight_factor") = a_k_vec, 
                            
                            Rcpp::Named("optimal_number_clusters") = opt_cl);
}



//------------------------------------------------------------------------------------
// dissimilarity and silhouette estimation for k-means based on the euclidean-distance
//------------------------------------------------------------------------------------

// the dissimilarities and silhouette widths of k-means differ from those of the k-medoids, because the latter are based on
// the dissimilarity matrix [ observations of the data are used as medoids ] whereas the centroids of the k-means algorithm
// depend on random initialization and the output-centroids are averages of multiple runs


// intra-cluster dissimilarity
//

// [[Rcpp::export]]
Rcpp::List INTRA_CLUSTER_DISS(arma::mat& data, Rcpp::List CLUSTERS) {
  
  Rcpp::List in_cluster_dist(CLUSTERS.size());
  
  for (int t = 0; t < CLUSTERS.size(); t++) {
    
    arma::uvec tmp_idx = Rcpp::as<arma::uvec>(CLUSTERS[t]);
    
    arma::mat sub_mat = data.rows(tmp_idx);
    
    arma::rowvec dist_VEC(sub_mat.n_rows);
    
    for (unsigned int d = 0; d < sub_mat.n_rows; d++) {
      
      double dist_ = 0.0;
      
      for (unsigned int w = 0; w < sub_mat.n_rows; w++) {
        
        dist_ += std::sqrt(arma::as_scalar(arma::accu(arma::square((arma::conv_to< arma::vec >::from(sub_mat.row(d)) - arma::conv_to< arma::vec >::from(sub_mat.row(w)))))));
      }
      
      dist_VEC(d) = dist_ / (sub_mat.n_rows - 1);                  // subtract 1 , as I do calculate the distance of the observation with itself
    }
    
    in_cluster_dist[t] = dist_VEC;
  }
  
  return in_cluster_dist;
}



// convert Rcpp::NumericMatrix to armadillo matrix
// [ http://stackoverflow.com/questions/31691130/conversion-of-r-matrices-to-armadillo-is-really-slow ]

// [[Rcpp::export]]
arma::mat Rcpp_2arma_mat(Rcpp::NumericMatrix x) {
  
  arma::mat res = Rcpp::as<arma::mat>(x);
  
  return(res);
}


// secondary function to calculate the silhouette metric
//

// [[Rcpp::export]]
Rcpp::List SILHOUETTE_metric(arma::mat& data, arma::vec CLUSTER, Rcpp::List tmp_clust, Rcpp::List in_cluster_dist) {

  arma::vec unq_values = arma::unique(CLUSTER);
  
  Rcpp::Environment gtools("package:gtools");
  
  Rcpp::Function permutations = gtools["permutations"];            // the permutations-function is necessary to calculate the minimum distance of each observation to the nearest cluster using as permutation pairs the clusters

  Rcpp::NumericMatrix idx = permutations(unq_values.n_elem, 2);
  
  arma::mat IDX = Rcpp_2arma_mat(idx) - 1;                         // subtract 1 from each cell, because gtools::permutations is an R function and indexing in cpp begins from 0
  
  Rcpp::List OUT_cluster_dist(IDX.n_rows);                         // calculate the average-outer-cluster-distances
  
  for (unsigned int t = 0; t < IDX.n_rows; t++) {
    
    arma::uvec tmp_idx1ST = Rcpp::as<arma::uvec>(tmp_clust[IDX(t,0)]);
    
    arma::uvec tmp_idx2ND = Rcpp::as<arma::uvec>(tmp_clust[IDX(t,1)]);
    
    arma::mat sub_mat1ST = data.rows(tmp_idx1ST);
    
    arma::mat sub_mat2ND = data.rows(tmp_idx2ND);
    
    arma::rowvec dist_VEC(sub_mat1ST.n_rows);
    
    for (unsigned int d = 0; d < sub_mat1ST.n_rows; d++) {
      
      double dist_out = 0.0;
      
      for (unsigned int w = 0; w < sub_mat2ND.n_rows; w++) {
        
        dist_out += std::sqrt(arma::as_scalar(arma::accu(arma::square(arma::conv_to< arma::vec >::from(sub_mat1ST.row(d)) - arma::conv_to< arma::vec >::from(sub_mat2ND.row(w))))));
      }
      
      dist_VEC(d) = dist_out / sub_mat2ND.n_rows;
    }
    
    OUT_cluster_dist[t] = dist_VEC;
  }
  
  // the 'OUT_cluster_dist' function returns the permutations of the clusters, I have to find and return the minimum distance of
  // each observation to the next nearest cluster using the permutation pairs
  // [ the output list 'OUT_cluster_dist1' should have the same dimensions as the 'INTRA_cluster_dissimilarity' list ]
  
  Rcpp::List OUT_cluster_dist1(in_cluster_dist.size());
  
  for (int iter = 0; iter < in_cluster_dist.size(); iter++) {
    
    arma::uvec tmp_find = arma::find(IDX.col(0) == iter);
    
    arma::mat TMP_mat;
    
    for (unsigned int n = 0; n < tmp_find.n_elem; n++) {
      
      arma::rowvec tmp_rowv = Rcpp::as<arma::rowvec>(OUT_cluster_dist[tmp_find(n)]);
      
      TMP_mat.set_size(tmp_find.n_elem, tmp_rowv.n_elem);
      
      TMP_mat.row(n) = tmp_rowv;
    }
    
    OUT_cluster_dist1[iter] = arma::min(TMP_mat);
  }
  
  // silhouette calculation
  
  Rcpp::List silhouet(OUT_cluster_dist1.size());
  
  for (int i = 0; i < OUT_cluster_dist1.size(); i++) {
    
    Rcpp::NumericVector midl_inter = in_cluster_dist[i];
    
    Rcpp::NumericVector midl_outer = OUT_cluster_dist1[i];
    
    Rcpp::NumericVector midl_silh(midl_inter.size());
    
    for (int j = 0; j < midl_inter.size(); j++) {
      
      midl_silh[j] = calc_silhouette(midl_inter[j], midl_outer[j]);
    }
    
    silhouet[i] = midl_silh;
  }
  
  return silhouet;
}



// dissimilarity or silhouette evaluation of clustering [ dissimilarity returns faster as it only requires intra-distances to be calculated ]
// http://blog.data-miners.com/2011/03/cluster-silhouettes.html
//

// [[Rcpp::export]]
Rcpp::List evaluation_rcpp(arma::mat& data, arma::vec CLUSTER, bool silhouette = false) {

  Rcpp::List tmp_clust = cluster_indices(CLUSTER);
  
  Rcpp::List in_cluster_dist = INTRA_CLUSTER_DISS(data, tmp_clust);
  
  if (!silhouette) {
    
    return(Rcpp::List::create(Rcpp::Named("clusters") = arma::conv_to< arma::rowvec >::from(CLUSTER),
                              
                              Rcpp::Named("cluster_indices") = tmp_clust,                                    // use the data indices, otherwise difficult to match clusters with silhouette or dissimilarity coefficients
                              
                              Rcpp::Named("INTRA_cluster_dissimilarity") = in_cluster_dist));                // the lower the better
  }
  
  else {
    
    Rcpp::List silhouet_out = SILHOUETTE_metric(data, CLUSTER, tmp_clust, in_cluster_dist);
    
    return(Rcpp::List::create(Rcpp::Named("clusters") = arma::conv_to< arma::rowvec >::from(CLUSTER),
                              
                              Rcpp::Named("cluster_indices") = tmp_clust,                                       // use the data indices, otherwise difficult to match clusters with silhouette or dissimilarity coefficients
                              
                              Rcpp::Named("INTRA_cluster_dissimilarity") = in_cluster_dist,                     // the lower the better
                              
                              Rcpp::Named("silhouette") = silhouet_out));                                       // the higher the better [ range between -1 and +1 ]
  }
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
Rcpp::List mini_batch_kmeans(arma::mat& data, unsigned int clusters, int batch_size, int max_iters, int num_init = 1, double init_fraction = 1.0, std::string initializer = "kmeans++", 
                             
                             int early_stop_iter = 10, bool verbose = false, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, double tol = 1e-4, double tol_optimal_init = 0.5, int seed = 1) {
  
  set_seed(seed);             // R's RNG
  
  if (clusters > data.n_rows - 2 || clusters < 1) { Rcpp::stop("the number of clusters should be at most equal to nrow(data) - 2 and never less than 1"); }
  
  bool flag = false;
  
  arma::mat CENTROIDS1;
  
  if (CENTROIDS.isNotNull()) {
    
    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);
    
    num_init = 1;
    
    flag = true;
  }
  
  arma::mat centers_out;
  
  arma::rowvec iter_before_stop(num_init, arma::fill::zeros);
  
  int end_init = 0;
  
  arma::rowvec bst_WCSS(clusters);
  
  bst_WCSS.fill(arma::datum::inf);           // initialize WCSS to Inf, so that in first iteration it can be compared with the minimum of the 'WSSE'
  
  if (verbose) { Rcpp::Rcout << " " << std::endl; }
  
  for (int init = 0; init < num_init; init++) {
    
    arma::mat update_centroids;
    
    if (!flag) {
      
      if (initializer == "kmeans++") {
        
        if (init_fraction == 1.0) {
          
          update_centroids = kmeans_pp_init(data, clusters, false);}
        
        if (init_fraction < 1.0 && init_fraction > 0.0) {
          
          int fract = std::ceil(data.n_rows * init_fraction);
          
          arma::uvec samp_init = arma::conv_to< arma::uvec >::from(sample_vec(fract, 0, data.n_rows - 1, false));
          
          arma::mat tmp_mat_init = data.rows(samp_init);
          
          update_centroids = kmeans_pp_init(tmp_mat_init, clusters, false);
        }
      }
      
      if (initializer == "random") {
        
        arma::uvec samp = arma::conv_to< arma::uvec >::from(sample_vec(clusters, 0, data.n_rows - 1, false));
        
        update_centroids = data.rows(samp);
      }
      
      if (initializer == "quantile_init") {
        
        int samp_ROWS = 0.0;
        
        if (init_fraction == 1.0) {
          
          samp_ROWS = data.n_rows / 10;}
        
        else {
          
          samp_ROWS = std::ceil(data.n_rows * init_fraction);
        }
        
        update_centroids = quantile_init_rcpp(data, samp_ROWS, clusters, false);
      }
      
      if (initializer == "optimal_init") {
        
        update_centroids = check_medoids(data, clusters, tol_optimal_init, false);       // tolerance parameter 'tol' here by default equals to 0.5 [ this parameter is important in case of duplicated rows ]
      }
    }
    
    else {
      
      update_centroids = CENTROIDS1;
    }
    
    arma::mat previous_centroids = update_centroids;
    
    arma::mat output_centroids;
    
    arma::rowvec output_SSE;
    
    double previous_cost = arma::datum::inf;
    
    int increment_early_stop = 0;
    
    int count = 0;
    
    for (int i = 0; i < max_iters; i++) {
      
      arma::uvec batch_idx = arma::conv_to< arma::uvec >::from(sample_vec(batch_size, 0, data.n_rows - 1, false));
      
      arma::mat batch_data = data.rows(batch_idx);
      
      arma::rowvec total_SSE(clusters, arma::fill::zeros);
      
      arma::rowvec CLUSTERS(batch_data.n_rows);
      
      for (unsigned int j = 0; j < batch_data.n_rows; j++) {
        
        arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(batch_data.row(j)), update_centroids);         // returns a rowvec with the SSE for each cluster
        
        int tmp_idx = MinMat(tmp_vec);                                                                              // returns the index of the tmp_vec with the lowest SSE
        
        total_SSE(tmp_idx) += tmp_vec(tmp_idx);                                                                     // assigns to total_SSE the minimum cost
        
        CLUSTERS(j) = tmp_idx;
      }
      
      arma::rowvec cluster_counts(clusters, arma::fill::zeros);
      
      double eta = 0.0;
      
      for (unsigned int k = 0; k < batch_data.n_rows; k++) {
        
        int idx = CLUSTERS(k);
        
        cluster_counts(idx) += 1;
        
        eta = 1.0 / cluster_counts(idx);
        
        arma::rowvec tmp_row = (1.0 - eta) * arma::conv_to< arma::rowvec >::from(update_centroids.row(idx)) + eta * arma::conv_to< arma::rowvec >::from(batch_data.row(k));
        
        update_centroids.row(idx) = tmp_row;
      }
      
      double tmp_norm = squared_norm(previous_centroids - update_centroids);
      
      double calc_cost = arma::accu(total_SSE);
      
      if (verbose) { Rcpp::Rcout << "iteration: " << i + 1 << "  --> total WCSS: " << calc_cost << "  -->  squared norm: " << tmp_norm << std::endl; }
      
      count = i;
      
      if (calc_cost < previous_cost) {
        
        previous_cost = calc_cost;
        
        increment_early_stop = 0;
        
        output_centroids = update_centroids;                // assign end-centroids and SSE when WCSS is minimal
        
        output_SSE = total_SSE;
      }
      
      if (calc_cost > previous_cost) {
        
        increment_early_stop += 1;
      }
      
      if (tmp_norm < tol || i == max_iters - 1 || increment_early_stop == early_stop_iter - 1) {
        
        // output_centroids = update_centroids;            // assign end-centroids and SSE when early_stop_iter == increment_early_stop [ repeated calculation of adjusted rand index shows slightly better results for the previous case ]
        // 
        // output_SSE = total_SSE;
        
        break;
      }
    }
    
    if (arma::accu(output_SSE) < arma::accu(bst_WCSS)) {
      
      end_init = init + 1;
      
      bst_WCSS = output_SSE;
      
      centers_out = output_centroids;
    }
    
    iter_before_stop(init) = count + 1;
    
    if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << "===================== end of initialization " << init + 1 << " =====================" << std::endl; Rcpp::Rcout << " " << std::endl; }
  }
  
  return Rcpp::List::create(Rcpp::Named("centroids") = centers_out, Rcpp::Named("WCSS_per_cluster") = bst_WCSS, 
                            
                            Rcpp::Named("best_initialization") = end_init, Rcpp::Named("iters_per_initialization") = iter_before_stop);
}




// predict function for mini-batch-kmeans, which takes a matrix of centroids
//

// [[Rcpp::export]]
Rcpp::List Predict_mini_batch_kmeans(arma::mat& data, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, bool fuzzy = false, double eps = 1.0e-6) {
  
  arma::mat CENTROIDS1;
  
  if (CENTROIDS.isNotNull()) {
    
    CENTROIDS1 = Rcpp::as<arma::mat>(CENTROIDS);
  }
  
  arma::rowvec CLUSTERS(data.n_rows);
  
  arma::mat soft_CLUSTERS(data.n_rows, CENTROIDS1.n_rows);
  
  for (unsigned int j = 0; j < data.n_rows; j++) {
    
    arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data.row(j)), CENTROIDS1);                  // returns a rowvec with the SSE for each cluster
    
    soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);
    
    int tmp_idx = MinMat(tmp_vec);                                                                        // returns the index of the tmp_vec with the lowest SSE
    
    CLUSTERS(j) = tmp_idx;
  }
  
  if (fuzzy) {
    
    arma::mat fuzzy_mat(soft_CLUSTERS.n_rows, soft_CLUSTERS.n_cols);
    
    for (unsigned int i = 0; i < soft_CLUSTERS.n_rows; i++) {
      
      fuzzy_mat.row(i) = norm_fuzzy(arma::conv_to< arma::rowvec >::from(soft_CLUSTERS.row(i)), eps);
    }
    
    return Rcpp::List::create(Rcpp::Named("clusters") = CLUSTERS, Rcpp::Named("fuzzy_clusters") = fuzzy_mat);}
  
  else {
    
    return Rcpp::List::create(Rcpp::Named("clusters") = CLUSTERS);
  }
}




//============================================ Gaussian Mixture Models (GMM) =========================================================================




// the GMM_arma returns the centroids, covariance matrix and the weights. If openmp exists then it uses all available threads
// distance is one of : "eucl_dist", "maha_dist"
// seed_mode is one of : "static_subset", "random_subset", "static_spread", "random_spread"  [ I excluded the "keep_existing" seed_mode ] -- user-defined parameter setting is not enabled
//

// [[Rcpp::export]]
Rcpp::List GMM_arma(arma::mat& data, unsigned int gaussian_comps, std::string dist_mode, std::string seed_mode, int km_iter, int em_iter, 
                    
                    bool verbose, double var_floor = 1e-10, int seed = 1) {
  
  arma::wall_clock timer;
  
  timer.tic();
  
  set_seed(seed);             // R's RNG
  
  arma::gmm_diag model;
  
  arma::mat means;
  
  bool status;
  
  if (seed_mode == "static_subset" && dist_mode == "eucl_dist") {
    
    status = model.learn(data.t(), gaussian_comps, arma::eucl_dist, arma::static_subset, km_iter, em_iter, var_floor, verbose);}
  
  else if (seed_mode == "random_subset" && dist_mode == "eucl_dist") {
    
    status = model.learn(data.t(), gaussian_comps, arma::eucl_dist, arma::random_subset, km_iter, em_iter, var_floor, verbose);}
  
  else if (seed_mode == "static_spread" && dist_mode == "eucl_dist") {
    
    status = model.learn(data.t(), gaussian_comps, arma::eucl_dist, arma::static_spread, km_iter, em_iter, var_floor, verbose);}
  
  else if (seed_mode == "random_spread" && dist_mode == "eucl_dist") {
    
    status = model.learn(data.t(), gaussian_comps, arma::eucl_dist, arma::random_spread, km_iter, em_iter, var_floor, verbose);}
  
  else if (seed_mode == "static_subset" && dist_mode == "maha_dist") {
    
    status = model.learn(data.t(), gaussian_comps, arma::maha_dist, arma::static_subset, km_iter, em_iter, var_floor, verbose);}
  
  else if (seed_mode == "random_subset" && dist_mode == "maha_dist") {
    
    status = model.learn(data.t(), gaussian_comps, arma::maha_dist, arma::random_subset, km_iter, em_iter, var_floor, verbose);}
  
  else if (seed_mode == "static_spread" && dist_mode == "maha_dist") {
    
    status = model.learn(data.t(), gaussian_comps, arma::maha_dist, arma::static_spread, km_iter, em_iter, var_floor, verbose);}
  
  else if (seed_mode == "random_spread" && dist_mode == "maha_dist") {
    
    status = model.learn(data.t(), gaussian_comps, arma::maha_dist, arma::random_spread, km_iter, em_iter, var_floor, verbose);}
  
  else {
    
    Rcpp::stop("Invalid seed_mode OR dist_mode. Valid 'seed_modes' are : 'static_subset', 'random_subset', 'static_spread' and 'random_spread'. Valid 'dist_modes' are : 'eucl_dist' and 'maha_dist'.");
  }
  
  // if(status == false) {
  //   
  //   Rcpp::Rcout << "learning failed" << std::endl;
  // }
  
  arma::mat loglik(data.n_rows, gaussian_comps, arma::fill::zeros);
  
  for (unsigned int j = 0; j < gaussian_comps; j++) {
    
    loglik.col(j) = arma::conv_to< arma::vec >::from(model.log_p(data.t(), j));
  }
  
  double n = timer.toc();
  
  if (verbose) { Rcpp::Rcout << "\ntime to complete : " << n << "\n" << std::endl; }
  
  return Rcpp::List::create( Rcpp::Named("centroids") = model.means.t(), Rcpp::Named("covariance_matrices") = model.dcovs.t(),           // each row of the 'covariance_matrices' is a different covariance matrix, use diag() to build each square diagonal matrix 
                             
                             Rcpp::Named("weights") = model.hefts.t(), Rcpp::Named("Log_likelihood_raw") = loglik, 
                             
                             Rcpp::Named("avg_Log_likelihood_DATA") = model.avg_log_p(data.t(), gaussian_comps - 1) );
}



// take a diagonal matrix in form of a vector and build a square diagonal matrix, then invert it
//

// [[Rcpp::export]]
arma::mat INV_COV(arma::vec COV_VEC) {
  
  arma::mat tmp_cov = arma::diagmat(COV_VEC);
  
  return arma::inv(tmp_cov);
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
// The predictions are based on centroids, covariance matrix and weights using the formulas and not the log-likelihoods.
//

// [[Rcpp::export]]
Rcpp::List predict_MGausDPDF(arma::mat data, arma::mat CENTROIDS, arma::mat COVARIANCE, arma::vec WEIGHTS, double eps = 1.0e-8) {
  
  arma::mat gaus_mat(data.n_rows, WEIGHTS.n_elem, arma::fill::zeros);
  
  arma::mat gaus_mat_log_lik(data.n_rows, WEIGHTS.n_elem, arma::fill::zeros);
  
  for (unsigned int i = 0; i < WEIGHTS.n_elem; i++) {
    
    arma::vec gaus_vec(data.n_rows, arma::fill::zeros);
    
    arma::vec gaus_vec_log(data.n_rows, arma::fill::zeros);
    
    for (unsigned int j = 0; j < data.n_rows; j++) {
      
      double n = data.n_cols;
      
      arma::vec tmp_vec = (arma::conv_to< arma::vec >::from(data.row(j)) - arma::conv_to< arma::vec >::from(CENTROIDS.row(i)));
      
      arma::mat tmp_cov_mt = arma::diagmat(arma::conv_to< arma::vec >::from(COVARIANCE.row(i)));
      
      double tmp_val = 1.0 / std::sqrt(2.0 * arma::datum::pi * arma::det(tmp_cov_mt));               // use determinant to get a single value
      
      double inner_likelih = 0.5 * (arma::as_scalar(tmp_vec.t() * INV_COV(arma::conv_to< arma::vec >::from(COVARIANCE.row(i))) * arma::conv_to< arma::mat >::from(tmp_vec)));
      
      gaus_vec_log(j) = -(n / 2.0) * std::log(2.0 * arma::datum::pi) - (1.0 / 2.0) * (std::log(arma::det(tmp_cov_mt))) - inner_likelih;
      
      gaus_vec(j) = tmp_val * std::exp(-inner_likelih);
    }
    
    gaus_mat.col(i) = arma::as_scalar(WEIGHTS(i)) * gaus_vec;
    
    gaus_mat_log_lik.col(i) = gaus_vec_log;
  }
  
  arma::mat loglik1(data.n_rows, WEIGHTS.n_elem, arma::fill::zeros);
  
  arma::rowvec loglik2(data.n_rows, arma::fill::zeros);
  
  for (unsigned int j = 0; j < loglik1.n_rows; j++) {
    
    arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(gaus_mat.row(j)) + eps;
    
    tmp_vec /= arma::sum(tmp_vec);                                                               // normalize row-data to get probabilities
    
    loglik1.row(j) = tmp_vec;                                                                    // assign probabilities
    
    arma::uvec log_lik_label = arma::find(tmp_vec == arma::max(tmp_vec));
    
    loglik2(j) = arma::as_scalar(log_lik_label(0));                                              // assign labels
  }
  
  return Rcpp::List::create( Rcpp::Named("Log_likelihood_raw") = gaus_mat_log_lik, Rcpp::Named("cluster_proba") = loglik1, 
                             
                             Rcpp::Named("cluster_labels") = loglik2 );
}




// function to calculate bic-aic
//

// [[Rcpp::export]]
arma::rowvec GMM_arma_AIC_BIC(arma::mat& data, unsigned int max_clusters, std::string dist_mode, std::string seed_mode,
                              
                              int km_iter, int em_iter, bool verbose, double var_floor = 1e-10, std::string criterion = "AIC", int seed = 1) {
  
  set_seed(seed);             // R's RNG
  
  arma::rowvec evaluate_comps(max_clusters, arma::fill::zeros), aic_avg_weights(max_clusters, arma::fill::zeros);
  
  for (unsigned int i = 0; i < max_clusters; i++) {
    
    if (verbose) { Rcpp::Rcout << "iteration: " << i + 1 << std::endl; }
    
    Rcpp::List gmm = GMM_arma(data, i + 1, dist_mode, seed_mode, km_iter, em_iter, false, var_floor = 1e-10);
    
    arma::mat loglik = Rcpp::as<arma::mat> (gmm[3]);
    
    arma::rowvec weights = Rcpp::as<arma::rowvec> (gmm[2]);
    
    arma::rowvec log_sum_exp(data.n_rows, arma::fill::zeros);                             
    
    for (unsigned int i = 0; i < loglik.n_rows; i++) {
      
      arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(loglik.row(i));
      
      tmp_vec += arma::log(weights);                                     // https://github.com/scikit-learn/scikit-learn/blob/51a765acfa4c5d1ec05fc4b406968ad233c75162/sklearn/utils/extmath.py
      
      double max_row = arma::max(tmp_vec);
      
      tmp_vec = arma::exp(tmp_vec - max_row);
      
      double tmp_log = arma::sum(tmp_vec);
      
      log_sum_exp(i) = max_row + std::log(tmp_log);
    }
    
    arma::mat centers = Rcpp::as<arma::mat> (gmm[0]);

    if (criterion == "AIC") {
      
      evaluate_comps(i) = -2.0 * arma::accu(log_sum_exp) + 2.0 * centers.n_rows * centers.n_cols;
    }
    
    if (criterion == "BIC") {
      
      evaluate_comps(i) = -2.0 * arma::accu(log_sum_exp) + std::log(data.n_rows) * centers.n_rows * centers.n_cols;
    }
  }
  
  return evaluate_comps;
}




//============================================ cluster Medoids ===================================================================================





// dissimilarity matrix using various distance metrics [ test cases for all methods as I changed some parameters ]
//

// [[Rcpp::export]]
arma::mat dissim_mat(arma::mat& data, std::string method, double minkowski_p = 1.0, bool upper = true, bool diagonal = true, int threads = 1, double eps = 1.0e-6) {
  
  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif
  
  arma::mat cov_mat;
  
  if (method == "mahalanobis") {
    
    cov_mat = arma::inv(arma::cov(data));
  }
  
  arma::mat mt(data.n_rows, data.n_rows);
  
  mt.fill(arma::datum::nan);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (unsigned int i = 0; i < data.n_rows - 1; i++) {
    
    int k = i;
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (unsigned int j = k + 1; j < data.n_rows; j++) {
      
      double tmp_idx;
      
      if (method == "euclidean") {
        
        tmp_idx = std::sqrt(arma::as_scalar(arma::accu(arma::square((data.row(i) - data.row(j))))));
      }
      
      else if (method == "manhattan") {
        
        tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - data.row(j)))));
      }
      
      else if (method == "chebyshev") {
        
        tmp_idx = arma::as_scalar(arma::max(arma::abs((data.row(i) - data.row(j)))));
      }
      
      else if (method == "canberra") {
        
        tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - data.row(j)) + eps)/(arma::abs(data.row(i)) + arma::abs(data.row(j)) + eps)));                 // added 1.0e-6 otherwise rstudio crashes
      }
      
      else if (method == "braycurtis") {
        
        tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - data.row(j))))/(arma::accu(arma::abs(data.row(i))) + arma::accu(arma::abs(data.row(j)))));
      }
      
      else if (method == "pearson_correlation") {
        
        tmp_idx = arma::as_scalar(1.0 - arma::cor(data.row(i), data.row(j)));
      }
      
      else if (method == "simple_matching_coefficient") {                                                    // for binary data
        
        double a = eps;
        double d = eps;
        
        for (unsigned int t = 0; t < data.row(i).n_elem; t++) {
          
          if (data.row(i)(t) == 1 && data.row(j)(t) == 1) {
            
            a += 1.0;}
          
          if (data.row(i)(t) == 0 && data.row(j)(t) == 0) {
            
            d += 1.0;
          }
        }
        
        tmp_idx = 1.0 - ((a + d) / data.row(i).n_elem);
      }
      
      else if (method == "minkowski") {                                                                                     // by default the order of the minkowski parameter equals k
        
        tmp_idx = std::pow(arma::as_scalar(arma::accu(arma::pow(arma::abs((data.row(i) - data.row(j))), minkowski_p))), 1.0 / minkowski_p);
      }
      
      else if (method == "hamming") {                                                                                     // for binary data
        
        tmp_idx = arma::as_scalar(arma::accu(data.row(i) != data.row(j))/(data.row(i).n_elem * 1.0));
      }
      
      else if (method == "mahalanobis") {                                                                                     // first create covariance matrix from data
        
        tmp_idx = arma::as_scalar(std::sqrt(arma::as_scalar(((data.row(i) - data.row(j)) * cov_mat) * (data.row(i) - data.row(j)).t())));
      }
      
      else if (method == "jaccard_coefficient") {                                                                                     // for binary data
        
        double a = eps;
        double b = eps;
        double c = eps;
        
        for (unsigned int t = 0; t < data.row(i).n_elem; t++) {
          
          if (data.row(i)(t) == 1 && data.row(j)(t) == 1) {
            
            a += 1.0;}
          
          if (data.row(i)(t) == 1 && data.row(j)(t) == 0) {
            
            b += 1.0;}
          
          if (data.row(i)(t) == 0 && data.row(j)(t) == 1) {
            
            c += 1.0;
          }
        }
        
        tmp_idx = 1.0 - (a / (a + b + c));
      }
      
      else if (method == "Rao_coefficient") {                                                                                     // for binary data
        
        double a = eps;
        
        for (unsigned int t = 0; t < data.row(i).n_elem; t++) {
          
          if (data.row(i)(t) == 1 && data.row(j)(t) == 1) {
            
            a += 1.0;
          }
        }
        
        tmp_idx = 1.0 - (a / data.row(i).n_elem);
      }
      
      else {
        
        tmp_idx = 0;                                             // default = 0; create exceptions in R, so that tmp_idx != 0;
      }
      
      if ( tmp_idx != tmp_idx ) {                                // handling of NAs, if NaN then distance 1.0 [  NaN will compare false to everything, including itself ], http://stackoverflow.com/questions/11569337/using-an-if-statement-to-switch-nan-values-in-an-array-to-0-0]
        
        tmp_idx = 1.0;
      }
      
      mt(j,i) = tmp_idx;
      
      if (upper) {
        
        mt(i,j) = tmp_idx;
      }
    }
  }
  
  if (diagonal) {
    
    mt.diag().zeros();
  }
  
  return(mt);
}



// secondary function to check if integer in vector of integers
//

// [[Rcpp::export]]
bool boolean_function(arma::rowvec x, int y) {
  
  bool flag = false;
  
  for (unsigned int i = 0; i < x.n_elem; i++) {
    
    if (x(i) == y) {
      
      flag = true;
      
      break;
    }
  }
  
  return flag;
}



// intra-outer-silhouette cluster values AND clustering-statistics
//

// [[Rcpp::export]]
Rcpp::List silhouette_matrix(arma::mat data, arma::rowvec end_indices_vec, arma::rowvec end_cost_vec, int threads = 1) {
  
  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif
  
  arma::rowvec sorted_medoids = arma::unique(end_indices_vec);
  
  arma::rowvec sorted_medoids_increment = arma::regspace<arma::rowvec> (1, 1, sorted_medoids.n_elem);
  
  arma::mat silhouette_matrix(end_indices_vec.n_elem, 7, arma::fill::zeros);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (unsigned int f = 0; f < end_indices_vec.n_elem; f++) {
    
    double intra_clust_diss = 0.0, diameter = 0.0, separation = arma::datum::inf;
    
    int intra_count_obs = 0;
    
    arma::mat outer_clust_diss(3, sorted_medoids.n_elem, arma::fill::zeros);              // temporary-outer-dissimilarities-matrix
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (unsigned int t = 0; t < end_indices_vec.n_elem; t++) {
      
      if (f != t) {
        
        if (end_indices_vec(f) == end_indices_vec(t)) {
          
          intra_clust_diss += data(f,t);
          
          if (diameter < data(f,t)) { diameter = data(f,t); }
          
          intra_count_obs += 1;
        }
        
        if (end_indices_vec(f) != end_indices_vec(t)) {
          
          for (unsigned int k = 0; k < sorted_medoids.n_elem; k++) {
            
            if (end_indices_vec(t) == sorted_medoids(k)) {
              
              outer_clust_diss(0,k) = sorted_medoids_increment(k);
              
              outer_clust_diss(1,k) += data(f,t);
              
              if (separation > data(f,t)) { separation = data(f,t); }
              
              outer_clust_diss(2,k) += 1;
            }
          }
        }
      }
    }
    
    outer_clust_diss.row(1) /= outer_clust_diss.row(2);                                                // returns Inf if the outer-cluster-count is 1 (a single value as cluster index)
    
    arma::uvec opt_idx = arma::find(outer_clust_diss.row(1) == arma::min(outer_clust_diss.row(1)));
    
    silhouette_matrix(f, 0) = end_indices_vec(f);                                                      // 1st column clusters
    
    silhouette_matrix(f, 1) = outer_clust_diss(0,opt_idx(0));                                          // 2nd column neighboring clusters [ using outer-cluster-dissimilarities ]
    
    silhouette_matrix(f, 2) = intra_clust_diss / intra_count_obs;                                      // 3rd column intra-cluster-dissimilarities
    
    silhouette_matrix(f, 3) = outer_clust_diss(1,opt_idx(0));                                          // 4th column outer-cluster-dissimilarities
    
    silhouette_matrix(f, 4) = calc_silhouette(silhouette_matrix(f,2), silhouette_matrix(f,3));         // 5th column silhouette widths
    
    silhouette_matrix(f, 5) = diameter;                                                                // diameter (maximal dissimilarity between two obs. of the same cluster)
    
    silhouette_matrix(f, 6) = separation;                                                              // separation (minimal dissimilarity between two obs. of different clusters)
    
  }
  
  arma::mat clustering_stats(6, sorted_medoids.n_elem, arma::fill::zeros);
  
  clustering_stats.row(5).fill(arma::datum::inf);
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2)
  #endif
  for (unsigned int s = 0; s < silhouette_matrix.n_rows; s++) {
    
    for (unsigned int g = 0; g < sorted_medoids.n_elem; g++) {
      
      if (sorted_medoids(g) == silhouette_matrix(s, 0)) {
        
        clustering_stats(0, g) = sorted_medoids(g);                                                                       // clustering labels
        
        clustering_stats(1, g) += 1;                                                                                      // number of obs in each cluster
        
        if (clustering_stats(2, g) < end_cost_vec(s)) { clustering_stats(2, g) = end_cost_vec(s); }                       // maximum cluster dissimilarity
        
        clustering_stats(3, g) += end_cost_vec(s);                                                                        // first sum dissimilarities, so that I can get the average cluster dissimilarity
        
        if (clustering_stats(4, g) < silhouette_matrix(s, 5)) { clustering_stats(4, g) = silhouette_matrix(s, 5); }       // diameter of cluster
        
        if (clustering_stats(5, g) > silhouette_matrix(s, 6)) { clustering_stats(5, g) = silhouette_matrix(s, 6); }       // separation of cluster
      }
    }
  }
  
  clustering_stats.row(3) /= clustering_stats.row(1);                                                                     // average intra cluster dissimilarity
  
  return Rcpp::List::create(Rcpp::Named("silhouette_matrix") = silhouette_matrix, Rcpp::Named("clustering_stats") = clustering_stats.t());
}



// remove items from a vector using item-values from a second vector
//

// [[Rcpp::export]]
arma::uvec subset_vec(arma::uvec x, arma::uvec y) {
  
  std::vector<double> vec = Rcpp::as<std::vector<double> >(Rcpp::wrap(x));
  
  for (unsigned int i = 0; i < y.n_elem; i++) {
    
    vec.erase(std::remove(vec.begin(), vec.end(), y(i)), vec.end());
  }
  
  return Rcpp::as<arma::uvec> (Rcpp::wrap(vec));
}




// the 'ClusterMedoids' function should accept either a matrix or a dissimilarity matrix
// the dissimilarity matrix in comparison to a single matrix will have nrows == ncols AND diagonal == 0.0 [ it will be also a matrix ]
//

// [[Rcpp::export]]
Rcpp::List ClusterMedoids(arma::mat& data, unsigned int clusters, std::string method, double minkowski_p = 1.0, int threads = 1, bool verbose = false, bool swap_phase = false, 
                          
                          bool fuzzy = false, int seed = 1) {
  
  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif
  
  set_seed(seed);             // R's RNG
  
  bool flag_dissim_mat = true;
  
  double sum_diag = arma::accu(data.diag());
  
  if (data.n_rows != data.n_cols && sum_diag != 0.0) {
    
    data = dissim_mat(data, method, minkowski_p, true, true, threads, 1.0e-6);
    
    flag_dissim_mat = false;
  }
  
  arma::vec sum_first = arma::conv_to< arma::vec >::from(arma::sum(data, 0));
  
  arma::uvec tmp_first_idx = arma::find(sum_first == arma::min(sum_first));
  
  int first_medoid = tmp_first_idx(0);
  
  double total_cost = arma::accu(data.col(first_medoid));
  
  arma::rowvec end_cost_vec = arma::conv_to< arma::rowvec >::from(data.row(first_medoid));
  
  arma::rowvec end_indices_vec(end_cost_vec.n_elem);
  
  end_indices_vec.fill(first_medoid);
  
  unsigned int count_clusters = 1;
  
  if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << "medoid " << first_medoid + 1 << " was added. Current dissimilarity of build phase  --> " << total_cost << std::endl; }
  
  while (true) {
    
    if (count_clusters > clusters - 1) { break; }
    
    double previous_cost = arma::datum::inf;
    
    int previous_idx = -1;
    
    arma::rowvec end_update_cost(end_cost_vec.n_elem), end_update_idx(end_indices_vec.n_elem);
    
    arma::rowvec tmp_unq_idxs = arma::unique(end_indices_vec);
    
    for (unsigned int i = 0; i < data.n_cols; i++) {
      
      arma::rowvec update_cost_vec = end_cost_vec;
      
      arma::rowvec update_indices_vec = end_indices_vec;
      
      bool tmp_flag = boolean_function(tmp_unq_idxs, i);
      
      if (!tmp_flag) {
        
        for (unsigned int j = 0; j < update_cost_vec.n_elem; j++) {
          
          if (update_cost_vec(j) != data(j, i)) {
            
            if (update_cost_vec(j) > data(j, i)) {
              
              update_cost_vec(j) = data(j, i);
              
              update_indices_vec(j) = i;
            }
          }
        }
      }
      
      double tmp_cost_val = arma::accu(update_cost_vec);
      
      if (tmp_cost_val < previous_cost) {
        
        previous_cost = tmp_cost_val;
        
        end_update_cost = update_cost_vec;
        
        end_update_idx = update_indices_vec;
        
        previous_idx = i;
      }
    }
    
    count_clusters += 1;
    
    end_cost_vec = end_update_cost;
    
    end_indices_vec = end_update_idx;
    
    if (verbose) { Rcpp::Rcout << "medoid " << previous_idx + 1 << " was added. Current dissimilarity of build phase  --> " << arma::accu(end_cost_vec) << std::endl; }
  }
  
  if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << " =================================================== " << std::endl; Rcpp::Rcout << " --- end of build phase ---> " << " dissimilarity: " << arma::accu(end_cost_vec) << std::endl; Rcpp::Rcout << " =================================================== " << std::endl;}
  
  Rcpp::List silh_lst;
  
  if (swap_phase) {
    
    arma::uvec swap_medoids = arma::conv_to< arma::uvec >::from(arma::unique(end_indices_vec));
    
    arma::uvec all_medoids = arma::regspace<arma::uvec>(0, 1, data.n_rows - 1);
    
    arma::uvec non_medoids = subset_vec(all_medoids, swap_medoids);
    
    double initial_cost = arma::accu(end_cost_vec);
    
    if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << "Initial dissimilarity of swap phase --> " << initial_cost << std::endl; }
    
    if (verbose && threads > 1) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << "NOTE: verbose is disabled in swap phase when threads > 1" << std::endl; }
    
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (unsigned int i = 0; i < non_medoids.n_elem; i++) {
      
      for (unsigned int j = 0; j < swap_medoids.n_elem; j++) {
        
        arma::uvec copy_medoids = swap_medoids;
        
        copy_medoids(j) = non_medoids(i);                // swap medoid with a non-medoid
        
        arma::mat new_cols = data.cols(copy_medoids);
        
        arma::rowvec tmp_min_cost(data.n_rows, arma::fill::zeros), tmp_clusters(data.n_rows, arma::fill::zeros);
        
        for (unsigned int i = 0; i < new_cols.n_rows; i++) {
          
          double tmp_cost = arma::min(new_cols.row(i));
          
          arma::uvec idx = arma::find(new_cols.row(i) == tmp_cost);
          
          tmp_min_cost(i) = tmp_cost;
          
          tmp_clusters(i) = copy_medoids(idx(0));
          
        }
        
        double new_cost = arma::accu(tmp_min_cost);
        
        if (new_cost < initial_cost) {
          
          initial_cost = new_cost;
          
          if (verbose && threads == 1) { Rcpp::Rcout << "swap of medoid " << swap_medoids(j) + 1 << " with the non-medoid " << non_medoids(i) + 1 << ". Current dissimilarity of the swap phase --> " << new_cost << std::endl; }
          
          swap_medoids(j) = non_medoids(i);
          
          end_cost_vec = tmp_min_cost;
          
          end_indices_vec = tmp_clusters;
        }
      }
    }
    
    if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << " ================================================== " << std::endl; Rcpp::Rcout << " --- end of swap phase ---> " << " dissimilarity: " << initial_cost << std::endl; Rcpp::Rcout << " ================================================== " << std::endl;}
  }
  
  arma::rowvec end_idxs = arma::unique(end_indices_vec);
  
  arma::rowvec true_idxs = arma::regspace<arma::rowvec>(0, 1, end_idxs.n_elem - 1);
  
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2)
  #endif
  for (unsigned int item = 0; item < end_indices_vec.n_elem; item++) {
    
    for (unsigned int item_1 = 0; item_1 < end_idxs.n_elem; item_1++) {
      
      if (end_indices_vec(item) == end_idxs(item_1)) {
        
        end_indices_vec(item) = true_idxs(item_1);
      }
    }
  }
  
  if (clusters > 1) { silh_lst = silhouette_matrix(data, end_indices_vec + 1, end_cost_vec, threads); }
  
  arma::mat fuz_out;
  
  if (fuzzy) {
    
    arma::uvec fuz_idx = arma::conv_to< arma::uvec >::from(end_idxs);
    
    arma::mat fuzz_mat = data.cols(fuz_idx);
    
    fuz_out.set_size(fuzz_mat.n_rows, fuzz_mat.n_cols);
    
    for (unsigned int i = 0; i < fuzz_mat.n_rows; i++) {
      
      arma::rowvec tmp_row = arma::abs(arma::conv_to< arma::rowvec >::from(fuzz_mat.row(i)));
      
      tmp_row = arma::abs(tmp_row);
      
      tmp_row += 1.0e-6;
      
      arma::rowvec d = arma::accu(tmp_row) / tmp_row;
      
      fuz_out.row(i) = d / arma::accu(d);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("medoids") = end_idxs, Rcpp::Named("cost") = arma::accu(end_cost_vec), Rcpp::Named("dissimilarity_matrix") = data, 
                            
                            Rcpp::Named("clusters") = end_indices_vec, Rcpp::Named("end_cost_vec") = end_cost_vec, Rcpp::Named("silhouette_matrix") = silh_lst[0], 
                                        
                                        Rcpp::Named("fuzzy_probs") = fuz_out, Rcpp::Named("clustering_stats") = silh_lst[1], Rcpp::Named("flag_dissim_mat") = flag_dissim_mat);
}


// calculate global dissimilarities for claraMedoids
//

// [[Rcpp::export]]
arma::mat dissim_MEDOIDS(arma::mat& data, std::string method, arma::mat MEDOIDS, double minkowski_p = 1.0, int threads = 1, double eps = 1.0e-6) {
  
  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif
  
  arma::mat cov_mat;
  
  if (method == "mahalanobis") {
    
    cov_mat = arma::inv(arma::cov(data));
  }
  
  arma::mat mt(data.n_rows, MEDOIDS.n_rows);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (unsigned int j = 0; j < MEDOIDS.n_rows; j++) {
    
    arma::vec tmp_vec(data.n_rows);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (unsigned int i = 0; i < data.n_rows; i++) {
      
      double tmp_idx;
      
      if (method == "euclidean") {
        
        tmp_idx = std::sqrt(arma::as_scalar(arma::accu(arma::square((data.row(i) - MEDOIDS.row(j))))));
      }
      
      else if (method == "manhattan") {
        
        tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - MEDOIDS.row(j)))));
      }
      
      else if (method == "chebyshev") {
        
        tmp_idx = arma::as_scalar(arma::max(arma::abs((data.row(i) - MEDOIDS.row(j)))));
      }
      
      else if (method == "canberra") {
        
        tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - MEDOIDS.row(j)) + eps)/(arma::abs(data.row(i)) + arma::abs(MEDOIDS.row(j)) + eps)));                 // added 1.0e-6 otherwise rstudio crashes
      }
      
      else if (method == "braycurtis") {
        
        tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - MEDOIDS.row(j))))/(arma::accu(arma::abs(data.row(i))) + arma::accu(arma::abs(MEDOIDS.row(j)))));
      }
      
      else if (method == "pearson_correlation") {
        
        tmp_idx = arma::as_scalar(1.0 - arma::cor(data.row(i), MEDOIDS.row(j)));
      }
      
      else if (method == "simple_matching_coefficient") {                                                    // for binary data
        
        double a = eps;
        double d = eps;
        
        for (unsigned int t = 0; t < data.row(i).n_elem; t++) {
          
          if (data.row(i)(t) == 1 && MEDOIDS.row(j)(t) == 1) {
            
            a += 1.0;}
          
          if (data.row(i)(t) == 0 && MEDOIDS.row(j)(t) == 0) {
            
            d += 1.0;
          }
        }
        
        tmp_idx = 1.0 - ((a + d) / data.row(i).n_elem);
      }
      
      else if (method == "minkowski") {                                                                                     // by default the order of the minkowski parameter equals k
        
        tmp_idx = std::pow(arma::as_scalar(arma::accu(arma::pow(arma::abs((data.row(i) - MEDOIDS.row(j))), minkowski_p))), 1.0 / minkowski_p);
      }
      
      else if (method == "hamming") {                                                                                     // for binary data
        
        tmp_idx = arma::as_scalar(arma::accu(data.row(i) != MEDOIDS.row(j))/(data.row(i).n_elem * 1.0));
      }
      
      else if (method == "mahalanobis") {                                                                                     // first create covariance matrix from data
        
        tmp_idx = arma::as_scalar(std::sqrt(arma::as_scalar(((data.row(i) - MEDOIDS.row(j)) * cov_mat) * (data.row(i) - MEDOIDS.row(j)).t())));
      }
      
      else if (method == "jaccard_coefficient") {                                                                                     // for binary data
        
        double a = eps;
        double b = eps;
        double c = eps;
        
        for (unsigned int t = 0; t < data.row(i).n_elem; t++) {
          
          if (data.row(i)(t) == 1 && MEDOIDS.row(j)(t) == 1) {
            
            a += 1.0;}
          
          if (data.row(i)(t) == 1 && MEDOIDS.row(j)(t) == 0) {
            
            b += 1.0;}
          
          if (data.row(i)(t) == 0 && MEDOIDS.row(j)(t) == 1) {
            
            c += 1.0;
          }
        }
        
        tmp_idx = 1.0 - (a / (a + b + c));
      }
      
      else if (method == "Rao_coefficient") {                                                                                     // for binary data
        
        double a = eps;
        
        for (unsigned int t = 0; t < data.row(i).n_elem; t++) {
          
          if (data.row(i)(t) == 1 && MEDOIDS.row(j)(t) == 1) {
            
            a += 1.0;
          }
        }
        
        tmp_idx = 1.0 - (a / data.row(i).n_elem);
      }
      
      else {
        
        tmp_idx = 0;                                             // default = 0; create exceptions in R, so that tmp_idx is never 0;
      }
      
      if ( tmp_idx != tmp_idx ) {                                // handling of NAs, if NaN then distance 1.0 [  NaN will compare false to everything, including itself ], http://stackoverflow.com/questions/11569337/using-an-if-statement-to-switch-nan-values-in-an-array-to-0-0]
        
        tmp_idx = 1.0;
      }
      
      tmp_vec(i) = tmp_idx;
    }
    
    mt.col(j) = tmp_vec;
  }
  
  return(mt);
}



// calculate fuzzy (soft) clusters and stats [ data is the subset of the dissimilarity matrix using the medoids as subset-indices ]
//

// [[Rcpp::export]]
Rcpp::List fuzzy_and_stats(arma::mat data, int threads = 1, double eps = 1.0e-6, bool fuzzy = false) {
  
  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif
  
  arma::mat fuz_out, stats_mat(5, data.n_cols, arma::fill::zeros);;
  
  arma::rowvec hard_out(data.n_rows, arma::fill::zeros);
  
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (unsigned int i = 0; i < data.n_rows; i++) {
    
    arma::rowvec tmp_row = arma::abs(arma::conv_to< arma::rowvec >::from(data.row(i)));
    
    arma::uvec tmp_hard = arma::find(tmp_row == arma::min(tmp_row));
    
    int idx = tmp_hard(0);
    
    stats_mat(1,idx) += 1;                                                                // 2nd column, number of observations in each cluster
    
    if (stats_mat(2,idx) < tmp_row(idx)) { stats_mat(2,idx) = tmp_row(idx); }             // 3rd column, maximum-intra-cluster-dissimilarity
    
    stats_mat(3,idx) += tmp_row(idx);                                                     
    
    hard_out(i) = idx;                                                                    // hard-labeling
    
    if (fuzzy) {
      
      fuz_out.set_size(data.n_rows, data.n_cols);
      
      tmp_row = arma::abs(tmp_row);
      
      tmp_row += eps;
      
      arma::rowvec d = arma::accu(tmp_row) / tmp_row;
      
      fuz_out.row(i) = d / arma::accu(d);
    }
  }
  
  stats_mat.row(3) /= stats_mat.row(1);                                                    // 4th column, average-intra-cluster-dissimilarity
  
  return Rcpp::List::create(Rcpp::Named("hard_labels") = hard_out, Rcpp::Named("fuzzy_labels") = fuz_out, Rcpp::Named("stats_mat") = stats_mat);
}



// isolation statistic for clara-medoids
// dissim_mat_subset is a subset of the dissimilarity matrix [ subset column-wise ]
//

// [[Rcpp::export]]
arma::rowvec isolation(arma::mat dissim_mat_subset, arma::uvec x) {
  
  arma::mat tmp_diss = dissim_mat_subset.rows(x);
  
  tmp_diss.diag().fill(arma::datum::inf);
  
  arma::rowvec max_diss(tmp_diss.n_cols, arma::fill::zeros);
  
  for (unsigned int i = 0; i < tmp_diss.n_cols; i++) {
    
    max_diss(i) = arma::min(tmp_diss.col(i));
  }
  
  return max_diss;
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
Rcpp::List ClaraMedoids(arma::mat& data, unsigned int clusters, std::string method, int samples, double sample_size, double minkowski_p = 1.0,
                        
                        int threads = 1, bool verbose = false, bool swap_phase = false, bool fuzzy = false, int seed = 1) {
  
  set_seed(seed);             // R's RNG
  
  int bst_sample = -1;
  
  double dism = arma::datum::inf;
  
  arma::mat probs, diss_glob_out;
  
  arma::uvec out_medoid, clr_split_out;
  
  Rcpp::List bst_lst;
  
  if (verbose) { Rcpp::Rcout << " " << std::endl; }
  
  for (int s = 0; s < samples; s++) {
    
    int samp_rows = std::ceil(data.n_rows * sample_size);
    
    arma::uvec clr_split = arma::conv_to< arma::uvec >::from(sample_vec(samp_rows, 0, data.n_rows - 1, false));
    
    if (s > 0) {                                                                                                            // include the bst medoids so far to the sample data set
      
      arma::uvec new_idx(clr_split.n_elem + out_medoid.n_elem);
      
      for (unsigned int i = 0; i < clr_split.n_elem; i++) {
        
        new_idx(i) = clr_split(i);
      }
      
      int count_clr = 0;
      
      for (unsigned int j = clr_split.n_elem; j < new_idx.n_elem; j++) {
        
        new_idx(j) = out_medoid(count_clr);
        
        count_clr += 1;
      }
      
      clr_split = arma::conv_to< arma::uvec >::from(arma::unique(new_idx));
    }
    
    arma::mat tmp_dat = data.rows(clr_split);                                                                                // use sample data for ClusterMedoids
    
    arma::mat copy_dat = tmp_dat;
    
    Rcpp::List clM_sblist = ClusterMedoids(tmp_dat, clusters, method, minkowski_p, threads, false, swap_phase, false);
    
    double local_dissim = Rcpp::as<double> (clM_sblist[1]);
    
    arma::uvec local_medoids = Rcpp::as<arma::uvec> (clM_sblist[0]);
    
    arma::mat tmp_glob = dissim_MEDOIDS(data, method, copy_dat.rows(local_medoids), minkowski_p , threads, 1.0e-6);          // use all data to calculate global dissimilarity   
    
    double global_dissimil = 0.0;
    
    for (unsigned int d = 0; d < tmp_glob.n_rows; d++) {
      
      global_dissimil += arma::min(tmp_glob.row(d));
    }
    
    if (verbose) { Rcpp::Rcout << "sample " << s + 1 << " -->  local (sample) dissimilarity : " << local_dissim << " -->  global dissimilarity : " << global_dissimil << std::endl; }
    
    if (global_dissimil < dism) {
      
      dism = global_dissimil;
      
      out_medoid = clr_split(local_medoids);
      
      clr_split_out = clr_split;
      
      bst_lst = clM_sblist;
      
      diss_glob_out = tmp_glob;
      
      bst_sample = s + 1;
    }
  }
  
  if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << " ====================================================== " << std::endl; Rcpp::Rcout << " --- best sample: " << bst_sample << " ---> " << " global dissimilarity: " << dism << std::endl; Rcpp::Rcout << " ====================================================== " << std::endl;}
  
  Rcpp::List fuz_and_stats = fuzzy_and_stats(diss_glob_out, threads, 1.0e-6, fuzzy);
  
  arma::rowvec hard_clust = Rcpp::as<arma::rowvec> (fuz_and_stats[0]);
  
  arma::rowvec isol_vec = isolation(diss_glob_out, out_medoid);
  
  arma::mat fuz_st_mat = Rcpp::as<arma::mat> (fuz_and_stats[2]);
  
  fuz_st_mat.row(4) = fuz_st_mat.row(2) / isol_vec;                  // divide the maximum intra-cluster-dissimilarity with the minimum distance between medoids to get isolation
  
  fuz_st_mat.row(0) = arma::regspace<arma::rowvec> (1, 1, out_medoid.n_elem);
  
  arma::mat subs_meds = data.rows(out_medoid);
  
  arma::mat bst_sample_dissm_mat = Rcpp::as<arma::mat> (bst_lst[2]);
  
  arma::mat bst_sample_silh_mat;
  
  if (clusters > 1) {
    
    bst_sample_silh_mat = Rcpp::as<arma::mat> (bst_lst[5]);
  }
  
  return Rcpp::List::create(Rcpp::Named("medoids") = subs_meds, Rcpp::Named("bst_dissimilarity") = dism, Rcpp::Named("medoid_indices") = out_medoid, 
                                        
                                        Rcpp::Named("sample_indices") = arma::conv_to< arma::rowvec >::from(clr_split_out), Rcpp::Named("clusters") = hard_clust, 
                                        
                                        Rcpp::Named("bst_sample_silhouette_matrix") = bst_sample_silh_mat, Rcpp::Named("fuzzy_probs") = fuz_and_stats[1], 
                                                                                                                                                     
                                                                                                                                                     Rcpp::Named("clustering_stats") = fuz_st_mat.t(), Rcpp::Named("bst_sample_dissimilarity_matrix") = bst_sample_dissm_mat);
}



// prediction function for the k-medoids
//

// [[Rcpp::export]]
Rcpp::List predict_medoids(arma::mat& data, std::string method, arma::mat MEDOIDS, double minkowski_p = 1.0, int threads = 1, bool fuzzy = false, double eps = 1.0e-6) {
  
  arma::mat tmp_dist = dissim_MEDOIDS(data, method, MEDOIDS, minkowski_p, threads, eps);
  
  arma::mat fuz_out;
  
  arma::rowvec hard_clust(tmp_dist.n_rows);
  
  double global_dissimil = 0.0;
  
  for (unsigned int i = 0; i < tmp_dist.n_rows; i++) {
    
    arma::rowvec tmp_row = arma::abs(arma::conv_to< arma::rowvec >::from(tmp_dist.row(i)));
    
    arma::uvec tmp_hard = arma::find(tmp_row == arma::min(tmp_row));
    
    int idx = tmp_hard(0);
    
    hard_clust(i) = idx;
    
    global_dissimil += tmp_dist(i,idx);
    
    if (fuzzy) {
      
      fuz_out.set_size(tmp_dist.n_rows, tmp_dist.n_cols);
      
      tmp_row = arma::abs(tmp_row);
      
      tmp_row += eps;
      
      arma::rowvec d = arma::accu(tmp_row) / tmp_row;
      
      fuz_out.row(i) = d / arma::accu(d);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("clusters") = hard_clust, Rcpp::Named("fuzzy_clusters") = fuz_out, Rcpp::Named("dissimilarity") = global_dissimil);
}




// function which takes an Rcpp List object of matrices from Cluster_Medoids OR Clara_Medoids and returns 
// those matrices split on the specific clusters


// [[Rcpp::export]]
Rcpp::List split_rcpp_lst(Rcpp::List lst) {
  
  arma::mat silh_mat = Rcpp::as<arma::mat> (lst[5]);
  
  arma::vec tmp_clust = arma::conv_to< arma::vec >::from(silh_mat.col(0));
  
  Rcpp::List idx_lst = cluster_indices(tmp_clust);
  
  Rcpp::List intra_dissm(idx_lst.size()), silhouet_lst(idx_lst.size());
  
  arma::rowvec intra_clust_disml(idx_lst.size(), arma::fill::zeros), silhouet(idx_lst.size(), arma::fill::zeros);
  
  for (int i = 0; i < idx_lst.size(); i++) {
    
    arma::uvec tmp_mat_idx = Rcpp::as<arma::uvec> (idx_lst[i]);
    
    arma::mat subs_mat = silh_mat.rows(tmp_mat_idx);
    
    arma::rowvec tmp_col_dism = arma::conv_to< arma::rowvec >::from(subs_mat.col(2));
    
    arma::rowvec tmp_col_silh = arma::conv_to< arma::rowvec >::from(subs_mat.col(4));
    
    intra_clust_disml(i) = arma::sum(tmp_col_dism);
    
    silhouet(i) = arma::sum(tmp_col_silh);
    
    intra_dissm[i] = tmp_col_dism;
    
    silhouet_lst[i] = tmp_col_silh;
  }
  
  double avg_intr_dis = arma::accu(intra_clust_disml);
  
  double avg_width_silh = arma::accu(silhouet);
  
  avg_intr_dis /= silh_mat.n_rows;
  
  avg_width_silh /= silh_mat.n_rows;
  
  return Rcpp::List::create(Rcpp::Named("avg_intra_clust_dissimilarity") = avg_intr_dis, 
                            
                            Rcpp::Named("sum_intra_dissim") = arma::accu(intra_clust_disml), Rcpp::Named("avg_width_silhouette") = avg_width_silh, 
                            
                            Rcpp::Named("list_intra_dissm") = intra_dissm, Rcpp::Named("list_silhouette") = silhouet_lst, Rcpp::Named("silhouette_plot") = true);
}




// optimal number of clusters for the k-medoids
//

// [[Rcpp::export]]
Rcpp::List OptClust(arma::mat& data, unsigned int iter_clust, std::string method, bool clara = false, int samples = 5, double sample_size = 0.001, double minkowski_p = 1.0, 
                    
                    std::string criterion = "dissimilarity", int threads = 1, bool swap_phase = false, bool verbose = false, int seed = 1) {
  
  set_seed(seed);             // R's RNG
  
  Rcpp::List medoids_object(iter_clust);
  
  if (verbose) { Rcpp::Rcout << " " << std::endl; }
  
  for (unsigned int iter = 0; iter < iter_clust; iter++) {
    
    if (iter == 0) { 
      
      if (verbose) { Rcpp::Rcout << "number of clusters: "<< iter + 1 << "  -->  " << criterion << ": " <<  arma::datum::inf << std::endl; }}
    
    else {
      
      if (!clara) {
        
        Rcpp::List cm_out = ClusterMedoids(data, iter + 1, method, minkowski_p, threads, false, swap_phase, false);
        
        medoids_object[iter] = split_rcpp_lst(cm_out);
        
        double dissiml = Rcpp::as<double> (cm_out[1]);
        
        if (verbose) { Rcpp::Rcout << "number of clusters: "<< iter + 1 << "  -->  " << criterion << ": " << dissiml << std::endl; }
      }
      
      else {
        
        Rcpp::List cl_out = ClaraMedoids(data, iter + 1, method, samples, sample_size, minkowski_p, threads, false, swap_phase,false);
        
        medoids_object[iter] = split_rcpp_lst(cl_out);
        
        double dissiml = Rcpp::as<double> (cl_out[1]);
        
        if (verbose) { Rcpp::Rcout << "number of clusters: "<< iter + 1 << "  -->  " << criterion << ": " <<  dissiml << std::endl; }
      }
    }
  }
  
  return medoids_object;
}


