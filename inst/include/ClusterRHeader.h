
#ifndef __ClusterRHeader__
#define __ClusterRHeader__

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>


namespace clustR {

  class ClustHeader {

    public:

      ClustHeader() { }

      //********************************************************************************************************************************************************************************************************
      //********************************************************************************************************************************************************************************************************

      //  the following correspond to "utils_rcpp.cpp" non-header ClusterR package functions  [ former ClusterR version ]

      //********************************************************************************************************************************************************************************************************
      //********************************************************************************************************************************************************************************************************


      // secondary function which takes a vector with the clusters and returns a list with the cluster indices
      //

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

      bool check_NaN_Inf(arma::mat x) {

        return x.is_finite();
      }


      // silhouette formula for the intra-, outer-cluster-dissimilarities
      // https://en.wikipedia.org/wiki/Silhouette_(clustering)
      //

      double calc_silhouette(double intra, double outer) {

        if (outer > intra) {

          return 1.0 - (intra / outer);}

        else if (outer < intra) {

          return (outer / intra) - 1.0;}

        else {                             // outer == intra

          return 0.0;
        }
      }



      // use base R's set.seed() in Rcpp for RNG
      // http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case
      //

      void set_seed(int seed) {

        if (Rcpp::all(Rcpp::is_na(Rcpp::IntegerVector{seed}))) return;
        Rcpp::Environment base_env("package:base");
        Rcpp::Function set_seed_r = base_env["set.seed"];
        set_seed_r(seed);
      }



      // sample with or without replacement [ 'num_elem' is an integer not a vector -- create a similar function for a vector ]
      //

      arma::rowvec sample_vec(int num_elem, int start, int end, bool replace) {

        if (replace) {
          return(arma::randi<arma::rowvec>( num_elem, arma::distr_param(start,end) ));
        }
        else {
          arma::rowvec tmp = arma::conv_to< arma::rowvec >::from(arma::shuffle(arma::regspace(start,end)));
          return(tmp.subvec(0,num_elem-1));
        }
      }



      // early-stopping criterium

      double squared_norm(arma::mat x) {
        return std::sqrt(arma::accu(arma::pow(x, 2)));
      }



      // this function gives the index of the minimum value in a vector
      //

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

      Rcpp::List validate_centroids(arma::mat& data,
                                    arma::mat init_centroids,
                                    int threads = 1,
                                    bool fuzzy = false,
                                    double eps = 1.0e-6) {

        #ifdef _OPENMP
        omp_set_num_threads(threads);
        #endif

        arma::rowvec tmp_idx(data.n_rows);
        arma::mat soft_CLUSTERS;

        if (fuzzy) {
          soft_CLUSTERS.set_size(data.n_rows, init_centroids.n_rows);
        }

        unsigned int k,j;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) shared(data, init_centroids, tmp_idx, soft_CLUSTERS, fuzzy) private(k,j)
        #endif
        for (k = 0; k < data.n_rows; k++) {

          arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data.row(k)), init_centroids);
          double iter_val = MinMat(tmp_vec);

          if (fuzzy) {
            for (j = 0; j < tmp_vec.n_elem; j++) {

              #ifdef _OPENMP
              #pragma omp atomic write
              #endif
              soft_CLUSTERS(k,j) = tmp_vec(j);
            }
          }

          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          tmp_idx(k) = iter_val;
        }

        if (fuzzy) {
          for (unsigned int i = 0; i < soft_CLUSTERS.n_rows; i++) {
            soft_CLUSTERS.row(i) = norm_fuzzy(arma::conv_to< arma::rowvec >::from(soft_CLUSTERS.row(i)), eps);
          }
        }

        return Rcpp::List::create(Rcpp::Named("clusters") = tmp_idx,
                                  Rcpp::Named("fuzzy_probs") = soft_CLUSTERS);
      }



      // kmeans++ distance
      //

      double kmeans_pp_dist(arma::rowvec vec, arma::rowvec centroid) {
        return arma::as_scalar(arma::accu(arma::pow(vec - centroid, 2)));
      }



      // kmeans++ algorithm to build initialization-centroids [ it works for 1:data.n_row-2 ]
      //
      // http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work
      // https://datasciencelab.wordpress.com/2014/01/15/improved-seeding-for-clustering-with-k-means/
      // https://en.wikipedia.org/wiki/K-means%2B%2B
      //

      arma::mat kmeans_pp_init(arma::mat& data, int clusters, bool medoids = false) {

        arma::rowvec centroid_data(data.n_rows);                                                 // begin with inf values as I'll pick in every iteration the minimum
        centroid_data.fill(arma::datum::inf);                                                    // initialize a vector of Inf [ max. values ]
        arma::mat centroids(clusters, data.n_cols);                                              // save the end-centroids
        arma::mat medoids_vec(1, clusters);                                                      // in case of medoids
        arma::rowvec first = sample_vec(1, 0, data.n_rows - 1, false);                           // choose the first observation at random
        arma::vec i_out = arma::regspace(0, data.n_rows - 1);

        int idx = first(0);                                                                      // update idx in every iteration [ see line 'idx = j' line ]

        for (int clust = 0; clust < clusters; clust++) {

          arma::rowvec choose_row = data.row(idx);
          arma::rowvec tmp_dist(data.n_rows);

          for (unsigned int i = 0; i < data.n_rows; i++) {                                       // iterate over the data-rows and calculate the distance between centroid and data

            double tmp_val = kmeans_pp_dist(arma::conv_to< arma::rowvec >::from(data.row(i)), choose_row);
            arma::vec tmp_cent = { centroid_data(i), tmp_val };                                  // then pick the minimum between the current and the previous iters values
            tmp_dist(i) = arma::min(tmp_cent);
          }

          centroid_data = tmp_dist;                                                              // at every iteration I overwrite the 'centroid_data' with the 'tmp_dist' values
          arma::rowvec prob = tmp_dist / arma::sum(tmp_dist);
          arma::rowvec cum_sum = arma::cumsum(prob);                                             // calculate the cumulative probability
          // cum_sum = arma::unique(cum_sum);                                                    // in my previous code version of the 'kmeans_pp_init()' function I continued with the unique values, however this can reduce the size and give an error in the next lines (I expect it to be the same size as the 'i_out_subs' vector). It doesn't make any difference because I iterate over the 'cum_sum_subs' vector and see if '(cum_sum_subs(j) > r)', thus duplicated values will be just skipped
          arma::rowvec rr = arma::conv_to< arma::rowvec >::from(arma::randu(cum_sum.n_elem));

          double r = rr(0);                                                                      // pick a random value using the randu() function of armadillo

          arma::vec i_out_subs(i_out);
          arma::vec cum_sum_subs(cum_sum.n_elem);

          for (unsigned int k = 0; k < cum_sum.n_elem; k++) {
            cum_sum_subs(k) = cum_sum(k);
          }

          arma::mat data_subs(data);

          if (clust > 0) {

            arma::uvec subs_idx = arma::conv_to< arma::uvec >::from(medoids_vec.row(0));           // use the medoids to subset the cum_sum and the data before picking a new centroid
            subs_idx = subs_idx.subvec(0, clust);                                                  // keep until the current iterations centroid-index
            i_out_subs.shed_rows(subs_idx);
            cum_sum_subs.shed_rows(subs_idx);
            data_subs.shed_rows(subs_idx);
          }

          for (unsigned int j = 0; j < cum_sum_subs.n_elem; j++) {

            if (cum_sum_subs(j) > r) {

              idx = i_out_subs(j);                                                                    // update (overwrite) the 'idx'
              centroids.row(clust) = arma::conv_to< arma::rowvec >::from(data_subs.row(j));           // add the centroids
              medoids_vec.col(clust) = i_out_subs(j);                                                 // the output 'medoids' are also the output row-data-indices  [ they return the index of the centroids ]

              break;
            }
          }
        }

        if (medoids) {
          return medoids_vec;
        }
        else {
          return centroids;
        }
      }


      //..................................................................................................... !! Previous version of the function which returned duplicated centroids [ see: https://github.com/mlampros/ClusterR/issues/25 ]
      // arma::mat kmeans_pp_init(arma::mat& data, int clusters, bool medoids = false) {
      //
      //   arma::rowvec centroid_data(data.n_rows);                                                 // begin with inf values as I'll pick in every iteration the minimum
      //
      //   centroid_data.fill(arma::datum::inf);
      //
      //   arma::mat centroids(clusters, data.n_cols);                                              // save the end-centroids
      //
      //   arma::mat medoids_vec(1, clusters);                                                      // in case of medoids
      //
      //   arma::rowvec first = sample_vec(1, 0, data.n_rows - 1, false);                           // choose the first observation at random
      //
      //   arma::rowvec indices(clusters);
      //
      //   int idx = first(0);                                                                      // update idx in every iteration
      //
      //   for (int clust = 0; clust < clusters; clust++) {
      //
      //     indices(clust) = idx;
      //
      //     arma::rowvec choose_row = data.row(idx);
      //
      //     arma::rowvec tmp_dist(data.n_rows);
      //
      //     for (unsigned int i = 0; i < data.n_rows; i++) {                                                // iterate over the data-rows and calculate the distance between centroid and data
      //
      //       double tmp_val = kmeans_pp_dist(arma::conv_to< arma::rowvec >::from(data.row(i)), choose_row);
      //
      //       arma::vec tmp_cent = { centroid_data(i), tmp_val };                                  // then pick the minimum between the current and the previous iters values
      //
      //       tmp_dist(i) = arma::min(tmp_cent);
      //     }
      //
      //     centroid_data = tmp_dist;
      //
      //     arma::rowvec prob = tmp_dist / arma::sum(tmp_dist);
      //
      //     arma::rowvec cum_sum = arma::cumsum(prob);                                             // calculate the cumulative probability
      //
      //     cum_sum = arma::unique(cum_sum);
      //
      //     arma::rowvec rr = arma::conv_to< arma::rowvec >::from(arma::randu(cum_sum.n_elem));
      //
      //     double r = rr(0);                                                                      // pick a random value using the randu() function of armadillo
      //
      //     for (unsigned int j = 0; j < cum_sum.n_elem; j++) {
      //
      //       if (cum_sum(j) > r) {
      //
      //         idx = j;                                                                           // update idx
      //
      //         centroids.row(clust) = arma::conv_to< arma::rowvec >::from(data.row(j));           // add the centroids
      //
      //         if (medoids) { medoids_vec.col(clust) = j; }
      //
      //         break;
      //       }
      //     }
      //   }
      //
      //   if (medoids) {
      //
      //     return medoids_vec;}
      //
      //   else {
      //
      //     return centroids;
      //   }
      // }
      //.....................................................................................................


      // calculate fuzzy (soft) clusters
      //

      arma::rowvec norm_fuzzy(arma::rowvec vec, double eps) {

        vec = arma::abs(vec);
        vec += eps;
        arma::rowvec d = arma::accu(vec) / vec;
        return d / arma::accu(d);
      }




      // secondary function for the 'quantile_init_rcpp' function
      //

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
            idx_out(j) = idx1_sorted(0);
          }
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


      //********************************************************************************************************************************************************************************************************
      //********************************************************************************************************************************************************************************************************

      //  the following correspond to "kmeans_miniBatchKmeans_GMM_Medoids.cpp" non-header ClusterR package functions  [ former ClusterR version ]

      //********************************************************************************************************************************************************************************************************
      //********************************************************************************************************************************************************************************************************


      //============================================ k-means ====================================================================================


      // center, square and take the sum of the data to calculate the total sum of squares
      // http://stackoverflow.com/questions/8637460/k-means-return-value-in-r
      //

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

      Rcpp::List KMEANS_rcpp(arma::mat& data,
                             int clusters,
                             int num_init = 1,
                             int max_iters = 200,
                             std::string initializer = "kmeans++",
                             bool fuzzy = false,
                             bool verbose = false,
                             Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue,
                             double tol = 1e-4,
                             double eps = 1.0e-6,
                             double tol_optimal_init = 0.5,
                             int seed = 1) {

        int dat_n_rows = data.n_rows;
        if (clusters > dat_n_rows - 2 || clusters < 1) {
          Rcpp::stop("the number of clusters should be at most equal to nrow(data) - 2 and never less than 1");
        }

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
        bst_WCSS.fill(arma::datum::inf);                                           // initialize WCSS to Inf, so that in first iteration it can be compared with the minimum of the 'WSSE'
        int end_init = 0;

        if (verbose) {
          Rcpp::Rcout << " " << std::endl;
        }

        arma::rowvec flag_exception(num_init, arma::fill::zeros);

        for (int init_count = 0; init_count < num_init; init_count++) {

          arma::mat conv;

          if (!flag) {

            if (initializer == "kmeans++") {
              conv = kmeans_pp_init(data, clusters, false);
            }
            if (initializer == "random") {
              arma::uvec samp = arma::conv_to< arma::uvec >::from(sample_vec(clusters, 0, data.n_rows - 1, false));
              conv = data.rows(samp);
            }
            if (initializer == "quantile_init") {
              int samp_ROWS = data.n_rows / 10;
              arma::uvec tmp_conv = quantile_init_rcpp(data, samp_ROWS, clusters);
              flag_exception(init_count) = duplicated_flag(tmp_conv);                           // print a warning in case of duplicated centroids
              conv = data.rows(tmp_conv);
            }
            if (initializer == "optimal_init") {
              arma::vec tmp_conv = check_medoids(data, clusters, tol_optimal_init);             // tolerance parameter 'tol' here by default equals to 0.5 [ this parameter is important in case of duplicated rows ]
              flag_exception(init_count) = !arma::is_finite(tmp_conv);                          // necessary in order to stop the function in case of UBSAN memory errors
              arma::uvec idx_out = arma::conv_to< arma::uvec >::from(tmp_conv);
              conv = data.rows(idx_out);
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

              arma::rowvec tmp_row = data.row(i);                                              // current row
              arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(tmp_row), conv);    // vector where each item corresponds to the WCSS betw. the row and the current centroids

              if (fuzzy) {
                soft_CLUSTERS.row(i) = arma::conv_to< arma::rowvec >::from(tmp_vec);
              }

              int tmp_idx = MinMat(tmp_vec);                         // returns the index of the tmp_vec with the lowest WSSE, that means this row (or observation) will be assigned to this cluster (centroid)
              total_WSSE(tmp_idx) += tmp_vec(tmp_idx);               // assigns to total_WSSE the minimum cost
              num_obs(tmp_idx) += 1;                                 // append this observation to the corresponding cluster
              new_centroids.row(tmp_idx) += tmp_row;                 // adds to the corresponding row of the tmp_clusters matrix the row of the data, so that the new centroids can be calculated
              CLUSTERS(i) = tmp_idx;
            }

            for (int j = 0; j < clusters; j++) {
              new_centroids.row(j) /= arma::as_scalar(num_obs(j));
            }

            double tmp_norm = squared_norm(conv - new_centroids);
            conv = new_centroids;

            if (verbose) {
              Rcpp::Rcout << "iteration: " << iter + 1 << " --> total WCSS: " << arma::accu(total_WSSE) << "  -->  squared norm: " << tmp_norm << std::endl;
            }

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

          double ACCUMUL_WCSS;
          ACCUMUL_WCSS = arma::as_scalar(arma::accu(bst_WCSS));

          if (arma::accu(SSE) < ACCUMUL_WCSS) {
            end_init = init_count + 1;
            bst_WCSS = SSE;
            bst_obs = OBS;
            centers_out = bst_centers;
            lst_out = CLUSTERS_OUT;

            if (fuzzy) {
              lst_fuzzy_out = fuzzy_OUT;
            }
          }

          if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << "===================== end of initialization " << init_count + 1 << " =====================" << std::endl; Rcpp::Rcout << " " << std::endl; }
        }

        int exc = arma::as_scalar(flag_exception(end_init - 1));                        // print warning OR stop function only if duplicates OR NA's are present in relevant output centroids [ end_init - 1 ]

        if (exc > 0) {

          if (initializer == "quantile_init") {
            std::string message = "the centroid matrix using 'quantile_init' as initializer contains duplicates for number of clusters equal to : " + std::to_string(clusters);
            Rcpp::warning(message);
          }

          if (initializer == "optimal_init") {
            std::string message = "The centroid matrix using 'optimal_init' as initializer contains NA's. Thus, the 'tol_optimal_init' parameter should be (probably) decreased for number of clusters equal to : " + std::to_string(clusters);
            Rcpp::stop(message);
          }
        }

        double tmp_sse = tot_ss_data(data);

        if (fuzzy) {

          arma::mat fuzzy_mat(lst_fuzzy_out.n_rows, lst_fuzzy_out.n_cols);

          for (unsigned int i = 0; i < lst_fuzzy_out.n_rows; i++) {
            fuzzy_mat.row(i) = norm_fuzzy(arma::conv_to< arma::rowvec >::from(lst_fuzzy_out.row(i)), eps);
          }

          return Rcpp::List::create(Rcpp::Named("clusters") = lst_out,
									Rcpp::Named("fuzzy_clusters") = fuzzy_mat,
									Rcpp::Named("centers") = centers_out,
									Rcpp::Named("total_SSE") = tmp_sse,
									Rcpp::Named("best_initialization") = end_init,
									Rcpp::Named("WCSS_per_cluster") = bst_WCSS,
									Rcpp::Named("obs_per_cluster") = bst_obs);
        }
        else {
          return Rcpp::List::create(Rcpp::Named("clusters") = lst_out,
                                    Rcpp::Named("centers") = centers_out,
                                    Rcpp::Named("total_SSE") = tmp_sse,
                                    Rcpp::Named("best_initialization") = end_init,
                                    Rcpp::Named("WCSS_per_cluster") = bst_WCSS,
                                    Rcpp::Named("obs_per_cluster") = bst_obs);
        }
      }



      // the KMEANS_arma returns only the centroids and if openmp exists then it uses all available threads
      // the number of columns of the data should be much larger than the number of clusters  (however it works for number_columns == clusters)
      // seed_mode is one of : "keep_existing" (I've to give the CENTROIDS), "static_subset", "random_subset", "static_spread", "random_spread"
      // in comparison to my 'KMEANS_rcpp' the 'kmeans_arma' does a single initialization and it suggests maximum 10 iterations
      //

      arma::mat KMEANS_arma(arma::mat& data, int clusters, int n_iter, bool verbose, std::string seed_mode = "random_subset",

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


      // 2-dimensional permutations [ returns a 2-column matrix ]
      //

      arma::mat PERMUTATIONS_2D(arma::uword UNQ_CLUSTERS) {
        arma::uword REPLICATE = UNQ_CLUSTERS - 1;
        arma::mat FIRST_COL(REPLICATE, 1);
        arma::mat SECOND_COL(UNQ_CLUSTERS, 1);
        for (arma::uword j = 0; j < UNQ_CLUSTERS; j++) {
          SECOND_COL(j, 0) = j;
        }

        arma::mat res_first_col, res_second_col;

        for (arma::uword i = 0; i < UNQ_CLUSTERS; i++) {
          FIRST_COL.fill(i);
          arma::mat SECOND_COL_ITER(SECOND_COL);                // copy the matrix in the current iteration
          SECOND_COL_ITER.shed_row(i);                          // remove row

          if (i == 0) {
            res_first_col = FIRST_COL;
            res_second_col = SECOND_COL_ITER;
          }
          else {
            res_first_col = arma::join_vert(res_first_col, FIRST_COL);
            res_second_col = arma::join_vert(res_second_col, SECOND_COL_ITER);
          }
        }

        arma::mat perm_mt = arma::join_horiz(res_first_col, res_second_col);
        return perm_mt;
      }


      // secondary function to calculate the silhouette metric
      //

      Rcpp::List SILHOUETTE_metric(arma::mat& data, arma::vec CLUSTER, Rcpp::List tmp_clust, Rcpp::List in_cluster_dist) {

        arma::vec unq_values = arma::unique(CLUSTER);
        arma::mat IDX = PERMUTATIONS_2D(unq_values.n_elem);
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

      Rcpp::List evaluation_rcpp(arma::mat& data, arma::vec CLUSTER, bool silhouette = false) {

        Rcpp::List tmp_clust = cluster_indices(CLUSTER);

        Rcpp::List in_cluster_dist = INTRA_CLUSTER_DISS(data, tmp_clust);

        if (!silhouette) {

          arma::rowvec befout_CLUSTER = arma::conv_to< arma::rowvec >::from(CLUSTER);

          return(Rcpp::List::create(Rcpp::Named("clusters") = befout_CLUSTER,

                                    Rcpp::Named("cluster_indices") = tmp_clust,                                    // use the data indices, otherwise difficult to match clusters with silhouette or dissimilarity coefficients

                                    Rcpp::Named("INTRA_cluster_dissimilarity") = in_cluster_dist));                // the lower the better
        }

        else {

          Rcpp::List silhouet_out = SILHOUETTE_metric(data, CLUSTER, tmp_clust, in_cluster_dist);

          arma::rowvec befout_CLUSTER = arma::conv_to< arma::rowvec >::from(CLUSTER);

          return(Rcpp::List::create(Rcpp::Named("clusters") = befout_CLUSTER,

                                    Rcpp::Named("cluster_indices") = tmp_clust,                                       // use the data indices, otherwise difficult to match clusters with silhouette or dissimilarity coefficients

                                    Rcpp::Named("INTRA_cluster_dissimilarity") = in_cluster_dist,                     // the lower the better

                                    Rcpp::Named("silhouette") = silhouet_out));                                       // the higher the better [ range between -1 and +1 ]
        }
      }


      // compute the silhouette width (including the clusters and intra cluster dissimilarity)
      // reference for the 'silhouette_global_average':
      //           https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html#sklearn.metrics.silhouette_score
      //

      Rcpp::List silhouette_clusters(arma::mat& data, arma::vec CLUSTER) {

        Rcpp::List obj_clust = evaluation_rcpp(data, CLUSTER, true);

        arma::rowvec clust_vec = Rcpp::as<arma::rowvec>(obj_clust["clusters"]);
        arma::rowvec unq_items = arma::unique(clust_vec) - 1;            // adjust indexing to C++ (the unique clusters are sorted in ascending order 0,1,2 etc.)
        arma::uword num_clusts = unq_items.n_elem;

        Rcpp::List clust_idx = Rcpp::as<Rcpp::List>(obj_clust["cluster_indices"]);
        Rcpp::List clust_disim = Rcpp::as<Rcpp::List>(obj_clust["INTRA_cluster_dissimilarity"]);
        Rcpp::List clust_silh = Rcpp::as<Rcpp::List>(obj_clust["silhouette"]);

        Rcpp::NumericVector unq_clusters_order(num_clusts), clusters_size(num_clusts), avg_disim(num_clusts), avg_silh(num_clusts);
        arma::mat clust_mt, disim_mt, silh_mt;

        Rcpp::List clusters_lst(num_clusts);
        for (arma::uword i = 0; i < num_clusts; i++) {
          arma::uvec idx_iter_clust = Rcpp::as<arma::uvec>(clust_idx[i]);

          arma::mat dat_iter_clust = clust_vec(idx_iter_clust);                           // clusters (based on indices), returns a matrix (1-row-matrix)
          arma::uword unq_clust = arma::as_scalar(arma::unique(dat_iter_clust));

          arma::mat dat_iter_disim = Rcpp::as<arma::mat>(clust_disim[i]);                 // cluster-dissimilarities (based on indices), returns a matrix
          double disim_avg_iter = arma::mean(dat_iter_disim.row(0));

          arma::rowvec dat_iter_silh = Rcpp::as<arma::rowvec>(clust_silh[i]);             // cluster-silhouette (based on indices), returns a vector
          double silh_avg_iter = arma::mean(dat_iter_silh);

          clusters_lst[i] = dat_iter_clust;
          unq_clusters_order[i] = unq_clust;
          clusters_size[i] = idx_iter_clust.n_elem;
          avg_disim[i] = disim_avg_iter;
          avg_silh[i] = silh_avg_iter;

          // rbind() the matrices
          arma::mat silh_vec2mat = arma::conv_to<arma::mat>::from(dat_iter_silh);     // convert first current iteration's silhouette vector to matrix

          if (i == 0) {                                                               // overwrite in iteration 0, then rbind()
            clust_mt = dat_iter_clust;
            disim_mt = dat_iter_disim;
            silh_mt = silh_vec2mat;
          }
          else {
            clust_mt = arma::join_vert(clust_mt, dat_iter_clust);                       // clusters         ( !! join_vert() )
            disim_mt = arma::join_horiz(disim_mt, dat_iter_disim);                      // dissimilarities  ( !! join_horiz() )
            silh_mt = arma::join_horiz(silh_mt, silh_vec2mat);                          // silhouette       ( !! join_horiz() )
          }
        }

        Rcpp::DataFrame df_summary = Rcpp::DataFrame::create(Rcpp::Named("cluster") = unq_clusters_order,
                                                             Rcpp::Named("size") = clusters_size,
                                                             Rcpp::Named("avg_intra_dissim") = avg_disim,
                                                             Rcpp::Named("avg_silhouette") = avg_silh);

        arma::mat clust_disim_silh = arma::join_horiz(clust_mt, disim_mt.t(), silh_mt.t());                 // transpose the matrices that I joint horizontal to come to the matrix
        double global_mean = arma::mean(silh_mt.row(0));

        Rcpp::List dat_out = Rcpp::List::create(Rcpp::Named("silhouette_matrix") = clust_disim_silh,
                                                Rcpp::Named("silhouette_summary") = df_summary,
                                                Rcpp::Named("silhouette_global_average") = global_mean);
        return dat_out;
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

      Rcpp::List mini_batch_kmeans(arma::mat& data, int clusters, int batch_size, int max_iters, int num_init = 1, double init_fraction = 1.0, std::string initializer = "kmeans++",

                                   int early_stop_iter = 10, bool verbose = false, Rcpp::Nullable<Rcpp::NumericMatrix> CENTROIDS = R_NilValue, double tol = 1e-4, double tol_optimal_init = 0.5, int seed = 1) {

        set_seed(seed);             // R's RNG

        int dat_n_rows = data.n_rows;

        if (clusters > dat_n_rows - 2 || clusters < 1) { Rcpp::stop("the number of clusters should be at most equal to nrow(data) - 2 and never less than 1"); }

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

        arma::rowvec flag_exception(num_init, arma::fill::zeros);

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

              arma::uvec tmp_update_centroids = quantile_init_rcpp(data, samp_ROWS, clusters);

              flag_exception(init) = duplicated_flag(tmp_update_centroids);                           // print a warning in case of duplicated centroids

              update_centroids = data.rows(tmp_update_centroids);
            }

            if (initializer == "optimal_init") {

              arma::vec tmp_update_centroids = check_medoids(data, clusters, tol_optimal_init);       // tolerance parameter 'tol' here by default equals to 0.5 [ this parameter is important in case of duplicated rows ]

              flag_exception(init) = !arma::is_finite(tmp_update_centroids);                          // necessary in order to stop the function in case of UBSAN memory errors

              arma::uvec idx_out = arma::conv_to< arma::uvec >::from(tmp_update_centroids);

              update_centroids = data.rows(idx_out);
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

        int exc = arma::as_scalar(flag_exception(end_init - 1));                        // print warning OR stop function only if duplicates OR NA's are present in relevant output centroids [ end_init - 1 ]

        if (exc > 0) {

          if (initializer == "quantile_init") {

            std::string message = "the centroid matrix using 'quantile_init' as initializer contains duplicates for number of clusters equal to : " + std::to_string(clusters);

            Rcpp::warning(message);}

          if (initializer == "optimal_init") {

            std::string message = "The centroid matrix using 'optimal_init' as initializer contains NA's. Thus, the 'tol_optimal_init' parameter should be (probably) decreased for number of clusters equal to : " + std::to_string(clusters);

            Rcpp::stop(message);
          }
        }

        return Rcpp::List::create(Rcpp::Named("centroids") = centers_out, Rcpp::Named("WCSS_per_cluster") = bst_WCSS,

                                  Rcpp::Named("best_initialization") = end_init, Rcpp::Named("iters_per_initialization") = iter_before_stop);
      }




      // predict function for mini-batch-kmeans, which takes a matrix of centroids
      //

      Rcpp::List Predict_mini_batch_kmeans(arma::mat& data, arma::mat& CENTROIDS, bool fuzzy = false, double eps = 1.0e-6) {

        arma::rowvec CLUSTERS(data.n_rows);

        arma::mat soft_CLUSTERS(data.n_rows, CENTROIDS.n_rows);

        for (unsigned int j = 0; j < data.n_rows; j++) {

          arma::vec tmp_vec = WCSS(arma::conv_to< arma::rowvec >::from(data.row(j)), CENTROIDS);                  // returns a rowvec with the SSE for each cluster

          soft_CLUSTERS.row(j) = arma::conv_to< arma::rowvec >::from(tmp_vec);

          int tmp_idx = MinMat(tmp_vec);                                                                          // returns the index of the tmp_vec with the lowest SSE

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


      // secondary function used in 'GMM_arma()'
      //

      template<class T>
      T GMM_arma_covariance_type(T model,
                                 arma::mat& data,
                                 int gaussian_comps,
                                 std::string dist_mode,
                                 std::string seed_mode,
                                 int km_iter,
                                 int em_iter,
                                 bool verbose,
                                 double var_floor = 1e-10,
                                 int seed = 1) {
        bool status;

        if (seed_mode == "static_subset" && dist_mode == "eucl_dist") {
          status = model.learn(data.t(), gaussian_comps, arma::eucl_dist, arma::static_subset, km_iter, em_iter, var_floor, verbose);
        }
        else if (seed_mode == "random_subset" && dist_mode == "eucl_dist") {
          status = model.learn(data.t(), gaussian_comps, arma::eucl_dist, arma::random_subset, km_iter, em_iter, var_floor, verbose);
        }
        else if (seed_mode == "static_spread" && dist_mode == "eucl_dist") {
          status = model.learn(data.t(), gaussian_comps, arma::eucl_dist, arma::static_spread, km_iter, em_iter, var_floor, verbose);
        }
        else if (seed_mode == "random_spread" && dist_mode == "eucl_dist") {
          status = model.learn(data.t(), gaussian_comps, arma::eucl_dist, arma::random_spread, km_iter, em_iter, var_floor, verbose);
        }
        else if (seed_mode == "static_subset" && dist_mode == "maha_dist") {
          status = model.learn(data.t(), gaussian_comps, arma::maha_dist, arma::static_subset, km_iter, em_iter, var_floor, verbose);
        }
        else if (seed_mode == "random_subset" && dist_mode == "maha_dist") {
          status = model.learn(data.t(), gaussian_comps, arma::maha_dist, arma::random_subset, km_iter, em_iter, var_floor, verbose);
        }
        else if (seed_mode == "static_spread" && dist_mode == "maha_dist") {
          status = model.learn(data.t(), gaussian_comps, arma::maha_dist, arma::static_spread, km_iter, em_iter, var_floor, verbose);
        }
        else if (seed_mode == "random_spread" && dist_mode == "maha_dist") {
          status = model.learn(data.t(), gaussian_comps, arma::maha_dist, arma::random_spread, km_iter, em_iter, var_floor, verbose);
        }
        else {
          Rcpp::stop("Invalid seed_mode OR dist_mode. Valid 'seed_modes' are : 'static_subset', 'random_subset', 'static_spread' and 'random_spread'. Valid 'dist_modes' are : 'eucl_dist' and 'maha_dist'.");
        }

        // if(status == false) {
        //   Rcpp::Rcout << "learning failed" << std::endl;
        // }

        return model;
      }


      // the GMM_arma returns the centroids, covariance matrix and the weights. If openmp exists then it uses all available threads
      // distance is one of : "eucl_dist", "maha_dist"
      // seed_mode is one of : "static_subset", "random_subset", "static_spread", "random_spread"  [ I excluded the "keep_existing" seed_mode ] -- user-defined parameter setting is not enabled
      //

      // [[Rcpp::export]]
      Rcpp::List GMM_arma(arma::mat& data,
                          int gaussian_comps,
                          std::string dist_mode,
                          std::string seed_mode,
                          int km_iter,
                          int em_iter,
                          bool verbose,
                          double var_floor = 1e-10,
                          int seed = 1,
                          bool full_covariance_matrices = false) {

        arma::wall_clock timer;
        timer.tic();
        set_seed(seed);             // R's RNG
        arma::mat means;

        if (full_covariance_matrices) {
          arma::gmm_full model_inner;
          model_inner = GMM_arma_covariance_type<arma::gmm_full>(model_inner,
                                                                 data,
                                                                 gaussian_comps,
                                                                 dist_mode,
                                                                 seed_mode,
                                                                 km_iter,
                                                                 em_iter,
                                                                 verbose,
                                                                 var_floor,
                                                                 seed);

          arma::mat loglik(data.n_rows, gaussian_comps, arma::fill::zeros);

          for (int j = 0; j < gaussian_comps; j++) {
            loglik.col(j) = arma::conv_to< arma::vec >::from(model_inner.log_p(data.t(), j));
          }

          arma::mat model_means = model_inner.means.t();
          arma::cube model_COVS = model_inner.fcovs;

          arma::rowvec model_hefts = arma::conv_to< arma::rowvec >::from(model_inner.hefts.t());
          double model_avg_log_p = model_inner.avg_log_p(data.t(), gaussian_comps - 1);

          double n = timer.toc();
          if (verbose) { Rcpp::Rcout << "\ntime to complete : " << n << "\n" << std::endl; }

          return Rcpp::List::create( Rcpp::Named("centroids") = model_means,
                                     Rcpp::Named("covariance_matrices") = model_COVS,
                                     Rcpp::Named("weights") = model_hefts,
                                     Rcpp::Named("Log_likelihood_raw") = loglik,
                                     Rcpp::Named("avg_Log_likelihood_DATA") = model_avg_log_p );
        }
        else {
          arma::gmm_diag model_inner;
          model_inner = GMM_arma_covariance_type<arma::gmm_diag>(model_inner,
                                                                 data,
                                                                 gaussian_comps,
                                                                 dist_mode,
                                                                 seed_mode,
                                                                 km_iter,
                                                                 em_iter,
                                                                 verbose,
                                                                 var_floor,
                                                                 seed);

          arma::mat loglik(data.n_rows, gaussian_comps, arma::fill::zeros);

          for (int j = 0; j < gaussian_comps; j++) {
            loglik.col(j) = arma::conv_to< arma::vec >::from(model_inner.log_p(data.t(), j));
          }

          arma::mat model_means = model_inner.means.t();
          arma::mat model_COVS = model_inner.dcovs.t();           // each row of the 'covariance_matrices' is a different covariance matrix, use diag() to build each square diagonal matrix

          arma::rowvec model_hefts = arma::conv_to< arma::rowvec >::from(model_inner.hefts.t());
          double model_avg_log_p = model_inner.avg_log_p(data.t(), gaussian_comps - 1);

          double n = timer.toc();
          if (verbose) { Rcpp::Rcout << "\ntime to complete : " << n << "\n" << std::endl; }

          return Rcpp::List::create( Rcpp::Named("centroids") = model_means,
                                     Rcpp::Named("covariance_matrices") = model_COVS,
                                     Rcpp::Named("weights") = model_hefts,
                                     Rcpp::Named("Log_likelihood_raw") = loglik,
                                     Rcpp::Named("avg_Log_likelihood_DATA") = model_avg_log_p );
        }
      }


      // take a diagonal matrix in form of a vector and build a square diagonal matrix, then invert it
      //

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
      // The predictions are based on centroids, covariance matrix and weights using the formulas and not on the log-likelihoods.
      //

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

            double tmp_determinant = arma::det(tmp_cov_mt);

            //if (tmp_determinant <= 1.0e-8) {
            if (tmp_determinant == 0.0) {

              Rcpp::stop("the determinant is zero or approximately zero. The data might include highly correlated variables or variables with low variance");
            }

            double tmp_val = 1.0 / std::sqrt(2.0 * arma::datum::pi * tmp_determinant);               // use determinant to get a single value

            double inner_likelih = 0.5 * (arma::as_scalar(tmp_vec.t() * INV_COV(arma::conv_to< arma::vec >::from(COVARIANCE.row(i))) * arma::conv_to< arma::mat >::from(tmp_vec)));

            gaus_vec_log(j) = -(n / 2.0) * std::log(2.0 * arma::datum::pi) - (1.0 / 2.0) * (std::log(tmp_determinant)) - inner_likelih;

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

        return Rcpp::List::create( Rcpp::Named("Log_likelihood_raw") = gaus_mat_log_lik,
                                   Rcpp::Named("cluster_proba") = loglik1,
                                   Rcpp::Named("cluster_labels") = loglik2 );
      }




      // function to calculate bic-aic
      //

      arma::rowvec GMM_arma_AIC_BIC(arma::mat& data, arma::rowvec max_clusters, std::string dist_mode, std::string seed_mode,

                                    int km_iter, int em_iter, bool verbose, double var_floor = 1e-10, std::string criterion = "AIC", int seed = 1) {

        int LEN_max_clust = max_clusters.n_elem;

        set_seed(seed);             // R's RNG

        arma::rowvec evaluate_comps(LEN_max_clust, arma::fill::zeros); //, aic_avg_weights(LEN_max_clust, arma::fill::zeros);

        for (int i = 0; i < LEN_max_clust; i++) {

          if (verbose) { Rcpp::Rcout << "iteration: " << i + 1 << "  num-clusters: " << max_clusters(i) << std::endl; }

          Rcpp::List gmm = GMM_arma(data, max_clusters(i), dist_mode, seed_mode, km_iter, em_iter, false, var_floor = 1e-10);

          arma::mat loglik = Rcpp::as<arma::mat> (gmm[3]);

          arma::rowvec weights = Rcpp::as<arma::rowvec> (gmm[2]);

          arma::rowvec log_sum_exp(data.n_rows, arma::fill::zeros);

          for (unsigned int j = 0; j < loglik.n_rows; j++) {

            arma::rowvec tmp_vec = arma::conv_to< arma::rowvec >::from(loglik.row(j));

            tmp_vec += arma::log(weights);                                     // https://github.com/scikit-learn/scikit-learn/blob/51a765acfa4c5d1ec05fc4b406968ad233c75162/sklearn/utils/extmath.py

            double max_row = arma::max(tmp_vec);

            tmp_vec = arma::exp(tmp_vec - max_row);

            double tmp_log = arma::sum(tmp_vec);

            log_sum_exp(j) = max_row + std::log(tmp_log);
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



      // the various distance metric methods for the 'dissim_mat' function
      //

      double METHODS(arma::mat& data,
                     arma::mat& data1,
                     std::string& method,
                     unsigned int i,
                     unsigned int j,
                     bool flag_isfinite,
                     arma::mat& cov_mat,
                     double minkowski_p = 1.0,
                     double eps = 1.0e-6,
                     bool exception_nan = true) {

        double tmp_idx;

        if (method == "euclidean") {

          if (flag_isfinite) {

            tmp_idx = std::sqrt(arma::as_scalar(arma::accu(arma::square((data.row(i) - data1.row(j))))));}

          else {

            arma::rowvec tmp_idx1;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f) - data1.row(j)(f);

                count += 1;
              }
            }

            tmp_idx = std::sqrt(arma::as_scalar(arma::accu(arma::square(tmp_idx1))));
          }
        }

        else if (method == "manhattan") {

          if (flag_isfinite) {

            tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - data1.row(j)))));}

          else {

            arma::rowvec tmp_idx1;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f) - data1.row(j)(f);

                count += 1;
              }
            }

            tmp_idx = arma::as_scalar(arma::accu(arma::abs(tmp_idx1)));
          }
        }

        else if (method == "chebyshev") {

          if (flag_isfinite) {

            tmp_idx = arma::as_scalar(arma::max(arma::abs((data.row(i) - data1.row(j)))));}

          else {

            arma::rowvec tmp_idx1;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f) - data1.row(j)(f);

                count += 1;
              }
            }

            tmp_idx = arma::as_scalar(arma::max(arma::abs(tmp_idx1)));
          }
        }

        else if (method == "canberra") {

          if (flag_isfinite) {

            tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - data1.row(j)) + eps)/(arma::abs(data.row(i)) + arma::abs(data1.row(j)) + eps)));}       // added 'eps' otherwise rstudio crashes

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              tmp_idx = arma::as_scalar(arma::accu(arma::abs((tmp_idx1 - tmp_idx2) + eps)/(arma::abs(tmp_idx1) + arma::abs(tmp_idx2) + eps)));}

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else if (method == "braycurtis") {

          if (flag_isfinite) {

            tmp_idx = arma::as_scalar(arma::accu(arma::abs((data.row(i) - data1.row(j))))/(arma::accu(arma::abs(data.row(i))) + arma::accu(arma::abs(data1.row(j)))));}

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              tmp_idx = arma::as_scalar(arma::accu(arma::abs((tmp_idx1 - tmp_idx2)))/(arma::accu(arma::abs(tmp_idx1)) + arma::accu(arma::abs(tmp_idx2))));}

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else if (method == "pearson_correlation") {

          if (flag_isfinite) {

            tmp_idx = arma::as_scalar(1.0 - arma::cor(data.row(i), data1.row(j)));}

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              tmp_idx = arma::as_scalar(1.0 - arma::cor(tmp_idx1, tmp_idx2));}

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else if (method == "cosine") {

          if (flag_isfinite) {

            tmp_idx = 1.0 - (arma::as_scalar(arma::accu(data.row(i) % data1.row(j))) / (std::sqrt(arma::accu(arma::pow(data.row(i), 2.0))) * std::sqrt(arma::accu(arma::pow(data1.row(j), 2.0)))));}

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              tmp_idx = 1.0 - (arma::as_scalar(arma::accu(tmp_idx1 % tmp_idx2)) / (std::sqrt(arma::accu(arma::pow(tmp_idx1, 2.0))) * std::sqrt(arma::accu(arma::pow(tmp_idx2, 2.0)))));}

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else if (method == "simple_matching_coefficient") {                                                    // for binary data

          if (flag_isfinite) {

            double a = eps;
            double d = eps;

            for (unsigned int t = 0; t < data.row(i).n_elem; t++) {

              if (data.row(i)(t) == 1 && data1.row(j)(t) == 1) {

                a += 1.0;}

              if (data.row(i)(t) == 0 && data1.row(j)(t) == 0) {

                d += 1.0;
              }
            }

            tmp_idx = 1.0 - ((a + d) / data.row(i).n_elem);
          }

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              double a = eps;
              double d = eps;

              for (unsigned int t = 0; t < tmp_idx1.n_elem; t++) {

                if (tmp_idx1(t) == 1 && tmp_idx2(t) == 1) {

                  a += 1.0;}

                if (tmp_idx1(t) == 0 && tmp_idx2(t) == 0) {

                  d += 1.0;
                }
              }

              tmp_idx = 1.0 - ((a + d) / tmp_idx1.n_elem);
            }

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else if (method == "minkowski") {                                                                                     // by default the order of the minkowski parameter equals k

          if (flag_isfinite) {

            tmp_idx = std::pow(arma::as_scalar(arma::accu(arma::pow(arma::abs((data.row(i) - data1.row(j))), minkowski_p))), 1.0 / minkowski_p);
          }

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              tmp_idx = std::pow(arma::as_scalar(arma::accu(arma::pow(arma::abs((tmp_idx1 -tmp_idx2)), minkowski_p))), 1.0 / minkowski_p);}

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else if (method == "hamming") {                                                                                     // for binary data

          if (flag_isfinite) {

            tmp_idx = arma::as_scalar(arma::accu(data.row(i) != data1.row(j))/(data.row(i).n_elem * 1.0));
          }

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              tmp_idx = arma::as_scalar(arma::accu(tmp_idx1 != tmp_idx2)/(tmp_idx1.n_elem * 1.0));}

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else if (method == "mahalanobis") {                                                                                     // first create covariance matrix from data

          tmp_idx = arma::as_scalar(std::sqrt(arma::as_scalar(((data.row(i) - data1.row(j)) * cov_mat) * (data.row(i) - data1.row(j)).t())));
        }

        else if (method == "jaccard_coefficient") {                                                                                     // for binary data

          if (flag_isfinite) {

            double a = eps;
            double b = eps;
            double c = eps;

            for (unsigned int t = 0; t < data.row(i).n_elem; t++) {

              if (data.row(i)(t) == 1 && data1.row(j)(t) == 1) {

                a += 1.0;}

              if (data.row(i)(t) == 1 && data1.row(j)(t) == 0) {

                b += 1.0;}

              if (data.row(i)(t) == 0 && data1.row(j)(t) == 1) {

                c += 1.0;
              }
            }

            tmp_idx = 1.0 - (a / (a + b + c));
          }

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              double a = eps;
              double b = eps;
              double c = eps;

              for (unsigned int t = 0; t < tmp_idx1.n_elem; t++) {

                if (tmp_idx1(t) == 1 && tmp_idx2(t) == 1) {

                  a += 1.0;}

                if (tmp_idx1(t) == 1 && tmp_idx2(t) == 0) {

                  b += 1.0;}

                if (tmp_idx1(t) == 0 && tmp_idx2(t) == 1) {

                  c += 1.0;
                }
              }

              tmp_idx = 1.0 - (a / (a + b + c));
            }

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else if (method == "Rao_coefficient") {                                                                                     // for binary data

          if (flag_isfinite) {

            double a = eps;

            for (unsigned int t = 0; t < data.row(i).n_elem; t++) {

              if (data.row(i)(t) == 1 && data1.row(j)(t) == 1) {

                a += 1.0;
              }
            }

            tmp_idx = 1.0 - (a / data.row(i).n_elem);
          }

          else {

            arma::rowvec tmp_idx1, tmp_idx2;

            int count = 1;

            for (unsigned int f = 0; f < data.row(i).n_elem; f++) {

              if (arma::is_finite(data.row(i)(f)) && arma::is_finite(data1.row(j)(f))) {

                tmp_idx1.set_size(count);

                tmp_idx2.set_size(count);

                tmp_idx1(count - 1) = data.row(i)(f);

                tmp_idx2(count - 1) = data1.row(j)(f);

                count += 1;
              }
            }

            if (!tmp_idx1.is_empty()) {

              double a = eps;

              for (unsigned int t = 0; t < tmp_idx1.n_elem; t++) {

                if (tmp_idx1(t) == 1 && tmp_idx2(t) == 1) {

                  a += 1.0;
                }
              }

              tmp_idx = 1.0 - (a / tmp_idx1.n_elem);
            }

            else {

              tmp_idx = arma::datum::nan;
            }
          }
        }

        else {

          if (exception_nan) {

            tmp_idx = arma::datum::nan;}                                             // default = 0; create exceptions in R, so that tmp_idx = arma::datum::nan; [ applies to "dissim_mat" ]

          else {

            tmp_idx = 1.0;                                                           // default = 1.0; create exceptions in R, so that tmp_idx is never 1.0;     [ applies to "dissim_MEDOIDS" ]
          }
        }

        return tmp_idx;
      }



      // calculate the 'inverse' AND in case of exception the 'Moore-Penrose pseudo-inverse' of the covariance matrix FOR the 'mahalanobis' distance
      // https://github.com/mlampros/KernelKnn/issues/1
      //

      arma::mat INV_EXC(arma::mat cov_data) {

        arma::mat inv_tmp;

        try {

          inv_tmp = arma::inv(arma::cov(cov_data));
        }

        catch(...) {

          Rcpp::warning("the input matrix seems singular. The Moore-Penrose pseudo-inverse of the covariance matrix will be calculated");
        }

        if (inv_tmp.empty()) {

          inv_tmp = arma::pinv(arma::cov(cov_data));
        }

        return inv_tmp;
      }



      // dissimilarity matrix using various distance metric methods [ the function can handle missing values by using pair-wise deletion ]
      //

      arma::mat dissim_mat(arma::mat& data, std::string& method, double minkowski_p = 1.0, bool upper = true, bool diagonal = true, int threads = 1, double eps = 1.0e-6) {

        #ifdef _OPENMP
        omp_set_num_threads(threads);
        #endif

        bool flag_isfinite = data.is_finite();

        if (!flag_isfinite && (method == "mahalanobis")) {

          Rcpp::stop("in case of missing values the mahalanobis distance calculation is not feasible");
        }

        unsigned int COLS = data.n_cols;

        arma::mat cov_mat(COLS, COLS);                            // initialize for the openmp clause [ cov_mat dimensions equals to (data_number_columns, data_number_columns) ]

        if (method == "mahalanobis") {

          cov_mat = INV_EXC(data);
        }

        unsigned int ITERS = data.n_rows;

        arma::mat mt(ITERS, ITERS);

        mt.fill(arma::datum::nan);

        unsigned int i,j;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) shared(ITERS, method, mt, data, upper, flag_isfinite, cov_mat, minkowski_p, eps) private(i,j)
        #endif
        for (i = 0; i < ITERS - 1; i++) {                                                                                   // rows

          for (j = i + 1; j < ITERS; j++) {                                                                                 // rows

            double val = METHODS(data, data, method, i, j, flag_isfinite, cov_mat, minkowski_p, eps, true);                 // various methods

            #ifdef _OPENMP
            #pragma omp atomic write
            #endif
            mt(j,i) = val;

            if (upper) {

              #ifdef _OPENMP
              #pragma omp atomic write
              #endif
              mt(i,j) = val;
            }
          }
        }

        if (diagonal) {

          mt.diag().zeros();
        }

        return mt;
      }




      // secondary function to check if integer in vector of integers
      //

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



      // secondary function for 'silhouette_matrix'
      //

      arma::field<arma::mat> inner_field_func(unsigned int f, unsigned int sorted_medoids_elem, unsigned int END_IDX_nelem, arma::rowvec& end_indices_vec, arma::mat& data,

                                              arma::rowvec& sorted_medoids, arma::rowvec& sorted_medoids_increment) {


        double intra_clust_diss = 0.0, diameter = 0.0, separation = arma::datum::inf;

        int intra_count_obs = 0;

        arma::mat outer_clust_diss(3, sorted_medoids_elem, arma::fill::zeros);

        for (unsigned int t = 0; t < END_IDX_nelem; t++) {

          if (f != t) {

            if (end_indices_vec(f) == end_indices_vec(t)) {

              intra_clust_diss += data(f,t);

              if (diameter < data(f,t)) { diameter = data(f,t); }

              intra_count_obs += 1;
            }

            if (end_indices_vec(f) != end_indices_vec(t)) {

              for (unsigned int k = 0; k < sorted_medoids_elem; k++) {

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

        //------------------
        // 1st field-matrix:
        //------------------

        arma::mat SCALARS(1, 5, arma::fill::zeros);

        SCALARS(0,0) = intra_clust_diss;
        SCALARS(0,1) = diameter;
        SCALARS(0,2) = separation;
        SCALARS(0,3) = intra_count_obs;

        //------------------
        // 2nd field-matrix:
        //------------------

        outer_clust_diss.row(1) /= outer_clust_diss.row(2);                                                // returns Inf if the outer-cluster-count is 1 (a single value as cluster index)

        //---------------------------------------
        // 1st field-matrix extend (add opt_idx):
        //---------------------------------------

        arma::uvec opt_idx = arma::find(outer_clust_diss.row(1) == arma::min(outer_clust_diss.row(1)));

        SCALARS(0,4) = arma::as_scalar(opt_idx(0));

        //---------------------
        // return field object:
        //---------------------

        arma::field<arma::mat> OUT_FIELD(2,1);

        OUT_FIELD(0,0) = SCALARS;
        OUT_FIELD(1,0) = outer_clust_diss;

        return OUT_FIELD;
      }



      // intra-outer-silhouette cluster values AND clustering-statistics
      //

      Rcpp::List silhouette_matrix(arma::mat& data, arma::rowvec end_indices_vec, arma::rowvec end_cost_vec, int threads = 1) {

        #ifdef _OPENMP
        omp_set_num_threads(threads);
        #endif

        unsigned int END_IDX_nelem = end_indices_vec.n_elem;

        arma::rowvec sorted_medoids = arma::unique(end_indices_vec);

        unsigned int sorted_medoids_elem = sorted_medoids.n_elem;

        arma::rowvec sorted_medoids_increment = arma::regspace<arma::rowvec> (1, 1, sorted_medoids_elem);

        arma::mat Silhouette_matrix(END_IDX_nelem, 7, arma::fill::zeros);

        unsigned int f;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) shared(END_IDX_nelem, sorted_medoids_elem, end_indices_vec, data, sorted_medoids, sorted_medoids_increment, Silhouette_matrix) private(f)
        #endif
        for (f = 0; f < END_IDX_nelem; f++) {

          //--------------
          // field object
          //--------------

          arma::field<arma::mat> RES = inner_field_func(f, sorted_medoids_elem, END_IDX_nelem, end_indices_vec, data, sorted_medoids, sorted_medoids_increment);

          arma::mat unl_scalars = RES(0,0);             // first matrix

          double intra_clust_diss = unl_scalars(0,0);
          double diameter = unl_scalars(0,1);
          double separation = unl_scalars(0,2);
          double intra_count_obs = unl_scalars(0,3);
          int opt_idx = unl_scalars(0,4);

          arma::mat outer_clust_diss = RES(1,0);        // second matrix

          //----------------------------
          // populate silhouette matrix
          //----------------------------

          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          Silhouette_matrix(f, 0) = end_indices_vec(f);                                                     // 1st column clusters

          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          Silhouette_matrix(f, 1) = outer_clust_diss(0,opt_idx);                                            // 2nd column neighboring clusters [ using outer-cluster-dissimilarities ]

          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          Silhouette_matrix(f, 2) = intra_clust_diss / intra_count_obs;                                     // 3rd column intra-cluster-dissimilarities

          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          Silhouette_matrix(f, 3) = outer_clust_diss(1,opt_idx);                                            // 4th column outer-cluster-dissimilarities

          double silh_in_first, silh_in_second;

          #ifdef _OPENMP
          #pragma omp atomic read
          #endif
          silh_in_first = Silhouette_matrix(f,2);

          #ifdef _OPENMP
          #pragma omp atomic read
          #endif
          silh_in_second = Silhouette_matrix(f,3);

          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          Silhouette_matrix(f, 4) = calc_silhouette(silh_in_first, silh_in_second);                          // 5th column silhouette widths

          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          Silhouette_matrix(f, 5) = diameter;                                                                // diameter (maximal dissimilarity between two obs. of the same cluster)

          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          Silhouette_matrix(f, 6) = separation;                                                              // separation (minimal dissimilarity between two obs. of different clusters)
        }

        arma::mat clustering_stats(6, sorted_medoids_elem, arma::fill::zeros);

        clustering_stats.row(5).fill(arma::datum::inf);

        for (unsigned int s = 0; s < Silhouette_matrix.n_rows; s++) {

          for (unsigned int g = 0; g < sorted_medoids_elem; g++) {

            if (sorted_medoids(g) == Silhouette_matrix(s, 0)) {

              clustering_stats(0, g) = sorted_medoids(g);                                                                       // clustering labels

              clustering_stats(1, g) += 1;                                                                                      // number of obs in each cluster

              if (clustering_stats(2, g) < end_cost_vec(s)) {

                clustering_stats(2, g) = end_cost_vec(s); }                                                                     // maximum cluster dissimilarity

              clustering_stats(3, g) += end_cost_vec(s);                                                                        // first sum dissimilarities, so that I can get the average cluster dissimilarity

              if (clustering_stats(4, g) < Silhouette_matrix(s, 5)) {

                clustering_stats(4, g) = Silhouette_matrix(s, 5); }                                                             // diameter of cluster

              if (clustering_stats(5, g) > Silhouette_matrix(s, 6)) {

                clustering_stats(5, g) = Silhouette_matrix(s, 6); }                                                             // separation of cluster
            }
          }
        }

        clustering_stats.row(3) /= clustering_stats.row(1);

        clustering_stats = clustering_stats.t();                                                                                // transpose before return the output

        return Rcpp::List::create(Rcpp::Named("silhouette_matrix") = Silhouette_matrix, Rcpp::Named("clustering_stats") = clustering_stats);
      }



      // remove items from a vector using item-values from a second vector
      //

      arma::uvec subset_vec(arma::uvec x, arma::uvec y) {

        std::vector<double> vec = Rcpp::as<std::vector<double> >(Rcpp::wrap(x));

        for (unsigned int i = 0; i < y.n_elem; i++) {

          vec.erase(std::remove(vec.begin(), vec.end(), y(i)), vec.end());
        }

        return Rcpp::as<arma::uvec> (Rcpp::wrap(vec));
      }




      // secondary function used in 'ClusterMedoids'
      //

      arma::field<arma::rowvec> field_cm_inner(arma::uvec& copy_medoids, arma::uvec& non_medoids, arma::mat& data, unsigned int i, unsigned int j) {

        copy_medoids(j) = non_medoids(i);                // swap medoid with a non-medoid

        arma::mat new_cols = data.cols(copy_medoids);

        arma::rowvec tmp_min_cost(data.n_rows, arma::fill::zeros), tmp_clusters(data.n_rows, arma::fill::zeros), new_cost_vec(1, arma::fill::zeros);

        for (unsigned int f = 0; f < new_cols.n_rows; f++) {

          double tmp_cost = arma::min(new_cols.row(f));

          arma::uvec idx = arma::find(new_cols.row(f) == tmp_cost);

          tmp_min_cost(f) = tmp_cost;

          tmp_clusters(f) = copy_medoids(idx(0));
        }

        double new_cost = arma::accu(tmp_min_cost);

        new_cost_vec(0) = new_cost;

        arma::field<arma::rowvec> tmp_field(3,1);

        tmp_field(0,0) = tmp_min_cost;
        tmp_field(1,0) = tmp_clusters;
        tmp_field(2,0) = new_cost_vec;

        return tmp_field;
      }



      // the 'ClusterMedoids' function should accept either a matrix or a dissimilarity matrix
      // the dissimilarity matrix in comparison to a single matrix will have nrows == ncols AND diagonal == 0.0 [ it should be also a matrix ]
      // Be aware this function is not the exact but the "approximate partition around medoids" algorithm
      //

      Rcpp::List ClusterMedoids(arma::mat& data, int clusters, std::string method, double minkowski_p = 1.0, int threads = 1, bool verbose = false, bool swap_phase = false,
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

        int count_clusters = 1;

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

          unsigned int i,j;

          #ifdef _OPENMP
          #pragma omp parallel for schedule(static) shared(end_indices_vec, end_cost_vec, initial_cost, data, non_medoids, swap_medoids, threads, verbose, Rcpp::Rcout) private(i,j)
          #endif
          for (i = 0; i < non_medoids.n_elem; i++) {

            for (j = 0; j < swap_medoids.n_elem; j++) {

              arma::uvec copy_medoids(swap_medoids.n_elem);

              for (unsigned int cm = 0; cm < swap_medoids.n_elem; cm++) {

                #ifdef _OPENMP
                #pragma omp atomic read
                #endif
                copy_medoids(cm) = swap_medoids(cm);
              }

              arma::field<arma::rowvec> unl_field = field_cm_inner(copy_medoids, non_medoids, data, i, j);

              arma::rowvec tmp_min_cost = unl_field(0,0);
              arma::rowvec tmp_clusters = unl_field(1,0);
              arma::rowvec new_cost_vec = unl_field(2,0);

              double new_cost = new_cost_vec(0);

              if (new_cost < initial_cost) {

                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                initial_cost = new_cost;

                if (verbose && threads == 1) {

                  int sw_md, no_md;

                  #ifdef _OPENMP
                  #pragma omp atomic read
                  #endif
                  sw_md = swap_medoids(j);

                  #ifdef _OPENMP
                  #pragma omp atomic read
                  #endif
                  no_md = non_medoids(i);

                  Rcpp::Rcout << "swap of medoid " << sw_md + 1 << " with the non-medoid " << no_md + 1 << ". Current dissimilarity of the swap phase --> " << new_cost << std::endl;
                }

                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                swap_medoids(j) = non_medoids(i);

                for (unsigned int a_VEC = 0; a_VEC < tmp_min_cost.n_elem; a_VEC++) {

                  #ifdef _OPENMP
                  #pragma omp atomic write
                  #endif
                  end_cost_vec(a_VEC) = tmp_min_cost(a_VEC);
                }

                for (unsigned int e_VEC = 0; e_VEC < tmp_clusters.n_elem; e_VEC++) {

                  #ifdef _OPENMP
                  #pragma omp atomic write
                  #endif
                  end_indices_vec(e_VEC) = tmp_clusters(e_VEC);
                }
              }
            }
          }

          if (verbose) { Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << " ================================================== " << std::endl; Rcpp::Rcout << " --- end of swap phase ---> " << " dissimilarity: " << initial_cost << std::endl; Rcpp::Rcout << " ================================================== " << std::endl;}
        }

        arma::rowvec end_idxs = arma::unique(end_indices_vec);

        arma::rowvec true_idxs = arma::regspace<arma::rowvec>(0, 1, end_idxs.n_elem - 1);

        unsigned int item, item_1;

        #ifdef _OPENMP
        #pragma omp parallel for collapse(2) shared(end_indices_vec, true_idxs, end_idxs) private(item, item_1)
        #endif
        for (item = 0; item < end_indices_vec.n_elem; item++) {

          for (item_1 = 0; item_1 < end_idxs.n_elem; item_1++) {

            if (end_indices_vec(item) == end_idxs(item_1)) {

              #ifdef _OPENMP
              #pragma omp atomic write
              #endif
              end_indices_vec(item) = true_idxs(item_1);
            }
          }
        }

        arma::mat befout_silhouette_matrix;

        arma::mat befout_clustering_stats;

        if (clusters > 1) {

          silh_lst = silhouette_matrix(data, end_indices_vec + 1, end_cost_vec, threads);

          befout_silhouette_matrix = Rcpp::as<arma::mat> (silh_lst[0]);

          befout_clustering_stats = Rcpp::as<arma::mat> (silh_lst[1]);
        }

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

        double end_cost_vec_scalar = arma::accu(end_cost_vec);

        return Rcpp::List::create(Rcpp::Named("medoids") = end_idxs,
								  Rcpp::Named("cost") = end_cost_vec_scalar,
								  Rcpp::Named("dissimilarity_matrix") = data,
                                  Rcpp::Named("clusters") = end_indices_vec,
								  Rcpp::Named("end_cost_vec") = end_cost_vec,
								  Rcpp::Named("silhouette_matrix") = befout_silhouette_matrix,
								  Rcpp::Named("fuzzy_probs") = fuz_out,
								  Rcpp::Named("clustering_stats") = befout_clustering_stats,
								  Rcpp::Named("flag_dissim_mat") = flag_dissim_mat);
      }


      // calculate global dissimilarities for claraMedoids
      // the function can handle missing values by using pair-wise deletion
      //

      arma::mat dissim_MEDOIDS(arma::mat& data, std::string& method, arma::mat MEDOIDS, double minkowski_p = 1.0, int threads = 1, double eps = 1.0e-6) {

        #ifdef _OPENMP
        omp_set_num_threads(threads);
        #endif

        bool flag_isfinite = data.is_finite();

        if (!flag_isfinite && (method == "mahalanobis")) {

          Rcpp::stop("in case of missing values the mahalanobis distance calculation is not feasible");
        }

        unsigned int COLS = data.n_cols;

        arma::mat cov_mat(COLS, COLS);                            // initialize for the openmp clause [ cov_mat dimensions equals to (data_number_columns, data_number_columns) ]

        if (method == "mahalanobis") {

          cov_mat = INV_EXC(data);
        }

        unsigned int ITERS = data.n_rows;

        arma::mat mt(ITERS, MEDOIDS.n_rows);

        unsigned int i,j;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) shared(ITERS, MEDOIDS, eps, minkowski_p, cov_mat, flag_isfinite, method, data, mt) private(i,j)
        #endif
        for (i = 0; i < ITERS; i++) {

          for (j = 0; j < MEDOIDS.n_rows; j++) {

            double val = METHODS(data, MEDOIDS, method, i, j, flag_isfinite, cov_mat, minkowski_p, eps, false);             // various methods (of type integer)

            #ifdef _OPENMP
            #pragma omp atomic write
            #endif
            mt(i,j) = val;
          }
        }

        return mt;
      }




      // calculate fuzzy (soft) clusters and stats [ data is the subset of the dissimilarity matrix using the medoids as subset-indices ]
      //

      Rcpp::List fuzzy_and_stats(arma::mat data, double eps = 1.0e-6, bool fuzzy = false) {

        arma::mat fuz_out, stats_mat(5, data.n_cols, arma::fill::zeros);;

        arma::rowvec hard_out(data.n_rows, arma::fill::zeros);

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
      // 2. Compute the PAM algorithm on each subset and choose the corresponding k representative objects (medoids). Assign each observation of the entire dataset to the nearest medoid.
      // 3. Calculate the mean (or the sum) of the dissimilarities (for instance 'euclidean' distance) of the observations to their closest medoid. This is used as a measure of the goodness of the clustering.
      // 4. Retain (keep) the sub-dataset for which the mean (or sum) is minimal. A further analysis is carried out on the final partition.

      //--------------------------------------------------------------------------------------------------------------------------------------------------
      // Important: if a medoid is significant --> add this medoid to the next sample by subsetting the data and appending the corresponding observations
      //--------------------------------------------------------------------------------------------------------------------------------------------------

      // Instead of finding representative objects for the entire data set, 'ClaraMedoids' draws a sample of the data set, applies the 'ClusterMedoids' function
      // on the sample and finds the medoids of the sample. The point is that if the sample is drawn in a sufficiently random way the medoids of the sample would
      // approximate the medoids of the entire data set. To come up with a better approximation, 'ClaraMedoids' draws mulitple samples and GIVES THE BEST CLUSTERING
      // as the OUTPUT. Here, for ACCURACY the QUALITY of a clustering is measured based on the AVERAGE DISSIMILARITY of all objects in the entire data set and NOT
      // ONLY of those objects IN THE SAMPLES. Experiments, indicate that 5 samples of size 40 + 2 * k (i.e. 40 + 2 * number_of_clusters) give satisfactory results.
      //

      Rcpp::List ClaraMedoids(arma::mat& data, int clusters, std::string method, int samples, double sample_size, double minkowski_p = 1.0,

                              int threads = 1, bool verbose = false, bool swap_phase = false, bool fuzzy = false, int seed = 1) {

        set_seed(seed);                                                              // R's RNG
        arma::rowvec set_seed_loop = sample_vec(samples, 0, 1000000, false);         // make sure that the same seed values do not appear in the inner-for-loop (samples) by picking random values based on a big integer (1000000)

        int bst_sample = -1;

        double dism = arma::datum::inf;

        arma::mat probs, diss_glob_out;

        arma::uvec out_medoid, clr_split_out;

        Rcpp::List bst_lst;

        if (verbose) { Rcpp::Rcout << " " << std::endl; }

        for (int s = 0; s < samples; s++) {

          int iter_seed = set_seed_loop(s);
          set_seed(iter_seed);                                    // set the seed also in the inner loop due to the sample

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

          arma::mat tmp_dat = data.rows(clr_split);                                                                                         // use sample data for ClusterMedoids

          arma::mat copy_dat = tmp_dat;

          Rcpp::List clM_sblist = ClusterMedoids(tmp_dat, clusters, method, minkowski_p, threads, false, swap_phase, false, iter_seed);    // the 'seed' here for 'ClusterMedoids' is 'iter_seed' but it will be removed in version 1.4.0 (deprecation warning)

          double local_dissim = Rcpp::as<double> (clM_sblist["cost"]);

          arma::uvec local_medoids = Rcpp::as<arma::uvec> (clM_sblist["medoids"]);

          arma::mat tmp_glob = dissim_MEDOIDS(data, method, copy_dat.rows(local_medoids), minkowski_p , threads, 1.0e-6);                  // use all data to calculate global dissimilarity

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

        Rcpp::List fuz_and_stats = fuzzy_and_stats(diss_glob_out, 1.0e-6, fuzzy);

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

        arma::rowvec clr_split_out_rowvec = arma::conv_to< arma::rowvec >::from(clr_split_out);

        arma::mat fuz_and_stats_mt = Rcpp::as<arma::mat> (fuz_and_stats[1]);

        fuz_st_mat = fuz_st_mat.t();

        return Rcpp::List::create(Rcpp::Named("medoids") = subs_meds,
                                  Rcpp::Named("best_dissimilarity") = dism,
                                  Rcpp::Named("medoid_indices") = out_medoid,
                                  Rcpp::Named("sample_indices") = clr_split_out_rowvec,
                                  Rcpp::Named("clusters") = hard_clust,
                                  Rcpp::Named("silhouette_matrix") = bst_sample_silh_mat,
                                  Rcpp::Named("fuzzy_probs") = fuz_and_stats_mt,
                                  Rcpp::Named("clustering_stats") = fuz_st_mat,
                                  Rcpp::Named("dissimilarity_matrix") = bst_sample_dissm_mat);
      }



      // prediction function for the k-medoids
      //

      Rcpp::List predict_medoids(arma::mat& data, std::string method, arma::mat MEDOIDS, double minkowski_p = 1.0, int threads = 1, bool fuzzy = false, double eps = 1.0e-6) {

        arma::mat tmp_dist = dissim_MEDOIDS(data, method, MEDOIDS, minkowski_p, threads, eps);

        arma::mat fuz_out;

        if (fuzzy) {
          fuz_out.set_size(tmp_dist.n_rows, tmp_dist.n_cols);
        }

        arma::rowvec hard_clust(tmp_dist.n_rows);

        double global_dissimil = 0.0;

        for (unsigned int i = 0; i < tmp_dist.n_rows; i++) {

          arma::rowvec tmp_row = arma::abs(arma::conv_to< arma::rowvec >::from(tmp_dist.row(i)));

          arma::uvec tmp_hard = arma::find(tmp_row == arma::min(tmp_row));

          int idx = tmp_hard(0);

          hard_clust(i) = idx;

          global_dissimil += tmp_dist(i,idx);

          if (fuzzy) {

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

      Rcpp::List split_rcpp_lst(Rcpp::List lst) {

        arma::mat silh_mat = Rcpp::as<arma::mat> (lst["silhouette_matrix"]);

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

        double intra_clust_disml_scalar = avg_intr_dis;            // return this value

        double avg_width_silh = arma::accu(silhouet);

        avg_intr_dis /= silh_mat.n_rows;                           // modify this value

        avg_width_silh /= silh_mat.n_rows;

        bool silh_plot_boolean = true;

        return Rcpp::List::create(Rcpp::Named("avg_intra_clust_dissimilarity") = avg_intr_dis,

                                  Rcpp::Named("sum_intra_dissim") = intra_clust_disml_scalar, Rcpp::Named("avg_width_silhouette") = avg_width_silh,

                                  Rcpp::Named("list_intra_dissm") = intra_dissm, Rcpp::Named("list_silhouette") = silhouet_lst, Rcpp::Named("silhouette_plot") = silh_plot_boolean);
      }




      // optimal number of clusters for the k-medoids
      //

      Rcpp::List OptClust(arma::mat& data, arma::rowvec iter_clust, std::string method, bool clara = false, int samples = 5, double sample_size = 0.001, double minkowski_p = 1.0,

                          std::string criterion = "dissimilarity", int threads = 1, bool swap_phase = false, bool verbose = false, int seed = 1) {

        set_seed(seed);             // R's RNG

        int LEN_max_clust = iter_clust.n_elem;

        Rcpp::List medoids_object(LEN_max_clust);

        if (verbose) { Rcpp::Rcout << " " << std::endl; }

        for (int iter = 0; iter < LEN_max_clust; iter++) {

          if (iter_clust(iter) == 1) {        // It could be the case that the user gives the cluster-vector unsorted ( for instances 1 might appear in index 2), thus check every time and pass inf to clusters = 1

            std::string tmp_c = criterion == "dissimilarity" ? "average dissimilarity" : "average silhouette";

            if (verbose) { Rcpp::Rcout << "number of clusters: "<< iter_clust(iter) << "  -->  " << tmp_c << ": " <<  arma::datum::inf << std::endl; }}

          else {

            if (!clara) {

              Rcpp::List cm_out = ClusterMedoids(data, iter_clust(iter), method, minkowski_p, threads, false, swap_phase, false);

              Rcpp::List tmp_split = split_rcpp_lst(cm_out);

              medoids_object[iter] = tmp_split;

              double tmp_val = 0.0;

              if (criterion == "dissimilarity" && verbose) {

                tmp_val =  Rcpp::as<double> (tmp_split[0]);

                Rcpp::Rcout << "number of clusters: "<< iter_clust(iter) << "  -->  " << "average dissimilarity: " << tmp_val << std::endl;
              }

              if (criterion == "silhouette" && verbose) {

                tmp_val =  Rcpp::as<double> (tmp_split[2]);

                Rcpp::Rcout << "number of clusters: "<< iter_clust(iter) << "  -->  " << "average silhouette: " << tmp_val << std::endl;
              }
            }

            else {

              Rcpp::List cl_out = ClaraMedoids(data, iter_clust(iter), method, samples, sample_size, minkowski_p, threads, false, swap_phase,false);

              Rcpp::List tmp_split = split_rcpp_lst(cl_out);

              medoids_object[iter] = tmp_split;

              double tmp_val = 0.0;

              if (criterion == "dissimilarity" && verbose) {

                tmp_val =  Rcpp::as<double> (tmp_split[0]);

                Rcpp::Rcout << "number of clusters: "<< iter_clust(iter) << "  -->  " << "average dissimilarity: " << tmp_val << std::endl;
              }

              if (criterion == "silhouette" && verbose) {

                tmp_val =  Rcpp::as<double> (tmp_split[2]);

                Rcpp::Rcout << "number of clusters: "<< iter_clust(iter) << "  -->  " << "average silhouette: " << tmp_val << std::endl;
              }
            }
          }
        }

        return medoids_object;
      }


      //==================== cluster Medoids (added code snippets to answer https://github.com/mlampros/ClusterR/issues/39) ====================


      // R's setdiff() in RcppArmadillo using std::set_difference() [see: https://stackoverflow.com/a/41237778/8302386]
      //

      arma::uvec std_setdiff(arma::uvec& x, arma::uvec& y) {

        std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
        std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
        std::vector<int> out;

        std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                            std::inserter(out, out.end()));

        return arma::conv_to<arma::uvec>::from(out);
      }


      // Function that computes:
      //     - the clusters
      //     - the total dissimilarity cost
      // using the distance matrix and the medoids
      //

      Rcpp::List dissimilarity_cost_clusters(arma::mat& dissim_mat,
                                             arma::uvec& medoids) {

        arma::mat subs_dissim = dissim_mat.cols(medoids);
        arma::ucolvec medoid_idxs = arma::index_min(subs_dissim, 1);

        double total_cost = 0.0;
        for (arma::uword t = 0; t < medoid_idxs.n_elem; t++) {
          total_cost += subs_dissim(t, medoid_idxs(t));
        }

        arma::uvec clusters_out = medoids(medoid_idxs);

        return Rcpp::List::create(Rcpp::Named("cost") = total_cost,
                                  Rcpp::Named("clusters") = clusters_out);
      }


      // Attempt to implement the "build" phase of the "partition around medoids" algorithm based on
      //        - the algorithm description of the "Finding Groups in Data" book, 1990 (Kaufman, Rousseeuw)
      //        - the algorithm description in https://www.cs.umb.edu/cs738/pam1.pdf (the PAM Clustering Algorithm), see also "inst/papers_references/the_pam_clustering_algorithm.pdf"
      // However, I observed that still it doesn't return the optimal results compared to the "build" phase inside the "ClusterMedoids()" function (of this file)
      // The only difference is that the "updated_BUILD()" function returns also the "Contribution" or 'Total Gain" as described in the algorithm
      //

      Rcpp::List updated_BUILD(arma::mat& dissim_mat,
                               arma::uword k,
                               bool verbose = false,
                               int seed = 1) {
        set_seed(seed);
        arma::vec init_meds(k);                          // I could use here 'arma::uvec' and avoid the conversion later using 'arma::conv_to<arma::uvec>::from()' however '.fill(arma::datum::nan)' does not work with integers (see: https://github.com/RcppCore/RcppArmadillo/issues/399)
        init_meds.fill(arma::datum::nan);
        arma::rowvec init_cost = arma::sum(dissim_mat, 0);
        arma::uword med1 = arma::index_min(init_cost);
        init_meds(0) = med1;
        arma::uword NCOLS = dissim_mat.n_cols;
        arma::uvec col_idxs = arma::regspace<arma::uvec>(0, 1, NCOLS - 1);
        arma::rowvec end_indices_vec(NCOLS);
        end_indices_vec.fill(med1);
        double total_gain = 0.0;

        if (verbose) {
          Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << "medoid " << med1 + 1 << " was added. Sum of Distances of the first Medoid to the other objects  --> " << init_cost(med1) << std::endl;
        }

        arma::uvec idx_meds_wo_nas, unq_meds_wo_nas, remaining_col_idxs, remaining_candidate_idxs, candidate_idx_uvec;
        arma::mat S, D_s_all;
        arma::uword nonselected_idx, candidate_idx, most_similar_idx, idx_gain;
        arma::vec diff_dissims;
        double D_s, D_i, CONTRIBUTE;

        for (unsigned int item_k = 1; item_k < k; item_k++) {                             // start with index 1 because the first medoid is already computed

          idx_meds_wo_nas = arma::find_finite(init_meds);                                 // unique medoids by removing the NA's (this vector will have at least 1 medoid which is the first medoid)
          arma::uvec init_meds_uvec = arma::conv_to<arma::uvec>::from(init_meds);
          unq_meds_wo_nas = init_meds_uvec(idx_meds_wo_nas);
          S = dissim_mat.cols(unq_meds_wo_nas);                                           // S corresponds to the selected medoids
          remaining_col_idxs = std_setdiff(col_idxs, unq_meds_wo_nas);                    // O - S = U
          arma::vec candidates_GAIN(NCOLS, arma::fill::zeros);

          for (unsigned int i = 0; i < remaining_col_idxs.n_elem; i++) {                  // I  [ i in U ]  ( consider i from unselected )

            candidate_idx = remaining_col_idxs(i);
            candidate_idx_uvec = { candidate_idx };
            remaining_candidate_idxs = std_setdiff(remaining_col_idxs, candidate_idx_uvec);

            for (unsigned int j = 0; j < remaining_candidate_idxs.n_elem; j++) {          // J  [ j in (U - i) ]  ( consider j from (unselected - i) )

              nonselected_idx = remaining_candidate_idxs(j);
              D_s_all = S.row(nonselected_idx);
              most_similar_idx = D_s_all.index_min();
              D_s = arma::as_scalar(D_s_all.col(most_similar_idx));                       // dissimilarity to the closest object in S
              D_i = dissim_mat(candidate_idx, nonselected_idx);                           // dissimilarity between object i and object j

              diff_dissims = {D_s - D_i, 0.0};
              CONTRIBUTE = arma::max(diff_dissims);

              if (D_s > D_i) {
                candidates_GAIN(candidate_idx) += CONTRIBUTE;
              }
            }
          }

          idx_gain = candidates_GAIN.index_max();
          init_meds(item_k) = col_idxs(idx_gain);
          total_gain += candidates_GAIN(idx_gain);

          if (verbose) {
            Rcpp::Rcout << "medoid " << col_idxs(idx_gain) + 1 << " was added. Contribution to Gain  --> " << candidates_GAIN(idx_gain) << "    Total Gain (build phase)  --> " << total_gain << std::endl;
          }
        }

        arma::uvec init_meds_uvec = arma::conv_to<arma::uvec>::from(init_meds);
        Rcpp::List cost_clusts_obj = dissimilarity_cost_clusters(dissim_mat, init_meds_uvec);
        double total_cost = Rcpp::as<double>(cost_clusts_obj["cost"]);
        arma::uvec clusters = Rcpp::as<arma::uvec>(cost_clusts_obj["clusters"]);

        if (verbose) {
          Rcpp::Rcout << " " << std::endl; Rcpp::Rcout << " ========================================================================= " << std::endl;
          Rcpp::Rcout << " --- end of build phase ---> " << " Total Gain: " << total_gain << "  dissimilarity: " << total_cost << std::endl;
          Rcpp::Rcout << " ========================================================================= " << std::endl;
          Rcpp::Rcout << "Appending objects to medoids (clustering)" << std::endl;
        }

        return Rcpp::List::create(Rcpp::Named("medoids") = init_meds,
                                  Rcpp::Named("cost") = total_cost,
                                  Rcpp::Named("gain") = total_gain,
                                  Rcpp::Named("clusters") = clusters);
      }


      ~ClustHeader() { }
  };

}

#endif
