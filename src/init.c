#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ClusterR_boolean_function(SEXP, SEXP);
extern SEXP ClusterR_calc_silhouette(SEXP, SEXP);
extern SEXP ClusterR_check_medoids(SEXP, SEXP, SEXP);
extern SEXP ClusterR_check_NaN_Inf(SEXP);
extern SEXP ClusterR_ClaraMedoids(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_cluster_indices(SEXP);
extern SEXP ClusterR_ClusterMedoids(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_dissim_mat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_dissim_MEDOIDS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_duplicated_flag(SEXP);
extern SEXP ClusterR_evaluation_rcpp(SEXP, SEXP, SEXP);
extern SEXP ClusterR_fuzzy_and_stats(SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_GMM_arma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_GMM_arma_AIC_BIC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_INTRA_CLUSTER_DISS(SEXP, SEXP);
extern SEXP ClusterR_INV_COV(SEXP);
extern SEXP ClusterR_isolation(SEXP, SEXP);
extern SEXP ClusterR_KMEANS_arma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_kmeans_pp_dist(SEXP, SEXP);
extern SEXP ClusterR_kmeans_pp_init(SEXP, SEXP, SEXP);
extern SEXP ClusterR_KMEANS_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_mini_batch_kmeans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_MinMat(SEXP);
extern SEXP ClusterR_norm_fuzzy(SEXP, SEXP);
extern SEXP ClusterR_OptClust(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_opt_clust_fK(SEXP, SEXP, SEXP);
extern SEXP ClusterR_predict_medoids(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_predict_MGausDPDF(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_Predict_mini_batch_kmeans(SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_quantile_init_rcpp(SEXP, SEXP, SEXP);
extern SEXP ClusterR_quantile_value(SEXP, SEXP);
extern SEXP ClusterR_Rcpp_2arma_mat(SEXP);
extern SEXP ClusterR_sample_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_SCALE(SEXP, SEXP, SEXP);
extern SEXP ClusterR_set_seed(SEXP);
extern SEXP ClusterR_silhouette_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_SILHOUETTE_metric(SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusterR_split_rcpp_lst(SEXP);
extern SEXP ClusterR_squared_norm(SEXP);
extern SEXP ClusterR_subset_vec(SEXP, SEXP);
extern SEXP ClusterR_tot_ss_data(SEXP);
extern SEXP ClusterR_validate_centroids(SEXP, SEXP);
extern SEXP ClusterR_WCSS(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ClusterR_boolean_function",          (DL_FUNC) &ClusterR_boolean_function,           2},
    {"ClusterR_calc_silhouette",           (DL_FUNC) &ClusterR_calc_silhouette,            2},
    {"ClusterR_check_medoids",             (DL_FUNC) &ClusterR_check_medoids,              3},
    {"ClusterR_check_NaN_Inf",             (DL_FUNC) &ClusterR_check_NaN_Inf,              1},
    {"ClusterR_ClaraMedoids",              (DL_FUNC) &ClusterR_ClaraMedoids,              11},
    {"ClusterR_cluster_indices",           (DL_FUNC) &ClusterR_cluster_indices,            1},
    {"ClusterR_ClusterMedoids",            (DL_FUNC) &ClusterR_ClusterMedoids,             9},
    {"ClusterR_dissim_mat",                (DL_FUNC) &ClusterR_dissim_mat,                 7},
    {"ClusterR_dissim_MEDOIDS",            (DL_FUNC) &ClusterR_dissim_MEDOIDS,             6},
    {"ClusterR_duplicated_flag",           (DL_FUNC) &ClusterR_duplicated_flag,            1},
    {"ClusterR_evaluation_rcpp",           (DL_FUNC) &ClusterR_evaluation_rcpp,            3},
    {"ClusterR_fuzzy_and_stats",           (DL_FUNC) &ClusterR_fuzzy_and_stats,            4},
    {"ClusterR_GMM_arma",                  (DL_FUNC) &ClusterR_GMM_arma,                   9},
    {"ClusterR_GMM_arma_AIC_BIC",          (DL_FUNC) &ClusterR_GMM_arma_AIC_BIC,          10},
    {"ClusterR_INTRA_CLUSTER_DISS",        (DL_FUNC) &ClusterR_INTRA_CLUSTER_DISS,         2},
    {"ClusterR_INV_COV",                   (DL_FUNC) &ClusterR_INV_COV,                    1},
    {"ClusterR_isolation",                 (DL_FUNC) &ClusterR_isolation,                  2},
    {"ClusterR_KMEANS_arma",               (DL_FUNC) &ClusterR_KMEANS_arma,                7},
    {"ClusterR_kmeans_pp_dist",            (DL_FUNC) &ClusterR_kmeans_pp_dist,             2},
    {"ClusterR_kmeans_pp_init",            (DL_FUNC) &ClusterR_kmeans_pp_init,             3},
    {"ClusterR_KMEANS_rcpp",               (DL_FUNC) &ClusterR_KMEANS_rcpp,               13},
    {"ClusterR_mini_batch_kmeans",         (DL_FUNC) &ClusterR_mini_batch_kmeans,         13},
    {"ClusterR_MinMat",                    (DL_FUNC) &ClusterR_MinMat,                     1},
    {"ClusterR_norm_fuzzy",                (DL_FUNC) &ClusterR_norm_fuzzy,                 2},
    {"ClusterR_OptClust",                  (DL_FUNC) &ClusterR_OptClust,                  12},
    {"ClusterR_opt_clust_fK",              (DL_FUNC) &ClusterR_opt_clust_fK,               3},
    {"ClusterR_predict_medoids",           (DL_FUNC) &ClusterR_predict_medoids,            7},
    {"ClusterR_predict_MGausDPDF",         (DL_FUNC) &ClusterR_predict_MGausDPDF,          5},
    {"ClusterR_Predict_mini_batch_kmeans", (DL_FUNC) &ClusterR_Predict_mini_batch_kmeans,  4},
    {"ClusterR_quantile_init_rcpp",        (DL_FUNC) &ClusterR_quantile_init_rcpp,         3},
    {"ClusterR_quantile_value",            (DL_FUNC) &ClusterR_quantile_value,             2},
    {"ClusterR_Rcpp_2arma_mat",            (DL_FUNC) &ClusterR_Rcpp_2arma_mat,             1},
    {"ClusterR_sample_vec",                (DL_FUNC) &ClusterR_sample_vec,                 4},
    {"ClusterR_SCALE",                     (DL_FUNC) &ClusterR_SCALE,                      3},
    {"ClusterR_set_seed",                  (DL_FUNC) &ClusterR_set_seed,                   1},
    {"ClusterR_silhouette_matrix",         (DL_FUNC) &ClusterR_silhouette_matrix,          4},
    {"ClusterR_SILHOUETTE_metric",         (DL_FUNC) &ClusterR_SILHOUETTE_metric,          4},
    {"ClusterR_split_rcpp_lst",            (DL_FUNC) &ClusterR_split_rcpp_lst,             1},
    {"ClusterR_squared_norm",              (DL_FUNC) &ClusterR_squared_norm,               1},
    {"ClusterR_subset_vec",                (DL_FUNC) &ClusterR_subset_vec,                 2},
    {"ClusterR_tot_ss_data",               (DL_FUNC) &ClusterR_tot_ss_data,                1},
    {"ClusterR_validate_centroids",        (DL_FUNC) &ClusterR_validate_centroids,         2},
    {"ClusterR_WCSS",                      (DL_FUNC) &ClusterR_WCSS,                       2},
    {NULL, NULL, 0}
};

void R_init_ClusterR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}