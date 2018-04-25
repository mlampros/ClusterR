#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ClusterR_boolean_function(SEXP, SEXP);
extern SEXP _ClusterR_calc_silhouette(SEXP, SEXP);
extern SEXP _ClusterR_check_medoids(SEXP, SEXP, SEXP);
extern SEXP _ClusterR_check_NaN_Inf(SEXP);
extern SEXP _ClusterR_ClaraMedoids(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_cluster_indices(SEXP);
extern SEXP _ClusterR_ClusterMedoids(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_dissim_mat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_dissim_MEDOIDS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_duplicated_flag(SEXP);
extern SEXP _ClusterR_evaluation_rcpp(SEXP, SEXP, SEXP);
extern SEXP _ClusterR_field_cm_inner(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_fuzzy_and_stats(SEXP, SEXP, SEXP);
extern SEXP _ClusterR_GMM_arma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_GMM_arma_AIC_BIC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_inner_field_func(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_INTRA_CLUSTER_DISS(SEXP, SEXP);
extern SEXP _ClusterR_INV_COV(SEXP);
extern SEXP _ClusterR_isolation(SEXP, SEXP);
extern SEXP _ClusterR_KMEANS_arma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_kmeans_pp_dist(SEXP, SEXP);
extern SEXP _ClusterR_kmeans_pp_init(SEXP, SEXP, SEXP);
extern SEXP _ClusterR_KMEANS_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_METHODS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_mini_batch_kmeans(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_MinMat(SEXP);
extern SEXP _ClusterR_norm_fuzzy(SEXP, SEXP);
extern SEXP _ClusterR_OptClust(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_opt_clust_fK(SEXP, SEXP, SEXP);
extern SEXP _ClusterR_predict_medoids(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_predict_MGausDPDF(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_Predict_mini_batch_kmeans(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_quantile_init_rcpp(SEXP, SEXP, SEXP);
extern SEXP _ClusterR_quantile_value(SEXP, SEXP);
extern SEXP _ClusterR_Rcpp_2arma_mat(SEXP);
extern SEXP _ClusterR_sample_vec(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_SCALE(SEXP, SEXP, SEXP);
extern SEXP _ClusterR_set_seed(SEXP);
extern SEXP _ClusterR_silhouette_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_SILHOUETTE_metric(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ClusterR_split_rcpp_lst(SEXP);
extern SEXP _ClusterR_squared_norm(SEXP);
extern SEXP _ClusterR_subset_vec(SEXP, SEXP);
extern SEXP _ClusterR_tot_ss_data(SEXP);
extern SEXP _ClusterR_validate_centroids(SEXP, SEXP);
extern SEXP _ClusterR_WCSS(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ClusterR_boolean_function",          (DL_FUNC) &_ClusterR_boolean_function,           2},
    {"_ClusterR_calc_silhouette",           (DL_FUNC) &_ClusterR_calc_silhouette,            2},
    {"_ClusterR_check_medoids",             (DL_FUNC) &_ClusterR_check_medoids,              3},
    {"_ClusterR_check_NaN_Inf",             (DL_FUNC) &_ClusterR_check_NaN_Inf,              1},
    {"_ClusterR_ClaraMedoids",              (DL_FUNC) &_ClusterR_ClaraMedoids,              11},
    {"_ClusterR_cluster_indices",           (DL_FUNC) &_ClusterR_cluster_indices,            1},
    {"_ClusterR_ClusterMedoids",            (DL_FUNC) &_ClusterR_ClusterMedoids,             9},
    {"_ClusterR_dissim_mat",                (DL_FUNC) &_ClusterR_dissim_mat,                 7},
    {"_ClusterR_dissim_MEDOIDS",            (DL_FUNC) &_ClusterR_dissim_MEDOIDS,             6},
    {"_ClusterR_duplicated_flag",           (DL_FUNC) &_ClusterR_duplicated_flag,            1},
    {"_ClusterR_evaluation_rcpp",           (DL_FUNC) &_ClusterR_evaluation_rcpp,            3},
    {"_ClusterR_field_cm_inner",            (DL_FUNC) &_ClusterR_field_cm_inner,             5},
    {"_ClusterR_fuzzy_and_stats",           (DL_FUNC) &_ClusterR_fuzzy_and_stats,            3},
    {"_ClusterR_GMM_arma",                  (DL_FUNC) &_ClusterR_GMM_arma,                   9},
    {"_ClusterR_GMM_arma_AIC_BIC",          (DL_FUNC) &_ClusterR_GMM_arma_AIC_BIC,          10},
    {"_ClusterR_inner_field_func",          (DL_FUNC) &_ClusterR_inner_field_func,           7},
    {"_ClusterR_INTRA_CLUSTER_DISS",        (DL_FUNC) &_ClusterR_INTRA_CLUSTER_DISS,         2},
    {"_ClusterR_INV_COV",                   (DL_FUNC) &_ClusterR_INV_COV,                    1},
    {"_ClusterR_isolation",                 (DL_FUNC) &_ClusterR_isolation,                  2},
    {"_ClusterR_KMEANS_arma",               (DL_FUNC) &_ClusterR_KMEANS_arma,                7},
    {"_ClusterR_kmeans_pp_dist",            (DL_FUNC) &_ClusterR_kmeans_pp_dist,             2},
    {"_ClusterR_kmeans_pp_init",            (DL_FUNC) &_ClusterR_kmeans_pp_init,             3},
    {"_ClusterR_KMEANS_rcpp",               (DL_FUNC) &_ClusterR_KMEANS_rcpp,               12},
    {"_ClusterR_METHODS",                   (DL_FUNC) &_ClusterR_METHODS,                   10},
    {"_ClusterR_mini_batch_kmeans",         (DL_FUNC) &_ClusterR_mini_batch_kmeans,         13},
    {"_ClusterR_MinMat",                    (DL_FUNC) &_ClusterR_MinMat,                     1},
    {"_ClusterR_norm_fuzzy",                (DL_FUNC) &_ClusterR_norm_fuzzy,                 2},
    {"_ClusterR_OptClust",                  (DL_FUNC) &_ClusterR_OptClust,                  12},
    {"_ClusterR_opt_clust_fK",              (DL_FUNC) &_ClusterR_opt_clust_fK,               3},
    {"_ClusterR_predict_medoids",           (DL_FUNC) &_ClusterR_predict_medoids,            7},
    {"_ClusterR_predict_MGausDPDF",         (DL_FUNC) &_ClusterR_predict_MGausDPDF,          5},
    {"_ClusterR_Predict_mini_batch_kmeans", (DL_FUNC) &_ClusterR_Predict_mini_batch_kmeans,  4},
    {"_ClusterR_quantile_init_rcpp",        (DL_FUNC) &_ClusterR_quantile_init_rcpp,         3},
    {"_ClusterR_quantile_value",            (DL_FUNC) &_ClusterR_quantile_value,             2},
    {"_ClusterR_Rcpp_2arma_mat",            (DL_FUNC) &_ClusterR_Rcpp_2arma_mat,             1},
    {"_ClusterR_sample_vec",                (DL_FUNC) &_ClusterR_sample_vec,                 4},
    {"_ClusterR_SCALE",                     (DL_FUNC) &_ClusterR_SCALE,                      3},
    {"_ClusterR_set_seed",                  (DL_FUNC) &_ClusterR_set_seed,                   1},
    {"_ClusterR_silhouette_matrix",         (DL_FUNC) &_ClusterR_silhouette_matrix,          4},
    {"_ClusterR_SILHOUETTE_metric",         (DL_FUNC) &_ClusterR_SILHOUETTE_metric,          4},
    {"_ClusterR_split_rcpp_lst",            (DL_FUNC) &_ClusterR_split_rcpp_lst,             1},
    {"_ClusterR_squared_norm",              (DL_FUNC) &_ClusterR_squared_norm,               1},
    {"_ClusterR_subset_vec",                (DL_FUNC) &_ClusterR_subset_vec,                 2},
    {"_ClusterR_tot_ss_data",               (DL_FUNC) &_ClusterR_tot_ss_data,                1},
    {"_ClusterR_validate_centroids",        (DL_FUNC) &_ClusterR_validate_centroids,         2},
    {"_ClusterR_WCSS",                      (DL_FUNC) &_ClusterR_WCSS,                       2},
    {NULL, NULL, 0}
};

void R_init_ClusterR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}