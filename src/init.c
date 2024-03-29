#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ClusterR_affinity_propagation(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_check_NaN_Inf(void *);
extern SEXP _ClusterR_ClaraMedoids(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_ClusterMedoids(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_cost_clusters_from_dis_meds(void *, void *);
extern SEXP _ClusterR_dissim_mat(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_dissim_MEDOIDS(void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_evaluation_rcpp(void *, void *, void *);
extern SEXP _ClusterR_GMM_arma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_GMM_arma_AIC_BIC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_KMEANS_arma(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_KMEANS_rcpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_mini_batch_kmeans(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_opt_clust_fK(void *, void *, void *);
extern SEXP _ClusterR_OptClust(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_predict_medoids(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _ClusterR_predict_MGausDPDF(void *, void *, void *, void *, void *);
extern SEXP _ClusterR_Predict_mini_batch_kmeans(void *, void *, void *, void *);
extern SEXP _ClusterR_preferenceRange(void *, void *, void *);
extern SEXP _ClusterR_SCALE(void *, void *, void *);
extern SEXP _ClusterR_silhouette_clusters(void *, void *);
extern SEXP _ClusterR_split_rcpp_lst(void *);
extern SEXP _ClusterR_validate_centroids(void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_ClusterR_affinity_propagation",        (DL_FUNC) &_ClusterR_affinity_propagation,         9},
    {"_ClusterR_check_NaN_Inf",               (DL_FUNC) &_ClusterR_check_NaN_Inf,                1},
    {"_ClusterR_ClaraMedoids",                (DL_FUNC) &_ClusterR_ClaraMedoids,                11},
    {"_ClusterR_ClusterMedoids",              (DL_FUNC) &_ClusterR_ClusterMedoids,               9},
    {"_ClusterR_cost_clusters_from_dis_meds", (DL_FUNC) &_ClusterR_cost_clusters_from_dis_meds,  2},
    {"_ClusterR_dissim_mat",                  (DL_FUNC) &_ClusterR_dissim_mat,                   7},
    {"_ClusterR_dissim_MEDOIDS",              (DL_FUNC) &_ClusterR_dissim_MEDOIDS,               6},
    {"_ClusterR_evaluation_rcpp",             (DL_FUNC) &_ClusterR_evaluation_rcpp,              3},
    {"_ClusterR_GMM_arma",                    (DL_FUNC) &_ClusterR_GMM_arma,                    10},
    {"_ClusterR_GMM_arma_AIC_BIC",            (DL_FUNC) &_ClusterR_GMM_arma_AIC_BIC,            10},
    {"_ClusterR_KMEANS_arma",                 (DL_FUNC) &_ClusterR_KMEANS_arma,                  7},
    {"_ClusterR_KMEANS_rcpp",                 (DL_FUNC) &_ClusterR_KMEANS_rcpp,                 12},
    {"_ClusterR_mini_batch_kmeans",           (DL_FUNC) &_ClusterR_mini_batch_kmeans,           13},
    {"_ClusterR_opt_clust_fK",                (DL_FUNC) &_ClusterR_opt_clust_fK,                 3},
    {"_ClusterR_OptClust",                    (DL_FUNC) &_ClusterR_OptClust,                    12},
    {"_ClusterR_predict_medoids",             (DL_FUNC) &_ClusterR_predict_medoids,              7},
    {"_ClusterR_predict_MGausDPDF",           (DL_FUNC) &_ClusterR_predict_MGausDPDF,            5},
    {"_ClusterR_Predict_mini_batch_kmeans",   (DL_FUNC) &_ClusterR_Predict_mini_batch_kmeans,    4},
    {"_ClusterR_preferenceRange",             (DL_FUNC) &_ClusterR_preferenceRange,              3},
    {"_ClusterR_SCALE",                       (DL_FUNC) &_ClusterR_SCALE,                        3},
    {"_ClusterR_silhouette_clusters",         (DL_FUNC) &_ClusterR_silhouette_clusters,          2},
    {"_ClusterR_split_rcpp_lst",              (DL_FUNC) &_ClusterR_split_rcpp_lst,               1},
    {"_ClusterR_validate_centroids",          (DL_FUNC) &_ClusterR_validate_centroids,           5},
    {NULL, NULL, 0}
};

void R_init_ClusterR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
