% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_functions.R
\name{predict_Medoids}
\alias{predict_Medoids}
\alias{predict.MedoidsCluster}
\title{Predictions for the Medoid functions}
\usage{
predict_Medoids(
  data,
  MEDOIDS = NULL,
  distance_metric = "euclidean",
  fuzzy = FALSE,
  minkowski_p = 1,
  threads = 1
)

\method{predict}{MedoidsCluster}(object, newdata, fuzzy = FALSE, threads = 1, ...)
}
\arguments{
\item{data}{matrix or data frame}

\item{MEDOIDS}{a matrix of initial cluster medoids (data observations). The rows of the MEDOIDS matrix should be equal to the number of clusters and the columns of the MEDOIDS matrix should be equal to the columns of the data.}

\item{distance_metric}{a string specifying the distance method. One of,  \emph{euclidean},  \emph{manhattan},  \emph{chebyshev},  \emph{canberra},  \emph{braycurtis},  \emph{pearson_correlation},  \emph{simple_matching_coefficient},  \emph{minkowski},  \emph{hamming},  \emph{jaccard_coefficient},  \emph{Rao_coefficient},  \emph{mahalanobis}, \emph{cosine}}

\item{fuzzy}{either TRUE or FALSE. If TRUE, then probabilities for each cluster will be returned based on the distance between observations and medoids.}

\item{minkowski_p}{a numeric value specifying the minkowski parameter in case that distance_metric = "minkowski"}

\item{threads}{an integer specifying the number of cores to run in parallel. Openmp will be utilized to parallelize the number of initializations (num_init)}

\item{object, newdata, ...}{arguments for the `predict` generic}
}
\value{
a list with the following attributes will be returned : clusters, fuzzy_clusters (if fuzzy = TRUE), dissimilarity.
}
\description{
Predictions for the Medoid functions
}
\examples{

data(dietary_survey_IBS)

dat = dietary_survey_IBS[, -ncol(dietary_survey_IBS)]

dat = center_scale(dat)

cm = Cluster_Medoids(dat, clusters = 3, distance_metric = 'euclidean', swap_phase = TRUE)

pm = predict_Medoids(dat, MEDOIDS = cm$medoids, 'euclidean', fuzzy = TRUE)
}
\author{
Lampros Mouselimis
}
