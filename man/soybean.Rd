\name{soybean}
\alias{soybean}
\docType{data}
\title{
The soybean (large) data set from the UCI repository
}
\description{
There are 19 classes, only the first 15 of which have been used in prior work. The folklore seems to be that the last four classes are unjustified by the data since they have so few examples. There are 35 categorical attributes, some nominal and some ordered. The value 'dna' means does not apply. The values for attributes are encoded numerically, with the first value encoded as '0', the second as '1', and so forth. Unknown values were imputated using the mice package.
}
\usage{data(soybean)}
\format{
  A data frame with 307 Instances and 36 attributes (including the class attribute, "class")
}
\details{
The column names of the data (including the class) appear in the following order: 

date, plant-stand, precip, temp, hail, crop-hist, area-damaged, severity, seed-tmt, germination, plant-growth, 
leaves, leafspots-halo, leafspots-marg, leafspot-size, leaf-shread, leaf-malf, leaf-mild, stem, lodging, 
stem-cankers, canker-lesion, fruiting-bodies, external decay, mycelium, int-discolor, sclerotia, fruit-pods, 
fruit spots, seed, mold-growth, seed-discolor, seed-size, shriveling, roots, class
}
\references{
R.S. Michalski and R.L. Chilausky, Learning by Being Told and Learning from Examples: An Experimental Comparison of the Two Methods of Knowledge Acquisition in the Context of Developing an Expert System for Soybean Disease Diagnosis, International Journal of Policy Analysis and Information Systems, Vol. 4, No. 2, 1980. 

Donor: Ming Tan & Jeff Schlimmer (Jeff.Schlimmer cs.cmu.edu)

download source: https://archive.ics.uci.edu/ml/datasets/Soybean+(Large)
}
\examples{

data(soybean)

X = soybean[, -ncol(soybean)]

y = soybean[, ncol(soybean)]
}
\keyword{datasets}
