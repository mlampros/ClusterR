FROM rocker/rstudio:devel

LABEL maintainer='Lampros Mouselimis'

RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update && \
 apt-get install -y libfftw3-dev libgmp-dev git-core pandoc pandoc-citeproc libpng-dev make libcurl4-openssl-dev libssl-dev && \
 apt-get install -y sudo && \
 apt-get install -y libarmadillo-dev && \
 apt-get install -y libblas-dev && \
 apt-get install -y liblapack-dev && \
 apt-get install -y libarpack++2-dev && \
 apt-get install -y gfortran && \
 apt-get install -y libgmp3-dev && \
 apt-get install -y libfftw3-dev && \
 apt-get install -y libtiff5-dev  && \
 apt-get install -y libxml2-dev && \
 apt-get install -y libssh2-1-dev && \
 apt-get install -y zlib1g-dev && \
 R -e "install.packages('devtools', dependencies = TRUE, repos = 'https://cloud.r-project.org/')" && \
 R -e "install.packages(c( 'Rcpp', 'graphics', 'grDevices', 'utils', 'stats', 'gmp', 'ggplot2', 'lifecycle', 'RcppArmadillo', 'OpenImageR', 'FD', 'testthat', 'covr', 'knitr', 'rmarkdown', 'remotes' ), repos =  'https://cloud.r-project.org/' )"

ADD http://www.random.org/strings/?num=10&len=8&digits=on&upperalpha=on&loweralpha=on&unique=on&format=plain&rnd=new uuid
ARG BUILD_DATE

RUN echo "$BUILD_DATE"
RUN R -e "remotes::install_github('mlampros/ClusterR', upgrade = 'never', dependencies = FALSE, repos = 'https://cloud.r-project.org/')" && \
 apt-get autoremove -y && \
 apt-get clean

ENV USER rstudio
