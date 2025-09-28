FROM rocker/rstudio:devel
LABEL maintainer='Lampros Mouselimis'

RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update && \
 apt-get install -y \
    libfftw3-dev \
    libgmp3-dev \
    git-core \
    pandoc \
    libpng-dev \
    make \
    libcurl4-openssl-dev \
    libssl-dev \
    sudo \
    libarmadillo-dev \
    libblas-dev \
    liblapack-dev \
    libarpack++2-dev \
    gfortran \
    libtiff5-dev \
    libxml2-dev \
    libssh2-1-dev \
    zlib1g-dev && \
 R -e "install.packages('devtools', dependencies = TRUE, repos = 'https://cloud.r-project.org/')" && \
 R -e "install.packages(c( 'Rcpp', 'graphics', 'grDevices', 'utils', 'stats', 'gmp', 'ggplot2', 'lifecycle', 'RcppArmadillo', 'OpenImageR', 'FD', 'testthat', 'covr', 'knitr', 'rmarkdown', 'remotes' ), repos =  'https://cloud.r-project.org/' )" && \
 apt-get autoremove -y && \
 apt-get clean && \
 rm -rf /var/lib/apt/lists/*

ADD http://www.random.org/strings/?num=10&len=8&digits=on&upperalpha=on&loweralpha=on&unique=on&format=plain&rnd=new uuid

ARG BUILD_DATE
RUN echo "$BUILD_DATE"

RUN R -e "remotes::install_github('mlampros/ClusterR', upgrade = 'never', dependencies = FALSE, repos = 'https://cloud.r-project.org/')"

ENV USER=rstudio
