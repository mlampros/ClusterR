FROM rocker/rstudio 

 
LABEL maintainer='Lampros Mouselimis' 

 
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update && \ 
 apt-get install -y libfftw3-dev libgmp-dev git-core pandoc pandoc-citeproc libpng-dev libssl-dev make libcurl4-openssl-dev && \ 
 apt-get install -y sudo && \ 
 apt-get install -y libgfortran3 && \ 
 apt-get install -y libarmadillo-dev && \ 
 apt-get install -y libblas-dev && \ 
 apt-get install -y liblapack-dev && \ 
 apt-get install -y libarpack++2-dev && \ 
 apt-get install -y gfortran && \ 
 apt-get install -y libgmp3-dev && \ 
 apt-get install -y libfftw3-dev && \ 
 apt-get install -y libtiff5-dev  && \ 
 sudo ln -s /usr/lib/x86_64-linux-gnu/libgfortran.so.3 /usr/lib/libgfortran.so && \ 
 apt-get -qq update && \ 
 R -e "install.packages(c( 'Rcpp', 'graphics', 'grDevices', 'utils', 'gmp', 'FD', 'stats', 'ggplot2', 'RcppArmadillo', 'OpenImageR', 'testthat', 'covr', 'knitr', 'rmarkdown', 'remotes' ), repos =  'https://cloud.r-project.org/' )" && \ 
 R -e "remotes::install_github('mlampros/ClusterR', upgrade = 'always', dependencies = TRUE, repos = 'https://cloud.r-project.org/')" && \ 
 apt-get autoremove -y && \ 
 apt-get clean 

 
ENV USER rstudio 

 
