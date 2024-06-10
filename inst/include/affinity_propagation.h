
/**
 *
 * Rcpp-implementation of the "affinity propagation clustering" algorithm and the "preferenceRange".
 *
 * See, "BJ Frey and D Dueck, Science 315, 972-976, Feb 16, 2007", for a description of the algorithm.
 * 
 * See, "https://www.psi.toronto.edu/index.php?q=affinity%20propagation" for more details [ especially for the required matlab files ]
 *
 * Copyright 2007, BJ Frey and Delbert Dueck. This software may be freely used and distributed for non-commercial purposes.
 *
 */


#ifndef __AffinityPropagation__
#define __AffinityPropagation__

#ifdef _OPENMP
#include <omp.h>
#endif

# include <limits>                                  // minimum, maximum values for double in c++
# include <unordered_map>
#include <iostream>
# include <cmath>



//****************************************************************************************************************************************** .h file


class Affinity_Propagation {

  public :

    Affinity_Propagation() {}

    void set_seed(int seed);
    arma::colvec max_min_col_items(arma::colvec &x, int y, bool maximum);
    int modulus (int a, int b);
    arma::uvec matlab_setdiff(arma::uvec x, arma::uvec y);
    Rcpp::List affinity_propagation(arma::mat &s, std::vector<double> p, int maxits, int convits, double dampfact,
                                    bool details, double nonoise, double eps, bool time);
    std::vector<double> preferenceRange(arma::mat &s, std::string method, int threads);

    ~Affinity_Propagation() {}
};


#endif              // __AffinityPropagation__

//****************************************************************************************************************************************** .cpp file


//---------------------------------------------------------------------------------
// use base R's set.seed() in Rcpp for RNG
// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case
//---------------------------------------------------------------------------------

void Affinity_Propagation::set_seed(int seed) {
  if (Rcpp::all(Rcpp::is_na(Rcpp::IntegerVector{seed}))) return;
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}


//-----------------------------------------------------------------------------
// this function if maximum = TRUE is equivalent to the matlab's max(vector, 0)
// and if maximum = FALSE is equivalent to the matlab's min(vector, 0)
//-----------------------------------------------------------------------------

arma::colvec Affinity_Propagation::max_min_col_items(arma::colvec &x, int y, bool maximum = false) {
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if (maximum) {
      if (y > x[i]) {      // maximum
        x[i] = y;
      }
    }
    else {                 // minimum
      if (y < x[i]) {
        x[i] = y;
      }
    }
  }
  return x;
}


//------------------
// remainder for int
//------------------

int Affinity_Propagation::modulus (int a, int b) {

  return(a % b);
}


//-----------------------------------------------------------------------------------------
// matlab's setdiff() function in rcpp-armadillo ( find the values in A that are not in B )
//-----------------------------------------------------------------------------------------

arma::uvec Affinity_Propagation::matlab_setdiff(arma::uvec x, arma::uvec y) {

  arma::uvec out;
  int increment = 0;

  for (unsigned int i = 0; i < x.n_elem; i++) {
    arma::uvec tmp = arma::find(x(i) == y);
    if (tmp.is_empty()) {
      increment++;
      out.resize(increment);
      out(increment - 1) = x(i);
    }
  }

  return arma::unique(out);
}


//-------------------------------
// affinity propagation algorithm
//-------------------------------

Rcpp::List Affinity_Propagation::affinity_propagation(arma::mat &s, std::vector<double> p, int maxits = 1000, int convits = 100, double dampfact = 0.9,
                                                      bool details = false, double nonoise = 0.0, double eps = 2.2204e-16, bool time = false) {

  arma::wall_clock timer;
  if (time || details) {
    timer.tic();
  }

  if (maxits <= 0) Rcpp::stop("the 'maxits' parameter must be a positive integer!\n");
  if (convits <= 0) Rcpp::stop("the 'convits' parameter must be a positive integer!\n");
  if ((dampfact < 0.5) || (dampfact >= 1.0)) Rcpp::stop("the 'dampfact' parameter must be >= 0.5 and < 1\n");
  if (dampfact > 0.9) Rcpp::Rcout << "\n*** Warning: Large damping factor in use. Turn on plotting\n" << "    to monitor the net similarity. The algorithm will\n" <<
    "    change decisions slowly, so consider using a larger value\n" << "    of convits.\n\n" << std::endl;

  //-----------------------------------------------------
  // Check that standard arguments are consistent in size
  //-----------------------------------------------------

  unsigned int N;
  unsigned int rows_3col = 0;                                                                       // initialize the 'rows_3col' to 0. By setting it to 'arma::datum::inf' it gives "runtime error: inf is outside the range of representable values of type 'unsigned int'" ( clang-UBSAN )

  if (s.n_cols == 3 && s.n_rows != 3) {                                                             // 3-column input matrix
    rows_3col = std::sqrt(s.n_rows);                                                                // rows should be equal to the square-root
    if (p.size() == 1) {
      N = rows_3col;
    }
    else {
      N = p.size();
    }
    if (rows_3col > N) {
      Rcpp::stop("data point index exceeds number of data points");
    }
  }
  else if (s.n_rows == s.n_cols) {                                                                  // similarity matrix
    N = s.n_rows;
    if ((p.size() != N) && (p.size() != 1)) {
      Rcpp::stop("p should be scalar or a vector of size N");
    }
  }
  else {
    Rcpp::stop("s must have 3 columns or be a square matrix!");
  }

  //----------------------------
  // Construct similarity matrix
  //----------------------------

  arma::mat S;
  bool symmetric = false;

  if (s.n_cols == 3 && s.n_rows != 3) {                   // create the similarity matrix from (i, j, similarity) 3-column data
    S.set_size(rows_3col, rows_3col);
    S.fill(-arma::datum::inf);
    for (unsigned int j = 0; j < s.n_rows; j++) {
      S(s(j,0),s(j,1)) = s(j,2);
    }
  }
  else {
    S = s;
  }
  if (S.is_symmetric()) {
    symmetric = true;
  }


  //---------------------------------------------------------------------------------
  // In case a user did not remove degeneracies from the input similarities, avoid
  // degenerate solutions by adding a small amount of noise to the input similarities
  //---------------------------------------------------------------------------------

  double realmin_ = std::numeric_limits<double>::min();
  double realmax_ = std::numeric_limits<double>::max();

  if (nonoise == 0.0) {
    set_seed(0);
    S = S + (eps * S + realmin_ * 100.0) % arma::randn(N,N);
  }

  //---------------------------------------
  // Place preferences on the diagonal of S
  //---------------------------------------

  if (p.size() == 1) {
    for (unsigned int i = 0; i < N; i++) {
      S(i,i) = p[0];                                                               // p is an std::vector<double> of length 1
    }
  }
  else {
    for (unsigned int i = 0; i < N; i++) {
      S(i,i) = p[i];                                                               // p is an std::vector<double> equal to the number of rows of the similarity matrix
    }
  }

  //--------------------------------------------------
  // Numerical stability -- replace -INF with -realmax
  //--------------------------------------------------

  arma::uvec n = arma::find(S < -realmax_);
  if (!n.is_empty()) {
    Rcpp::Rcout << "Warning: -INF similarities detected; changing to -REALMAX to ensure numerical stability!\n" << std::endl;
    S(n).fill(-realmax_);
  }

  arma::uvec nn = arma::find( S > realmax_, 1);
  if (!nn.is_empty()) {
    Rcpp::stop("+INF similarities detected; change to a large positive value (but smaller than +REALMAX)");
  }

  //---------------------------------
  // Allocate space for messages, etc
  //---------------------------------

  arma::colvec dS = S.diag();
  arma::mat A = arma::zeros(N,N);
  arma::mat R = arma::zeros(N,N);
  arma::mat netsim, idx, dpsim, expref;

  if (details) {
    idx.set_size(N,maxits + 1);
    idx.fill(0.0);
    netsim.set_size(1, maxits + 1);
    netsim.fill(0.0);
    dpsim.set_size(1, maxits + 1);
    dpsim.fill(0.0);
    expref.set_size(1, maxits + 1);
    expref.fill(0.0);
  }

  //----------------------------------------------
  // Execute parallel affinity propagation updates
  //----------------------------------------------

  arma::mat e = arma::zeros(N, convits);
  bool dn = false;
  int i = 0;
  arma::mat ST;

  if (symmetric) {                                                                 // saves memory if it's symmetric
    ST = S;
  }
  else {
    ST = S.t();
  }

  int tmpdpsim = 0;                                                                // initialize 'tmpdpsim' to 0 . By setting this to 'arma::datum::nan' it gives "runtime error: nan is outside the range of representable values of type 'int'" ( clang-UBSAN )
  int tmpnetsim,tmpexpref;
  arma::uvec tmpidx;
  bool unconverged = true;                                                         // unconverged set initially to 'true' [ see line 478 ]

  while (!dn) {

    i = i + 1;

    //-------------------------
    // Compute responsibilities
    //-------------------------

    A = A.t();
    R = R.t();

    for (unsigned int ii = 0; ii < N; ii++) {

      arma::colvec old = R.col(ii);
      arma::colvec AS = A.col(ii) + ST.col(ii);

      double Y = arma::max(AS);
      int I = arma::index_max(AS);

      AS(I) = -arma::datum::inf;                                    // set the 1st max. value to -inf

      double Y2 = arma::max(AS);

      R.col(ii) = ST.col(ii) - Y;                                   // add the first max. value to ST
      R(I,ii) = ST(I,ii) - Y2;
      R.col(ii) = (1.0 - dampfact) * R.col(ii) + dampfact * old;    // Damping
      arma::uvec idx_in = arma::find(R.col(ii) > realmax_);

      if (!idx_in.is_empty()) {                                     // continue only if uvec is not-empty
        for (unsigned int f = 0; f < idx_in.n_elem; f++) {
          R(f,ii) = realmax_;
        }
      }
    }

    //-----------------------
    // Compute availabilities
    //-----------------------

    A = A.t();
    R = R.t();

    for (unsigned int jj = 0; jj < N; jj++) {
      arma::colvec old = A.col(jj);
      arma::colvec R_jj_in = R.col(jj);
      arma::colvec Rp = max_min_col_items(R_jj_in, 0.0, true);              // max(vector, 0)

      Rp(jj) = R(jj,jj);
      A.col(jj) = arma::accu(Rp) - Rp;
      double dA = A(jj,jj);
      arma::colvec A_jj_in = A.col(jj);
      A.col(jj) = max_min_col_items(A_jj_in, 0.0, false);                   // min(vector, 0)
      A(jj,jj) = dA;
      A.col(jj) = (1.0 - dampfact) * A.col(jj) + dampfact * old;            // Damping
    }


    //----------------------
    // Check for convergence
    //----------------------

    arma::colvec E = arma::conv_to<arma::colvec>::from((A.diag() + R.diag()) > 0);
    e.col(modulus(i-1, convits)) = E;                                                  // I removed the + 1     [ dif in indexing to matlab ]
    int K = arma::accu(E);                                                             // K parameter is specific to the while inner loop

    arma::colvec se;
    if ((i >= convits) || (i >= maxits)) {
      se = arma::sum(e,1);
      unconverged = (arma::accu((se == convits) + (se == 0)) != N);
      if ((!unconverged && (K > 0)) || (i == maxits)) {
        dn = true;
      }
    }


    //-----------------------------------------------------
    // Handle plotting and storage of details, if requested   [ plotting not supported in Rcpp ]
    //-----------------------------------------------------

    if (details) {

      if (K==0) {
        tmpnetsim = 0;                                               // initialize this to 0. By setting it to 'arma::datum::nan' it gives "runtime error: nan is outside the range of representable values of type 'int'" ( clang-UBSAN )
        tmpdpsim = 0;                                                // same as line 352
        tmpexpref = 0;                                               // same as line 352
        tmpidx.set_size(N);
        tmpidx.fill(arma::datum::nan);                               // 'tmpidx' can take 'arma::datum::nan' because it is initialized as 'arma::uvec' ( of type double )
      }
      else {
        arma::uvec I = arma::find(E == 1);                           // 'I' can be empty or having 1 or more items
        arma::uvec notI = arma::find(E != 1);
        arma::uvec c;
        if (!I.is_empty()) {
          c = arma::index_max(S.cols(I), 1);

          arma::uvec tmp_c = arma::unique(c);
          arma::uvec tmp_Ic = I(c);
          c(I) = arma::regspace<arma::uvec>( 0, 1, K-1 );
          tmpidx = I(c);
          arma::umat tmp_um(2, notI.n_elem);
          arma::uvec tmp_v_um = tmpidx(notI);
          for (unsigned int f = 0; f < notI.n_elem; f++) {
            tmp_um(0,f) = notI(f);
            tmp_um(1,f) = tmp_v_um(f);
          }
          arma::uvec sub2_idx;
          if (tmp_um.n_cols == 1) {                                                          // case where notI & tmpidx(notI) consist of 1 element
            arma::uword tmp_idx_um = arma::sub2ind(size(S), tmp_um(0,0), tmp_um(1,0));       // size(S) is (N,N) ; USE THE 'matrix_of_subscripts' to store the indices to resemble the matlab 'sub2ind' function
            sub2_idx.set_size(1);
            sub2_idx(0) = tmp_idx_um;
          }
          else {
            sub2_idx = arma::sub2ind(size(S), tmp_um);
          }
          tmpdpsim = arma::accu(S(sub2_idx));
        }
        tmpexpref = arma::accu(dS(I));
        tmpnetsim = tmpdpsim + tmpexpref;
      }
    }

    if (details) {
      netsim.col(i-1) = tmpnetsim;                                   // subtract 1 because i now is used for indexing [ applies to the next 2 lines too ]
      dpsim.col(i-1) = tmpdpsim;
      expref.col(i-1) = tmpexpref;

      arma::uvec tmp_unq = arma::unique(tmpidx);
      idx.col(i-1) = arma::conv_to<arma::colvec>::from(tmpidx);
    }
  }                                                                 // end of while loop

  arma::uvec I = arma::find((A.diag() + R.diag()) > 0);
  int K = I.n_elem;                                                 // Identify exemplars  [ the uvec 'I' contains the exemplars ]

  if (K > 0) {
    arma::uvec c = arma::index_max(S.cols(I), 1);
    c(I) = arma::regspace<arma::uvec>( 0, 1, K-1 );                 // Identify clusters [ converted 'K' to 'K-1' due to dif in indexing compared to matlab ]


    //------------------------------------------------------------------
    // Refine the final set of exemplars and clusters and return results
    //------------------------------------------------------------------

    for (int k = 0; k < K; k++) {
      arma::uvec ii = arma::find(c == k);
      arma::rowvec tmp_j = arma::sum(S(ii,ii), 0);
      int j = arma::index_max(tmp_j);
      I(k) = ii(j);
    }
    arma::uvec notI = matlab_setdiff(arma::regspace<arma::uvec>( 0, 1, N-1 ), I);   // converted 'N' to 'N-1' due to dif in indexing compared to matlab && I didn't use reshape() compared to initial code
    c = arma::index_max(S.cols(I), 1);
    c(I) = arma::regspace<arma::uvec>( 0, 1, K-1 );                                 // Identify clusters [ converted 'K' to 'K-1' due to dif in indexing compared to matlab ]
    tmpidx = I(c);
    arma::umat tmp_um(2, notI.n_elem);
    arma::uvec tmp_v_um = tmpidx(notI);
    for (unsigned int f = 0; f < notI.n_elem; f++) {
      tmp_um(0,f) = notI(f);
      tmp_um(1,f) = tmp_v_um(f);
    }
    arma::uvec sub2_idx;
    if (tmp_um.n_rows == 1) {                                                        // case where notI & tmpidx(notI) consist of 1 element
      arma::uword tmp_idx_um = arma::sub2ind(size(S), tmp_um(0,0), tmp_um(1,0));     // size(S) is (N,N) ; Use the 'matrix_of_subscripts' to store the indices to resemble the matlab 'sub2ind' function
      sub2_idx.set_size(1);
      sub2_idx(0) = tmp_idx_um;
    }
    else {
      sub2_idx = arma::sub2ind(size(S), tmp_um);
    }
    tmpdpsim = arma::accu(S(sub2_idx));
    tmpexpref = arma::accu(dS(I));
    tmpnetsim = tmpdpsim + tmpexpref;
  }
  else {
    tmpidx.set_size(N);
    tmpidx.fill(arma::datum::nan);                          // 'tmpidx' can take 'arma::datum::nan' because it is initialized as 'arma::uvec' ( of type double )        
    tmpnetsim = 0;                                          // initialize this to 0. By setting it to 'arma::datum::nan' it gives "runtime error: nan is outside the range of representable values of type 'int'" ( clang-UBSAN )
    tmpexpref = 0;                                          // same as line 445
  }

  if (details) {
    netsim(i) = tmpnetsim;                                  // here rather than i+1 I use i [ applies to the next lines too ]
    netsim = netsim.submat(0, 0, 0, i);
    dpsim(i) = tmpdpsim;
    dpsim = dpsim.submat(0, 0, 0, i);
    expref(i) = tmpexpref;
    expref=expref.submat(0, 0, 0, i);
    idx.col(i) = arma::conv_to<arma::colvec>::from(tmpidx);
    idx = idx.submat(0, 0, N-1, i);
  }
  else {
    netsim=tmpnetsim;
    dpsim=tmpdpsim;
    expref=tmpexpref;
    idx = arma::conv_to<arma::mat>::from(tmpidx);
  }

  double n_time = timer.toc();

  if (details) {
    Rprintf("\nNumber of exemplars identified: %d  (for %d data points)\n", K, N);
    Rprintf("Net similarity: %d\n", tmpnetsim);
    Rprintf("  Similarities of data points to exemplars: %g\n", dpsim(0,i));
    Rprintf("  Preferences of selected exemplars: %d\n", tmpexpref);
    Rprintf("Number of iterations: %d\n\n", i);
  }

  if (time || details) {
    Rprintf("Elapsed time: %g sec\n", n_time);
  }

  if (unconverged) {
    Rprintf("\n*** Warning: Algorithm did not converge. Activate plotting\n");
    Rprintf("    so that you can monitor the net similarity. Consider\n");
    Rprintf("    increasing maxits and convits, and, if oscillations occur\n");
    Rprintf("    also increasing dampfact.\n\n");
  }

  //----------------------------------
  // return clusters for each exemplar
  //----------------------------------

  std::unordered_map<int, std::vector<int> > clusters;

  for (unsigned int t = 0; t < N; t++) {
    if (details) {
      clusters[idx(t, i)].push_back(t);
    }
    else {
      clusters[idx(t,0)].push_back(t);
    }
  }

  //----------------------------------------------------------------------------------------
  // return the netsim, dpsim, expref as doubles && transpose the idx's (if details = FALSE)
  //----------------------------------------------------------------------------------------

  double netsim_, dpsim_, expref_;

  if (details) {
    netsim_ = netsim(0,i);
    dpsim_ = dpsim(0,i);
    expref_ = expref(0,i);
  }
  else {
    netsim_ = netsim(0,0);
    dpsim_ = dpsim(0,0);
    expref_ = expref(0,0);
    idx = idx.t();               // so that I get a row-vector
  }

  //------------------------------------
  // return the exemplars as std::vector
  //------------------------------------

  std::vector<int> exempl(K);   // K is the number of exemplars

  for (int s = 0; s < K; s++) {
    exempl[s] = I(s);
  }

  return Rcpp::List::create( Rcpp::Named("K") = K, Rcpp::Named("N") = N, Rcpp::Named("netsim") = netsim_, Rcpp::Named("dpsim") = dpsim_, Rcpp::Named("expref") = expref_,
                                         Rcpp::Named("iterations") = i, Rcpp::Named("exemplars") = exempl, Rcpp::Named("idx") = idx, Rcpp::Named("clusters") = clusters);
}



//----------------------------------------------------------------------
// inner function for the exact method of the 'preferenceRange' function
//----------------------------------------------------------------------

double inner_exact(int j21, int j22, arma::mat &S) {
  
  arma::uvec j_in(2);
  j_in(0) = j21;
  j_in(1) = j22;
  arma::mat mt_in = S.cols(j_in);
  arma::colvec tmp_j_in = arma::max(mt_in, 1);
  return arma::accu(tmp_j_in);
}



//---------------------------
// 'p' preference-range value
//---------------------------


std::vector<double> Affinity_Propagation::preferenceRange(arma::mat &s, 
                                                          std::string method = "bound", 
                                                          int threads = 1) {

  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif
  
  int N;
  arma::mat S;
  
  if (s.n_cols == 3 && s.n_rows != 3) {                                                             // 3-column input matrix
    N = std::sqrt(s.n_rows);                                                                        // rows should be equal to the square-root
  }
  else if (s.n_rows == s.n_cols) {
    N = s.n_rows;
  }
  else {
    Rcpp::stop("s must have 3 columns or be a square matrix!");
  }
  
  if (s.n_cols == 3 && s.n_rows != 3) {                   // create the similarity matrix from (i, j, similarity) data
    S.set_size(N, N);
    S.diag().zeros();
    S.fill(-arma::datum::inf);
    for (unsigned int j = 0; j < s.n_rows; j++) {
      S(s(j,0),s(j,1)) = s(j,2);
    }
  }
  else {
    S = s;
  }
  
  double pmax = 0.0;
  double pmin, tmp;                                       // initialize pmin, pmax, tmp
  arma::colvec m;
  
  arma::rowvec tmp_dpsim1 = arma::sum(S, 0);
  double dpsim1 = arma::max(tmp_dpsim1);
  if (dpsim1 == -arma::datum::inf) {
    pmin = arma::datum::nan;                              // 'pmin' can be set to 'arma::datum::nan' as it is of type double
  }
  else if (method == "bound") {
    for (int k = 0; k < N; k++) {
      S(k,k) = -arma::datum::inf;
    }
    m = arma::max(S, 1);
    tmp = arma::accu(m);
    double yy = arma::min(m);
    int ii = arma::index_min(m);                          // differences might appear if more than one vector-items have the same max. value
    double TMP = arma::datum::inf;
    for (int k = 0; k < ii - 1; k++) {
      if (m(k) < TMP) {
        TMP = m(k);
      }
    }
    for (int kk = (ii + 1); kk < N; kk++) {
      if (m(kk) < TMP) {
        TMP = m(kk);
      }
    }
    tmp = tmp - yy - TMP;
    pmin = dpsim1 - tmp;
  }
  else {
    
    double dpsim2 = -arma::datum::inf;
    
    int j21, j22;
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) shared(S, dpsim2, N) private(j21, j22)
    #endif
    for (j21 = 0; j21 < N-1; j21++) {
      for (j22 = j21+1; j22 < N; j22++) {
        
        double tmp = inner_exact(j21, j22, S);
        if (tmp > dpsim2) {
          
          #ifdef _OPENMP
          #pragma omp atomic write
          #endif
          dpsim2 = tmp;
        }
      }
    }
    pmin = dpsim1 - dpsim2;
  }
  for (int ff = 0; ff < N; ff++) {
    S(ff,ff) = -arma::datum::inf;
  }
  
  pmax = S.max();
  
  std::vector<double> min_max;
  min_max.push_back(pmin);
  min_max.push_back(pmax);
  
  return min_max;
}

