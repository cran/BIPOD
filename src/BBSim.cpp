#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
extern "C" SEXP BBSim(SEXP start,
		      SEXP end,
		      SEXP n,
		      SEXP Sigma,
		      SEXP T,
		      SEXP t0,
		      SEXP seed,
		      SEXP p
		      ) {
  try{
    arma::mat start_  = Rcpp::as<arma::mat>(start);
    arma::mat end_    = Rcpp::as<arma::mat>(end);
    arma::mat Sigma_  = Rcpp::as<arma::mat>(Sigma);
    double n_         = Rcpp::as<double>(n);
    double T_         = Rcpp::as<double>(T);
    double t0_        = Rcpp::as<double>(t0);
    int seed_         = Rcpp::as<int>(seed);
    int p_            = Rcpp::as<int>(p);

    double dt=(T_-t0_)/(n_-1);
    mat Mean(1,p_);
    mat SD(p_,p_);
    mat BB(n_,p_);
    BB.row(0) = start_;
    //srand(seed_);
    if(seed_>0){arma_rng::set_seed(seed_);}
    if(seed_<=0){arma_rng::set_seed_random();}
    mat WWW = randn<mat>(n_-1,p_);

    for(int i=1; i<n_; ++i){
      Mean = BB.row(i-1)+1/(n_-i)*(end_-BB.row(i-1));
      SD   = sqrt((n_-i-1)*dt/(n_-i))*Sigma_;
      BB.row(i) = (WWW.row(i-1)*SD)+Mean;
    }

    //    return Rcpp::List::create(Rcpp::Named("BB") = BB);

			      return Rcpp::wrap(BB);


    // Below possible exceptions are taken care off
  }
  catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)");
  }
  // return to R
  return R_NilValue;
}

