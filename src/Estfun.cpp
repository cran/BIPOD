#include <RcppArmadillo.h>

using namespace std;
using namespace arma;

// Declare the BIPOD class
class BIPOD {
public:
  mat Drift;
  mat Diff;
  Mat<int> PathAcc;
  Mat<int> PathAccPoints;
  Mat<int> DiffAcc;
  Mat<int> DriftAcc;
  mat LatentPath;
  bool flag;
};

#ifndef __X_H_INCLUDED__   // if x.h hasn't been included yet...
#define __X_H_INCLUDED__   //   #define this so the compiler knows it has been included
#include "OUdrift.h"
#include "CIRdrift.h"
#include "FHNdrift.h"
#include "FHN5drift.h"
#include "driftfun.h"
#include "Radon.h"
//#include "DiffSimInt.h"
#include "BBSimInt.h"
#include "PostDistDrift.h"
#include "MvrNorm.h"
#include "maincode.h"
#endif

extern "C" SEXP Estfun(SEXP Delta,
		       SEXP data,
		       SEXP ImputeN,
		       SEXP GibbsN,
		       SEXP seed,
		       SEXP Start,
		       SEXP parMat,
		       SEXP diffPriorMean,
		       SEXP diffPriorCovar,
		       SEXP diffRW,
		       SEXP LatentPathStart,
		       SEXP Model,
		       SEXP driftPriorMean,
		       SEXP driftPriorCovar,
		       SEXP driftRW,
		       SEXP MHdrift,
		       SEXP LatentMeanY0,
		       SEXP LatentVarY0,
		       SEXP RWrhoPaths,
		       SEXP RWrho2PathPoints){
  try{
    int M   = Rcpp::as<int>(ImputeN);
    int seed_   = Rcpp::as<int>(seed);
    int GibbsN_ = Rcpp::as<int>(GibbsN);
    double RWrhoPaths_ = Rcpp::as<double>(RWrhoPaths);
    double RWrho2PathPoints_ = Rcpp::as<double>(RWrho2PathPoints);
    double Delta_   = Rcpp::as<double>(Delta);
    double LatentMeanY0_ = Rcpp::as<double>(LatentMeanY0);
    double LatentVarY0_ = Rcpp::as<double>(LatentVarY0);
    mat V      = Rcpp::as<mat>(data);
    mat Start_  = Rcpp::as<mat>(Start);
    mat parMat_  = Rcpp::as<mat>(parMat);
    mat diffPriorMean_  = Rcpp::as<mat>(diffPriorMean);
    mat diffPriorCovar_  = Rcpp::as<mat>(diffPriorCovar);
    mat diffRW_  = Rcpp::as<mat>(diffRW);
    mat LatentPathStart_  = Rcpp::as<mat>(LatentPathStart);
    string Model_  = Rcpp::as<string>(Model);
    int MHdrift_  = Rcpp::as<int>(MHdrift);
    mat driftPriorMean_  = Rcpp::as<mat>(driftPriorMean);
    mat driftPriorCovar_  = Rcpp::as<mat>(driftPriorCovar);
    mat driftRW_  = Rcpp::as<mat>(driftRW);

    //srand(seed_);
    if(seed_>0){arma_rng::set_seed(seed_);}
    if(seed_<=0){arma_rng::set_seed_random();}

    BIPOD driftEsts;
    driftEsts = maincode(V,
		        Delta_,
			 Start_,
			 M,
			 GibbsN_,
			 parMat_,
			 diffPriorMean_,
			 diffPriorCovar_,
			 diffRW_,
			 LatentPathStart_,
			 Model_,
			 driftPriorMean_,
			 driftPriorCovar_,
			 driftRW_,
			 MHdrift_,
			 LatentMeanY0_,
			 LatentVarY0_,
			 RWrhoPaths_,
			 RWrho2PathPoints_);
    return Rcpp::List::create(
			      Rcpp::Named("Drift") = driftEsts.Drift,
			      Rcpp::Named("Diff") = driftEsts.Diff,
			      Rcpp::Named("PathAcc") = driftEsts.PathAcc,
			      Rcpp::Named("PathAccPoints") = driftEsts.PathAccPoints,
			      Rcpp::Named("flag") = driftEsts.flag,
			      Rcpp::Named("LatentPath") = driftEsts.LatentPath,
			      Rcpp::Named("DiffAcc") = driftEsts.DiffAcc,
			      Rcpp::Named("DriftAcc") = driftEsts.DriftAcc);
  }
  catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)");
  }
  // return to R
  return R_NilValue;
}

