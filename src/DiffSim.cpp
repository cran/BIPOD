#include <RcppArmadillo.h>
using namespace std;
using namespace arma;

mat OUdrift2(double x, double y){
  mat fAll(7,2);
  fAll << 0 << 0 << endr
       << -x << 0 << endr
       << -y << 0 << endr
       << 0 << -x << endr
       << 0 << -y << endr
       << 1 << 0 << endr
       << 0 << 1 << endr;
  return(fAll);
}

mat CIRdrift2(double x, double y){
  mat fAll(7,2);
  fAll << 0 << 0 << endr
       << -x << 0 << endr
       << -y << 0 << endr
       << 0 << -x << endr
       << 0 << -y << endr
       << 1 << 0 << endr
       << 0 << 1 << endr;
  return(fAll);
}

mat FHNdrift2(double x, double y){
  mat fAll(5,2);
  fAll << 0 << -y << endr
       << x-pow(x,3)-y << 0 << endr
       << 1 << 0 << endr
       << 0 << x << endr
       << 0 << 1 << endr;
  return(fAll);
}


mat FHN5drift2(double x, double y){
  mat fAll(6,2);
  fAll << 0 << -y << endr
       << -pow(x,3) << 0 << endr
       << x-y << 0 << endr
       << 1 << 0 << endr
       << 0 << x << endr
       << 0 << 1 << endr;
  return(fAll);
}

mat driftfun2(double x,
	      double y,
	      mat driftpar,
	      string Model)
{
  mat BB(1,2);
  mat one(1,1);
  one.fill(1);
  if(Model=="OU"){
    BB = join_rows(one,driftpar)*OUdrift2(x,y);
  }
  if(Model=="CIR"){
    BB = join_rows(one,driftpar)*CIRdrift2(x,y);
  }
  else if(Model=="FHN"){
    BB = join_rows(one,driftpar)*FHNdrift2(x,y);
  }
  else if(Model=="FHN5"){
    BB = join_rows(one,driftpar)*FHN5drift2(x,y);
  }
  return trans(BB);}

cx_mat invComp(cx_mat A){
  /// Compute the inverse of two by two complex matrix.
  cx_double DetA(1);
  cx_mat B(2,2);
  DetA = A(0,0)*A(1,1)-A(0,1)*A(1,0);

  if((real(DetA)==0) && imag(DetA)==0){
    return(A);
  }
  else{
    B(0,0) = (cx_double) A(1,1)/DetA;
    B(0,1) = (cx_double) -A(0,1)/DetA;
    B(1,0) = (cx_double) -A(1,0)/DetA;
    B(1,1) = (cx_double) A(0,0)/DetA;
  return(B);
  }
  // if(1==3){
  //   mat Re_A = real(A);
  //   mat Im_A = imag(A);
  //   mat ResRe;
  //   mat ResIm;
  //   if(1==2){
  //     ResRe = inv(Re_A);
  //     ResIm = Re_A*0;
  //   }
  //   else{
  //     double detRe = det(Re_A);
  //     double detIm = det(Im_A);
  //     double tmp = accu(trans(Re_A)%inv(Im_A));
  //     double denum = pow(detRe-detIm,2)+pow(detIm*tmp,2);
  //     ResRe = (inv(Re_A)*(pow(detRe,2)-detRe*detIm)+pow(detIm,2)*inv(Im_A)*tmp)/denum;
  //     ResIm =
  // 	(inv(Im_A)*(detRe*detIm-pow(detIm,2))-detRe*detIm*tmp*inv(Re_A))/denum;
  //   }
  //   cx_mat B(ResRe,ResIm);
}

mat MvrNorm2(mat mean, mat covar){
  int p = covar.n_cols;
  mat A = chol(covar);
  mat Z = randn<mat>(p,1);
  return(mean+trans(A)*Z);
}

extern "C" SEXP DiffSim(SEXP n,
			SEXP start,
			SEXP Delta,
			SEXP driftpar,
			SEXP Sigma,
			SEXP seed,
			SEXP thin,
			SEXP Model){
  try{
    int seed_       = Rcpp::as<int>(seed);
    mat start_      = Rcpp::as<mat>(start);
    mat driftpar_   = Rcpp::as<mat>(driftpar);
    cx_mat Sigma_   = Rcpp::as<cx_mat>(Sigma);
    double n_       = Rcpp::as<double>(n);
    double thin_    = Rcpp::as<double>(thin);
    double Delta_   = Rcpp::as<double>(Delta);
    string Model_   = Rcpp::as<string>(Model);

    int nInt_       = Rcpp::as<int>(n);


    mat BBB(2,2);
    cx_vec BBBeigv;
    cx_mat BBBeigm;
    cx_mat CC2;
    cx_mat CC;
    cx_mat Mexp(2,2);
    cx_mat condVar(2,2);
    cx_mat expB(2,2);
    cx_mat MeanTerm(2,2);



    //    srand(seed_);
    if(seed_>0){arma_rng::set_seed(seed_);}
    if(seed_<=0){arma_rng::set_seed_random();}


    /*mean=0 and var=Delta. */
    mat B  = (0+sqrt(Delta_)*randn<mat>(nInt_-1,2));
    /* Automatically inserts 0 in Y*/
    mat Y(nInt_,2);
    mat Ythin( (int) (floor(n_/thin_)),2);
    mat Sigmatmp(2,2);

    Y.row(0) = start_;
    int exactOU=0;
    /* Drift matrix for OU */
    BBB << driftpar_(0,0) << driftpar_(0,1) << endr
	<< driftpar_(0,2) << driftpar_(0,3) << endr;
    if(Model_=="OU" && arma::rank(BBB) == 2){
      exactOU = 1;}

    if(Model_=="CIR"){
      for(int i=1; i<nInt_;++i){
	Sigmatmp << sqrt(Y(i-1,0)) << 0 << endr
		 << 0 << sqrt(Y(i-1,1)) << endr;
    	Y.row(i) = Y.row(i-1)+trans(driftfun2(Y(i-1,0),Y(i-1,1),driftpar_,Model_))*Delta_+B.row(i-1)*real(Sigma_)*Sigmatmp;
      }
    }
    else if(exactOU == 1){
      mat A(2,1);
      A << driftpar_(0,4) << endr
	<< driftpar_(0,5) << endr;

      eig_gen(BBBeigv,BBBeigm,BBB);

      /// Conditional variance: Var(Y_i+1|Y_i),
      CC = invComp(BBBeigm)*Sigma_;
      CC2 =CC*trans(CC);
      Mexp(0,0) = CC2(0,0)*(1.0-exp(-2.0*BBBeigv(0)*Delta_))/(2.0*BBBeigv(0));
      Mexp(0,1) =
	CC2(0,1)*(1.0-exp(-(BBBeigv(0)+BBBeigv(1))*Delta_))/(BBBeigv(0)+BBBeigv(1));
      Mexp(1,0) = Mexp(0,1);
      Mexp(1,1) = CC2(1,1)*(1.0-exp(-2.0*BBBeigv(1)*Delta_))/(2.0*BBBeigv(1));

      condVar = BBBeigm*Mexp*trans(BBBeigm);

      /// Matrix exponential exp(-B*Delta)
      expB = BBBeigm*diagmat(exp(-Delta_*BBBeigv))*invComp(BBBeigm);

      MeanTerm = BBBeigm*diagmat((1-exp(-Delta_*BBBeigv))/BBBeigv)*invComp(BBBeigm)*A;
      //      MeanTerm = BBBeigm*tmpmat*invComp(BBBeigm)*A;

      for(int i=1; i<nInt_;++i){
	Y.row(i) = trans(MvrNorm2(real(expB)*trans(Y.row(i-1))+real(MeanTerm),real(condVar)));
      }
      //      return Rcpp::wrap(Y);
    }
    else{
      if(Model_ == "OUa"){Model_ = "OU";}
      for(int i=1; i<nInt_;++i){
	Y.row(i) = Y.row(i-1)+trans(driftfun2(Y(i-1,0),Y(i-1,1),driftpar_,Model_))*Delta_+B.row(i-1)*real(Sigma_);
      }
    }

    if( ((int) thin_) == 1 || exactOU == 1){
      return Rcpp::wrap(Y);
    }
    else{
      for(int i=0; i< floor(n_/thin_); ++i){
	Ythin.row(i) = Y.row(i* (int) thin_);
      }
      return Rcpp::wrap(Ythin);
    }


    //    return Rcpp::List::create(Rcpp::Named("DiffSim") = FH);

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


