mat Radon(mat Z,
	  mat driftpar,
	  double sigma1,
	  double sigma2,
	  double Delta,
	  string Model){
  int n=Z.n_rows;
  double DeltaM = Delta/(n-1);
  /* mat SigmaSol(2,2); */
  /* SigmaSol << 1/sigma1 << 0 << endr */
  /* 	   << 0 << 1/sigma2 << endr; */
  mat Sigma(2,2);
  Sigma << sigma1 << 0 << endr
	<< 0 << sigma2 << endr;
  mat integ1(1,1);
  integ1.fill(0);
  mat integ2(1,1);
  integ2.fill(0);
  mat tmp(1,2);
  mat tmpZ(2,2);
  if(Model=="CIR"){
    for(int i=1; i<n; ++i){
      tmpZ << 1.0/Z(i-1,0) << 0 << endr
	   << 0 << 1.0/Z(i-1,1) << endr;

      tmp = 0.5*tmpZ*inv(Sigma*Sigma)*(4.0*driftfun(pow(sigma1*Z(i-1,0),2)/4.0,
						   pow(sigma2*Z(i-1,1),2)/4.0,
						   driftpar,
						   Model)-diagvec(Sigma*Sigma));
      integ1 = integ1 + (Z.row(i)-Z.row(i-1))*tmp;
      integ2 = integ2 + (trans(tmp)*tmp)*DeltaM;
    }
  }
  else{
    for(int i=1; i<n; ++i){
      tmp = inv(Sigma)*driftfun(sigma1*Z(i-1,0),
			      sigma2*Z(i-1,1),
			      driftpar,
			      Model);
      integ1 = integ1 + (Z.row(i)-Z.row(i-1))*tmp;
      integ2 = integ2 + (trans(tmp)*tmp)*DeltaM;
    }
  }
  return(integ1-.5*integ2);
}
