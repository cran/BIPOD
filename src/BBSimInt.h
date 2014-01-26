mat BBSimInt(mat start,
	     mat end,
	     double n,
	     double sigma1,
	     double sigma2,
	     double T,
	     double t0){
  mat Sigma;
  Sigma << sigma1 << 0 << endr
	<< 0 << sigma2 << endr;

  double dt =  (T-t0)/(n-1);
  mat Mean(1,2);
  mat SD(2,2);
  mat BB(n,2);
  BB.row(0) = start;
  mat WWW = randn<mat>(n-1,2);
  for(int i=1; i<n; ++i){
    Mean = BB.row(i-1)+1/(n-i)*(end-BB.row(i-1));
    SD   = sqrt((n-i-1)*dt/(n-i))*Sigma;
    BB.row(i) = (WWW.row(i-1)*SD)+Mean;
  }

  return(BB);
}
