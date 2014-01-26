mat PostDistDrift(mat Z,
		  double Delta,
		  double sigma1,
		  double sigma2,
		  int parlen,
		  string Model){
  int n = Z.n_rows;
  mat R(parlen+1,parlen+1);
  mat one(1,1);
  one.fill(1.0);
  R.fill(0);
  mat I(parlen+1,1);
  I.fill(0);
  mat F(parlen,1);
  F.fill(0);
  mat Gamma(2,2);
  Gamma << 1.0/pow(sigma1,2) << 0 << endr
	<< 0 << 1.0/pow(sigma2,2) << endr;
  cube f (n-1,2,parlen+1);
  cube Gam2(n-1,2,1);
  Gam2.fill(0);
  if(Model=="OU"){
    for(int i=0; i<(n-1); ++i){
      f.subcube(i,0,0,i,1,parlen) = trans(OUdrift(Z(i,0),Z(i,1)));
    }
  }
  else if(Model=="CIR"){
    for(int i=0; i<(n-1); ++i){
      f.subcube(i,0,0,i,1,parlen) = trans(CIRdrift(Z(i,0),Z(i,1)));
      Gam2.subcube(i,0,0,i,0,0) = one*1.0/(pow(sigma1,2)*Z(i,0));
      Gam2.subcube(i,1,0,i,1,0) = one*1.0/(pow(sigma2,2)*Z(i,1));
    }
  }
  else if(Model=="FHN"){
    for(int i=0; i<(n-1); ++i){
      f.subcube(i,0,0,i,1,parlen) = trans(FHNdrift(Z(i,0),Z(i,1)));
    }
  }
  else if(Model=="FHN5"){
    for(int i=0; i<(n-1); ++i){
      f.subcube(i,0,0,i,1,parlen) = trans(FHN5drift(Z(i,0),Z(i,1)));
    }
  }
  //  else cout << "ERROR: Model wrong!";

  if(Model=="CIR"){
    for(int i=0; i<=parlen; ++i){
      for(int j=0; j<=i; ++j){
	R(i,j) = Delta*accu(f.slice(i)%Gam2.slice(0)%f.slice(j));
      }
      if(i>0){
	I(i) = accu(f.slice(i)%Gam2.slice(0)%(Z.rows(1,n-1)-Z.rows(0,n-2)));
	F(i-1)=I(i)-R(i,0);
      }
    }
    //    R.fill(1);
    //    I.fill(1);
    //	F.fill(1);
  }
  else{
    for(int i=0; i<=parlen; ++i){
      for(int j=0; j<=i; ++j){
	R(i,j) = Delta*accu((f.slice(i)*Gamma)%f.slice(j));
      }

      if(i>0){
	I(i) = accu((f.slice(i)*Gamma)%(Z.rows(1,n-1)-Z.rows(0,n-2)));
	F(i-1)=I(i)-R(i,0);
      }
    }
  }
  mat Var = inv(symmatl(R.submat(1,1,parlen,parlen)));
  mat Mean = Var*F;
  mat out(parlen,parlen+1);
  out.col(0) = Mean;
  out.submat(0,1,parlen-1,parlen) = Var;
  return(out);
}
