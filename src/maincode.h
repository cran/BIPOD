BIPOD maincode(mat V,
	       double Delta,
	       mat Start,
	       int M,
	       int GibbsN,
	       mat ParMat,
	       mat diffPriorMean,
	       mat diffPriorCovar,
	       mat diffRW,
	       mat LatentPathStart,
	       string Model,
	       mat driftPriorMean,
	       mat driftPriorCovar,
	       mat driftRW,
	       int MHdrift,
	       double LatentMeanY0,
	       double LatentVarY0,
	       double RWrhoPaths,
	       double RWrho2PathPoints){

  //// Defining variables BEGIN /////

  mat DeltaMat(1,1); // Used as matrix arg MvrNorm
  DeltaMat << Delta << endr;
  mat ZZZ(1,2);

  int J = GibbsN; // Number of iterations of the Gibbs sampler
  int parlen=0;
  // Number of parameters in the drift
  if(Model=="OU"){
    parlen = OUdrift(1,1).n_rows-1;
  }
  if(Model=="CIR"){
    parlen = CIRdrift(1,1).n_rows-1;
  }
  if(Model=="FHN"){
    parlen = FHNdrift(1,1).n_rows-1;
  }
  if(Model=="FHN5"){
    parlen = FHN5drift(1,1).n_rows-1;
  }


  mat likelihood(1,1);
  mat likelihoodProp(1,1);
  double driftPrior;
  double driftPriorProp;
  mat driftSample(1,parlen);

  mat TmpVar(4,4);
  mat TmpMean(4,1);

  int N = V.n_rows; // Number of observed records
  int NN = N-1;  // For easier notation
  double MM=M-1; // For easier notation

  mat JDet2(N,1);
  mat JDet2Prop(N,1);

  mat weight(NN,1);
  mat weight2(NN,1);
  mat weightSigma(1,1);
  mat weightSigmaProp(NN,1);

  mat weightTmp(NN,1);
  mat weightTmp2(NN,1);
  mat Tmp(parlen,parlen+1);
  mat Z(N,2);

  bool flag = (V.n_cols==2);
  int drift;
  int diffusion;

  if(flag==0){
    mat And(NN+1,1);
    if(LatentPathStart.n_rows!=V.n_rows){
      And.fill(LatentPathStart(0,0));
      V = join_rows(V,And);
    }
    if(LatentPathStart.n_rows==V.n_rows){
      V = join_rows(V,LatentPathStart);
    }
  }

  //// Defining variables END /////


  //// Specifying initial value for drift and diffusion parameters BEGIN ////
  mat driftEsts(J+1,parlen);
  for(int i=0;i<parlen;++i){
    driftEsts(0,i) = Start(0,i);
  }

  mat diffEsts(J+1,2);
  diffEsts(0,0) = Start(0,parlen);
  diffEsts(0,1) = Start(0,parlen+1);

  // Setting initial values for drift and diffusion parameters and
  // checking if all drift or diffusion parameters are known.

  int pmnc = ParMat.n_cols;

  drift = 1;
  diffusion = 1;

  if(ParMat.n_rows==3){
    for(int i=0;i<pmnc;i++){
      if(ParMat(2,i)==1.0){
	driftEsts(0,ParMat(0,i)-1)=ParMat(1,i);
      }
      if(ParMat(2,i)==2){
	diffEsts(0,ParMat(0,i)-1-parlen)=ParMat(1,i);
      }
    }
    if(pmnc>=parlen && ParMat(2,parlen-1)==1.0){
      drift=0;
    }
    if(pmnc>=2 && ParMat(2,pmnc-2)==2){
      diffusion=0;
    }
  }

  //// Specifying initial value for drift and diffusion parameters END ////


  int TopN;
  int TopN2;
  TopN = floor((double)(J)/10)-1; /*FFF*/
  TopN2 = min(10000-1,TopN);
  mat driftEsts2(TopN2,parlen);
  mat diffEsts2(TopN2,2);

    Mat<int> PathAcc1(NN,J);
    Mat<int> PathAcc12(NN,TopN2);

    mat Vupdate(N,J+1);
    mat Vupdate2(N,TopN2);
    Mat<int> PathAcc2(N,J);
    Mat<int> PathAcc22(N,TopN2);

    Vupdate.col(0) = V.col(1);
    mat Vtilde(NN*MM+1,2);
    mat VtildeProp(NN*MM+1,2);
    mat VtildeProp22(NN*MM+1,2);
    mat Ztilde(NN*MM+1,2);
    mat ZtildeProp(NN*MM+1,2);
    mat ZtildeProp2(NN*MM+1,2);

    mat start1(1,2);
    start1.fill(0);

    mat NUMVEC(M,1);
    for(int i=0; i<M; ++i){
      NUMVEC(i,0)=i;
    }
    mat NUMVEC2(2*M-1,1);
    for(int i=0; i<(2*M-1); ++i){
      NUMVEC2(i,0)=i;
    }
    mat JDet(1,1);
    mat JDetProp(1,1);
    mat sigmaprior(1,1);
    mat sigmapriorProp(1,1);
    mat tmpG(NN,1);
    mat logG(1,1);
    mat logGProp(1,1);
    mat pathWeight1(1,1);
    mat pathWeight2(1,1);
    Mat<int> diffAcc(J,1);
    Mat<int> diffAcc2(TopN2,1);
    Mat<int> driftAcc(J,1);
    Mat<int> driftAcc2(TopN2,1);
    mat sigmaProp(1,2);
    mat ZZTmp(N,2);

    mat ZZ(N,2);

    if(Model=="CIR"){
      for(int i=0; i<N; ++i){
	ZZ(i,0) = 2.0*sqrt(V(i,0))/diffEsts(0,0);
	ZZ(i,1) = 2.0*sqrt(V(i,1))/diffEsts(0,1);
      }
    }
    else{
      ZZ = V*inv(diagmat(diffEsts.row(0)));
    }
    for(int i=0; i<NN; ++i){
      Vtilde.rows(i*MM,(i+1)*MM) = BBSimInt(start1,
					    start1,
					    M,
					    1,
					    1,
					    Delta,
					    0);

      Ztilde.rows(i*MM,(i+1)*MM) =
	Vtilde.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZ.row(i)+(NUMVEC/MM)*ZZ.row(i+1);
    }

    for(int j=1; j<=J; ++j){
      /// Updating path except latent coordinate at observed time points ///
      if(Model=="CIR"){
	for(int i=0; i<N; ++i){
	  ZZ(i,0) = 2*sqrt(V(i,0))/diffEsts(j-1,0);
	  ZZ(i,1) = 2*sqrt(V(i,1))/diffEsts(j-1,1);
	}
      }
      else{
	ZZ = V*inv(diagmat(diffEsts.row(j-1)));
      }
      for(int i=0; i<NN; ++i){
	Ztilde.rows(i*MM,(i+1)*MM) =
	  Vtilde.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZ.row(i)+(NUMVEC/MM)*ZZ.row(i+1);

	weight.row(i) = Radon(Ztilde.rows(i*MM,(i+1)*MM),
			      driftEsts.row(j-1),
			      diffEsts(j-1,0),
			      diffEsts(j-1,1),
			      Delta,
			      Model);

	////// Proposal computations //////
	VtildeProp.rows(i*MM,(i+1)*MM) =
	  ((1.0-NUMVEC/MM)*Vtilde.row(i)+(NUMVEC/MM)*Vtilde.row(i+1))*(1-RWrhoPaths)+
	  RWrhoPaths*Vtilde.rows(i*MM,(i+1)*MM)+sqrt(1.0-pow(RWrhoPaths,2))*
	  BBSimInt(start1,
		   start1,
		   M,
		   1,
		   1,
		   Delta,
		   0);

	ZtildeProp.rows(i*MM,(i+1)*MM) =
	  VtildeProp.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZ.row(i)+(NUMVEC/MM)*ZZ.row(i+1);

	weightTmp.row(i)= Radon(ZtildeProp.rows(i*MM,(i+1)*MM),
				driftEsts.row(j-1),
				diffEsts(j-1,0),
				diffEsts(j-1,1),
				Delta,
				Model);

      	if(log(randu())< as_scalar(weightTmp.row(i)-weight.row(i))){
	  Vtilde.rows(i*MM,(i+1)*MM) = VtildeProp.rows(i*MM,(i+1)*MM);
	  Ztilde.rows(i*MM,(i+1)*MM) = ZtildeProp.rows(i*MM,(i+1)*MM);
	  PathAcc1(i,j-1) = 1;}
	else{PathAcc1(i,j-1) = 0;}
      }

      /// Updating latent coordinate at observed time points ///

      if(flag==0){
	for(int i=0; i<(NN-1); ++i){
	  if(Model=="CIR"){
	    ZZ(i,0)   = 2.0*sqrt(V(i,0))/diffEsts(j-1,0);
	    ZZ(i+1,0) = 2.0*sqrt(V(i+1,0))/diffEsts(j-1,0);
	    ZZ(i+2,0) = 2.0*sqrt(V(i+2,0))/diffEsts(j-1,0);
	    ZZ(i,1)   = 2.0*sqrt(V(i,1))/diffEsts(j-1,1);
	    ZZ(i+1,1) = 2.0*sqrt(V(i+1,1))/diffEsts(j-1,1);
	    ZZ(i+2,1) = 2.0*sqrt(V(i+2,1))/diffEsts(j-1,1);
	  }
	  else{
	    ZZ.submat(i,1,i,1)     = V.submat(i,1,i,1)/diffEsts(j-1,1);
	    ZZ.submat(i+1,1,i+1,1) = V.submat(i+1,1,i+1,1)/diffEsts(j-1,1);
	    ZZ.submat(i+2,1,i+2,1) = V.submat(i+2,1,i+2,1)/diffEsts(j-1,1);
	  }
	  Ztilde.rows(i*MM,(i+1)*MM) =
	    Vtilde.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZ.row(i)+(NUMVEC/MM)*ZZ.row(i+1);
	  Ztilde.rows((i+1)*MM,(i+2)*MM) =
	    Vtilde.rows((i+1)*MM,(i+2)*MM)+(1.0-NUMVEC/MM)*ZZ.row(i+1)+(NUMVEC/MM)*ZZ.row(i+2);

	  ZZZ.submat(0,0,0,0) = ZZ.submat(i+1,0,i+1,0);

	  // Proposal:
	  ZZZ.submat(0,1,0,1) =
	    .5*(ZZ.submat(i,1,i,1)+ZZ.submat(i+2,1,i+2,1))*(1-RWrho2PathPoints)+
	    RWrho2PathPoints*ZZ.submat(i+1,1,i+1,1)+
	    sqrt(1.0-pow(RWrho2PathPoints,2))*
	    MvrNorm(0*DeltaMat,.5*DeltaMat);

	  ZtildeProp2.rows(i*MM,(i+1)*MM) =
	    Vtilde.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZ.row(i)+(NUMVEC/MM)*ZZZ.row(0);
	  ZtildeProp2.rows((i+1)*MM,(i+2)*MM) =
	    Vtilde.rows((i+1)*MM,(i+2)*MM)+(1.0-NUMVEC/MM)*ZZZ.row(0)+(NUMVEC/MM)*ZZ.row(i+2);


	  /// Recomputing the log weight ///
	  weight2.row(i) = Radon(Ztilde.rows(i*MM,(i+2)*MM),
				 driftEsts.row(j-1),
				 diffEsts(j-1,0),
				 diffEsts(j-1,1),
				 2*Delta,
				 Model);

	  /// Computing the log weight for the proposal path ///
	  weightTmp2.row(i) = Radon(ZtildeProp2.rows(i*MM,(i+2)*MM),
				    driftEsts.row(j-1),
				    diffEsts(j-1,0),
				    diffEsts(j-1,1),
				    2*Delta,
				    Model);

	  if(log(randu())< as_scalar(weightTmp2.row(i)-weight2.row(i))){
	    Ztilde.submat((i+1)*MM,1,MM*(1+i),1) = ZZZ.submat(0,1,0,1);
	    /// Transforming proposal, 2. coord ///
	    if(Model=="CIR"){
	      Ztilde.submat((i+1)*MM,0,MM*(1+i),0) = ZZZ.submat(0,0,0,0);
	      V(i+1,1) = pow(Ztilde((i+1)*MM,1)*diffEsts(j-1,1),2)/4.0;
	    }
	    else{
	      V(i+1,1) = Ztilde((i+1)*MM,1)*diffEsts(j-1,1);
	    }
	    PathAcc2(i+1,j-1) = 1;
	  }
	  else{
	    PathAcc2(i+1,j-1) = 0;
	  }
	  Vupdate(i+1,j) = V(i+1,1);
	}//i
      }//flag==0

      //FIRST ENDPOINT
      if(flag==0){//0
	  if(Model=="CIR"){
	    ZZ(0,0) = 2.0*sqrt(V(0,0))/diffEsts(j-1,0);
	    ZZ(1,0) = 2.0*sqrt(V(1,0))/diffEsts(j-1,0);

	    ZZ(0,1) = 2.0*sqrt(V(0,1))/diffEsts(j-1,1);
	    ZZ(1,1) = 2.0*sqrt(V(1,1))/diffEsts(j-1,1);
	  }
	  else{
	  ZZ.submat(0,1,0,1) = V.submat(0,1,0,1)/diffEsts(j-1,1);
	  ZZ.submat(1,1,1,1) = V.submat(1,1,1,1)/diffEsts(j-1,1);
	  }
	  Ztilde.rows(0,MM) =
	    Vtilde.rows(0,MM)+(1.0-NUMVEC/MM)*ZZ.row(0)+(NUMVEC/MM)*ZZ.row(1);

	  // Proposal
	  ZZZ.submat(0,0,0,0) = ZZ.submat(0,0,0,0);
	  ZZZ.submat(0,1,0,1) =
	    //	  MvrNorm(ZZ.submat(1,1,1,1),0.5*DeltaMat);
	    MvrNorm(ZZ.submat(1,1,1,1),1.0*DeltaMat);

	  ZtildeProp2.rows(0,MM) =
	    Vtilde.rows(0,MM)+(1.0-NUMVEC/MM)*ZZZ.row(0)+(NUMVEC/MM)*ZZ.row(1);

	  /// Recomputing the log weight ///
	  weight2.row(0) = Radon(Ztilde.rows(0,MM),
				 driftEsts.row(j-1),
				 diffEsts(j-1,0),
				 diffEsts(j-1,1),
				 Delta,
				 Model);

	  /// Computing the log weight for the proposal path ///
	  weightTmp2.row(0) = Radon(ZtildeProp2.rows(0,MM),
				    driftEsts.row(j-1),
				    diffEsts(j-1,0),
				    diffEsts(j-1,1),
				    Delta,
				    Model);

//	  if(log(randu())< as_scalar(weightTmp2.row(0)-weight2.row(0))){
	  if(log(randu())<
	     as_scalar(weightTmp2.row(0)-
		       0.5/LatentVarY0*pow(ZtildeProp2(0,1)-LatentMeanY0,2)-
		       weight2.row(0)+
		       0.5/LatentVarY0*pow(Ztilde(0,1)-LatentMeanY0,2))){
	    Ztilde.submat(0,1,0,1) = ZZZ.submat(0,1,0,1);
	    // Transforming proposal, 2. coord ///
	    if(Model=="CIR"){
	      //    Ztilde.submat(0,0,0,0) = ZZZ.submat(0,0,0,0);
	      V(0,1) = pow(Ztilde(0,1)*diffEsts(j-1,1),2)/4.0;
	    }
	    else{
	      V(0,1) = Ztilde(0,1)*diffEsts(j-1,1);
	    }
	    PathAcc2(0,j-1) = 1;
	  }
	  else{
	    PathAcc2(0,j-1) = 0;
	  }
	  Vupdate(0,j) = V(0,1);
      }

      //LAST ENDPOINT
      if(flag==0){//0
	  if(Model=="CIR"){
	    ZZ(NN-1,0) = 2.0*sqrt(V(NN-1,0))/diffEsts(j-1,0);
	    ZZ(NN,0)   = 2.0*sqrt(V(NN,0))/diffEsts(j-1,0);

	    ZZ(NN-1,1) = 2.0*sqrt(V(NN-1,1))/diffEsts(j-1,1);
	    ZZ(NN,1)   = 2.0*sqrt(V(NN,1))/diffEsts(j-1,1);
	  }
	  else{
	    ZZ.submat(NN-1,1,NN-1,1) = V.submat(NN-1,1,NN-1,1)/diffEsts(j-1,1);
	    ZZ.submat(NN,1,NN,1)     = V.submat(NN,1,NN,1)/diffEsts(j-1,1);
	  }
	  Ztilde.rows((NN-1)*MM,NN*MM) =
	    Vtilde.rows((NN-1)*MM,NN*MM)+(1.0-NUMVEC/MM)*ZZ.row(NN-1)+(NUMVEC/MM)*ZZ.row(NN);

	  ZZZ.submat(0,0,0,0) = ZZ.submat(NN,0,NN,0);
	  // Proposal
	  //	  ZZZ.submat(0,1,0,1) =	  MvrNorm(ZZ.submat(NN-1,1,NN-1,1),0.5*DeltaMat);
	  ZZZ.submat(0,1,0,1) =  MvrNorm(ZZ.submat(NN-1,1,NN-1,1),1.0*DeltaMat);


	  ZtildeProp2.rows((NN-1)*MM,NN*MM) =
	    Vtilde.rows((NN-1)*MM,NN*MM)+(1.0-NUMVEC/MM)*ZZ.row(NN-1)+(NUMVEC/MM)*ZZZ.row(0);

	  /// Recomputing the log weight ///
	  weight2.row(NN-1) = Radon(Ztilde.rows((NN-1)*MM,NN*MM),
				    driftEsts.row(j-1),
				    diffEsts(j-1,0),
				    diffEsts(j-1,1),
				    Delta,
				    Model);

	  /// Computing the log weight for the proposal path ///
	  weightTmp2.row(NN-1) = Radon(ZtildeProp2.rows((NN-1)*MM,NN*MM),
				       driftEsts.row(j-1),
				       diffEsts(j-1,0),
				       diffEsts(j-1,1),
				       Delta,
				       Model);

	  if(log(randu())< as_scalar(weightTmp2.row(NN-1)-weight2.row(NN-1))){
	    Ztilde.submat(NN*MM,1,NN*MM,1) = ZZZ.submat(0,1,0,1);
	    // Transforming proposal, 2. coord ///
	    if(Model=="CIR"){
	      Ztilde.submat(NN*MM,0,NN*MM,0) = ZZZ.submat(0,0,0,0);
	      V(NN,1) = pow(Ztilde(NN*MM,1)*diffEsts(j-1,1),2)/4.0;
	    }
	    else{
	      V(NN,1) = Ztilde(NN*MM,1)*diffEsts(j-1,1);
	    }
	    PathAcc2(NN,j-1) = 1;
	  }
	  else{
	    PathAcc2(NN,j-1) = 0;
	  }
	  Vupdate(NN,j) = V(NN,1);
      }

      /// Updating drift parameters ///
      //Test BEGIN. CIR not included AAA MH step
      if(drift==1 && MHdrift==2){

      driftPrior =
	-0.5*as_scalar((driftEsts.row(j-1)-driftPriorMean)*inv(driftPriorCovar)*trans(driftEsts.row(j-1)-driftPriorMean));
      likelihood = Radon(Ztilde,
			 driftEsts.row(j-1),
			 diffEsts(j-1,0),
			 diffEsts(j-1,1),
			 Delta*NN,
			 Model);

      driftSample = trans(MvrNorm(trans(driftEsts.row(j-1)),driftRW));

      if(ParMat.n_rows==3){
	for(int i=0;i<pmnc;i++){
	  if(ParMat(2,i)==1.0){
	    driftSample(ParMat(0,i)-1)=ParMat(1,i);
	  }
	}
      }
      driftPriorProp =
      	-0.5*as_scalar((driftSample-driftPriorMean)*inv(driftPriorCovar)*trans(driftSample-driftPriorMean));
      likelihoodProp = Radon(Ztilde,
      			     driftSample,
      			     diffEsts(j-1,0),
      			     diffEsts(j-1,1),
			     Delta*NN,
      			     Model);

      if(log(randu())< as_scalar(driftPriorProp+likelihoodProp-driftPrior-likelihood)){
      	driftEsts.row(j) = driftSample;
	driftAcc(j-1,0) = 1;
      }
      else{
	driftEsts.row(j) = driftEsts.row(j-1);
	driftAcc(j-1,0) = 0;
      }

      }
      //TEST END
      //TEST2 BEGIN Sampling from posterior using prior
      if(drift==1 && MHdrift==1){
	Tmp =
	  PostDistDrift(Ztilde*diagmat(diffEsts.row(j-1)),Delta/MM,diffEsts(j-1,0),diffEsts(j-1,1),parlen,Model);
	TmpVar = inv(inv(Tmp.cols(1,parlen))+inv(driftPriorCovar));
	TmpMean = TmpVar*(inv(Tmp.cols(1,parlen))*Tmp.col(0)+inv(driftPriorCovar)*trans(driftPriorMean));

	driftEsts.row(j) =
	  trans(MvrNorm(TmpMean,TmpVar));
	}
      //TEST2 END
      //Sampling from posterior without using prior
      if(drift==1 && MHdrift==0){
	if(Model=="CIR"){
	  Tmp =
	    PostDistDrift(pow(Ztilde,2)*diagmat(pow(diffEsts.row(j-1),2))/4.0,Delta/MM,diffEsts(j-1,0),diffEsts(j-1,1),parlen,Model);
	  driftEsts.row(j) =
	    trans(MvrNorm(Tmp.col(0),Tmp.cols(1,parlen)));
	}
	else{
	  Tmp =
	    PostDistDrift(Ztilde*diagmat(diffEsts.row(j-1)),Delta/MM,diffEsts(j-1,0),diffEsts(j-1,1),parlen,Model);
	  driftEsts.row(j) =
	    trans(MvrNorm(Tmp.col(0),Tmp.cols(1,parlen)));
	}
      }

      if(drift==0){
	driftEsts.row(j) = driftEsts.row(j-1);
      }

      if(ParMat.n_rows==3){
	for(int i=0;i<pmnc;i++){
	  if(ParMat(2,i)==1.0){
	    driftEsts(j,ParMat(0,i)-1)=ParMat(1,i);
	  }
	}
      }

      /// Updating diffusion parameters ///
      if(diffusion==0){
	diffEsts.row(j) = diffEsts.row(j-1);
      }
      if(diffusion==1){

	if(Model=="CIR"){
	  for(int i=0; i<N; ++i){
	    ZZ(i,0) = 2.0*sqrt(V(i,0))/diffEsts(j-1,0);
	    ZZ(i,1) = 2.0*sqrt(V(i,1))/diffEsts(j-1,1);
	    /* JDet2(i) = */
	    /*   1.0/(diffEsts(i,0)*diffEsts(i,1)*sqrt(V(i,0)*V(i,1))); */
	    JDet2(i) =
	      1.0/(pow(diffEsts(i,0)*diffEsts(i,1),2)*ZZ(i,0)*ZZ(i,1)/4.0);
	  }
	  for(int i=0;i<NN;++i){
	    Ztilde.rows(i*MM,(i+1)*MM) =
	      Vtilde.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZ.row(i)+(NUMVEC/MM)*ZZ.row(i+1);
	  }
	}


	else{
	  ZZ = V*inv(diagmat(diffEsts.row(j-1)));
	  JDet = det(inv(diagmat(diffEsts.row(j-1))));
	  for(int i=0;i<NN;++i){
	    Ztilde.rows(i*MM,(i+1)*MM) =
	      Vtilde.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZ.row(i)+(NUMVEC/MM)*ZZ.row(i+1);
	  }
	}

	//////////////////////////////
	sigmaprior =
	  -0.5*(log(diffEsts.row(j-1))-diffPriorMean)*inv(diffPriorCovar)*trans(log(diffEsts.row(j-1))-diffPriorMean);
	for(int ii=0;ii<NN;++ii){
	  tmpG.row(ii) = (ZZ.row(ii+1)-ZZ.row(ii))*trans(ZZ.row(ii+1)-ZZ.row(ii));
	}
        logG = -1.0/(2.0*Delta)*accu(tmpG);
	weightSigma = Radon(Ztilde,
			    driftEsts.row(j),
			    diffEsts(j-1,0),
			    diffEsts(j-1,1),
			    Delta*NN,
			    Model);

	///////////////////
	if(Model=="CIR"){
	  //	  pathWeight1 =	sigmaprior+log(accu(JDet2))+logG+weightSigma;
	  pathWeight1 = sigmaprior+accu(log(JDet2))+logG+weightSigma;
	}
	else{
	  pathWeight1 = sigmaprior+N*log(JDet)+logG+weightSigma;
	}

	/// Proposal

	sigmaProp = trans(exp(MvrNorm(log(trans(diffEsts.row(j-1))),diffRW)));
	if(ParMat.n_rows==3){
	  for(int i=0;i<pmnc;i++){
	    if(ParMat(2,i)==2){
	      sigmaProp(0,ParMat(0,i)-1-parlen)=ParMat(1,i);
	    }
	  }
	}
	///////////////
	if(Model=="CIR"){
	  for(int i=0; i<N; ++i){
	    ZZTmp(i,0) = 2.0*sqrt(V(i,0))/sigmaProp(0,0);
	    ZZTmp(i,1) = 2.0*sqrt(V(i,1))/sigmaProp(0,1);
	    /* JDet2Prop(i) =
	  1.0/(sigmaProp(0,0)*sigmaProp(0,1)*sqrt(V(i,0)*V(i,1))); */
	    JDet2Prop(i) =
	      1.0/(pow(sigmaProp(0,0)*sigmaProp(0,1),2)*ZZTmp(i,0)*ZZTmp(i,1)/4.0);
	  }
	  for(int i=0;i<NN;++i){
	    ZtildeProp.rows(i*MM,(i+1)*MM) =
	      Vtilde.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZTmp.row(i)+(NUMVEC/MM)*ZZTmp.row(i+1);
	  }
	}
	///////////

	else{
	  ZZTmp = V*inv(diagmat(sigmaProp));
	  JDetProp = det(inv(diagmat(sigmaProp)));
	  for(int i=0;i<NN;++i){
	    ZtildeProp.rows(i*MM,(i+1)*MM) =
	      Vtilde.rows(i*MM,(i+1)*MM)+(1.0-NUMVEC/MM)*ZZTmp.row(i)+(NUMVEC/MM)*ZZTmp.row(i+1);
	  }
	}

	///////////
	sigmapriorProp = //-0.5*log(det(diffPriorCovar))
	  -0.5*(log(sigmaProp)-diffPriorMean)*inv(diffPriorCovar)*trans(log(sigmaProp)-diffPriorMean);
	for(int ii=0;ii<NN;++ii){
	  tmpG.row(ii) = (ZZTmp.row(ii+1)-ZZTmp.row(ii))*trans(ZZTmp.row(ii+1)-ZZTmp.row(ii));
	}
        logGProp = -1.0/(2.0*Delta)*accu(tmpG);
	weightSigmaProp = Radon(ZtildeProp,
				driftEsts.row(j),
				sigmaProp(0,0),
				sigmaProp(0,1),
				Delta*NN,
				Model);

	if(Model=="CIR"){
	  //	  pathWeight2 =	sigmapriorProp+log(accu(JDet2Prop))+logGProp+weightSigmaProp;
	  pathWeight2 = sigmapriorProp+accu(log(JDet2Prop))+logGProp+weightSigmaProp;
	}
	else{
	  pathWeight2 =
	    sigmapriorProp+N*log(JDetProp)+logGProp+weightSigmaProp;
	}

	if(log(randu())< as_scalar(pathWeight2-pathWeight1)){
	  diffEsts.row(j) = sigmaProp;
          ZZ = ZZTmp;
          diffAcc(j-1,0) = 1;
	}
	else{
	  diffEsts.row(j) = diffEsts.row(j-1);
	  diffAcc(j-1,0) = 0;
	}
      }//diffusion==1
    }//j

    BIPOD BIPODESTS;

    if(J>10000){
      for(int k=0;k<TopN2;++k){
	driftEsts2.row(k)=driftEsts.row(10*k);
	diffEsts2.row(k)=diffEsts.row(10*k);
	PathAcc12.col(k)=PathAcc1.col(10*k);
	diffAcc2.row(k)=diffAcc.row(10*k);
	driftAcc2.row(k)=driftAcc.row(10*k);
	if(flag==0){
	  Vupdate2.col(k)=Vupdate.col(10*k);
	  PathAcc22.col(k)=PathAcc2.col(10*k);
	}
      }

      if(drift==1){
	BIPODESTS.Drift = driftEsts2;
      }
      if(diffusion==1){
	BIPODESTS.Diff       = diffEsts2;
	BIPODESTS.DiffAcc    = diffAcc2;
      }
      if(drift==1 && MHdrift==2){
	BIPODESTS.DriftAcc = driftAcc2;
      }
      if(flag==0){
	BIPODESTS.PathAccPoints   = PathAcc22;
	BIPODESTS.LatentPath = Vupdate2;
      }
      BIPODESTS.PathAcc   = PathAcc12;
    }
    if(J<=10000){
      if(drift==1){
	BIPODESTS.Drift      = driftEsts;
      }
      if(diffusion==1){
	BIPODESTS.Diff       = diffEsts;
	BIPODESTS.DiffAcc    = diffAcc;
      }
      if(drift==1 && MHdrift==2){
	BIPODESTS.DriftAcc = driftAcc;
      }

      if(flag==0){
	BIPODESTS.PathAccPoints   = PathAcc2;
	BIPODESTS.LatentPath = Vupdate;
      }
      BIPODESTS.PathAcc   = PathAcc1;
    }

    if(drift==0){
	BIPODESTS.Drift << 0 << endr;
    }
    if(diffusion==0){
      BIPODESTS.Diff << 0 << endr;
      BIPODESTS.DiffAcc << 0 << endr;
    }
    if(flag==1){
      BIPODESTS.PathAccPoints << 0 << endr;
      BIPODESTS.LatentPath << 0 << endr;
    }
    BIPODESTS.flag       = flag;
    return(BIPODESTS);
}
