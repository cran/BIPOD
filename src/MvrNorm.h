mat MvrNorm(mat mean, mat covar){
  int p = covar.n_cols;
  mat A = chol(covar);
  mat Z = randn<mat>(p,1);
  return(mean+trans(A)*Z);
}
