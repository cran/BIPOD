mat FHNdrift(double x, double y){
  mat fAll(5,2);
  fAll << 0 << -y << endr
       << x-pow(x,3)-y << 0 << endr
       << 1 << 0 << endr
       << 0 << x << endr
       << 0 << 1 << endr;
  return(fAll);
}
