mat OUdrift(double x, double y){
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
