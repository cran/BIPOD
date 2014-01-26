mat driftfun(double x,
	     double y,
	     mat driftpar,
	     string Model)
{
  mat BB(1,2);
  mat one(1,1);
  one.fill(1);
  if(Model=="OU"){
    BB = join_rows(one,driftpar)*OUdrift(x,y);
  }
  if(Model=="CIR"){
    BB = join_rows(one,driftpar)*CIRdrift(x,y);
  }
  if(Model=="FHN"){
    BB = join_rows(one,driftpar)*FHNdrift(x,y);
  }
  if(Model=="FHN5"){
    BB = join_rows(one,driftpar)*FHN5drift(x,y);
  }
  return trans(BB);
}
