///////////////////////////////////////////////////// 
// This code is part of the diagfrob MAGMA library //
//                                                 //
// copyright (c) 2017 Jan Tuitman                  //
/////////////////////////////////////////////////////


load "diagfrob.m";


ZT<T>:=PolynomialRing(RationalField());

// diagonal hypersurface of degree d in P^n over Z with coefficients a_0,..,a_n:

n:=3;
d:=4;
alist:=[1,2,3,4];


for i:=5 to 30 do
  tijd:=Cputime();
  p:=NextPrime(2^i); 
  // chi:=Lfactor(p,n,d,alist:bgs:=false);   // don't use sqrt(p) trick 
  chi:=Lfactor(p,n,d,alist:bgs:=true);    // do use sqrt(p) trick
  print p, Cputime(tijd);
end for;


