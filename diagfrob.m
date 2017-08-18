///////////////////////////////////////////////////// 
// This code is part of the diagfrob MAGMA library //
//                                                 //
// copyright (c) 2017 Jan Tuitman                  //
/////////////////////////////////////////////////////


load "linearrecurrence.m"; // Implementation of Bostan-Gaudry-Schost by Minzlaff and Fontein


print "";
print "*******************************";
print "* diagfrob Magma library v1.1 *";
print "*                             *";
print "* by Jan Tuitman              *";
print "*******************************";
print "";


basis_sets:=function(S, d)

  // Constructing a monomial basis for the cohomology

  n := Rank(S)-1;
  B := [* *];
  for k:=1 to n do
    if k*d lt n+1 then
      Append(~B, [* *]);
      continue;
    end if;
    f := S!0;
    for i:=1 to n+1 do
      f := f + S.i;
    end for;
    Lk := Monomials(f^(k*d - (n+1)));
    Bk := [* *];
    for f in Lk do
      for i:=1 to n+1 do
        if Degree(f, i) ge d-1 then
          continue f;
        end if;
      end for;
      Append(~Bk, f);
    end for;
    Append(~B, Bk);
  end for;
	
  return B;

end function;


find_uv:=function(p,n,d)

  // Finding the basis vectors of the cohomology and
  // their images under p^(-1)*Frob_p.

  S:=PolynomialRing(IntegerRing(),n+1);
  I:=basis_sets(S,d);

  ulist:=[];
  for i:=1 to #I do
    for j:=1 to #I[i] do
      ulist:=Append(ulist,Exponents(I[i][j]));
    end for;
  end for;

  Zd:=IntegerRing(d);
  vlist:=[];

  for i:=1 to #ulist do
    u:=ulist[i];
    v:=[];
    for j:=1 to #u do
      v[j]:=IntegerRing()!Zd!(p*(u[j]+1)-1);
    end for;
    vlist:=Append(vlist,v);
  end for;

return ulist,vlist;

end function;


precisions_chi:=function(p,n,d)

  // Compute the precisions to which the coefficients of the L-factor 
  // chi have to be known, for it to be determined exactly (for a 
  // hypersurface of degree d in P^n over F_p with p prime).

  ulist,vlist:=find_uv(p,n,d);
  b:=#ulist;
  N:=[];
  for i:=1 to b do
    N:=Append(N,Floor(Log(p,2*(b/i)*p^(i*(n-1)/2)))+1);
  end for;

  S:=PolynomialRing(IntegerRing(),n+1);
  I:=basis_sets(S,d);
  
  Gamma:=[0];
  h:=0;
  for i:=1 to #I do
    for j:=1 to (#I[i]) do
      h:=h+(i-1);
      Gamma:=Append(Gamma,h);
    end for;
  end for;

  Nprime:=[];
  if p ge n then
    for i:=1 to #N do
      Nprime[i]:=N[i]-Gamma[i];
    end for;
  end if;

  return N,Maximum(Nprime); 

end function;


precisions_Phi0:=function(p,n,NPhi0)

  // Compute the precision parameters for computing p^(-1)*Frob_p
  // on H^n of the complement of a hypersurface in P^n over F_p with
  // p prime to p-adic precision NPhi0.

  delta:=Valuation(Factorial(n-1),p);
  for i:=1 to n-1 do
    delta:=delta+Floor(Log(p,i));
  end for;

  NPhi0prime:=NPhi0+(n-1)+Valuation(Factorial(n-1),p)+2*delta;

  M:=Ceiling((p^2/(p-1))*(NPhi0prime+Log(p,NPhi0prime+3)+4));
  R:=Floor(M/p);

  return NPhi0prime,M,R;

end function;


rising_fac:=function(w,r)

  prod:=Parent(w)!1;
  if r eq 0 then
    return prod;
  end if;
  for i:=1 to r do
    prod:=prod*w;
    w:=w+1;
  end for;
  return prod;

end function;


factorials:=function(p,N,R,m0list)

  // Compute the factorials needed in compute_Phi0 naively

  Qp:=pAdicField(p:Precision:=N);

  facsm:=[];
  for i:=1 to #m0list do
    m0:=m0list[i];
    facsm0:=[];
    for j:=0 to R do
      facsm0[j+1]:=Qp!Factorial(m0+p*j);
    end for;
    facsm:=Append(facsm,facsm0);
  end for;

  return facsm;

end function;


factorials_bgs:=function(p,N,R,m0list);

  // Compute the factorials needed in compute_Phi0 using the
  // Baby Step Biant Step following Bostan-Gaudry-Schost. 

  N:=N+Floor((Maximum(m0list)+p*R)/(p-1));
  Qp:=pAdicRing(p:Precision:=N);
  Qpx<x>:=PolynomialRing(Qp);

  M:=Matrix(Qpx,1,1,[x]);

  mlist:=[];
  for i:=1 to #m0list do
    for j:=0 to R do
      mlist:=Append(mlist,m0list[i]+p*j);
    end for;
  end for;
  mlist:=Sort(mlist);

  L_:=[];
  R_:=[];
  for i:=1 to #mlist do
    L_[i]:=0;
    R_[i]:=mlist[i];
  end for;
  s := Ilog(4,R_[#R_]);
  DD := UpperCaseDD(Qp!1,2^s,2^s);
  DDi := DD^(-1);

  seq:=LinearRecurrence(M,L_,R_,DDi,s); // Call to code of Minzlaff and Fontein

  facs:=[];
  for i:=1 to #m0list do
    facsm0:=[];
    for j:=0 to R do
      m:=m0list[i]+p*j;
      for k:=1 to #mlist do
        if mlist[k] eq m then
          facsm0[j+1]:=seq[k][1,1];
        end if;
      end for;
    end for; 
    facs:=Append(facs,facsm0);
  end for;

  return facs;

end function;


mu_list:=function(p,N,R,m0,a,facsm0) // for given m0 and a find mu for all r    
    
  Qp:=pAdicField(p:Precision:=N);

  facs:=[];
  for j:=0 to R do
    facs[j+1]:=Qp!Factorial(j);
  end for;

  aipow:=[];
  aipow[1]:=Qp!1;
  aipow[2]:=(Qp!a)^(p-1);
  for j:=3 to R+1 do
    aipow[j]:=aipow[j-1]*aipow[2];
  end for;

  mulist:=[];
  mulist[1]:=1/(facsm0[1]);
  for r:=1 to R do
    sum:=Qp!0;
    for j:=0 to r do
      sum:=sum+p^(r-j)/(facsm0[r-j+1]*facs[j+1]*aipow[j+1]);
    end for;
    mulist:=Append(mulist,sum);
  end for;

  return mulist;

end function;


alpha:=function(p,n,d,N,M,u,v,mu_matrix,m0alist,alist)

  Qp:=pAdicField(p:Precision:=N);
  
  prod:=Qp!1;
  for i:=0 to n do
    
    apow:=(Qp!alist[i+1])^(p-1);
    m0:=IntegerRing()!((p*(u[i+1]+1)-(v[i+1]+1))/d);
    for j:=1 to #m0alist do
      if [m0,alist[i+1]] eq m0alist[j] then
        m0arank:=j;
      end if;
    end for;
    R:=Floor(M/p);
    
    sum:=Qp!0;
    for r:=0 to R do
      sum:=sum+apow^(r)*rising_fac(Qp!(u[i+1]+1)/d,r)*mu_matrix[m0arank][r+1];
    end for;
    
    prod:=prod*alist[i+1]^(m0)*sum;
  end for;

  return prod;

end function;


compute_Phi0:=function(p,n,d,alist,NPhi0:bgs:=true)

  // Computes the action of p^(-1)*Frob_p on H^n of the complement of the 
  // diagonal hypersurface defined by a_0 x_0^d + .. + a_n x_n^d in P^n 
  // over the finite field F_p with p prime to p-adic precision NPhi0.

  ulist,vlist:=find_uv(p,n,d);

  m0list:=[];
  m0alist:=[];
  for i:=1 to #ulist do
    u:=ulist[i];
    v:=vlist[i];
    for j:=0 to n do
      m0:=IntegerRing()!((p*(u[j+1]+1)-(v[j+1]+1))/d);
      a:=alist[j+1];
      if not m0 in m0list then
        m0list:=Append(m0list,m0);
      end if;
      if not [m0,a] in m0alist then 
        m0alist:=Append(m0alist,[m0,a]);
      end if;
    end for;
  end for;

  NPhi0prime,M,R:=precisions_Phi0(p,n,NPhi0);

  if bgs then
    facsm:=factorials_bgs(p,NPhi0prime,R,m0list); 
  else
    facsm:=factorials(p,NPhi0prime,R,m0list);
  end if;

  mu_matrix:=[];
  for i:=1 to #m0alist do
    m0:=m0alist[i][1];
    a:=m0alist[i][2];
    for j:=1 to #m0list do
      if m0 eq m0list[j] then
        m0index:=j;
      end if;
    end for;
    mu_matrix:=Append(mu_matrix,mu_list(p,NPhi0prime,R,m0,a,facsm[m0index]));
  end for;

  sum_u:=0;
  for i:=0 to n do
    sum_u:=sum_u+u[i+1];
  end for;
  ku:=IntegerRing()!((n+1+sum_u)/d);

  Qp:=pAdicField(p:Precision:=NPhi0);
  Phi0:=ZeroMatrix(Qp,#ulist,#ulist);

  for i:=1 to #ulist do
    u:=ulist[i];
    v:=vlist[i];
    
    sum_u:=0;
    for i:=0 to n do
      sum_u:=sum_u+u[i+1];
    end for;
    ku:=IntegerRing()!((n+1+sum_u)/d);

    sum_v:=0;
    for i:=0 to n do
      sum_v:=sum_v+v[i+1];
    end for;
    kv:=IntegerRing()!((n+1+sum_v)/d);

    for j:=1 to #vlist do
      if vlist[i] eq ulist[j] then
        Phi0[i,j]:=(-1)^kv*p^(n-ku)*(Factorial(kv-1)/Factorial(ku-1))*alpha(p,n,d,NPhi0prime,M,u,v,mu_matrix,m0alist,alist)^(-1);
      end if;
    end for;
  end for;

  return Phi0;

end function; 


revcharpoly:=function(p,n,d,Phi) 

  // Compute the reverse characteristic polynomial of the Frobenius matrix
  // Phi0 using the Newton-Girard identities.

  b:=NumberOfRows(Phi);

  N:=precisions_chi(p,n,d);

  PhiQ:=ZeroMatrix(RationalField(),NumberOfRows(Phi),NumberOfColumns(Phi));
  for i:=1 to NumberOfRows(Phi) do
    for j:=1 to NumberOfColumns(Phi) do
      PhiQ[i,j]:=(RationalField()!Phi[i,j]); 
    end for;
  end for;

  ZT<T>:=PolynomialRing(IntegerRing()); 
  chi:=Reverse(Coefficients(ZT!CharacteristicPolynomial(PhiQ))); 
  
  chi[1]:=1;
  s:=[];
  s[1]:=b;

  for i:=1 to b do
    sum:=0;
    for j:=2 to i do
      sum:=sum-s[i-j+2]*chi[j];
    end for;
    lbound:=Ceiling(sum/i-p^(N[i])/2); 
    chi[i+1]:=lbound+((chi[i+1]-lbound) mod p^(N[i]));
    s[i+1]:=sum-i*chi[i+1];
  end for;

  return ZT!chi;

end function;


test_Lfactor:=function(chi,p,n)

  // Test whether the integer polynomial chi has roots of absolute value p^(-(n-1)/2)

  QT:=PolynomialRing(RationalField());
 
  f:=Evaluate(chi,QT.1)*Evaluate(chi,-QT.1);
 
  g:=QT!0;
  for i:=0 to Degree(chi) do
    g:=g+Coefficient(f,2*i)*QT.1^i;
  end for;

  g:=Evaluate(g,(QT.1)/p^(n-1));

  return HasAllRootsOnUnitCircle(g);

end function;



Lfactor:=function(p,n,d,alist:bgs:=true,test:=true)

  // Computes the polynomial chi(T) such that the diagonal hypersurface
  // defined by a_0 x_0^d + .. + a_n x_n^d in P^n over the finite field
  // F_p with p prime has zeta function given by:
  //
  // chi^{(-1)^n} / (1-T)(1-pT)...(1-p^{n-1}T)
  //
  // optional paramters:
  //
  // bgs  - if true use Baby Step Giant Step following Bostan-Gaudry-Schost (default).
  // test - if true test whether the output is a Weil polynomial of the right weight (default). 

  N,NPhi0:=precisions_chi(p,n,d);
  Phi0:=compute_Phi0(p,n,d,alist,NPhi0:bgs:=bgs);
  chi:=revcharpoly(p,n,d,Phi0);

  if not test_Lfactor(chi,p,n) then
    print "p =", p, "n =", n, "d =", d, "alist = ", alist, "bgs =", bgs;
    print "Output seems to be wrong. Please report this example to jan.tuitman@kuleuven.be.";
  end if;

  return chi;

end function;



