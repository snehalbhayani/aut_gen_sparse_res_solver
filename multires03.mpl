##
##    Title: 	multires
##    Created:	Mon Jun 27 15:23:33 1998
##    Authors: 	Laurent Buse <lbuse@sophia.inria.fr>
##              Bernard Mourrain <mourrain@sophia.inria.fr>
##    Contributions of:
##              S. Tonelli (geometric decomposition 1996)
##		I. Emiris, J. Canny, P. Pedersen (sparse resultant 2000)   
##              O. Ruatta  (duality)
##
## Description: Elimination of variables, multivariate bezoutian, 
##              multivariate resultant, sparse resultant, ...
##
## Version last modified 08/2002
## ----------------------------------------------------------------------
##
##	The multires Maple file gathers some functions for dealing with
##      algebraic problems. You can use it with a simple : 
##      read "multires.mpl";
##      This files is divided in the following parts:
##
##	1) Polynomial manipulations: Basic functions to work with polynomials
##         are developped. There is functions to get coefficents of
##         polynomials as well as functions to put polynomials in 
##         coefficients matrices.
##	2) Resultant matrix construction: In this part some matrices of 
##         different kinds of resultants are developped (Bezoutians, 
##         Macaulay matrices, Jouanolou matrices, residual resultant, ...). 
##	3) Matrix operations: some functions dealing with eigenvalues, 
##         eigenvectors and Schur complement are given. These functions 
##         are generally used with resultant matrices in order
##         to solve polynomial systems.
##	4) Geometric decomposition: this part consists in an algorithm 
##         to compute the geometric decomposition of a variety: 
##         function "decomp".
##	5) Elimination procedure: melim, mesolve, ...
##	6) Duality: tools for dealing with dual forms on polynomials.
##	7) Toric resultant: in this big part of this file a function is 
##         developped to compute matrices of sparse resultants.
##	8) Geometry of curves and surfaces: this part regroups some usefull
##         functions to deal with implicit curves and surfaces in 2D or 3D.
##
## NB:The "warnings" obtained when loading multires are due to the Maple 
## packages linalg and simplex which are used here.
## 
## ----------------------------------------------------------------------
with(linalg): with(combinat,choose):
__ := NULL: addtohelp := proc() global __; __ := __, args: end:


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# POLYNOMIAL MANIPULATIONS  (Coefficients and matrices formulations)     #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

#List the monomials of a polynomial p in the variables var:
lstmonof := proc(p,var)
local lm,r,c;
  lm := NULL;  
  r := sort(p,var,tdeg);
  while r <> 0 do
     c := tcoeff(r,var,'m');
     lm := lm, m;
     r := expand(r-c*m);
  od;
  lm; 
end:

listofallmon := proc(var,n)
  lstmonof(expand((1+convert(var, `+`))^n),var); 
end:

#----------------------------------------------------------------------#
#return of monomials of degree nu in the variable var:

listmonhomdegnu := proc(var,nu)
	{lstmonof(expand((convert(var, `+`))^nu),var)};
end:

#----------------------------------------------------------------------#
#return the coeff. of p in variables var of the monomial m:
coeffof := proc(m,p,var)
local lm,lc,k;
  lc := [coeffs(p,var,'lm')];
  lm := [lm];
  if member(m,lm,'k') then lc[k] else 0 fi;
end:

#----------------------------------------------------------------------
addtohelp(coefofexp):
#----------------------------------------------------------------------
## HELP coefofexp
##
## coefofexp - monomial coefficient in a polynomial.
##
## CALLING SEQUENCE:
##     coefofexp(p, v, d)
##
## PARAMETERS:
##     p  - polynomial
##     v  - array of variables
##     d  - array of integers
## DESCRIPTION:
## - Compute the coefficient of p corresponding to v^d.  
##   Arrays v and d must be 1-dimensional and of the same size.
##
## EXAMPLES:
##> coefofexp(x^2-y*z+y*x-1,array(1..3,[x,y,z]),array(1..3,[0,1,1]));
##
##                              -1
##
##> coefofexp(x^2-y*z+y*x-1,array(1..3,[x,y,z]),array(1..3,[1,0,2]));
##
##                              0
## SEE ALSO:
## coeffs
##
coefofexp := proc(p::polynom, v::array, d::array)
  local i, ans, sz;
  if (not(evalb(`toric/type1Array`(v,name)))) then
    ERROR(`Coeffofmon requires second argument`,v,`be an array of names`);
  fi;
  if (not(evalb(`toric/type1Array`(d,integer)))) then
    ERROR(`Coeffofmon requires third argument`, d,`be an array of ints`);
  fi;

  sz := nops(convert(d,list));
  if (not(evalb(nops(convert(v,list))=sz))) then
    ERROR(`Coeffofmon requires size(variables)=size(degrees)`);
  fi;

  ans := p;
  for i from 1 to sz do
    ans := coeff(ans,v[i],d[i]);
  od;

  RETURN(eval(ans));

end:    

#----------------------------------------------------------------------#
addtohelp(coeffsof):
## ----------------------------------------------------------------------
## HELP coeffsof
##
## coeffsof - list of coefficient of a polynomial.
##
## CALLING SEQUENCE:
##     coeffsof(p)
##     coeffsof(p, var)
##     coeffsof(p, lm)
##     coeffsof(p, 'l')
##
## PARAMETERS:
##     p  - polynomial
##     var - list of a list of variables
##     lm  - list of monomials
##     l   - the name of a variable 
## DESCRIPTION:
## - Compute the coefficients  of the polynomial p.
##
## - If the second argument is a list of list var, it is the coefficient
##   coeff with respect to all the monomials in the variables var.
## 
## - If the second argument is a simple list of monomials, it is the
##   coefficient  coeff with respect to this list of monomials. The other
##   coefficients are ignored.
## 
## - If the second argument is the name of variable, the monomials
##   indexing the coefficients.
## 
## - The default value corresponds to the case var.
##
## EXAMPLES:
##> readlib(multires);
##> coeffsof(a^2-b*c+b-1); 
## [a^2, b*c, b, 1]
##
##                              [1, -1, 1, -1]
##
##> coeffsof(a^2-b*c+b-1,[a,b],'l'); l;
##
##                              [1, -c + 1, -1]
##
##                                  2
##                                [a , b, 1]
##
##> coeffsof(a^2-b*c+b-1,[[1,b,a*b]]); 
## 
##                              [-1, -c + 1, 0]
## SEE ALSO:
## coeffs
## 
coeffsof := proc(p,lx,l)
local var,lm,i,j,ll,pe;
  if nargs >1 then 
     if type(args[2],listlist) then 
        var := sort([op(indets(op(args[2])))]); 
        lm  := op(args[2]);
     elif type(args[2],list) then 
	var := args[2]; 
	lm := lstmonof(p,var);
	lm := [op(sort(convert([op({lm})],`+`),var))];
	if nargs>2 then l:=lm; else lprint(lm); fi;
     else 
	var := sort([op(indets(p))]); 
	lm := lstmonof(p,var);
	lm := [op(sort(convert([op({lm})],`+`),var,tdeg))];
	lx :=lm;
     fi;
  else 
    var := sort([op(indets(p))]);  
    lm  := NULL;
    lm  := lstmonof(p,var,tdeg); 
    lm  := [op(sort(convert([op({lm})],`+`),var))]; 
    lprint(lm);
  fi;
  pe := expand(p);
  ll := NULL;
  for i to nops(lm) do
     ll := ll, coeffof(lm[i],pe,var);
  od;
  [ll];
end:


# ----------------------------------------------------------------------
addtohelp(matrixof):
## ----------------------------------------------------------------------
## HELP matrixof
##
## matrixof - Coefficient matrix of a list of polynomials
##
## CALLING SEQUENCE:
##     matrixof(lp)
##     matrixof(lp, var)
##     matrixof(lp, lm)
##     matrixof(lp, 'l')
##
## PARAMETERS:
##     lp  - list of polynomials
##     var - list of a list of variables
##     lm  - list of monomials
##     l   - the name of a variable 
## DESCRIPTION:
## - Compute the coefficient matrix of lp.
##
## - If the second argument is a list of list var, it is the coefficient
##   matrix with respect to all the monomials in the variables var.
## 
## - If the second argument is a simple list of monomials, it is the
##   coefficient  matrix with respect to this list of monomials. The other
##   coefficients are ignored.
## 
## - If the second argument is the name of variable, the monomials
##   indexing the columns are stored in this variable.
## 
## - The default value corresponds to the case var.
##
## EXAMPLES:
##> readlib(multires);
##> matrixof([a^2,b^2-1,a*b-1]); 
##
##                              [1     0     0]
##                              [             ]
##                              [0     0     1]
##                              [             ]
##                              [0     1     0]
##                              [             ]
##                              [0    -1    -1]
##
##> matrixof([a^2,b^2-c,a*b-c^2],[a,b]); 
##
##                             [1    0      0 ]
##                             [              ]
##                             [0    0      1 ]
##                             [              ]
##                             [0    1      0 ]
##                             [              ]
##                             [             2]
##                             [0    -c    -c ]
##
##> matrixof([a^2,b^2-1,a*b-1],[[1,a^2,a*b]]); 
##
##                              [0    -1    -1]
##                              [             ]
##                              [1     0     0]
##                              [             ]
##                              [0     0     1]
##
##> matrixof([a^2,b^2-1,a*b-1],'lm'); lm;	
##
##                              [1     0     0]
##                              [             ]
##                              [0     0     1]
##                              [             ]
##                              [0     1     0]
##                              [             ]
##                              [0    -1    -1]
##
##                               2        2
##                             [a , a b, b , 1]
##
## SEE ALSO:
## deltaof
## 
matrixof := proc(lp ::list,lx)
local var,lm,i,j,ll;
  if nargs >1 then 
     if type(args[2],listlist) then 
        var := sort([op(indets(op(args[2])))]); 
        lm  := op(args[2]);
     elif type(args[2],list) then 
	var := args[2]; 
	lm := NULL;
	for i to nops(lp) do lm := lm, lstmonof(lp[i],var);  od;
	lm := [op(sort(convert([op({lm})],`+`),var))];
     else 
	var := sort([op(indets(lp))]); 
	lm := NULL;
	for i to nops(lp) do lm := lm, lstmonof(lp[i],var);  od;
	lm := [op(sort(convert([op({lm})],`+`),var,tdeg))];
	lx :=lm;
     fi;
  else 
    var := sort([op(indets(lp))]);  
    lm := NULL;
    for i to nops(lp) do lm := lm, lstmonof(lp[i],var,tdeg);  od;
    lm :=[op(sort(convert([op({lm})],`+`),var))]; 
    lprint(lm);
  fi;
  ll := NULL;
  for i to nops(lm) do
    for j to nops(lp) do
     ll := ll, coeffof(lm[i],expand(lp[j]),var);
    od;
  od;
  matrix(nops(lm),nops(lp),[ll]);
end:
# ----------------------------------------------------------------------
matrixtof := proc() transpose(matrixof(args)) end:
# ----------------------------------------------------------------------
exponentsof := proc(p,v,l)
local lm,lc,ll,var,i,j,M;
  if nargs >1 then 
     if type(args[2],listlist) then
	var := sort([op(indets(p))]);  
	lc := coeffsof(p,args[2],lm);
        if nargs>2 then l := lc fi;
     elif type(args[2],list) then 
        var := args[2];
	lc := coeffsof(p,args[2],lm);
        if nargs>2 then l := lc fi;
     else
        var:= sort([op(indets(p))]);
	lc := coeffsof(p,var,lm);
        v := lc;
     fi;
  else 
     var := sort([op(indets(p))]);  
     lc := coeffsof(p,var,lm);
  fi;
  M := matrix(nops(var),nops(lm));
  for i to nops(var) do
    for j to nops(lm) do
       M[i,j] := degree(lm[j],var[i]);
    od;
  od;
  evalm(M);
end:

# ----------------------------------------------------------------------
addtohelp(Delta):
## ----------------------------------------------------------------------
## HELP Delta
##
## Delta - Compute the matrix associated with a polynomial in x and y.
##
## CALLING SEQUENCE:
##     Delta(p)
##     Delta(p, X)
##     Delta(p, X, l1, l2)
##
## PARAMETERS:
##     p     - polynomial in var and _y[i]
##     X     - list of variables
##     l1,l2 - variables used to get the list of monomials which index the 
##             rows and columns of the output matrix.
##
## DESCRIPTION:
## - Compute the matrix of the coefficients t[alpha,beta] in the decomposition
##   p= sum t[alpha,beta] x^alpha y^beta
##
## - It is assumed that the variables _y[i] are indexed from 1 to i = nops(X).
##
## EXAMPLES:
##> Delta(a^2+2*a*b+b^2);
##[a*b, a^2, b^2]
##[1]
##                                     [2]
##                                     [ ]
##                                     [1]
##                                     [ ]
##                                     [1]
## SEE ALSO:
##    matrixof
##
Delta := proc(P,var,lx,lxi)
# dev. of P in the basis of mon. in X and _y[i];
local i,j,vx,vxi,c,nbase,lbase,lcx,nx,mxi,nxi,mtr,mr,lambda,v,lr,dt;

  if nargs <2 then 
    vx  := sort([op(indets(P))]);
  else 
    vx := args[2];
  fi;
  for i from 1 to nops(vx) do vx := subs(_y[i]=NULL,vx) od;
  vxi := [_y[j]$j=1..nops(vx)];

  dt :=collect(expand(P),vx,distributed);
  lcx:= [coeffs(sort(dt),vx,'ll')];

  if nargs>2 then lx:=[ll] else lprint([ll]) fi;

  nx := nops(lcx);
  nxi := 0;
  mr := [];
  mxi := NULL;

  for j from 1 to nops(lcx) do 
    c := lcx[j];
    lambda := [coeffs(c,vxi,'v')];
    v := [v];
    if  nxi=0 
        then lr :=[]
        else lr := convert(vector(nxi,0),list)
    fi;

    for i from 1 to nops(v) do 
       if assigned(lbase[v[i]]) 
          then 
	       lr  := subsop(lbase[v[i]] = lambda[i], lr);
          else 
	       lr  := [op(lr),lambda[i]];
	       nxi := nxi+1; 
	       lbase[v[i]] := nxi;
	       mxi := mxi,v[i];
       fi;
    od;
    mr := [op(mr),lr]; 
  od; 

  if nargs >3 then lxi := [mxi] else lprint([mxi]); fi;

  for  i from 1 to nops(mr) do
    if nops(mr[i]) <> nxi then 
       mr := subsop(i= [op(mr[i]), 
        op(convert(vector(nxi-nops(mr[i]), 0),list))],mr);
    fi;
  od;
  matrix(nx,nxi,mr);
end:

#-------------------------------------------------------
deltaof := proc(Th,lx ::list,lxi::list)
local var,lm,i,j,ll,c,vx,vxi,p;
  vx  := sort([op(indets(lx))]);
  vxi := sort([op(indets(lxi))]);
  ll := NULL;
  p:= expand(Th);
  for i to nops(lx) do
    c := coeffof(lx[i],p,vx);
    for j to nops(lxi) do
     ll := ll, coeffof(lxi[j],c,vxi);
    od;
  od;
  matrix(nops(lx),nops(lxi),[ll]);
end:


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#   RESULTANT MATRIX CONSTRUCTION                                      #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

## --------------------------------------------------------------------
addtohelp(Jd):
## --------------------------------------------------------------------
## HELP Jd
##
## Jd - discrete jacobian matrix 
##
## CALLING SEQUENCE:
##     Jd(lp)
##     Jd(lp, var)
##
## PARAMETERS:
##     lp  - list of polynomials
##     var - list of variables
##
## DESCRIPTION:
## - Compute the discrete Jacobian matrix with respect to the variables var.
##
## - The default for var is the indeterminates of lp.
## 
## - The other set of variables is {_y[1], _y[2], ...}
## 
## EXAMPLES:
## > Jd([a^2,b^2-1,a*b-1]);
##                     [   2                             ]
##                     [  a        a + _y[1]        0    ]
##                     [                                 ]
##                     [ 2                               ]
##                     [b  - 1         0        b + _y[2]]
##                     [                                 ]
##                     [a b - 1        b          _y[1]  ]
##
## SEE ALSO:
## Theta
##
Jd := proc(lp::list) 
local N,i,j,k,vx,vxi,mtr,s;
  if nargs <2 then 
    vx := sort(convert(indets(lp),list));
  else
    vx := args[2]
  fi;
  N := nops(lp);
  if nops(vx) <> N-1 then 
    ERROR(`the number of polynomials must be the number of variables +1 `) 
  fi; 
  mtr := array(1..N,1..N);
  for i from 1 to N do mtr[i,1] := lp[i] od;
  for j from 2 to N do 
     for i from 1 to N do
       mtr[i,j]:= subs(vx[j-1]=_y[j-1],mtr[i,j-1])
     od;
  od;
  for j from N to 2 by -1 do 
     for i from 1 to N do
       mtr[i,j] := quo(mtr[i,j]-mtr[i,j-1], _y[j-1] - vx[j-1], _y[j-1]);
     od;
  od;
  #dt :=det(convert(mtr,matrix));
  #collect(expand(dt),vx,distributed);
  evalm(mtr);
end:

## ----------------------------------------------------------------------
addtohelp(Theta):
## ----------------------------------------------------------------------
## HELP Theta
##
## Theta - Compute the discrete jacobian
##
## CALLING SEQUENCE:
##     Theta(lp)
##     Theta(lp, var)
##
## PARAMETERS:
##     lp  - list of polynomials
##     var - list of variables
##
## DESCRIPTION:
## - Compute the discrete Jacobian with respect to the variables var.
##
## - The default for var is the indeterminates of lp.
## 
## - Uses the function Jd.
## 
## EXAMPLES:
##> Theta([a^2,b^2-1,a*b-1]);
##  2      2                  2
##-b  _y[1]  + _y[1] a + _y[1]  + a b _y[1] _y[2] - a b - a _y[2] - _y[1] b
##
##   - _y[1] _y[2]
##
## SEE ALSO:
## Jd
Theta := proc(lp::list)
 det(Jd(args)):
end:

## ----------------------------------------------------------------
addtohelp(mbezout):
## ----------------------------------------------------------------
## HELP mbezout
##
## mbezout - Compute the Bezoutian matrix of a list of polynomials
##
## CALLING SEQUENCE:
##     mbezout(lp)
##     mbezout(lp, var)
##     mbezout(lp, var, l1)
##     mbezout(lp, var, l1, l2)
##
## PARAMETERS:
##     lp  - list of polynomials.
##     var - list of variables (optional).
##     l1  - unassigned variable where are stored the monomials indexing 
##           the rows of the matrix (optional).
##     l2  - unassigned variable where are stored the monomials indexing 
##           the columns of the matrix (optional).
##
## DESCRIPTION:
## - Compute the bezoutian matrix of lp with respect to the variables var. 
##
## - The default value for var is the set of indeterminates of lp.
## 
## - Uses the function Theta.
##
## EXAMPLES:
##> readlib(multires);
##> mbezout([a^2,b^2-1,a*b-1]); 
## [1, b, a, a*b, b^2]
## [_y[1]*_y[2], _y[1]^2, _y[1], _y[2], 1]
## 
##                          [-1     1     0     0     0]
##                          [                          ]
##                          [ 0     0    -1     0     0]
##                          [                          ]
##                          [ 0     0     1    -1     0]
##                          [                          ]
##                          [ 1     0     0     0    -1]
##                          [                          ]
##                          [ 0    -1     0     0     0]
## 
##> mbezout([a^2-x,b^2-y,a*b-z],[a,b]); 
## [1, b, a, a*b, b^2]
## [_y[1]*_y[2], _y[1]^2, _y[1], _y[2], 1]
## 
##                          [-z    y     0     0     0 ]
##                          [                          ]
##                          [0     0     -z    x     0 ]
##                          [                          ]
##                          [0     0     y     -z    0 ]
##                          [                          ]
##                          [1     0     0     0     -z]
##                          [                          ]
##                          [0     -1    0     0     x ]
## 
## SEE ALSO:
## Jd, Theta
##
mbezout := proc(lpol::list)
# lpol := [p1, ..., pn];
local var,M,t;
  if nargs <2 then 
    var  := sort([op(indets(lpol))]);
  else 
    var := args[2];
  fi;
  t := Theta(lpol,var);
  if nargs>3 then   
     Delta(t,var,args[3],args[4]);
  elif nargs>2 then
     Delta(t,var,args[3]);
  else
     Delta(t,var);
  fi;
end:

## ----------------------------------------------------------------------
addtohelp(gbezout):
## ----------------------------------------------------------------
## HELP gbezout
##
## gbezout - Compute the Bezoutian matrix, of generic polynomials with
##           support in the list of monomial lm.
##
## CALLING SEQUENCE:
##     gbezout(lm)
##     gbezout(lm, l1)
##     gbezout(lm, l1, l2)
##
## PARAMETERS:
##     lm  - list of monomials
##     l1  - unassigned variable where are stored the monomials indexing 
##           the rows of the matrix (optional).
##     l2  - unassigned variable where are stored the monomials indexing 
##           the columns of the matrix (optional).
##
## DESCRIPTION:
## - Compute the bezoutian matrix of generic polynomials with support in lm.
## 
## - The output bezoutian matrix involves the n+1 by n+1 determinants 
##   Xi(i0, ...,in) of the generic coefficients. Xi(i0, ...,in) is the
##   determinant  of the coefficients of the monomials of index i0,...,in of
##   the list lm in the generic polynomials.
## 
## - The number of variables is the number of indeterminates, which appear in
##   the list lm.
## 
## - It uses the function Delta.
##
## EXAMPLES:
##> readlib(multires);
##> gbezout([1,a,b,a^2,a*b,b^2]); 
##
## SEE ALSO:
##  mbezout
gbezout := proc(lm::list)
local i,j,r,x,var,la,l,n;
  var := sort([op(indets(lm))]);
  n := nops(var);
  l := choose([i$i=1..nops(lm)],n+1); 
  r := 0;
  for x in l do 
   la := NULL;
   for j to nops(x) do la := la, lm[x[j]] od;
   r := r + Xi(op(x))*Theta([la],var);
  od:
   Delta(r,var,args[2..nargs]);
end:

## ----------------------------------------------------------------------
addtohelp(mresultant):
## ----------------------------------------------------------------------
## HELP mresultant
##
## mresultant - Macaulay resultant matrix of a list of polynomials
##
## CALLING SEQUENCE:
##     mresultant(lp)
##     mresultant(lp, var)
##
## PARAMETERS:
##     lp  - list of polynomials
##     var - list of variables
##
## DESCRIPTION:
## - Compute the Macaulay resultant matrix of lp.
##
## - The second argument can be a list of variables. The number of 
##   polynomials in lp should one more than in var.
## 
## - The default value corresponds to the case where var is the set of 
##   indeterminates in lp.
## 
## - The determinant of this matrix is a multiple of the resultant of the
##   polynomials lp OVER P^n. If the determinant is zero, either there is 
##   common root in PP^n, or there are obvious relations between the 
##   polynomials.
##
## - The size of the matrix is the number of monomials of degree less than
##   (d1-1)+ ...+(dn+1-1) +1 where di is the degree of lp[i] and n is the
##   number of polynomials in lp.
##
## EXAMPLES:
##> readlib(multires);
##> mresultant([a+1,2*b^2+2,3*a^2-3]); 
## [1, b, a, a*b, b^3, a^3, a*b^2, a^2*b, a^2, b^2]
##
##            [1    0    0    0    2    0    0    -3     0     0]
##            [                                                 ]
##            [0    1    0    0    0    2    0     0     0    -3]
##            [                                                 ]
##            [1    0    1    0    0    0    2     0    -3     0]
##            [                                                 ]
##            [0    1    0    1    0    0    0     0     0     0]
##            [                                                 ]
##            [0    0    0    0    0    2    0     0     0     0]
##            [                                                 ]
##            [0    0    0    0    0    0    0     0     3     0]
##            [                                                 ]
##            [0    0    0    0    0    0    2     0     0     0]
##            [                                                 ]
##            [0    0    0    1    0    0    0     0     0     3]
##            [                                                 ]
##            [0    0    1    0    0    0    0     3     0     0]
##            [                                                 ]
##            [0    0    0    0    2    0    0     0     0     0]
##
##> mresultant([u*a+v,b^2-1,a^2+1],[a,b]); 
## [1, b, a, a*b, b^3, a^3, a*b^2, a^2*b, a^2, b^2]
##
##            [v    0    0    0    -1     0     0    1    0    0]
##            [                                                 ]
##            [0    v    0    0     0    -1     0    0    0    1]
##            [                                                 ]
##            [u    0    v    0     0     0    -1    0    1    0]
##            [                                                 ]
##            [0    u    0    v     0     0     0    0    0    0]
##            [                                                 ]
##            [0    0    0    0     0     1     0    0    0    0]
##            [                                                 ]
##            [0    0    0    0     0     0     0    0    1    0]
##            [                                                 ]
##            [0    0    0    0     0     0     1    0    0    0]
##            [                                                 ]
##            [0    0    0    u     0     0     0    0    0    1]
##            [                                                 ]
##            [0    0    u    0     0     0     0    1    0    0]
##            [                                                 ]
##            [0    0    0    0     1     0     0    0    0    0]
##
## SEE ALSO:
##    mbezout, spresultant
##
mresultant := proc(l ::list,vv,lll)
local n,nu,ld,i,j,var,svar, Si,S,lp,lm,m,ll,ltmp;
  n := nops(l)-1;
  if nargs >1 then var := args[2] else  var := sort([op(indets(l))]); fi;
  svar := {op(var)};
  if(nops(var)<>(nops(l)-1)) then
     ERROR(`the number of polynomials in the list must be the number of variables +1 `);
  fi;
  ld := NULL;
  for i to nops(l) do
   ld := ld,degree(l[i],svar);
  od;
  ld := [ld]; 
  nu := convert(ld,`+`)-nops(l)+1;
  lp := NULL;
  S := {listofallmon(var,nu)};
  lm := [op(S)];
  for i from nops(l) to 2 by -1 do
    Si := {lstmonof(expand(var[i-1]^ld[i]*(1+convert(var, `+`))^(nu-ld[i])),var)};
    Si := Si intersect S;
    for m in Si do
	 lp :=expand(m/var[i-1]^ld[i]*l[i]),lp;
     od;
     S := S minus Si; 
  od;
  S := [lstmonof(convert(S,`+`),var)];  
  ltmp := NULL;
  for m in S do 
      ltmp := ltmp,expand(m*l[1]);
  od;

  lp := [ltmp,lp]; 
  lm := [op(S), op({op(lm)} minus {op(S)})]; 
  if nargs>2 then lll:=lm else lprint(lm); fi;
  ll := NULL;
  for i to nops(lm) do
    for j to nops(lp) do
     ll := ll, coeffof(lm[i],lp[j],var);
    od;
  od;
  matrix(nops(lm),nops(lp),[ll]);
end:

#-------------------------------------------------------
mresultantlist := proc(l ::list,vv,lll)
local n,nu,ld,i,j,var,svar, Si,S,lp,lm,m,ll,ltmp;
  n := nops(l)-1;
  if nargs >1 then var := args[2] else  var := sort([op(indets(l))]); fi;
  ld := NULL; 
  if(nops(var)<>(nops(l)-1)) then
     ERROR(`the number of polynomials in the list must be the number of variables +1 `);
  fi;
  svar := {op(var)};
  for i to nops(l) do
   ld := ld,degree(l[i],svar);
  od;
  ld := [ld];
  if nargs>2 then lll:=ld else lprint(ld); fi;
  nu := convert(ld,`+`)-nops(l)+1;lprint(nu);
  lp := NULL;
  S := {listofallmon(var,nu)};
  lm := [op(S)];
  for i from nops(l) to 2 by -1 do
    Si := {lstmonof(expand(var[i-1]^ld[i]*(1+convert(var, `+`))^(nu-ld[i])),var)};
    Si := Si intersect S;
    S := S minus Si; 
  od;
  S := [lstmonof(convert(S,`+`),var)];  
  [op(S), op({op(lm)} minus {op(S)})];
end:

## ----------------------------------------------------------------------
addtohelp(jresultant):
## ----------------------------------------------------------------------
## HELP jresultant
##
## jresultant -  An inertia form of a list of polynomials
##
## CALLING SEQUENCE:
##     jresultant(lp,var)
##     jresultant(lp,var,mu)
##
## PARAMETERS:
##     lp  - list of polynomials
##     var - list of variables
##     mu  - integer
##
## DESCRIPTION:
## - Compute the matrix introduces by J.P. Jouanolou in "Formes d'inertie
##   et resultant : un formulaire" of lp.
##
## - The second argument is a list of variables. The number of 
##   polynomials in lp should one more than in var. The third argument 
##   is an integer with value between 0 and nu (= the sum of the degree minus
##   one of all polynomials).
## 
## - The determinant of this matrix is a multiple of the resultant of the
##   polynomials lp OVER P^n. If the determinant is zero, either there is 
##   a common root in P^n, or there are obvious relations between the 
##   polynomials.
##
## - The size of the matrix is the sum of the number of monomials of degree
##   nu and the number of monomials "d-repu" of degree d-mu.
##  
##
## EXAMPLES:
##> readlib(multires);
##> jresultant([a*b*u+v,b^2-1,a^2],[a,b]); 
##
##                     [0    -u    -v    0     0    0 ]
##                     [                              ]
##                     [0    -v    0     0     0    0 ]
##                     [                              ]
##                     [0    0     0     -v    0    0 ]
##                     [                              ]
##                     [0    0     0     u     0    v ]
##                     [                              ]
##                     [0    0     0     0     1    -1]
##                     [                              ]
##                     [1    0     0     0     0    0 ]
##
## SEE ALSO:
##    mbezout, mresultant, spresultant 
## 
jresultant := proc (lp::list,var::list,mu::integer) 
local  i,j,k,n,nu,varx,Var,Varx,vary,varxy,debut,dindex,
       d,delta,X,repnu,Repnu,repdnu,Repdnu,monnu,Monnu,
       mondnu,Mondnu,N,morl,J,f,s,t,L,C;   
	
	n:=nops(lp);	               
	Var:={seq(var[i],i=1..n-1)};
	varx:=[seq(_x[i],i=1..n)];
	Varx:={seq(_x[i],i=1..n-1)};
	vary:=[seq(_y[i],i=1..n)];
	varxy:=[seq(_x[i],i=1..n),seq(_y[i],i=1..n)];
	d:=[seq(degree(lp[i],Var),i=1..n)]; 
	
	for i from 1 to n do
	       f[i]:=expand(subs(seq(var[j]=_x[j],j=1..n-1),lp[i]));   
	       s:={lstmonof(f[i],varx)};
	       t:=nops(s);
	       f[i]:=sum('coeffof(s[j],f[i],varx)*s[j]*_x[n]^(d[i]-degree(s[j],Varx))','j'=1..t);
	od;
		
	delta:=sum('d[i]-1','i'=1..n);
	if (nargs < 3) then nu:=iquo(delta,2) else nu:=mu fi;
		         
	X:=convert(varx,`+`);
	for i from 1 to n do 
		if (nu<d[i]) then repnu[i]:={};
			    else repnu[i]:={lstmonof(expand((_x[i]^d[i])*(X^(nu-d[i]))),varx)};
		fi;
	od;
	for i from 2 to n do 
		for j from 1 to i-1 do repnu[i]:=repnu[i] minus repnu[j] od;
	od;
	Repnu:=sum('nops(repnu[i])','i'=1..n);
	for i from 1 to n do 
		if (delta-nu<d[i]) then repdnu[i]:={};
			    else repdnu[i]:={lstmonof(expand((_x[i]^d[i])*(X^(delta-nu-d[i]))),varx)};
		fi;
	od;
	for i from 2 to n do 
		for j from 1 to i-1 do repdnu[i]:=repdnu[i] minus repdnu[j] od;
	od; 
	Repdnu:=sum('nops(repdnu[i])','i'=1..n);
	
	monnu:=listmonhomdegnu(varx,nu);
	Monnu:=nops(monnu);
	mondnu:=listmonhomdegnu(varx,delta-nu);
	Mondnu:=nops(mondnu);
	
	N:=Monnu+Repdnu;
	
	J:=array(1..N,1..N);L:=[];C:=[];

	#haut gauche
	morl:=morley([seq(f[i],i=1..n)],varx);
	for i from 1 to Monnu do 
		for j from 1 to Mondnu do 
			J[i,j]:=coeffof(subs(seq(_x[l]=_y[l],l=1..n),monnu[i])*mondnu[j],morl,varxy);
		od; 
	od;
	#low left
	for dindex from 1 to n do 
		debut:=Monnu+sum('nops(repdnu[i])','i'=1..(dindex-1));
		for i from 1 to nops(repdnu[dindex]) do
			for j from 1 to Mondnu do
				J[debut+i,j]:=coeffof(mondnu[j],expand(repdnu[dindex][i]/(_x[dindex]^d[dindex])*f[dindex]),varxy);
 				
			od;
		if dindex < n then
			for k from dindex+1 to n do
				if degree(repdnu[dindex][i],_x[k]) >= d[k] then L:=[op(L),debut+i]; break; fi;
			od;
		fi;
		od; 		
	od;
	#top right
	for dindex from 1 to n do
		debut:=Mondnu+sum('nops(repnu[i])','i'=1..(dindex-1));
		for i from 1 to nops(repnu[dindex]) do
			for j from 1 to Monnu do
				J[j,debut+i]:=coeffof(monnu[j],expand(repnu[dindex][i]/(_x[dindex]^d[dindex])*f[dindex]),varxy);
			od;
		if dindex < n then
			for k from dindex+1 to n do
				if degree(repnu[dindex][i],_x[k]) >= d[k] then C:=[op(C),debut+i]; break; fi;
			od;
		fi;
		od; 		
	od;
      	
	for i from Monnu+1 to N do
		for j from Mondnu+1 to N do J[i,j]:=0; od;
	od;

	RETURN(evalm(J));

end:	


## ----------------------------------------------------------------------
addtohelp(morley):
## ----------------------------------------------------------------------
## HELP morley
##
## morley -  The matrix of a sort of bezoutian for homogeneous polynomials
##
## CALLING SEQUENCE:
##     morley(lp,var)
##
## PARAMETERS:
##     lp  - list of homogeneous polynomials
##     var - list of variables
##
## DESCRIPTION:
## - Compute the matrix of the morley's form introduces by J.P. Jouanolou in
##   "Formes d'inertie et resultant : un formulaire" of lp.
##   The new set of variables are _y[i].
##
## - The second argument is a list of variables. The number of 
##   polynomials in lp should be the same than the length of var.
## 
## EXAMPLES:
##> readlib(multires);
##> morley([a*b-c^2,b^2-c^2,a^2],[a,b,c]);
##       2              2
## -_y[1]  _y[3] - _y[1]  c + _y[1] _y[3] _y[2] + _y[1] _y[3] b + _y[1] c _y[2]
## 
##      + _y[1] b c - a _y[1] _y[3] - a _y[1] c + a _y[3] _y[2] + a _y[3] b
## 
##      + a c _y[2] + a b c
## 
## SEE ALSO:
## mbezout, mresultant
##
morley := proc(lp::list,varx::list) 
local  i,j,l,s,n,p,f,mx,d,Mondx,H,aj,pd,morl;      
        
	n:=nops(lp);
	d:=[seq(degree(lp[i],varx),i=1..n)];
	Mondx:=[seq(listmonhomdegnu(varx,d[i]),i=1..n)]; 
	H:=array(1..n,1..n);
	for i from 1 to n do 
		for j from 1 to n do
			H[i,j]:=0;
			for s from 1 to nops(Mondx[i]) do  
				mx:=Mondx[i][s];  
				f:=0;	
				aj:=subs(seq(varx[l]=1,l=1..n),diff(mx,varx[j]));  
				if (aj>0) then pd:=sum('varx[j]^p*_y[j]^(aj-1-p)','p'=0..aj-1);
                                               f:=subs(seq(varx[l]=1,l=j..n),seq(varx[l]=_y[l],l=1..j),mx)*subs(seq(varx[l]=1,l=1..j),mx)*pd;
			        fi;				
				H[i,j]:=H[i,j]+f*coeffof(mx,lp[i],varx);	
                        od;  
               od;
	od;
	RETURN(expand(det(H)));
end:

## ----------------------------------------------------------------------
addtohelp(bkmresultant):
## ----------------------------------------------------------------------
## HELP bkmresultant
##
## bkmresultant - Compute a resultant matrix for a residual intersection
##
## CALLING SEQUENCE:
##     bkmresultant(lp,M,var,reg)
##
## PARAMETERS:
##     lp  - list of homogeneous polynomials
##     M   - matrix
##     var - list of variables
##     reg - integer
##
## DESCRIPTION:
##
## - Compute the matrix of the first application in the resolution of (I:J) 
##   given in the article of Bruns, Kustin and Miller ( bkm ) :
##   "The resolution of the generic residual intersection of a 
##     complete intersection", Journal of Algebra 128.
##
## - The first argument is a list of homogeneous polynomials I=(f1,..,fm). 
##   Given a homogeneous complete intersection J=(g1,..,gn), such that I 
##   is inclued in J and (I:J) is a residual intersection,
##   we want to compute a sort of resultant of (I:J). 
##   The matrix M is the matrix such that I=J.M.
##   The integer reg must be superior or equal to the regularity of the
##   ideal (I:J).
##
## - The result of bkmresultant is a surjective nxm matrix such that the 
##   determinant of a nxn minor is a multiple of the resultant of I on 
##   the closure of V(I)\V(J). This minor can be obtain with "hmaxminor".
## 
##
## EXAMPLES:
##> bkmresultant([a[0]*z^2+a[1]*x*z+a[2]*y*z+a[3]*(x^2+y^2),b[0]*z^2+b[1]*x*z+b[2]*y*z+b[3]*(x^2+y^2),c[0]*z^2+c[1]*x*z+c[2]*y*z+c[3]*(x^2+y^2)],matrix([[a[0]*z+a[1]*x+a[2]*y,b[0]*z+b[1]*x+b[2]*y,c[0]*z+c[1]*x+c[2]*y],[a[3],b[3],c[3]]]),[x,y,z],2); 
##>hmaxminor(%);
##
## SEE ALSO:
##    mbezout, mresultant, spresultant, jresultant, bkmdegree
## 
bkmresultant:= proc ( lp::list, M::matrix , var::list , reg::integer)
local i,j,k,m,n,combinaison,BKM,bkm,Base,base,Var,N,P,nu,monomes;

	m:=coldim(M); n:=rowdim(M);
	Var:={seq(var[i],i=1..m)};
	combinaison:=choose([seq(i,i=1..m)],n); #liste pour les determinants de Phi.
	Base:=listmonhomdegnu(var,reg); base:=nops(Base);
	BKM:=[];
	
	#Macaulay part

	for i from 1 to m do
		nu:=reg-degree(lp[i],Var);
		monomes:=listmonhomdegnu(var,nu);
		bkm:=matrix(base,nops(monomes));
		for j from 1 to nops(monomes) do
			for k from 1 to base do	
				bkm[k,j]:=coeffof(Base[k],expand(lp[i]*monomes[j]),var);
			od;
		od;
		BKM:=[op(BKM),eval(bkm)];
	od;
	BKM:=concat(op(BKM));
	
	#Partie des det de Phi

	for i from 1 to nops(combinaison) do
		N:=submatrix(M,[seq(j,j=1..n)],combinaison[i]); #La sous-matrice de M pour i1,..,in
		P:=det(N); #le polynome a decomposer
		nu:=reg-sum('degree(N[j,j],Var)','j'=1..n);
		monomes:=listmonhomdegnu(var,nu);
		bkm:=matrix(base,nops(monomes)); #sous matrice de BKM
		for j from 1 to nops(monomes) do
			for k from 1 to base do	
				bkm[k,j]:=coeffof(Base[k],expand(P*monomes[j]),var);
			od;
		od;
		BKM:=[op(BKM),eval(bkm)];
	od;
	RETURN(concat(op(BKM)));
end:

## ----------------------------------------------------------------------
addtohelp(bkmdegree):
## ----------------------------------------------------------------------
## HELP bkmdegree
##
## bkmdegree - Compute the degree of the residual resultant bkmresultant 
##
## CALLING SEQUENCE:
##     bkmdegree(ld,lk)
##
## PARAMETERS:
##     ld  - list of integers
##     lk  - list of integers
##
## DESCRIPTION:
##
## - Calling with ld:=[d1,...,dm] and lk:=[k1,...,kn], bkmdegree gives
##   the degree of the bkmresultant in the coefficient of the polynomial
##   f_0.
##
## EXAMPLES:
##> readlib(multires);
##>bkmdegree([2,2,2],[1,1]);
##
## SEE ALSO:
##    bkmresultant
## 

bkmdegree:=proc(ld,lk)
local A,B,N,i,j,m,n,t;

         m:=nops(ld);
         n:=nops(lk);
         t:=m+1-n;
         A:=array(1..t,1..t);
         if t >= n+1 then
                for i from 1 to n+1 do
                     for j from 1 to i do
                         A[i,j]:=(-1)^(i+j)*sympol(lk,i-j);
                     od;
                     for j from i+1 to t do
                         A[i,j]:=0;
                     od;
                 od;    
                 for i from n+2 to t do
                     for j from 1 to i-n-1 do A[i,j]:=0; od;
                     for j from i-n to i do A[i,j]:=(-1)^(i-j)*sympol(lk,i-j); od;
                     for j from i+1 to t do A[i,j]:=0; od;
                 od;        
         else for i from 1 to t do
                     for j from 1 to i do
                         A[i,j]:=(-1)^(i+j)*sympol(lk,i-j);
                     od;
                     for j from i+1 to t do
                         A[i,j]:=0;
                     od;
                 od;  

         fi;
         if m > n then
                 B:=delrows(A,1..1);
                 N:=0;
                 for i from n to m do N:=N-sympol(ld,m-i)*det(delcols(B ,i+1-n..i+1-n)); od;
                 N:=N*sympol(lk,n)+sympol(ld,m);
                  else 
                 N:=sympol(ld,m)-sympol(lk,n);
         fi;

RETURN(N);
end:

#used in bkmdegree:
sympol:=proc(l,n)
   RETURN(coeff(expand(x*product(x-l[i],i=1..nops(l))), x^(nops(l)-n+1))*(-1)^n)
end:

## ----------------------------------------------------------------------
addtohelp(cm2resultant):
## ----------------------------------------------------------------------
## HELP cm2resultant
##
## bkmresultant - Compute a resultant matrix for a residual intersection
##
## CALLING SEQUENCE:
##     cm2resultant(H,R,var,reg)
##
## PARAMETERS:
##     H  - matrix
##     R   - matrix
##     var - list of variables
##     reg - integer
##
## DESCRIPTION:
##
## - Compute the first map of the complex which computes the residual
##   resultant of a local complete intersection Cohen-Macaulay of 
##   codimension two (see PhD of Laurent Buse).
##
##  - Given a homogeneous ideal l.c.i. CM codim 2 J=(g1,..,gn), such that 
##   I=(f1,..,fm) is inclued in J and (I:J) is a residual intersection, the
##   matrix H is such that I=J.H. The matrix R is the matrix of the first
##   syzygies of J. var denotes the variables and reg the critical degree.
##
## - The result of cm2resultant is a surjective  matrix such that the 
##   determinant of a maximal minor is a multiple of the resultant of I on 
##   the closure of V(I)\V(J). This minor can be obtain with "hmaxminor".
## 
##
## EXAMPLE:
##>H:=matrix([[b_0*z,c_1*z,a_2*z],[f_0*x,f_1*y,d_2*x]]);
##>Phi:=matrix([[-y*x],[x^2+y^2]]);
##>cm2resultant(H,Phi,[x,y,z],3);
##>factor(det(hmaxminor(%)));
##                        4                          2    2
##                -a_2 f_1  f_0  (-b_0 d_2 + a_2 f_0)  b_0
##
## SEE ALSO:
##   bkmresultant;
##
cm2resultant:= proc (H::matrix, R::matrix, var::list, reg::integer)
local m,n,Var,combinaison,Base,base,B,i,N,P,nu,monomes,Bt,j,k,M;

        M:=concat(H,R);
        m:=coldim(M); n:=rowdim(M);
        Var:={seq(var[i],i=1..nops(var))};
        combinaison:=choose([seq(i,i=1..m)],n); #liste pour les det.
        Base:=listmonhomdegnu(var,reg); base:=nops(Base);
        B:=[];
        for i from 1 to nops(combinaison) do
                N:=submatrix(M,[seq(j,j=1..n)],combinaison[i]);
                P:=expand(det(N));
                nu:=reg-degree(P,Var);
                monomes:=listmonhomdegnu(var,nu);
                Bt:=matrix(base,nops(monomes));
                for j from 1 to nops(monomes) do
                        for k from 1 to base do
                                Bt[k,j]:=coeffof(Base[k],expand(P*monomes[j]),var);
                        od;
                od;
                B:=[op(B),eval(Bt)];
        od;
        RETURN(concat(op(B)));
end:

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# MATRIX OPERATIONS                                                    #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

# return the sequence of non-zero terms of the last non zero row of
# a matrix.
lasts := proc(T)
# derniers coefficients non nul de la matrice.
local i,j,r,l;
 i := rowdim(T);
 r := 0;
 while i > 0 and r =0 do 
   j:=1;
   while T[i,j] = 0 and j < coldim(T) do j := j+1; od;
   r := T[i,j];
   i := i-1;
 od;
 if i>=0 then    lprint(i); RETURN(T[i+1,l]$l=j..coldim(T)) else RETURN(0); fi;
end:

## ----------------------------------------------------------------------
addtohelp(uschur, lschur, triang, shape):
## ----------------------------------------------------------------------
## HELP uschur
##
## uschur - compute the upper schur complement 
##
## CALLING SEQUENCE:
##     lschur(M,n)
##
## PARAMETERS:
##     M  - matrix
##     n  - integer
## DESCRIPTION:
## - Compute the upper triangular schur complement of size, inverting the
##   lower diagonal block of size coldim(M)-n.
##
## - The resulting matrix is of size n.
##
## EXAMPLES:
##> readlib(multires);
##> A := matrix(3,3,[0,1,0,0,0,1,1,1,1]);
##                             [0    1    0]
##                             [           ]
##                        A := [0    0    1]
##                             [           ]
##                             [1    1    1]
##> uschur(A,2);
##                                  [ 0     1]
##                                  [        ]
##                                  [-1    -1]
## SEE ALSO:
##  invert, lschur
##
uschur := proc(M,n) 
  evalm( submatrix(M,1..n,1..n) -  submatrix(M,1..n,n+1..coldim(M)) &*  inverse(submatrix(M,n+1..rowdim(M),n+1..coldim(M)))&* submatrix(M,n+1..rowdim(M),1..n));
end:

## ----------------------------------------------------------------------
## HELP lschur
##
## lschur - Compute the lower Schur complement.
##
## CALLING SEQUENCE:
##     lschur(M,n)
##
## PARAMETERS:
##     M  - matrix
##     n  - integer
## DESCRIPTION:
## - Compute the lower triangular schur complement of size n, inverting the
##   uppper diagonal block of size n
##
## - The resulting matrix is of size coldim(M)-n.
##
## EXAMPLES:
##> readlib(multires);
##> A := matrix(3,3,[0,1,0,1,0,1,1,1,1]);
##                             [0    1    0]
##                             [           ]
##                        A := [1    0    1]
##                             [           ]
##                             [1    1    1]
##> lschur(A,2);
##                                      [0]
## SEE ALSO:
##  invert, uschur
##
lschur := proc(M,n)
  evalm( submatrix(M,n+1..rowdim(M),n+1..coldim(M)) 
-  submatrix(M,n+1..rowdim(M),1..n) &*  inverse(submatrix(M,1..n,1..n))
  &* submatrix(M,1..n,n+1..coldim(M)));
end:

#----------------------------------------------------------------------#
nzt := proc(V)
local r,i,L;
 L := convert(V,list);
 r :=0;
 for i to nops(L) do
  if(L[i]<>0) then r := r+1 fi;
 od;
 r;
end:

idx := proc(V)
local r,i,L;
 L := convert(V,list);
 i := 1;
 while(i <= nops(L) and L[i]=0) do i := i+1 od;
 nops(L)-i;
end:

lexo := proc(l1,l2)
 if (l1[1]<l2[1])    then true 
 elif (l1[1]> l2[1]) then false else
     evalb(l1[2]<l2[2]) 
 fi
end:

nlexo := proc(l1,l2) evalb(not(lexo(l1,l2))) end:

## ----------------------------------------------------------------------
## HELP triang
##
## triang - Pseudo-triangulation of the matrix by rows and columns 
##          permutations.
##
## CALLING SEQUENCE:
##     triang(m);
##
## PARAMETERS:
##     m  - matrix
##
## DESCRIPTION:
## - Put the matrix in a form which is closed to a triangular form
##   by rows and columns permutations.
##
## - Count the number of non-zero elements per column and sort the columns
##   according to this number. Then for each row, compute the index of the
##   first non-zero coefficient and sort the rows according to this index.
## 
##
## EXAMPLES:
##> readlib(multires);
##> A := matrix(3,3,[0,1,0,0,0,1,1,1,1]);
##                             [0    1    0]
##                             [           ]
##                        A := [0    0    1]
##                             [           ]
##                             [1    1    1]
##> triang(A);
##                             [1    1    1]
##                             [           ]
##                             [0    1    0]
##                             [           ]
##                             [0    0    1]
## SEE ALSO:
## gausselim, ffgausselim 
##
triang := proc(M)
local Lc,Lr, i,N:
 Lc := NULL;
 for i to coldim(M) do
   Lc := Lc, [nzt(col(M,i)),i];
 od;
 Lc:= map(u->op(2,u),sort([Lc],lexo));
 N := submatrix(M,1..rowdim(M),Lc);

 Lr := NULL;
 for i to rowdim(N) do
   Lr := Lr, [idx(row(N,i)),i];
 od;
 Lr:= map(u->op(2,u),sort([Lr],nlexo));
 submatrix(N,Lr,1..coldim(N));
end:
## ----------------------------------------------------------------------
upper := proc(N)
local Lr,i;
 Lr := NULL;
 for i to rowdim(N) do
   Lr := Lr, [idx(row(N,i)),i];
 od;
 Lr:= map(u->op(2,u),sort([Lr],nlexo));
 submatrix(N,Lr,1..coldim(N));
end:

## ----------------------------------------------------------------------
## HELP shape
##
## shape - Output the shape of the matrix.
##
## CALLING SEQUENCE:
##     shape(m);
##
## PARAMETERS:
##     m  - matrix
##
## DESCRIPTION:
## - Output a x for a non-zero coefficient and white space for zero 
##   coefficients 
##
## EXAMPLES:
##> readlib(multires);
##> A := matrix(3,3,[0,1,0,0,0,1,1,1,1]);
##                             [0    1    0]
##                             [           ]
##                        A := [0    0    1]
##                             [           ]
##                             [1    1    1]
##> shape(A);
## 3   3
## [ x ] 1
## [  x] 2
## [xxx] 3
##
shape := proc(M)
local i,j, s;
 lprint(rowdim(M), coldim(M));
 for i to rowdim(M) do
   s := ``;
   for j to coldim(M) do
     if(evalb(M[i,j]=0)) then s:=s,` ` else s := s,`x` fi
   od;
   lprint(cat(`[`,s,`] `,i));
 od;
end:

## ----------------------------------------------------------------------
addtohelp(compagnon):
## ----------------------------------------------------------------------
## HELP compagnon
##
## compagnon - Compute the compagnon matrix of a polynomial matrix
##
## CALLING SEQUENCE:
##     compagnon(M,x);
##
## PARAMETERS:
##     M  - matrix with entries which are polynomial in u
##     u  - variable 
## DESCRIPTION:
## - return a list [A,B] of matrices such that det(A-x* B) is equal 
##   (up to a sign) to det(M).
## 
## - The matrix B may not be inevrtible.
## 
## - Reduce the problem of solving M v=0 to an eigenvector problem.
##
## - The default value for u is 1.
##
## EXAMPLES:
##> readlib(multires);
##> M := matrix(3,3,[a+u^2,1-u,u^2, 2*u^2-b,u-2,u^2-2,u-1,u^2,u^2+c]);
##
##                         [      2                 2  ]
##                         [ a + u      1 - u      u   ]
##                         [                           ]
##                    M := [   2                  2    ]
##                         [2 u  - b    u - 2    u  - 2]
##                         [                           ]
##                         [              2       2    ]
##                         [ u - 1       u       u  + c]
## 
##> A:= compagnon(M,u);
##
##        [0      0    0      1     0    0]  [1    0    0    0    0    0]
##        [                               ]  [                          ]
##        [0      0    0      0     1    0]  [0    1    0    0    0    0]
##        [                               ]  [                          ]
##        [0      0    0      0     0    1]  [0    0    1    0    0    0]
##  A := [[                               ], [                          ]]
##        [-a    -1    0      0     1    0]  [0    0    0    1    0    1]
##        [                               ]  [                          ]
##        [b      2    2      0    -1    0]  [0    0    0    2    0    1]
##        [                               ]  [                          ]
##        [1      0    -c    -1     0    0]  [0    0    0    0    1    1]
##
##> det(A[1]-u*A[2])-det(M);
##                               0
## SEE ALSO:
## 
##
compagnon :=proc(M,x)
local i,j,k,d,n,A,B;
  d := -1;
  n := rowdim(M);
  for i to rowdim(M) do
    for j to coldim(M) do
      d :=  max(d,degree(M[i,j],x));
    od;
  od:
  d;
  A := matrix(n*d,n*d,0);
  for i to n*(d-1) do A[i,n+i] := 1 od;
  B := matrix(n*d,n*d,(i,j) -> if i=j then 1 else 0 fi);

  for i to rowdim(M) do
    for j to coldim(M) do
       for k from 0 to d-1 do
         A[n*(d-1)+i,n*k+j]:= -coeff(M[i,j],x,k);
       od;
       B[n*(d-1)+i,n*(d-1)+j] := coeff(M[i,j],x,d);
    od;
  od:
  [evalm(A),evalm(B)];
end:

## --------------------------------------------------------------------
addtohelp(neigenvects):
## --------------------------------------------------------------------
## HELP neigenvects
##
## neigenvects - Compute the normalized eigenvectors.
##
## CALLING SEQUENCE:
##     neigenvects(m,u);
## PARAMETERS:
##     m  - matrix
##     u  - (optional) index of the which will be 1 in the eigenvector.
## DESCRIPTION:
## - Compute the eigenvectors of the matrix m and divide by the ith 
##   coefficient. It should not be zero.
## - The default value for u is 1.
##
## EXAMPLES:
##> readlib(multires);
##> A := matrix(3,3,[0,1.,0,0,0,1.,1,1.,1]);
##                             [0    1    0]
##                             [           ]
##                        A := [0    0    1]
##                             [           ]
##                             [1    1    1]
##> neigenvects(A);
## 
## [1.839286757, 1, {[1.000000000, 1.839286756, 3.382975769]}], [
##
##                                      [                            -10
##    -.4196433793 + .6062907299 I, 1, {[1.000000000 + .5503714261 10    I,
##
##                                                              ]
##    -.4196433774 + .6062907281 I, -.1914878844 - .5088517781 I]}], [
##
##                                      [                            -10
##    -.4196433793 - .6062907299 I, 1, {[1.000000000 - .5503714261 10    I,
##
##                                                              ]
##    -.4196433774 - .6062907281 I, -.1914878844 + .5088517781 I]}]
##
##
## SEE ALSO:
## eigenvects, eigenvals
##
neigenvects := proc(M)
local u, i, r, ev;
 if nargs >1 then u:=args[2] else u:=1 fi;
 ev := [eigenvects(M)];
 r := NULL;
 for i from 1 to nops(ev) do 
  r := r, [ev[i][1],ev[i][2], {evalm(ev[i][3][1]/ev[i][3][1][u]),
		               ev[i][3][j]$j=2..nops(ev[i][3])}];
 od; 
 r;
end:

## --------------------------------------------------------------------
addtohelp(eigensolve):
## --------------------------------------------------------------------
## HELP eigensolve
##
## eigensolve - Compute subvectors of the eigenvectors.
##
## CALLING SEQUENCE:
##     eigensolve(M,li);
##
## PARAMETERS:
##     M  - matrix
##     li - list of index.
## DESCRIPTION:
## - Compute the subvectors obtained by normalisation by the index li[1]
##   and by returning the entries li[2], ...,li[nops(li)] of the normalized
##   eigenvector.
##
## - If M is the matrix of multiplication by an element modulo an ideal
##   and if li are the index of the monomials 1,x1,...,xn in the basis
##   of the quotient, then  this procedure output the coordinates of the
##  roots of the system (if they are simple).
##
## EXAMPLES:
##> readlib(multires);
##> A := matrix(4,4,[0,1.,0,0,0,1.,1,1.,1,2,1,0,0,0,2,-1]);
##> eigensolve(A,[1,2,3]);
## [-.6307079644, .1603082442], [2.885337042, 3.591227439],
##
##    [-.6273145348 + 1.120269926 I, .7492321665 - .8610490231 I],
##
##    [-.6273145348 - 1.120269926 I, .7492321665 + .8610490231 I]
##
##
## SEE ALSO:
## eigenvects, eigenvals
##		
eigensolve := proc(A,li::list)
local ev, r,i,j,l;
 ev := [neigenvects(A,li[1])];
 r :=  NULL;
 for i from 1 to nops(ev) do
   l := NULL;
   for j from 2 to nops(li) do l := l,ev[i][3][1][li[j]]; od;
   r := r, [l];
 od;
 r;
end:

#----------------------------------------------------------------------#
subgausselim := proc(M,l::list)
local N,lr,i,j;
N := gausselim(submatrix(M,l,1..coldim(M)));
lr := NULL;
j := 1;
for i to rowdim(M) do
   if member(i,{op(l)}) then lr := lr,convert(row(N,j),list) ; j := j+1 
   else
    lr := lr, convert(row(M,i),list);
   fi;
od;
matrix([lr]);
end:

#----------------------------------------------------------------------#
ured := proc(M,n) 
  concat(array(identity, 1..n,1..n),
         evalm(-submatrix(M,1..n,n+1..coldim(M))
	 &*inverse(submatrix(M,n+1..rowdim(M),n+1..coldim(M))))); 
end:

#----------------------------------------------------------------------#
Kernel := proc(M)
local U,M1,M2,N,d,W,i,V;
   if type(evalm(M), matnum) then 
      evalf(Svd(M,U,`left`));
      convert(linalg[col](U,coldim(U)),list)
   else
      M1 := map(Re,evalm(M));
      M2 := map(Im,evalm(M));
      N := stack(concat(M1,-M2), concat(M2,M1));
      evalf(Svd(N,U,`left`));
      V := linalg[col](U,coldim(U));
      d := iquo(vectdim(V),2);
      W := NULL;
      for i to d do 
      	 W := W, V[i] + I* V[d+i]
      od;
      [W];
   fi;
end: 

## ----------------------------------------------------------------------
addtohelp(hrank):
## ----------------------------------------------------------------------
_random := rand(-100..100):
## ----------------------------------------------------------------------
## HELP hrank
##
## hrand - Compute the << generic >> rank of a matrix depending on parameters
##
## CALLING SEQUENCE:
##     hrank(M);
##
## PARAMETERS:
##     M  - matrix
## DESCRIPTION:
## - Compute the << generic >> rank of a matrix, whose coefficients depends
##   on some parameters. These parameters are substituted by random  numbers
##   using the function _random() (default: rand(-100..100)/100)
##   and the rank is computed with the function rank.
##
## - The matrix M should contains exact coefficients.
##
## EXAMPLES:
##> readlib(multires);
##> A:=matrix(4,4,
##	[0,1+u,u^2,0,0,1,1,u*v-1,a^2-1,2*a^2-2,a^2*b-b,0,w,2*w,w*b,0]);
##> hrank(A);
##
##                                  3
## 
## SEE ALSO:
## rand, rank
##
hrank := proc(M)
local x, Ms;
  x := indets(convert(evalm(M),listlist));
  lprint(coldim(M),rowdim(M));
  Ms := subs(map(u->(u= _random()/100),x),evalm(M));
  rank(Ms);
end:

#----------------------------------------------------------------------#
hrankdcmp := proc(M)
local h, x, Ms;
  h := rand(-100..100);
  x := indets(convert(evalm(M),listlist));
  Ms := subs(map(u->(u= h()/100),x),evalm(M));
  LUdecomp(Ms,L='l',U='u',U1='u1',R='r',P='p',det='d',rank='ran');
  submatrix(evalm(inverse(p)&* M), 1..ran, 1..coldim(M));
end:

#----------------------------------------------------------------------#
addtohelp(hmaxminor):
## ----------------------------------------------------------------------
## HELP hmaxminor
##
## hrand - Compute a maximal minor of a matrix depending on parameters
##
## CALLING SEQUENCE:
##     hmaxminor(M);
##
## PARAMETERS:
##     M  - matrix
## DESCRIPTION:
## - Compute a maximal minor of a matrix, whose coefficients depends
##   on some parameters. These parameters are substituted by random  numbers
##   using the function _random() (default: rand(-100..100)/100)).
##
## - The matrix M should contains exact coefficients.
##
## EXAMPLES:
##> readlib(multires);
##> A:=matrix(4,4,
##	[0,1+u,u^2,0,0,1,1,u*v-1,a^2-1,2*a^2-2,a^2*b-b,0,w,2*w,w*b,0]);
##> hmaxminor(A);
##                      [ 2           2         2      ]
##                      [a  - 1    2 a  - 2    a  b - b]
##                      [                              ]
##                      [  0          1           1    ]
##                      [                              ]
##                      [                          2   ]
##                      [  0        1 + u         u    ]
## 
## SEE ALSO:
## rand, rank
##
hmaxminor := proc(M)
local T, d, N;
  N := hrankdcmp(transpose(M));
  triang(hrankdcmp(transpose(N)));
end:

## ----------------------------------------------------------------------
addtohelp(berkowitz):
## ----------------------------------------------------------------------
## HELP berkowitz
##
## berkowitz  - determinant by Berkowitz formula
##
## CALLING SEQUENCE:
##     berkowitz(M);
##
## PARAMETERS:
##     M - a matrix
##    
## DESCRIPTION:
## - Compute the determinant of M by Berkowitz formula.
## 
## - The determinant is not expanded.
## 
## EXAMPLES:
##> readlib(multires);
##> A := matrix(4,4,[0,u,1,0,1,v,2,1,1,2,1,-1,0,1,-1,u]);
##> expand(berkowitz(A));
## 
##                                               2
##                              2 - u v + 4 u + u
##
## SEE ALSO:
##  det
##
berkowitz := proc(M)
local U,c,i;
    U := evalm(M);
    for i from 2 to coldim(M) do
        U := evalm((trace(U)*M - (i - 1)*U &* M )/i)
    od;
    trace(U);
end:

## ----------------------------------------------------------------------
addtohelp(signature):
## ----------------------------------------------------------------------
## HELP signature
##
## signature  - Signature of quadratic form associated with a symmetric 
##              numeric matrix Signature of a quadratic form 
##
## CALLING SEQUENCE:
##     signature(Q);
##
## PARAMETERS:
##     Q  - matrix
##    
## DESCRIPTION:
## - Compute the signature [p,q] of the quadratic form defined by the 
##   matrix Q, where p is the number of positive components and q the
##   number of negative terms. The rank of the matrix is p+q.
##
## - Q should be a real symmetric matrix.
##
## - Count the sign of the eigenvalues (computed by the function eigenvals). 
##
## EXAMPLES:
##> readlib(multires);
##> Q := matrix(4,4,[0,1,1,0,1,1,2,1,1,2,1,-1,0,1,-1,-2]);
##> signature(Q);
## 
##                         [1, 2]
##
## SEE ALSO:
## eigenvals
##
signature := proc(Q)
local lv,p,q,i;
  lv := [evalf(eigenvals(Q))];
  p :=0; q := 0;
  for i to nops(lv) do
    if lv[i] >0 then  p := p+1 
      elif lv[i] < 0 then q := q+1;
    fi;
  od;
  [p,q];
end:

## ----------------------------------------------------------------------
addtohelp(melim):
## ----------------------------------------------------------------------
## HELP melim
##
## melim - Pseudo-melimulation of the matrix by rows and columns 
##          permutations.
##
## CALLING SEQUENCE:
##     melim(l,v);
##
## PARAMETERS:
##     l  -  list of polynomials
##     v  - list of variables that you want to eliminate.
##
## DESCRIPTION:
## - Eliminate the variables v amoung the polynomials in l.
##
## - The number of polynomials should be one more than the number 
##   of variables.
## 
## - Compute a maximal non-zero minor (Bareiss method, implemented in 
##   ffgausselim) of the corresponding Bezoutian matrix.
##
## EXAMPLES:
##> readlib(multires);
##> den:= r^3+s*r^2;
##> melim([r*(s^2-r^2+1)-x*den,(s^2+r*s+1)-y*den,(s^3+r*s^2)-z*den],[r,s]):
##> factor(%);
##  2                2         3       2    2      2  2       4        5    6
## y  (-1 - 6 x - 4 z  x - 20 x  - 15 x  - z  - 6 x  z  + 16 x  z - 6 x  - x
## 
##           2         4      2    2             3        2          2      4
##      - 2 y  z - 15 x  - 2 y  z x  + 2 z + 24 x  z - 4 y  z x + 6 y  x - y
## 
##           2  3      2  2      2             6      4  2         2      5
##      + 2 y  x  + 6 y  x  + 2 y  + 10 z x + x  z - x  z  + 21 z x  + 6 x  z
## 
##           3  2
##      - 4 x  z )
##> melim([u[0]+u[1]*a+u[2]*b,a^2-1,b^2-4],[a,b]):
##> factor(%);
##     -(u[0] + u[1] + 2 u[2]) (u[0] - u[1] - 2 u[2]) 
##
##                  (u[0] + u[1] - 2 u[2]) (u[0] -u[1] + 2 u[2])
##
## SEE ALSO:
## mbezout, ffgausselim 
##
melim := proc(lp::list, var::list)
local M, f,lmm,r,i,T;   
   if nargs<=2 then
       M := mbezout(lp,var);
   else
       f := args[3];
       M := f(lp,var);
   fi;	
   T := linalg[ffgausselim](M):
   lmm := [lasts(T)];
   if nops(lmm)>0 then 
      r := lmm[1]; 
   else
      ERROR(`Null Matrix `);
   fi;
   for i from 2 to nops(lmm) do 
     r := gcd(r,lmm[i]);
   od;
   RETURN(r);
end:

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# GEOMETRIC DECOMPOSITION OF A VARIETY                                 #
# This part has been written by S. Tonelli (DEA '98, Univ. de Nice).   #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

readlib(factors) : 

mineur:=proc(lp::list,var::list)
  local b,i,n,m;
  n :=nops(var);
  b :=melim(lp,var);
  m :=collect(b,[seq(u[i],i=1...n)],distributed);
  RETURN(m);
end : 


rep:=proc(lp,mn,var)
  local i,j,k,h,l1,l2,l3,l4,l5,l,r,c,n,m,d,s,f,p,g,lg,lr;

  n:=nops(var); 
  m:=nops(lp);
  h:=rand(-10...10)/10;

  d:=simplify(mn/gcd(mn,diff(mn,u[0])));
  
lprint(`Squarefree part of the maximal minor:`);
print(factor(d));

 # Calcul des di
   l1:=d;
   for i from 1 to n do 
    l1:=subs(u[i]=h()+t*u[i],l1);
   od;
   l2:=subs(u[0]=u[0]+h(),l1);
   l3:=collect(l2,t,distributed);
   l4:=coeff(l3,t);
   l5:=collect(l4,[seq(u[i],i=1...n)]);
   d[0]:=subs(t=0,l3);

   for i from 1 to n do
    d[i]:=coeff(l5, u[i]);
lprint(cat(`Numerator d`,i,`:`)); print(factor(d[i]));
   od;

 for i from 1 to n do
  s[i]:=d[i]/diff(d[0],u[0]); 
 od;

 # Factorisation de d(u0) et recuperation des facteurs irreductibles. 
  f:=factors(d[0]); lprint(`Factorisation of d0: `); print(f);
  for i from 1 to nops(op(2,f)) do
   p[i]:=op(1,op(i,op(2,f)));
  od;

 #Calcul des gi(u0)
  for i from 1 to m do
   l:=lp[i];
   for j from 1 to n do
    l:=subs(var[j]=s[j],l);
   od;
   g[i]:=simplify(l);
  od;

 #Liste des numerateurs des gi. 
  lg:=NULL;
  for i from 1 to m do
   lg:=lg,numer(g[i]);
  od;

 # On garde les pi qui divisent tous les numerateurs. 
  lr:=NULL;
  for i from 1 to nops(op(2,  f)) do
   c:=0;
   for j from 1 to m do
    if rem(lg[j],p[i],u[0])=0 then
     c:=c+1;
    fi;
   od;
   if c=m then
    r :=p[i];
    for k from 1 to n do
     r:=r, simplify(rem(d[k],p[i],u[0])/rem(diff(d[0],u[0]),p[i],u[0]));
    od;
    lr:=lr,[r]; 
   fi;
  od;

lprint(cat(`Rational minimal representation of the components of dimension `,m-n));
RETURN(lr);
end : 

## ----------------------------------------------------------------------
addtohelp(decomp):
## ----------------------------------------------------------------------
## HELP decomp
##
## decomp - Decompose a variety into irreducible components.
##
## CALLING SEQUENCE:
##     decomp(lp, var)
##
## PARAMETERS:
##     lp  - list of polynomials
##     var - list of variables
##
## DESCRIPTION:
## - Compute a decomposition of the variety defined by lp, into
##   irreducible components, described by a rational representation.
##
## EXAMPLES:
## SEE ALSO:
## 
decomp:=proc(lp,var)
 
  local i,j,n,lv,m,f,h,r,delta,s,d;
  n := nops(var);
  if nargs>2 then d := n-args[3] else   d := 0 fi;
  lv:=var;
  m:=nops(lp); 
 
  for i from 1 to m do
   f[i]:=lp[i];
  od; 

  h:=rand(-10...10)/10;

  while evalb(n<>0 and n>=d) do
   if m>n then 
    for i from 1 to n do
     f[i]:=sum('h()*f[j]', 'j'=1...m);
    od; 
    m:=n;
    elif m<n then 
    lv:=[seq(var[i],i=1...m)];
    n:=nops(lv);
   fi; 
   delta:=mineur([u[0]+sum('u[i]*lv[i]', 'i'=1...n),seq(f[i],i=1.. n)],lv);
   r:=rep(lp,delta,lv); print(r);
   lv:=[seq(lv[i],i=2...n)];
   n:=nops(lv); 
   r;
  od;
end : 

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# ELIMINATION PROCEDURES                                               #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

mesolve := proc(lp::list)
local var, mthd, M, p, lr, x, s, lm, ln, V, W, k;
  if nargs =1 then
     var := sort([op(indets(lp))]);
     var := [var[1..nops(var)-1]];
     mthd := mbezout;
  elif nargs=2 then 
     var := args[2];
     mthd := mbezout;
  else 
     var := args[2];
     mthd := args[3];
  fi;
   M  := mthd(lp, var, 'lm', 'ln');
   k := 1;
   member(1,lm,'k');
   p  := melim(lp,var, mthd);
   x  := op(indets(p));
   lm := [op(lm),x];
   lr := [fsolve(p,x,complex)];
   V := NULL;
   for s in lr do
      W := Kernel(subs(x=s,evalm(M)));
      W := convert(evalm(W/W[k]),list);
      V := V, [op(W),s];
   od;
   lm,V;
end:  

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# DUALITY                                                              #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# Product of a linear form and a polynomial.

`&.` := proc(a,b)
local v,w,ai,i;
     v := [op(indets(a))];
     if(op(0,v[1]) = 'd') then 
          w := map(u_->1/u_,subs(d=x,v));
     else
          w := map(u_->1/u_,subs(x=d,v));
     fi;	  
     ai := a;
     for i to nops(v) do
        ai := subs(v[i]=w[i],ai);
     od;
     select(u_->evalb(type(u_,monomial)),expand(ai*b))
end:

dualbasis := proc(lp)
local l,i,n,R,var,nu;
  if nargs =1 then
     var := sort([op(indets(lp))]);
  else 
     var := args[2];
  fi;
  R := mresultant(lp,var,l);
  n := 1; nu := degree(lp[1])+1;
  for i from 2 to nops(lp) do 
     n :=n *degree(lp[i],var);
     nu := nu+ degree(lp[i],var)-1;
  od;
  lprint(l[1..n]);
  evalm(ured(R,n) &* subs(x=d,l) +O(d^(nu)));
end:  

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#  TORIC RESULTANTS
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# Maple8 code for sparse (toric) resultant matrices.
# Author: Ioannis Z Emiris, Univ. of Athens, Greece
# emiris@di.uoa.gr, http://www.di.uoa.gr/~emiris
# based on code by John Canny and Paul Pedersen in 1993.
# This is software under development, and is free for distribution.
# The author assumes no responsibility, but comments are welcome.

# Created: 5/96
# Last update: 06/2003

# Main function call:
# spresultant(polynomial_list, eliminated_variable_list,
#	[,geometric_perturbation_vector [,lifting_matrix]] ) 
#	RETURNS sparse (toric) resultant matrix in input coefficients

# References:
# For the subdivision-based algorithm: J.F. Canny and I.Z. Emiris,
# "An efficient algorithm for the sparse mixed resultant", in Proc.
# AAECC-93, LNCS 263, pp. 89--104.  Final version in J. ACM, 2000.
 
with(simplex):
with(linalg):
readlib(unassign):

############# utilities #####################

RandSign := proc() eval(eval(rand(0..1)())*2-1) end:

# Return true/false depending on whether 'tarray' is an array.
#   For 1-dimensional arrays (vectors) check whether its elements
#   are of type 'ttype' and return true/false accordingly.
#
`toric/type1Array` := proc(tarray::array,ttype::name)
  local i, size;
  if not type(tarray,vector) then RETURN(true); fi;	
  # else 1-dim array ie. vector and will check elements
  size := nops(convert(tarray,list));
  for i from 1 to size do
    if (not(evalb(type(tarray[i],ttype)))) then RETURN(false); fi;
  od;
  RETURN(true);
end:	# `toric/type1Array`

# computes the 2-norm of a vector expressed as (1Xn) matrix
#
`toric/twonorm` := proc(x::matrix)
  local i,r;
    r := x[1,1]^2;
    for i from 2 to coldim(x) do r := r + x[1,i]^2 od;
    convert(sqrt(r), float)
end: # `toric/twonorm`

# map a function 'fnc' on all coefficients of polynomial 'poly'
# return the new polynomial
#
`toric/map_coefs` := proc (fnc,poly::polynom)
  local cvec, mvec, mons;
  cvec := vector(map(fnc,[coeffs(expand(poly),indets(expand(poly)),'mons')]));
  mvec := vector([mons]);
  # lprint(`map_coef_returns`,dotprod (cvec, mvec, 'orthogonal'));
  dotprod (cvec, mvec, 'orthogonal');
end: 	# `toric/map_coefs`

#################### format conversion #############################
#----------------------------------------------------------------------
addtohelp(`toric/polytope`):
#---------------------------------------------------------------------- 
## HELP `toric/polytope`
##
## `toric/polytope` - compute supports and Newton polytopes of polynomial system
##
## CALLING SEQUENCE:
##     `toric/polytope`(pol, var)
##
## PARAMETERS:
##     pol  - list or vector of polynomials
##     var  - list or vector of variables
## DESCRIPTION:
## - Compute the supports and Newton polytopes of polynomial system
##   pol as function of variables var.
##
## - The function returns a sequence of two arrays for the supports
##   and the Newton polytopes
##
## EXAMPLES:
##> `toric/polytope`([x^2-x+y*x-1,x*y-3,2*x-x*y-1],[x,y]);
##
##   [[0   1   2   1]  [0   1]  [0   1   1]]
##   [[             ], [     ], [         ]], [{1, 3, 4}, {1, 2}, {1, 2, 3}]
##   [[0   0   0   1]  [0   1]  [0   0   1]]
##
## SEE ALSO:
## `toric/mixed`, spresultant
##
`toric/polytope` := proc (pol::list, var::list)
local polylist, varlist,
      apol, npols,
      supps, hulls;
  polylist := convert(pol,list);
  varlist := convert(var,list);

  npols := nops(polylist);
  # printf("working on %d polynomials in %d variables\n",npols,nops(varlist));

  supps := array (1..npols);
  hulls := array (1..npols);

  # printf("Newton polytope vertex cardinalities are");
  for apol from 1 to npols do
    supps[apol] := `toric/supp_mat` (polylist[apol], varlist);		# mat[dim,num]
    hulls[apol] := `toric/comp_hull` (supps[apol], 2*nops(varlist));	# 2 = arbitrary
    # printf(", %d",nops(hulls[apol]));
  od;
  # printf(".\n");

  RETURN (eval(supps), eval(hulls));
end:	# `toric/polytope`

# in_poly = input polynomial
# vars = list of variables to be eliminated (after hiding)
# return the support of polyn as num_vars x num_monoms matrix
#
`toric/supp_mat` := proc (in_poly::polynom, vars::list)
  local v, nvars,
	poly, supp,
	mons, monl,
	m, nmons;
  # lprint(`got polyn`,in_poly,`in vars`,vars);
  poly := collect(simplify (expand(in_poly)), vars, distributed);
  coeffs( eval(poly), vars, 'mons');
  # lprint(`put monoms in`,mons);

  monl := sort([mons]);
  nmons := nops(monl);
  # lprint(`collected polyn`,poly,`w/monoms`,monl);

  nvars := nops (vars);
  supp := matrix (nvars, nmons);

  for m from 1 to nmons do
    for v from 1 to nvars do supp[v,m] := degree (monl[m], vars[v]); od;
  od;	# for monomial

  evalm (supp);
end:	# `toric/supp_mat`

# polyn = polynomial in any variables
# return polynomial multiplied by min. int to get rid of denominators
#   so that max abs.value of new coeffs doesnt exceed local bound
#
`toric/scaled` := proc (polyn::polynom)
  local coefs, c,			# all Numeric coefficients
        mx, mult,			# max coef, denominator multiplier
        too_big_constant;
  too_big_constant := 2^30;

  coefs := coeffs (collect(expand(polyn),indets(polyn),distributed));
  mx := max(op(map(abs,[coefs])));

  mult := 1;
  for c in coefs do if not type(c,integer) then
    c := convert(c,rational,exact);
    mult := ilcm (mult, denom(c));
  fi; od;
  if evalf(abs(mult*mx)) >= too_big_constant then
    mult := max(1, floor( too_big_constant/mx ));
  else # lprint(`multiplying by _lcm_ of denoms=`,mult);
  fi;
  eval(mult * polyn);
end:  # `toric/scaled`

################## matrix construction #####################

# checks if "el" lies in the set or list "vset", using equal rather than =
#
`toric/vmem` := proc(el, vset)
  local found, el2;
  found := false;
  for el2 in vset do
    if equal(el, el2) then found := true; break fi;
  od;
  found;
end: # `toric/vmem`

# define some random affine functionals.  they shall determine
# a mixed subdivision by lifting then taking the lower envelope
#
`toric/randaff` := proc(d::integer)
  local i,j,l,die;
  die := rand(1..1000);
  l := matrix(d,d);
  for i from 1 to d do
    for j from 1 to d do l[i,j] := die() od;
  od;
  eval(l);
end: # `toric/randaff`

# random displacement (perturbation) vector of length n
#
`toric/randvect` := proc(n::integer)
  local i, rproc, sgn, r;
  rproc := rand(1..1000);
  sgn := rand(1..2);
  r := array(1..n);

  for i to n do
    r[i] := convert(rproc()/10000, float);
    if sgn()=1 then r[i] := -r[i]; fi
  od;
  eval(r);
end: # `toric/randvect`

# A = array of matrices, each contains vertices as columns.
# alift = random affine functional expressed as matrix.
# d = length of A = dimension - 1.
#
# Constructs the following input equations for LP package:
#   objective             = \sum_i l_i(\sum_j \lambda_{i,j} A_{i,j}) + l_i[const],
#   vertex equations lhs  = \sum_{ij} \lam_{ij} A_{ij},
#   convexity constraints = \sum \lambda_{i,j} = 1;
# where
#   l_i = alift_i = random affine functional, l_i[const] = 0 iff lin.functional
#   \lambda_{i,j} = LP variables,
#
`toric/make_eqns` := proc(A::array, alift::matrix, d::integer)
  local i, j,
      C, cd, conveq, convsum, ipoly, lambda, nverts, objective, vertlhs,
      linpart;

  objective := 0;				# objective function
  vertlhs := vector(d-1, 0);		# left hand sides vertex equations
  conveq := {};				# convexity constraints
  nverts := dotprod(vector(d,1), map(coldim, A));# total number of vertices
  linpart := submatrix (alift, 1..d, 1..(d-1));

  for i from 1 to d do
        ipoly := A[i];
        convsum := 0;
        cd := coldim(ipoly);
        lambda := array(1..cd);			# array for new variables

        for j from 1 to cd do
            lambda[j] := `lam`||i||`x`||j;
            convsum := convsum + lambda[j];
        od;
 
        conveq := conveq union {convsum = 1};
        C := multiply(ipoly, lambda);
        vertlhs := matadd(vertlhs, C);
 
        objective := objective + multiply (subvector(linpart,i,1..(d-1)),C); 
    	if coldim(alift)=d then objective := objective + alift[i,d] fi;
  od;

  # print([objective, conveq, vertlhs, nverts]);
  [objective, conveq, vertlhs, nverts];
end: # `toric/make_eqns`

# computes Center Of Gravity of input polynomial
# hull=arr[d,#monoms.in.hull], exponents are columns
# returns matrix[#monoms,1]
#
`toric/cog` := proc(hull::array)
  local numonom, v;
  numonom := coldim(hull);
  v := multiply(hull, convert(vector(numonom,1), matrix));
  scalarmul(v, 1/numonom);
end: # `toric/cog`

# oldvec = vector w/entries in {-1,0,1}
# expresses assignment, will return next higher assignment
# 
`toric/incr3asg` := proc (oldvec::vector)
  local i, k, dim, lis;
  lis := convert(oldvec,list);
  dim := nops(lis);
  k := dim;
  while lis[k] = 1 do k:=k-1 od;
  lis[k] := lis[k]+1;
  for i from k+1 to dim do lis[i] := -1 od;
  eval(convert(lis,vector));
end:	# `toric/incr3asg`

# hullArr expresses Newton polytopes
# dim = dimension of lifted space = number of polytopes
# edata = list of objects for linear optimizationS
# delta = vector of geometric perturbation
# return Integer point guaranteed to lie inside perturbed Mink.sum
#   ie. it is vertex in perturbed lifted M.sum.  The last optimization
#   returns info that can be used later but for now it is wasted.
#
`toric/start_vertex` := proc (hullArr::array, dim::integer, edata::list, delta::vector)
  local i,
	asg, count_asg,
	imax, overts,
	init_v, 			# sum of cog's
	ivt;				# rounded point

    init_v := `toric/cog`(hullArr[1]);		# start at the center of gravity..
    for i from 2 to dim do		# ..of Minkowski sum
      init_v := matadd(init_v,`toric/cog`(hullArr[i]))  # =mat[rowdim,1]=[(d-1),1]
    od;
    init_v := map(round,init_v);

    ivt := eval(init_v);			# (dim-1)x1 matrix
    imax := 0;
    asg := vector(dim-1,-1);
    count_asg := 1;

    while imax = 0 do
      # print(`Start point candidate`,eval(ivt));
      overts := `toric/opt_verts`(edata, matadd(ivt,delta,1,-1), dim);	# lower hull vertex
      imax := dim;
      while imax>0 and overts[imax]=0 do imax:=imax-1 od;
      if imax=0 then
	ivt := matadd(init_v,asg,1,1);
	count_asg := count_asg + 1;				# asg to be comp'd
	if count_asg > 3^(dim-1) then ERROR(`cant find valid starting point`) fi;
	asg := `toric/incr3asg` (asg);
      fi;
    od; # while
    eval(ivt);
end:	# `toric/start_vertex`

# must turn `lamMxN` into list [M, N]
#
`toric/recode` := proc(z0)
  local i, L, M, N, s, state, z,
	char0, char1, char2, char3, char4,
	char5, char6, char7, char8, char9;
  for i from 0 to 9 do `char`||`i` := i; od;
  z := convert(z0,string);
  L := length(z);
  M := 0;
  N := 0;
  state := 1;
  for i from 4 to L do
    s := substring(z,i..i);
    if s = "x" then state := 2 fi;
    if state = 1 then M := 10*M + `char`||s;		# for decimals
    elif state = 2 then state := 3;
    elif state = 3 then N := 10*N + `char`||s; fi;
  od;
  [M, N];
end:		# `toric/recode`

# Check which polytope points contribute to optimum sum for vert.
# edata = list of objects including equations from `toric/make_eqns`(),
# vert  = vertex in Minkowski sum,
# d     = dimension + 1 = number of polynomials.
# return optverts = d-dim vector, where d=n+1
#   optverts[i] = summand vertex of support i, 0 if no vertex from polytope i
#
`toric/opt_verts` := proc(edata::list, vert, d::integer)
  local	assigns, eqns, ivert, j, nverts, objective,
	optverts, optvtx,
	verteq, vertlhs, a,
 	epsilon,
	start_timer;
  global tot_time_lp;
     epsilon := 10^(-round(2*Digits/3));

   # Pick apart input from `toric/make_eqns`()
     objective := op(1, edata);		# objective function
     eqns := op(2, edata);		# convexity constraints
     vertlhs := op(3, edata);		# vertex eqns left hand sides
     nverts := op(4, edata);		# total number of input vertices

   # Run LP package to find optimal sum (on lower envelope) for vert
    for j to d-1 do
        eqns := eqns union {vertlhs[j] = vert[j,1]}
    od;
    # printf("Calling minimize() with: d=%d, Opt_verts=%a\n",d,eval(vert));
    # printf("   obj.fnc=%a, eqns=%a\n",eval(objective),eval(eqns));

    start_timer := time();
    assigns := minimize(objective, eqns, NONNEGATIVE);
    tot_time_lp := tot_time_lp + time() - start_timer;

    if assigns=NULL then ERROR(`unbounded`);
    # elif assigns={} then ERROR(`infeasible`); 
    fi;

   # find vertices appearing in the representation
    optverts := vector(d,0);			# unset

    for a in assigns do
        if (abs(rhs(a)-1.0) < epsilon) then
	    optvtx := `toric/recode`(lhs(a));		# list [pol,mon]
	    if optverts[optvtx[1]]>0 then ERROR(`optimal vertex already assigned`) fi;
            optverts[optvtx[1]] := optvtx[2];
        fi
    od;   
    optverts;
end: 			# `toric/opt_verts`

# greedy algorithm for finding rows of sparse (toric) resultant matrix
#  H 	 = c. hulls: arrays of arrays ((d-1) X num vertices) of input vertices
#  S     = supports: arrays of arrays of input vertices plus interior points
#  edata = equations defining lower hull suitable for simplex package
#  v     = vertex in the Minkowski sum of the inputs H[i], coldim x 1 matrix
#  delta = random small displacement vector
#  d     = dimension + 1 = lifted dimension = Number of polynomials
# tracing= debug/message level during execution.  0 for silence
#
# returns sorted list "guide" indexing the matrix rows as polyn*monom products
# each row coded as list [ pt in Mink.sum, polynomial indx, monom in its support ]
# so that the support pt (or monom) w/max polyn.index appears in the optimal sum
# expressing the Mink.sum point 
#
`toric/compute_rows` := proc(H::array, S::array, edata::list, init_v,\
	delta::vector, d::integer, tracing::integer)
  local cols, guide,
 	i, imax, ipoly, jmax, overts, todo, vt1, vt2, counter,
	vt;
    counter := 0;
    guide := {};			# init set of rows
    todo := {eval(init_v)};
    cols := {eval(init_v)};
    # printf("Constructing matrix rows:");

    while not (todo = {}) do		# i.e. continue while vt is defined

      vt := op(1, todo);
      todo := todo minus {eval(vt)};
	if tracing>1 then
	    counter := counter+1;
	    print(`compute rows: Row monomial`, counter, `equals`, vt);
	fi;
   	# if modp(counter,10)=0 then printf("%d ",counter); fi;

      overts := `toric/opt_verts`(edata, matadd(vt,delta,1,-1), d);	# lower hull vertex

 	# find largest index of input polytope contributing a vertex to optimal sum
      imax := d;
	while overts[imax]=0 do imax:=imax-1 od;
      # for vt1 in overts do if
	# vt1[1] > imax then imax := vt1[1]; jmax := vt1[2] fi od;
        jmax := overts[imax];

        guide := guide union {[eval(vt),imax,jmax]};

	# lprint(`values of imax, vt, overts equal`,imax,v,overts);

        vt2 := matadd(vt, col(H[imax], jmax), 1, -1);
        ipoly := S[imax];
        for i from 1 to coldim(ipoly) do	# find exponents occurring as..
            vt1 := matadd(vt2, col(ipoly, i));	# ..multiples of monomials in S[imax]
            if not `toric/vmem`(vt1, cols) then	
                cols := cols union {eval(vt1)};
                todo := todo union {eval(vt1)};
	    fi
        od;
    od;	# while
    # printf("\n");
    sort(convert(guide, list), `toric/guide_order`);
end: # `toric/compute_rows`

# same as `toric/compute_rows` for dimension = #vars = 2 so d = #pols = 3
# row_guide, todo, allcols: include column monoms ie. in Mink.sum
# invariance: allcols = row_guide union todo.  row_guide, todo distinct
# add to todo if allcols unassigned, then assign
# hence can use lists and check allcols for adding, instead of set union
#
`toric/comp_rows_2` := proc(H::array, S::array, edata::list, init_v,\
	delta::vector, tracing::integer)
local	allcols,			# array of column monoms
	row_guide, 			# matrix rows
	maxcoor,			# max coor in Mink.sum of d polytopes
	i, coor, 			# indices
	d, vt,				# #pols=#vars+1: for generalization
	imax, ipoly, jmax, overts,
	todo,				# (column) monoms to be examined
	vt1, vt2, monoms, counter;

    d := 3;
    counter := 0;			# only if tracing>1

    row_guide := [];			# init rows
    todo := { eval(init_v) };

    maxcoor := vector(d-1);
    for coor from 1 to d-1 do		# init
	maxcoor[coor] := max(op(convert(row(H[1],coor),list)));
    od;
    for i from 2 to d do for coor from 1 to d-1 do
	maxcoor[coor] := maxcoor[coor] + max(op(convert(row(H[i],coor),list)));
    od; od;
    if tracing>0 then lprint(`  maximum coordinates`,eval(maxcoor)); fi;
	
    allcols := array(0..maxcoor[1],0..maxcoor[2]); # not assigned: cheaper test
    allcols[init_v[1,1],init_v[2,1]] := 1;		# assigned
    # printf("Constructing matrix rows:");

    while todo <> {} do 			# loop while vt defined
        vt := op(1, todo);
        todo := todo minus {eval(vt)};

	counter := counter+1;
        # if modp(counter,10)=0 then printf("%d ",counter); fi;
        overts := `toric/opt_verts`(edata, matadd(vt,delta,1,-1), 3);	# lower hull vertex
	if tracing>1 then 
	  print(`toric/comp_rows_2: Row monomial`,counter,`equals`,vt,`expressed`,overts);
	fi;

 	# find max index `imax` of polytope contributing vertex to optimal sum
        imax := 3;
	while overts[imax]=0 do imax:=imax-1; if imax=0 then ERROR(`0-imax`) fi od;
        jmax := overts[imax];
        row_guide := [ op(row_guide), [eval(vt),imax,jmax] ];

	if tracing>1 then
	  lprint(`values of imax,vt,overts equal`,imax,eval(vt),eval(overts));
        fi;

        vt2 := matadd(vt, col(H[imax], jmax), 1, -1);
        ipoly := S[imax];
	monoms := coldim(ipoly);

        for i from 1 to monoms do		# find exponents occurring as..
            vt1 := matadd(vt2, col(ipoly, i));	# ..multiples of monomials in S[imax]
            if not assigned( allcols[vt1[1,1],vt1[2,1]] ) then	
		allcols[vt1[1,1],vt1[2,1]] := 1;
                todo := todo union {eval(vt1)};
	    fi;
        od;
    od;
    
    RETURN( sort (row_guide, `toric/guide_order`));

end: # `toric/comp_rows_2`

# criterion for sorting "guide" list of matrix rows
#
`toric/guide_order` := proc(a, b) `toric/v_order`(op(1, a), op(1, b)) end:

# true iff a<b lexicographically
#
`toric/v_order` := proc(a::matrix, b::matrix)
  local bigger, d, i;
    d := rowdim(a);
    bigger := false;
    for i to d do
       if a[i,1] < b[i,1] then
           bigger := true;
           break;
       elif a[i,1] > b[i,1] then
           bigger := false;
           break;
       fi;
    od;
    bigger;
end:		# `toric/v_order`

# suppMat[d-1,#monoms] = matrix of support points
# hullSet = set of indices in suppMat defining c.hull
# sort columns of suppMat on `toric/v_order`, affecting suppMat
# return hull2supp[h]=s indexing: hull-col h = supp-col s
#   submatrix of suppMat comprised of hull vertices
#   must be in same order as supportMatrices to detect diag.monoms
#
`toric/sort_vecs` := proc(hullSet::set,suppMat::matrix)
  local rows, cols, i, j, k,
	d, n, noflips, temp,
	isVertex,
	hull2supp;
  # print('input_hullset,suppmat',hullSet,suppMat);
  d := rowdim(suppMat);
  n := coldim(suppMat);
  isVertex := vector(n,0);			# init
  for i in hullSet do isVertex[i] := 1 od;

  for i from n-1 by -1 to 1 do
    noflips := true;
    for j from 1 to i do
      if `toric/v_order`(submatrix(suppMat,1..d,[j+1]), submatrix(suppMat,1..d,[j])) then
        noflips := false;
        if isVertex[j] <> isVertex[j+1] then
          temp:=isVertex[j]; isVertex[j]:=isVertex[j+1]; isVertex[j+1]:=temp;
        fi;
        for k from 1 to d do
          temp := suppMat[k, j+1];
          suppMat[k, j+1] := suppMat[k, j];
          suppMat[k, j] := temp;
        od;#k
      fi;
    od;#j
    if noflips then break fi;
  od;
  # lprint(`This is a sorted support`, eval(suppMat));

  hull2supp := vector(nops(hullSet),0);		# init
  j:=1;
  for i from 1 to n do
    if isVertex[i]=1 then
      hull2supp[j] := i;			# j-th hull col = i-th supp col
      j := j+1;
    fi;
  od;
  # print('ret hull2set',hullSet); 
  [eval(hull2supp), evalm(submatrix(suppMat,1..d,convert(hull2supp,list)))];
end:	# `toric/sort_vecs`
# ----------------------------------------------------------------------
addtohelp(`toric/mixed`):
# ----------------------------------------------------------------------
## HELP `toric/mixed`
##
## `toric/mixed` - sparse (toric) resultant matrix of supports
##
## CALLING SEQUENCE:
##     `toric/mixed`(hulls, supps, trace);
##     `toric/mixed`(hulls, supps, trace, delta);
##     `toric/mixed`(hulls, supps, trace, delta, lifting);
##
## PARAMETERS:
##     hulls - array of sets defining Newton polytopes
##     supps - 3-dimensional array exrepssing supports
##     trace - non-negative integer
##     delta - vector or list of real numbers
##     lifting - 2-dimensional array or matrix of lifting
## DESCRIPTION:
## - Compute the sparse (toric) resultant matrix of the system defined by
##   the given supports, containing indeterminate coefficients.
##   Coefficients are labeled with respect to the sorted supports.
##   The number of rows in the coefficients of the first polynomial
##   is optimal.
##   
## - Letting n be the number of variables, n+1 must be the number of
##   polynomials.  Then supps must be an array of n+1 arrays of size
##   n x m, where m = the number of points in that support.
##   
## - Integer trace determines the amount of printed information
##   during execution: 0 for minimal information, 2 for maximal.
## 
## - If delta is given, it is used as the geometric perturbation
##   applied to the Minkowski sum of the Newton polytopes.
##   The n entries must be sufficiently small in order not to perturb
##   points in the sum outside adjacent cells in the sum's subdivision.
##   If delta is not given, it is chosen randomly.
##
## - If lifting is given, it specifies the affine lifting functions
##   of the Newton polytopes, thus defining the mixed subdivision of
##   their Minkowski sum which will specify the matrix.
##   It must be a square matrix of dimension n+1, each row corresponding
##   to a Newton polytope.  If lifting is not given, it is chosen randomly.
##
## - The function returns a sequence of 4 objects.  The first is a list
##   of two matrices: the resultant matrix and the denominator matrix,
##   if RATIONAl_FORMULA=1, otherwise an empty matrix.  The second one
##   is a list of the rows with information about how each was obtained:
##   the corrresponding Minkowski sum point, the chosen polynomial and
##   its monomial in the row content function of the point.  Third is the
##   array of the sorted supports, and last is an array of the monomials
##   appearing on the matrix diagonal.
##
## - The function is based on an implementation by Canny and Pedersen.
##   The algorithm is the greedy version, by Canny and Pedersen 
##   (Tech. Report 1394, C.S. Dept, Cornell University, 1993),
##   of the algorithm by Canny and Emiris (Proc. AAECC-1993, LNCS 263,
##   pp. 89.  Final version J. ACM, 2000).  See these references for a
##   definition of the row content function.
## 
## EXAMPLES:
##> supp_hull := `toric/polytope`([x^2-x+y*x-1,x*y-3,2*x-x*y-1],[x,y]):\
##> `toric/mixed`(supp_hull[2], supp_hull[1], 0);
##
## [[c2x1     0      c2x2     0       0  ]   ]
## [[                                    ]   ]
## [[c3x1    c3x2    c3x3     0       0  ]   ]
## [[                                    ]   ]
## [[c1x1    c1x2    c1x3    c1x4     0  ],[]],
## [[                                    ]   ]
## [[ 0      c3x1     0      c3x2    c3x3]   ]
## [[                                    ]   ]
## [[ 0      c2x1     0       0      c2x2]   ]
##
##     [2]          [3]          [3]          [4]          [4]
##   [[[ ], 2, 1], [[ ], 3, 2], [[ ], 1, 2], [[ ], 3, 2], [[ ], 2, 2]],
##     [1]          [1]          [2]          [1]          [2]
##
##   [[0    1    1    2]  [0    1]  [0    1    1]]
##   [[                ], [      ], [           ]], [{3}, {1, 2}, {2}]
##   [[0    0    1    0]  [0    1]  [0    0    1]]
##
## SEE ALSO:
## spresultant, `toric/polytope`
##
`toric/mixed` := proc(hulls::array, supps::array, trace::integer )
#	   						, delta , lifting
local	d,				# number of polytopes = dimension + 1
	delta,				# pertubration vector
	hullArr,			# hullArr[p]=[d-1,#hull.verts] 
	hull2supp,			# hull2supp[i]=vector of #hull monoms
	edata, i, ip, ipoly, ivert,
	j, k, zeroes,			# indices
	msize,				# sparse resultant (Newton) matrix M size
	mons,				# number of monoms in one polynomial
	v, vt1, vt2,
	lift, alift,			# linear lifting, affine : matrix d rows
	rowsperpol,			# number of matrix M rows per polynomial
	guide_rows, M,			# outputs: row coding, resultant matrix
	diagMons,			# array of sets w/diag.monoms per poly
	total_time,
	rows_time, constr_time,
	outList;
global  tot_time_lp;

    total_time := time();			# init
    d := op(op(2,eval(hulls)))[2];		# #polytopes = dimension + 1
    # d := rowdim(convert(hulls,matrix));

    delta := `toric/randvect` (d-1);	 		# fix type (array)
    if nargs>3 then
	for i to d-1 do delta[i] := args[4][i]; od;
    	if trace>0 then printf("toric/mixed(): Given Delta =") fi;
    else 
	if trace>0 then printf("toric/mixed(): Random Delta =") fi;
    fi;
    if trace>0 then 
	for i from 1 to d-1 do printf(" %.3f",delta[i]); od; printf("\n");
    fi;

    if nargs>4 then
	alift := args[5];
    	if trace>0 then print(`Given lifting =`,eval(alift)) fi;
    else
    	alift := `toric/randaff`(d);		# random affine functional, matrix dxd
	if trace>0 then printf("Random affine lifting =%a\n",eval(alift)) fi;
    fi;

    tot_time_lp := 0;				# initialize global timer
    if trace>0 then printf("toric/mixed(): sorting support points\n"); fi;

    hullArr := array(1..d);			# explicit subsets of supp
    for i from 1 to d do			# also sort arguments
      outList := `toric/sort_vecs` (hulls[i],supps[i]);
      hull2supp[i] := outList[1];		# hull2supp[hull.col]=supp.col
      hullArr[i] := outList[2];			# submatrix of new supps[i]
      # print('hull_to_supp@hull_arr',eval(hull2supp[i]),eval(hullArr[i]));
    od;

    if trace>0 then
      lprint(`The coeff labels correspond to the sorted supports:`);
      if d=2 then print(eval(supps[1]),eval(supps[2]));
      elif d=3 then print(eval(supps[1]),eval(supps[2]),eval(supps[3]));
      else for i from 1 to d do lprint(eval(supps[i])); od;
      fi;
    fi;

    if trace>0 then
      printf("toric/mixed(): computing data common to all vertex equations\n");
    fi;
    edata := `toric/make_eqns` (hullArr,alift,d);  # =[objective,conveq,vertlhs] 

    v := `toric/start_vertex` (hullArr,d,edata,delta);
    if trace>0 then printf("start vertex=%a.\n",eval(v)); fi;

    if trace>0 then printf("greedy search for row indices\n"); fi;
    rows_time := time();
    if d=3 then
	if trace>0 then printf("computing rows by special 2-dim routine\n"); fi;
	guide_rows := `toric/comp_rows_2` (hullArr, supps, edata, v, delta, trace);
    else
	if trace>0 then printf("computing rows by general-dim routine\n"); fi;
    	guide_rows := `toric/compute_rows` (hullArr, supps, edata, v, delta, d, trace);
    fi;
    rows_time := time() - rows_time;

    # if trace>0 then printf("toric/mixed(): got guide_rows="); lprint(guide_rows) fi;

    constr_time := time();
    msize := nops(guide_rows);
    M := array(1..msize, 1..msize);
    rowsperpol := vector(d,0); 

    # diagonal monomials

    diagMons := array(1..d);
    for i from 1 to d do diagMons[i] := {}; od;
    if trace>0 then printf("toric/mixed: init'd diagMons\n") fi;

    for i from 1 to msize do
      ivert := guide_rows[i];	# one row as list of 3 elts:
			      	# [Mink.sum pt, poly.indx, monom in hullArr]
      ip := ivert[2];		# 1 <= ip <= numpolys=d
      rowsperpol[ip] := rowsperpol[ip] + 1;
      vt1 := matadd(ivert[1], col(hullArr[ip], ivert[3]), 1, -1);
      ipoly := supps[ip];
      mons := coldim(ipoly);

      j := 1;
      for k from 1 to mons do
            vt2 := matadd(col(ipoly,k), vt1);	# vt2 = col monom
	    while not equal(eval(vt2), guide_rows[j][1]) do
		M[i, j] := 0;
		j := j + 1;
	    od;
	    M[i, j] := `c`||ip||`x`||k;
            j := j + 1;
      od;	# k-th monomial

      zeroes := j;				# rest of row is zero
      for j from zeroes to msize do M[i, j] := 0; od;

      # to deal w/degeneracy: indices of monomials on diagonal
      # ivert[3] indexes in hullArr, only supps passed on w/diagMons
      diagMons[ip] := diagMons[ip] union { hull2supp[ip][ivert[3]] };
      # print(eval(ivert[1]), eval(col(hullArr[ip],ivert[3])),eval(vt1),
      # ivert[3],eval(hull2supp[ip][ivert[3]]));

    od;		# i-th row

    constr_time := time() - constr_time; 
    if trace>0 then
	printf("toric/mixed: Delta = %a\n", eval(delta));
	printf("Greedy subdivision resultant matrix of dim %d x %d.\n",rowdim(M),coldim(M));
    fi;
    if rowdim(M)=1 then printf("SPARSE RESULTANT MAY BE POWER OF MATRIX DETERMINANT\n") fi;
    if trace>0 then
      printf("Number of rows per polynomial=%a\n",eval(rowsperpol));
      printf("Monomials on the diagonal are:\n",eval(diagMons));
    fi;

    if trace>0 then
      total_time := time() - total_time;
      printf("timings[secs] and percents:\n");
      printf("optimize=%as (%2.0f%%), find rows=%as (%2.0f%%),",
	tot_time_lp,100*tot_time_lp/total_time,rows_time,100*rows_time/total_time);
      printf("construct=%as (%2.0f%%), total=%as\n",
	constr_time,100*constr_time/total_time,total_time);
    fi;

    RETURN(evalm(M), eval(guide_rows), eval(supps),eval(diagMons));
end: # toric/mixed
#----------------------------------------------------------------------
addtohelp(spresultant):
#----------------------------------------------------------------------
## HELP spresultant
##
## spresultant - sparse (toric) resultant matrix of polynomials
##
## CALLING SEQUENCE:
##     spresultant(pols, vars);
##     spresultant(pols, vars, delta);
##     spresultant(pols, vars, delta, lifting);
##
## PARAMETERS:
##     pols - list of polynomials
##     vars - list of variables
##     delta - vector or list of real numbers
##     lifting - 2-dimensional array of lifting
## DESCRIPTION:
## - Compute the sparse (toric) resultant matrix of the given polynomials
##   by eliminating variables vars.  The number of rows in the
##   coefficients of the first polynomial is optimal (equal to the
##   respective degree of the sparse resultant).
##
## - Letting n be the number of variables, n+1 must be the number of
##   polynomials.
##
## - If delta is given, it is used as the geometric perturbation
##   applied to the Minkowski sum of the Newton polytopes.
##   It must contain n entries, sufficiently small in order not to perturb
##   points in the sum outside adjacent cells in the sum's subdivision.
##   If delta is not given, it is chosen randomly.
##
## - If lifting is given, it specifies the affine lifting functions
##   of the Newton polytopes, thus defining the mixed subdivision of
##   their Minkowski sum which will specify the matrix.
##   It must be a square matrix of dimension n+1, each row corresponding
##   to a Newton polytope.  If lifting is not given, it is chosen randomly.
##
## - The function returns a square matrix, whose determinant is a multiple
##   of the sparse resultant for generic coefficients.  For degenerate
##   coefficients, a perturbation can be defined by setting local variable
##   PERT_DEGEN_COEFS, as described by D'Andrea and Emiris (Computing Sparse
##   Projection Operators, In "Symbolic Computation: Solving Equations in
##   Algebra, Geometry, and Engineering, AMS, 2001).  The 2nd returned item
##   is a list of the monomials indexing the columns.
##
## - The matrix construction implements the greedy version, by Canny and
##   Pedersen (Tech. Report 1394, C.S. Dept, Cornell University, 1993),
##   of the algorithm by Canny and Emiris (Proc. AAECC-1993, LNCS 263,
##   pp. 89.  Final version J. ACM, 2000).
##
## EXAMPLES:
##> spresultant([x^2-x*z^2+y*x+z-1,x*y*z-3,2*x-x*y-z],[x,y]);
##
##   [-3   0   z  0  0]
##   [ 0  -3   0  0  z]
##   [-z   2  -1  0  0], [[1],[2],[2],[3],[3]]
##   [      2         ]   [1] [1] [2] [1] [2]
##   [z-1 -z   1  1  0]
##   [ 0  -1   0  2 -1]
##
##> spresultant([x+3*y,2*x-y+5*x^2*y,u0+u1*x+u2*y],[x,y],
## [.028,.078],matrix([[15,149,392],[222,381,684],[35,350,57]]));
##
##   [ 3  0 1  0  0  0]
##   [ 0  3 0  1  0  0]
##   [-1  0 2  0  0  5], [[1],[1],[2],[2],[2],[3]]
##   [u0 u2 0 u1  0  0]   [2] [3] [1] [2] [3] [2]
##   [ 0  0 0  0  3  1]
##   [ 0  0 0 u0 u2 u1]
##
## SEE ALSO:
## `toric/mixed`
##
spresultant := proc (in_polyns::list, elvars::list )
#			      		, delta, lifting
  local suppArr, outList,
	symbMat, specMat, 		# matrices
	elvarsArr, specSet,
	diagMon, colMons, pertCoef, coefval,
	polyns,
	n, i, j, pol,
	flag,
	req_args,
	PERT_DEGEN_COEFS,
	INTEGER_COEFFS;

  unassign('x'); 		# strings to be used later
  unassign('lam');

  # 0: no pert, 1: overcons, random pert, 2: u-res/random pert
  PERT_DEGEN_COEFS := 0;
  if PERT_DEGEN_COEFS=1 then
    printf("WILL PERTURB MATRIX BY RANDOM PERTURBATION OF OVERCONSTRAINED SYSTEM\n");
  elif PERT_DEGEN_COEFS=2 then
    printf("WILL PERTURB WELLCONSTRAINED SYSTEM, LEAVE 1st POLYNOMIAL AS IS\n");
  fi;

  INTEGER_COEFFS := true;
  if INTEGER_COEFFS then polyns := map(`toric/scaled`,convert(in_polyns,list));
  else polyns := convert(in_polyns,list); fi;

  req_args := 2;
  if nargs<req_args then ERROR(`small number of arguments to spresultant()`) fi;
  n := nops(elvars);
  if nops(polyns) <> (n+1) then ERROR(`bad #polynomials=`,nops(polyns)) fi;

  polyns := map(expand,polyns);
  for i from 1 to n+1 do for j from 1 to n do
    polyns[i] := simplify (polyns[i] * elvars[j]^(-ldegree(polyns[i],elvars[j])));
  od; od;

  # no output to file; outList=[suppArr, hullSet]
  # suppArr[p]=mat[dim,num], hullSet[p]=set of columns in suppArr[p]

  outList := `toric/polytope` (polyns, elvars, 0);

  # no detailed trace info printed by toric/mixed()
  # outList <- symb.matrix, row guide list, sorted supps, diag.monoms

  if nargs=req_args then
    outList := `toric/mixed` (eval(outList[2]), eval(outList[1]), 0);
  elif nargs=req_args+1 then
    outList := `toric/mixed` (eval(outList[2]), eval(outList[1]), 0, args[3]);
  elif nargs=req_args+2 then
    outList := `toric/mixed` (eval(outList[2]), eval(outList[1]), 0, args[3], args[4]);
  fi;

  # lprint(`spresultant() got`,outList,`now to de-assemble`);
  symbMat := eval(outList[1]);		# array/list

  # printf("The monomials indexing the matrix columns (2nd returned item):\n");
  colMons := map(el->el[1],outList[2]);
  # print (colMons);		 	# row guide list

  suppArr := eval(outList[3]);		# sorted supports used for cixj in symbMat
  if PERT_DEGEN_COEFS>0 then 
    diagMon := eval(outList[4]);	# array(1..d) of sets
  fi;

  elvarsArr := array(1..n,elvars);
  specSet := {};

  # treat special/degenerate coefficients by adding t-system

  if PERT_DEGEN_COEFS>0 then
    pertCoef := array(1..(n+1));
    # lprint(`spresultant() fills`,eval(pertCoef),`from`,eval(diagMon)); 
    for i from 1 to (n+1) do
      pertCoef[i] := vector(coldim(suppArr[i]),0);	# some pertCoef=0
      flag := false;					# exists pertCoef=1
      for j in diagMon[i] do				# index into support
	if not flag then pertCoef[i][j] := 1; flag := true;
	else
	  pertCoef[i][j] := RandSign()*rand(1..(10^ceil(evalf(log[2](n)))))();
 	  if not type(pertCoef[i][j],numeric) then ERROR(`no numeric pert.coeff`) fi;
	fi;
      od;
    od;	#for i
  fi;

  # i=1 must be u-polynomial in u-resultant
  for j from 1 to coldim(suppArr[1]) do
    coefval := coefofexp(polyns[1],elvarsArr,col(suppArr[1],j));
    if INTEGER_COEFFS then coefval := `toric/map_coefs`(round,coefval); fi;
    if PERT_DEGEN_COEFS>0 and type(PERT_DEGEN_COEFS,odd) then		
      coefval := coefval+t*pertCoef[1][j]; 	# overconstrained
    fi;
    specSet := eval(specSet) union {`c1x`||j = eval(coefval) }
  od;
  for i from 2 to (n+1) do for j from 1 to coldim(suppArr[i]) do
    coefval := coefofexp(polyns[i],elvarsArr,col(suppArr[i],j));
    if INTEGER_COEFFS then coefval := `toric/map_coefs`(round,coefval); fi;
    if PERT_DEGEN_COEFS>0 then coefval := coefval+t*pertCoef[i][j]; fi;
    specSet := eval(specSet) union {`c`||i||`x`||j = eval(coefval) }
  od; od;

  # printf("using specialized coefficients:%a\n",eval(specSet));

  specMat := matrix (rowdim(symbMat),coldim(symbMat),0); 
  for i from 1 to rowdim(symbMat) do for j from 1 to coldim(symbMat) do
      if symbMat[i,j]<>0 then
	specMat[i,j] := subs (specSet, symbMat[i,j]);
      fi;
  od; od;
  
  RETURN (evalm(specMat),eval(colMons));
end:  	# spresultant

# inputs:
# pointArr = dxn array point coords: d=dimension, n=#points
# limit = iteration limit
# output
# subset of pointArr column indices, expressing vertices
#
# Algorithm for (approximating) vertices of convex hull by
# sorting input points in random directions to find some vertices
#   but certain vertex may not be seen
#   in lo-dim hull, int.pt may be both min & max: dont accept it
#   tie in certain direction may give boundary point but not vertex
# Canny-Pedersen used to obtain randomized superset of vertices
# Now used to obtain certified subset of vertices
# deterministically use LinProg to check remaining points
#
`toric/comp_hull` := proc(pointArr::matrix, limit::integer)

  local	dim, i, j,
	mx, maxindex, maxtie, mn, minindex, mintie,
	num, rv,
      	vertset, unprocessed, interior,
	vertlhs, conveq, convsum, lambda, sub, temp;

    dim := rowdim(pointArr);		# dim = dimension
    num := coldim(pointArr);		# num = numbers of points

    # randomized computation of subset of verts if hull full dimensional

    vertset := {};
    unassign('i');
    unprocessed := {i$i=1..num};

    for i from 1 to limit do
        rv := randmatrix(1,dim);
        rv := scalarmul(rv, 1/ `toric/twonorm`(rv));
	# lprint(`randvector`,eval(rv));

        minindex := 1;
        maxindex := 1;
        mn := multiply(rv, col(pointArr, 1));
        mn := mn[1];
        mx := eval(mn);
	mintie := false; maxtie := false;

        for j from 2 to num do
            temp := multiply(rv, col(pointArr, j));
            temp := temp[1];
            if temp < mn then
                minindex := j;
                mn := temp;
		mintie := false;
            elif temp > mx then
                maxindex := j;
                mx := temp;
		maxtie := false;
	    else
		if temp=mn then mintie := true; fi;
 	    	if temp=mx then maxtie := true; fi;
            fi;
        od; #j

	# mn=mx iff degenerate hull, then dont decide
	# mintie or maxtie may have extreme pt, not vertex: dont decide
	if minindex<>maxindex then
	  if not mintie then
            vertset := vertset union {minindex};
            unprocessed := unprocessed minus {minindex};
  	    # lprint(`vertex w/min inn.prod`,minindex,eval(mn));
	  fi;
	  if not maxtie then
            vertset := vertset union {maxindex};
            unprocessed := unprocessed minus {maxindex};
  	    # lprint(`vertex w/max inn.prod`,maxindex,eval(mx));
	  fi;
	fi;

    od; #i: ends `limit` random iterations

    # lprint (`undecided`,unprocessed); 

    # Now use LP to check other vertices

    lambda := matrix(1, num);	 		# new variables
    for i from 1 to num do lambda[1,i] := cat(`lam`,i) od;
    interior := {};
    unassign('j');
    temp := { j$j=1..num };			# check pt against hull(temp)
    
    # lprint (`subset of verts`,vertset);
    # lprint(`check undecided:`,eval(unprocessed),`by LP using:`,eval(temp));

    # while (work minus work1 <> {}) or (work1 minus work <> {}) do work:=work1;

    for i in unprocessed do
    	    temp := temp minus {i};
            sub := convert(temp,list);		# superset of vertices
            convsum := 0;
            for j in temp do
	      convsum := convsum + lambda[1,j]
	    od;
            conveq := {convsum = 1}; 			# convexity
            vertlhs := multiply(submatrix(pointArr, 1..dim, sub),
                                transpose(submatrix(lambda, [1], sub)));
            for j from 1 to dim do
                conveq := conveq union {vertlhs[j, 1] = pointArr[j, i]}
            od;
            # lprint (`running LP on`,col(pointArr,i));
            if feasible(conveq, NONNEGATIVE) then 
                # lprint (`LP:interior point`,col(pointArr,i));
		interior := interior union {i};
		if member(i,vertset) then
		  ERROR(`interior pt in vtx subset`)
		fi;
	    else
		vertset := vertset union {i};
    	    	temp := temp union {i};
            fi;
    od; #i
    eval(vertset);
    # submatrix(pointArr, 1..dim, convert(vertset, list));
end: # `toric/comp_hull`

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
# GEOMETRY OF CURVES AND SURFACES                                      #
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
## Functions for dealing with implicit curves and surfaces in 2D or 3D.
## ---------------------------------------------------------------------
addtohelp(develop2d):
## ---------------------------------------------------------------------
## HELP develop2d
##
## develop2d - developped curve of an implicite planar curve
##
## CALLING SEQUENCE:
##     develop2d(p)
##     develop2d(p,var)
##
## PARAMETERS:
##     p   - implicit equation of the curve
##     var - list of two variables (optional)
## 
## DESCRIPTION:
## - Compute the curve of centers of curvature.
##
## - The equation may contain extra factors corresponding to 
##   circles centred at the singular points of the curve p(x,y)=0.
## 
## EXAMPLES:
## >develop2d(x^2-y^3);
## [x, y^4, x*y^3, y^2*x, x*y, y^3, y^2, y^5]
## [_y[1]^3, _y[1]*_y[2]^2, _y[1]*_y[2]^3, _y[1]*_y[2], _y[1]^2, 
##  _y[1]^2*_y[2], _y[2]^2, _y[2]^3, _y[1]]
## 7
##           4         3         2                  2           4        2
## -6613488 y  (18432 y  + 6144 y  + 512 y + 15552 x  y + 6561 x  + 288 x )
##
## SEE ALSO:
## melim
##
develop2d := proc(f)
local L1,L2,a,b,c,var,mthd;
  if nargs <2 then var := [op(sort(indets(f)))] else var:=args[3] fi;
  if nargs <4 then mthd := op(1,melim) else mthd := args[4] fi;
  a := diff(f,var[2]); b := -diff(f,var[1]); 
  c := diff(f,var[2])*(-var[1])+diff(f,var[1])*var[2];
  L1 := a*u +b*v +c;
  L2 := (diff(a,var[1])*diff(f,var[2])-diff(a,var[2])*diff(f,var[1]))*u
       +(diff(b,var[1])*diff(f,var[2])-diff(b,var[2])*diff(f,var[1]))*v
       +(diff(c,var[1])*diff(f,var[2])-diff(c,var[2])*diff(f,var[1]));
  subs(u=var[1],v=var[2], melim([L1,L2,f],[var[1],var[2]]));
  factor(%);
end:

## ----------------------------------------------------------------------
addtohelp(offset2d):
## ----------------------------------------------------------------------
## HELP offset2d
##
## offset2d - offset of an implicite planar curve
##
## CALLING SEQUENCE:
##     offset2d(p,d)
##     offset2d(p,d,var);
##     offset2d(p,d,var,mthd);
##
## PARAMETERS:
##     p    - implicit equation of the curve
##     d    - distance between the offset and the initial curves.
##     var  - the list of variables (optional)
##     mthd - the method of elimination of variables (optional)
##
## DESCRIPTION:
##
## - Compute a mutiple of the equation of the set of points at distance d
##   from the implicite curve p(x,y)=0.
##
## - The equation may contain extra factors corresponding to 
##   circles centred at the singular points of the curve p(x,y)=0.
## 
## EXAMPLES:
## >offset2d(x^2-y^3,1);
## [1, x, y, y^2*x, x*y, y^4, y^3, y^2]
## [_y[1]^2*_y[2], _y[1]*_y[2], _y[1]^2, _y[2]^2,_y[1]*_y[2]^2,_y[1],_y[2],1]
## 7
##                3         2           3  2         4           4        2  2
## 9 (529 + 4072 y  + 5814 x  y - 4892 y  x  - 4158 x  y - 1685 x  - 297 x  y
## 
##              4        2        2  4         4  2        6         3  4
##      + 3870 y  + 427 x  + 729 y  x  - 2619 y  x  + 729 x  - 1458 y  x
## 
##              5  2         5         6        7        8                  2
##      - 1458 y  x  - 2376 y  - 2900 y  + 216 y  + 729 y  - 1656 y - 1188 y
## 
##             6  2    2    2     3
##      + 729 y  x ) (y  + x  - 1)
## SEE ALSO:
## melim
##
offset2d := proc (f,d)
local var,mthd;
 if nargs <3 then var := [op(sort(indets(f)))] else var:=args[3] fi;
 if nargs <4 then mthd := op(1,melim) else mthd := args[4] fi;
 factor(mthd([f,diff(f,var[2])*(var[1]-_u)-diff(f,var[1])*(var[2]-_v),(var[1]-_u)^2+(var[2]-_v)^2-d^2],[var[1],var[2]]));
 subs(_u=var[1],_v=var[2],%);
end:
## ----------------------------------------------------------------------
addtohelp(median2d):
## ----------------------------------------------------------------------
## HELP median2d
##
## median2d - the mediatrice of two implicite planar curves.
##
## CALLING SEQUENCE:
##     median2d(f,g)
##     median2d(f,g,var);
##     median2d(f,g,var,mthd);
##
## PARAMETERS:
##     f,g  - implicit equations of the curves.
##     var  - the list of variables (optional)
##     mthd - the method of elimination of variables (optional)
##
## DESCRIPTION:
##
## - Compute a mutiple of the equation of the set of points at equidistance 
##   from the implicite curves  f(x,y)=0 and g(x,y)=0.
##
## - The equation may contain extra factors.
## 
## EXAMPLES:
## >median2d(x^2-y^3,x*y-1);
##                3         2           3  2         4           4        2  2
##                           3        8  3        7          2  6        3  5
## -1073741824 (-64 + 768 x y  - 112 x  y  - 192 y  x + 144 x  y  + 384 x  y
## 
##             4  4        5  3       7         8        7  2        3  6
##      - 396 x  y  + 288 x  y  + 48 x  y + 36 y  + 192 x  y  + 192 x  y
## 
##            4  5        6  3        5  4       8         8       6
##      + 96 x  y  + 496 x  y  - 768 x  y  + 12 x  y - 49 x  + 32 y
## 
##             2  3        4         10       9  2       4  7        6  5
##      + 128 x  y  + 384 x  y - 12 y   + 12 y  x  - 96 x  y  + 168 x  y
## 
##            9         8  2       4  6       7  3       5  5       8  2
##      + 24 y  x + 48 y  x  - 76 x  y  - 96 y  x  + 96 x  y  - 12 x  y
## 
##             6  4                  2        4    12        2  2       4
##      - 162 x  y  + 384 x y - 192 y  - 144 y  + y   - 864 x  y  - 16 x
## 
##             2        10  2       4  8      6  6      10        6
##      - 192 x  y - 6 y   x  + 12 x  y  - 8 x  y  + 8 x   - 120 x  y
## 
##             4  3        5  2        2  5        2  4       4  2
##      + 576 x  y  - 768 x  y  - 192 x  y  - 960 x  y  - 96 x  y
## 
##             3  3       5          5         6
##      + 896 x  y  - 96 x  y + 192 y  x + 72 x )
## 
##       4      2      2  2            4          2 3
##     (x  - 4 x  - 2 x  y  + 8 x y + y  - 4 - 4 y )         
## 
## SEE ALSO:
## melim
##
median2d := proc(f,g)
local fu, gv, dfu, dgv, dst,var,mthd;  
global melim;
 if nargs <3 then var := [op(sort(indets(f)))] else var:=args[3] fi;
 if nargs <4 then mthd := op(1,melim) else mthd := args[4] fi;
 fu:= subs(var[1]=_u[1],var[2]=_u[2],f);
 gv:= subs(var[1]=_v[1],var[2]=_v[2],g);
 dfu := diff(fu,_u[1])*(x-_u[1]) +  diff(fu,_y[1])*(y-_u[2]);
 dgv := diff(gv,_v[1])*(x-_v[1]) +  diff(gv,_v[1])*(y-_v[2]);
 dst := expand((x-_u[1])^2+(y-_u[2])^2-(x-_v[1])^2-(y-_v[2])^2);
 factor(mthd([fu,gv,dfu,dgv,dst],[_u[1],_u[2],_v[1],_v[2]]));
end:

# ----------------------------------------------------------------------
# UNDOCUMENTED
# ----------------------------------------------------------------------
matrixout := proc(lp ::list, var, lom)
local lm,i,j,ll;
  lm := NULL;
  for i to nops(lp) do lm := lm, lstmonof(lp[i],var);  od;
  lm := {lm} minus {op(lom)};
  lm := sort([op(lm)]);
  lprint(lm);
  ll := NULL;
  for i to nops(lm) do
     for j to nops(lp) do
      ll := ll, coeffof(lm[i],expand(lp[j]),var);
     od;
  od;
  matrix(nops(lm),nops(lp),[ll]);
end:

#----------------------------------------------------------------------#
`type/matnum` := proc(x)
local u,r,i,j;
  if type(x,matrix) then
   r := true;
   for i to rowdim(x) while r do
     for j to coldim(x) while r do
        r := r and type(x[i,j],numeric);
     od;
   od;
   r;
  else
   false;
  fi;
end:
#----------------------------------------------------------------------#
[__];
