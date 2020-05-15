# mixvol.mpl
#
# Ioannis Z Emiris (c) INRIA Sophia-Antipolis
# Ioannis.Emiris@inria.fr, http://www-sop.inria.fr/galaad/emiris
#
# created 02/01
# last update 12/01
#
# MapleVrel5 routines
# calling (stable) mixed volume C program

with(linalg):
#----------------------------------------------------------------------
addtohelp(mixvol):
#----------------------------------------------------------------------
## HELP mixvol
##
## mixvol -- (stable) mixed volumes of square Laurent polynomial systems
##
## CALLING SEQUENCE:
##     mixvol ( polynomials );
##     mixvol ( polynomials, sv );
##     mixvol ( polynomials, sv, filename );
##
## PARAMETERS:
##     polynomials - list of Laurent polynomials (integer exponents)
##     sv - 1 in order for stable mixed volume to be computed
##     filename - prefix of files used
##
## DESCRIPTION:
## - The number of given polynomials must equal the total number
##   of variables; exponents must be integer, including negative ints.
##
## - Writes filename.sites containing supersets of the polynomial
##   supports; the C program shall compute Newton polytopes.
##   If filename not given, writes MVtemp.sites.
## 
## - Calls C program computing (stable) mixed volume.
##
## - Reads in the (stable) mixed volume output in filename.out
##   or MVtemp.out.  Both input/output files are removed.  
##   Return the integer (stable) mixed volume from the output file.
##
## - Possible to include further options of the C program such as
##   stable volumes with respect to selected coordinates, cells,
##   tracing info.
##
## EXAMPLES:
##> mixvol([3*x^3+x^2-y*x,2-y]);      
## C program took 0sec, computed mixed volume = 2
##                                             2
##
##>  mixvol([3*x^3+x^2-y*x,2-y],1);
## C program took 0sec, computed stable volume = 3
##                                             3
##
mixvol := proc (in_polys::list)
#				, sv, fname)
  local varslist, polylist, numvars, numpols,
	reqargs,
	infile, outfile, out,		# files used by C program
	supps, mons,			# array of supports; monomials
	result,				# list
	cmd, sv,			# command string, flag
        mostneg,			# integer
	p, v, m;			# indices

  reqargs := 1;

  varslist := NULL;
  for p in in_polys do varslist := eval(varslist),op(indets(p)) od;
  varslist := sort([op({varslist})]);
  numvars := nops(varslist);

  polylist := map(p->collect(simplify(expand(p)), varslist, distributed), in_polys);
  numpols := nops(polylist);
  if numpols < numvars then
    print(polylist,varslist);
    ERROR(`fewer polyns than vars`,numpols,numvars);
  elif numpols > numvars then
    lprint(`fewer variables than polynomials, hence MV = 0`);
    RETURN(0);
  fi;

  sv := false;
  if nargs>reqargs and not type(args[reqargs+1],integer) then
    ERROR(`1st optional argument must be integer`);
    elif nargs>reqargs and args[reqargs+1]=1 then sv:= true;
  fi;

  supps := array (1..numpols);
  for p from 1 to numpols do
    if type(polylist[p],`+`) then mons := op(polylist[p]);
      else mons := polylist[p];
    fi;
    supps[p] := sort([mons]);
  od;

  if nargs = 3 then
    infile := cat(args[3], ".sites");
    outfile := cat(args[3], ".out");
  else
    infile := "MVtemp.sites"; outfile := "MVtemp.out";
  fi;

  writeto(infile);
  printf("%d\n",numvars);
  printf("%d\n",numpols);
  printf("\n");
  for p from 1 to numpols do printf("%d\n",nops(supps[p])); od;
  printf("\n");
  
  mostneg := 0;
  for p from 1 to numpols do
    for m from 1 to nops(supps[p]) do
      for v from 1 to numvars do
	printf("%d ", degree(supps[p][m],varslist[v]));
	if degree(supps[p][m],varslist[v])<mostneg then
	  mostneg:=degree(supps[p][m],varslist[v]) fi;
      od;
      printf("\n");
    od;
    printf("\n");
  od;

  writeto(terminal);
  if sv then cmd := "S" else cmd := "" fi;
  cmd := cat(eval(cmd),"Mixvol -n ");
  if sv then cmd := cat(eval(cmd),"-G "); fi;
  if mostneg<0 then cmd := cat(eval(cmd),"-a ",eval(-mostneg)); fi;
  cmd := cat(eval(cmd), " ", infile, " ", outfile);
  cmd := cat(eval(cmd), "; grep Total..v.= MVtemp.out |awk '{print $4}' | sed 's/[.]//g' > MVtemp.mv");
  printf("calling %s\n",cmd);
  system(eval(cmd));
  # system(cat("rm ",infile)); 
  print(11111);
  out := fopen(`MVtemp.mv`,READ,TEXT);
  result := fscanf(out,"%d\n");
if false then # broken on Maple7
  out := fopen(outfile,READ,TEXT);
  result := [];
  # print(outfile,result);
  if sv then while nops(result)=0 do
    result := fscanf(out,"\nTotal sv = %d. Total time = %dh %dm %ds.");
  od; else while nops(result)=0 do
    # result := fscanf(out,"\nTotal mv = %d. Total"); # Broken on maple6
    result := fscanf(out,"\nTotal mv = %d. Total time = %dh %dm %ds.");
  od; fi;
fi:
  close(out);
  # print(result);

if false then
  system(cat("rm ",outfile)); 
  printf("C program took %dsec, computed ",op(2,result)*3600+op(3,result)*60+op(4,result));
  if sv then
    printf("stable volume = %d\n",op(1,result));
    else printf("mixed volume = %d\n",op(1,result));
  fi;
fi;

  RETURN(op(1,result));

end:# mixvol
