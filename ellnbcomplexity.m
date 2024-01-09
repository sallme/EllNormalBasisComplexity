/*************************** Reference ******************************************
** The algorithms described in this package come from the following paper
** "The complexity of elliptic normal bases", D. Panario, M. Sall, and Q. Wang.
** To appear 
*********************************************************************************/

/**************************** Summary  *******************************************
The package contains two main functions. The function 'ENBparamsComputation' 
which takes as input a finite field Fq and an extension degree n and then returns 
(if possible) all the necessary parameters to the construction of an elliptic normal 
basis of F_{q^n} over Fq. The function 'ENBcomplexityBounds' which computes the 
weight of the special vectors and gives a lower and an upper bound on the complexity 
of the basis. These functions make use of several preliminary functions. Among them
we have the following

ComponentWiseProduct    ---> Computes the component wise product of two vectors.
KCyclicShift            ---> Computes the k-cyclic shift of a vector V.
ConvolutionProductNaive ---> Computes the convolution product of two vectors.
InvConvolutionProduct   ---> Computes the convolutional inverse of a vector.

EllCurvePrescribedOrderPossibility --> Tests if a curve with some order exists.
EllCurvePrescribedOrder            --> Define an elliptic curve of given order.
EllNormalBasisPossibility          --> Tests if an elliptic normal basis exists. 
**********************************************************************************/



/********************* Functions on vectors of length n ************************/

ComponentWiseProduct := function(V, W)
   U := 0*V;
   n := Ncols(V);
   for i:=1 to n do
	U[i] := V[i]*W[i];
   end for;
   return U;
end function;

KCyclicShift := function(V, k)
   W := 0*V; n := Ncols(V);
   for i:=0 to n-1 do
	W[i+1] := V[1+((i-k) mod n)];
   end for;
   return W;
end function;

ConvolutionProductNaive := function(U, V)
   W := 0*U; n := Ncols(U);
   for j := 0 to n-1 do
       for i:= 0 to n-1 do
	  W[j+1] +:= U[i+1]*V[1+((j-i) mod n)];
       end for;
   end for;
   return W;
end function;

InvConvolutionProduct := function(V, Fq) 
   Vi := 0*V; n := Ncols(V);
   Px := PolynomialRing(Fq); x := Px.1;
   _, vx, _ := ExtendedGreatestCommonDivisor(Px!Eltseq(V), x^n-1);
   VX := Eltseq(vx); VX cat:= [Fq!0 : i in [#VX+1..n]]; 
   Vi := VectorSpace(Fq, n)!VX;
   return Vi;
end function;


/************************ Functions on elliptic curve group *********************/

EllCurvePrescribedOrderPossibility := function(Fq, d)

    q := #Fq; p:= Characteristic(Fq); n := Degree(Fq); m := q+1-d;

    /* d can't be too small nor too large */
    if m^2 gt 4*q then return false; end if;

    /* When p divides m, only the (very few) supersingular curves exists  */
    if (m mod p) eq 0 then
	if IsOdd(n) then
	    if m ne 0 and
	       m ne Floor(Sqrt(p*q)) and
	       m ne -Floor(Sqrt(p*q)) then
	       return false;
	    end if;
	else
	    if (m ne 2*Floor(Sqrt(q)) and
	        m ne -2*Floor(Sqrt(q))) and
	       ((m ne Floor(Sqrt(q)) and
	         m ne -Floor(Sqrt(q))) or
		(p ne 3 and (p mod 3) ne 2)) and
	       ((m ne 0) or
	       (((p mod 4) ne 2) and ((p mod 4) ne 3))) then
	       return false;
	    end if;
	end if;
    end if;

    return true;
    
end function;



EllCurvePrescribedOrder := function(Fp, order)

    if EllCurvePrescribedOrderPossibility(Fp, order) eq false then return false; end if;

    found := false;
    repeat
	for i := 1 to 10 do
	    repeat
		repeat
		    a1 := Random(Fp); a2 := Random(Fp); a3 := Random(Fp);
		    a4 := Random(Fp); a6 := Random(Fp);

		    Disc := -a1^6*a6 + a1^5*a3*a4 - a1^4*a2*a3^2 - 12*a1^4*a2*a6 + a1^4*a4^2 + 
			8*a1^3*a2*a3*a4 + a1^3*a3^3 + 36*a1^3*a3*a6 - 8*a1^2*a2^2*a3^2 - 
			48*a1^2*a2^2*a6 + 8*a1^2*a2*a4^2 - 30*a1^2*a3^2*a4 + 72*a1^2*a4*a6 + 
			16*a1*a2^2*a3*a4 + 36*a1*a2*a3^3 + 144*a1*a2*a3*a6 - 96*a1*a3*a4^2 - 
			16*a2^3*a3^2 - 64*a2^3*a6 + 16*a2^2*a4^2 + 72*a2*a3^2*a4 + 288*a2*a4*a6 - 
			27*a3^4 - 216*a3^2*a6 - 64*a4^3 - 432*a6^2;

		until Disc ne 0;
		Ec := EllipticCurve([a1, a2, a3, a4, a6]);
	    until (#Ec mod order) eq 0;

	    P := Random(Ec);
	    if (Order(P) mod order) eq 0 and  (Order(P) mod order^2) ne 0 then
		found := true; break;
	    end if;
	end for;
    until found eq true;

    return Ec, (Order(P) div order)*P;
	
end function;

/* Checking if there exists an elliptic normal basis of degre d over Fq */
EllNormalBasisPossibility := function(Fq, d)

    q := #Fq; p:= Characteristic(Fq); n := Degree(Fq); m := q+1-d;

    if d le 1 then return false; end if;
    if d eq q-1 and not q in [2, 2^2, 2^3, 3, 5, 7] then return false; end if;

    k := Ceiling(q+1-2*Sqrt(q)) div d;
    dp := k*d; found := false;
    while dp le Floor(q+1+2*Sqrt(q)) and found eq false do
	found := EllCurvePrescribedOrderPossibility(Fq, dp);
	dp +:= d;
    end while;

    return found;
	
end function;

EllPointPrescribedOrder := function(Ec, m)
	
	if (#Ec mod m) ne 0 then return "Impossible to find this point"; end if; 

	repeat t := Random(Ec); until IsOrder(t, m) eq true;

    return t;

end function;

/* The function f_OA evaluated at P */
ENBfOA := function(A, P, Ec)

    if A eq Ec!0 then return 1; end if;

    a1, a2, a3, a4, a6 := Explode(aInvariants(Ec));
    x1, y1, _ := Explode(Eltseq(A));
    x, y,   _ := Explode(Eltseq(P));

    return (y+y1+a1*x1+a3)/(x-x1);
end function;


/* The function f_AB evaluated at P */
ENBfAB := function(A, B, P, Ec)

    a1, a2, a3, a4, a6 := Explode(aInvariants(Ec));

    if A eq Ec!0 then
	return ENBfOA(B, P, Ec);
    end if;

    if B eq Ec!0 then
	return -ENBfOA(A, P, Ec)-a1;
    end if;

    x1, y1, _ := Explode(Eltseq(A));
    x2, y2, _ := Explode(Eltseq(B));
    x,  y,  _ := Explode(Eltseq(P));

    if A eq -B then
	return -(3*x1^2 + 2*a2*x1 + a4 -a1*y1) / (2*y1+a1*x1+a3)-(x*a1+a3+2*y1)/(x-x1);
    end if;
 
    return (y1+y2+a3+a1*x1)/(x2-x1)-((-x2+x1)*y+(-y2+y1+a1*(x1-x2))*x+(-(y1+a1*x1+a3)*x2+x1*(y2+a1*x2+a3)))/(x-x1)/(x-x2);

end function;


/************************* Elliptic Normal Basis Parameters ********************/

ENBparamsComputation := function(Fq, n)    
	if Type(EllCurvePrescribedOrder(Fq, n)) eq BoolElt then
	    return "Sorry, we can't define an elliptic normal basis over this finite field extension";
	else 
        Ec, t := EllCurvePrescribedOrder(Fq, n);
	end if;
    q := #Fq; Px := PolynomialRing(Fq); x := Px.1;	    
    printf "(Confirmation) The order of the Kernel is %o \n", Order(t);
    /* Computation of an isogeny with kernel order n */
	DX := &*{(x-(k*t)[1]) : k in [1..Order(t)-1]};
    Ei, In := IsogenyFromKernel(Ec, DX);
    nx := IsogenyMapPhi(In); dx := IsogenyMapPsiSquared(In);
    GCDx := Gcd(nx, dx); nx := nx / GCDx; dx := dx / GCDx;

    /****** Constructing an extension from the isogeny ********/
    found := false;
    repeat pt := Random(Ei); until pt ne Ei!0;
    printf "(CheckProcess) We try the point %o of order %o.\n", pt, Order(pt);
    eqP := Px!(nx-pt[1]*dx);
    printf "(CheckProcess) Is the rational point generating an irreducible preimage? %o! \n", IsIrreducible(eqP);
    if not IsIrreducible(eqP) then 
        printf "(DeclineTrial) So, we can't use this point. Let's try again!\n";
    end if;    
    while IsIrreducible(eqP) eq false do
	repeat pt := Random(Ei); until pt ne Ei!0;
        eqP := Px!(nx-pt[1]*dx);
        printf "(CheckProcess) We have tried the point %o of order %o \n", pt, Order(pt);
    end while;     
    printf "(Confirmation) Now the point %o of order %o works. It gives:\n", pt, Order(pt);
    printf "(FoundElement) %o.\n", eqP;
    printf "(Verification) Is the preimage of the point irreducible? %o!\n", IsIrreducible(eqP);
    Fqn := ext<Fq|eqP/Coefficient(eqP, Degree(eqP))>;   

    /********** Looking for a generic point Pb ****************/
    En := BaseExtend(Ec, Fqn);
    y  := Fqn!1;
    Pb := Points(En, y)[1];
    repeat
      if n*Pb eq En!0 then
 	        printf "(CheckProcess) The point %o is a %o torsion point. Let's try again.\n", Pb, n;
      end if;
      Pb := Random(En);// Points(En, y)[1];
    until n*Pb ne En!0;
    printf "(Confirmation) We found a generic point:\n";
    printf "(FoundElement) %o\n", Pb;
        
    /****** Looking for the constant scalares ********/
    printf "(CheckProcess) Now let us look for the constant scalars.\n";
    p := Characteristic(Fq);    
    if p eq 2 or p eq 3 then
       Fq2 := GF(#Fq^2);
    else
       Fq2 := GF(#Fq^2);
    end if;
    Ec2 := BaseExtend(Ec, Fq2);
    repeat Pp := Random(Ec2); until Evaluate(DX, Pp[1]) ne 0;
    Qp := Ec!0; Q := Ec!t;
    cstc := Fq2!ENBfAB(Qp, Q, Pp, Ec);
    for i:=1 to n-1 do
	    Qp := Q; Q +:= t; cstc +:= Fq2!ENBfAB(Qp, Q, Pp, Ec);
    end for;
    cstc := Fq!cstc;    
    if cstc eq 0 then
	    csta := Fq!1; cstb := (Fq!1) / n;
    else
	    csta := 1 / cstc; cstb := Fq!0;
    end if;
    printf "(Confirmation) We found the following three constants: (%o, %o, %o).\n", cstc, csta, cstb; 

    /**************** Elliptic Normal Basis Construction ******************/
    printf "(Confirmation) Everything is settled, we can construct a elliptic normal basis.\n";
    ENB := []; //VectorSpace(Fqn, n)!0;
    ENB[1] := csta*ENBfOA(En!t, Pb, En);
    printf "(CheckProcess) We work with %o\n", ENB[1];
    ENBtest := IsNormal(ENB[1]);
    printf "(Verification) Is this found element an elliptic normal element? %o\n", ENBtest; 
    for i:=2 to n do 
        ENB[i] := ENB[1]^(q^(i-1));
        printf "(Verification) Are the elements normal? %o \n", IsNormal(ENB[i]);
    end for;
	
    /******* Looking for parameters for complexity computation *******/
    Gp, Mp := AbelianGroup(Ec); R := Ec!0;
    for g in Setseq(Generators(Gp)) do
	if n mod Order(Mp(g)) ne 0 and Order(Mp(g)) ne 1 then
	   R := Mp(g);
           break;
	end if;
    end for;    
    if R eq Ec!0 then 
	    printf "(DeclineTrial) But it does not allow fast multiplication. We give up this trial!\n";
        return "(Work-arounds) If the field extension is appropriate, you can try another curve!";
    else
	    Xr := VectorSpace(Fq, n)!0; Ur := VectorSpace(Fq, n)!0;
	    Q := R; Xr[1] := Eltseq(Q)[1]; Ur[1] := csta * ENBfOA(t, Q, Ec) + cstb;
	    for i := 1 to n-1 do
	        Q +:= t; Xr[i+1] := Eltseq(Q)[1]; Ur[i+1] := csta*ENBfOA(t, Q, Ec)+cstb; 
	    end for;
        printf "(Confirmation) A point %o of order %o outside the kernel is found.\n", R, Order(R);
	    printf "(Confirmation) Further, we found the following principal special vectors:\n";
        printf "(SpecialsVect) %o\n", Ur;
        Uri := InvConvolutionProduct(Ur, Fq);
        printf "(InvSpeciVect) %o\n", Uri;
        printf "(Verification) %o\n", ConvolutionProductNaive(Ur, Uri);
        printf "(SpecialsVect) %o\n", Xr;
        printf "(CheckProcess) Special vectors from component-wise product cyclic shift:\n";
	        for k:=1 to n do
                    Vck := KCyclicShift(Ur, k);
                    Vrk := ComponentWiseProduct(Ur, Vck); 
		    printf "(SpecialsVect) Case k = %o: %o \n", k, Vrk;
	        end for;
    end if;

    /***************** Constructing a structure that store the parameters *******************/
    ENBformat := recformat<Ec, En, n, In, Ei, Pb, Pt, ca, cb, cc, ENB, R, Vr, Vrx>;
    ENBparams := rec<ENBformat | 
                 Ec := Ec, En:= En, In := In, Ei := Ei, n:=n, Pb := Pb, Pt := t, 
                 ca := csta, cb := cstb, cc := cstc, ENB := ENB, R := R, Vr := Ur, Vrx := Xr>;
    printf "(Confirmation) All the necessary parameters are returned!!! \n";

    return ENBparams;

end function;


/************** Bounds on the complexity of elliptic normal basis **************/
ENBcomplexityBounds := function(Params)
   Fq  := BaseField(Params`Ec); n := Params`n;
   Vr  := Params`Vr; Vrx := Params`Vrx;
   Vri := InvConvolutionProduct(Vr, Fq);
   sum := 0;
   printf "(CheckProcess) Convolution product inverse vector R and k-cyclic shift:\n";
   for k:=2 to n-2 do
       Vrk := ComponentWiseProduct(Vr,  KCyclicShift(Vr, k));   
       VriVrk := ConvolutionProductNaive(Vri, Vrk); 
       wt := Weight(VriVrk);
       printf "(SpecialsVect) Case k = %o:  %o. The weight is %o.\n", k, VriVrk, wt;
       sum +:= wt;
   end for;
   LowB := sum + 3;
   UpB  := sum + 3*n; UpBnb := n^2 - n + 1;
   if UpB ge UpBnb then UpB := UpBnb; end if; 
   printf "The lower bound of the elliptic normal basis complexity is %o.\n", LowB;
   printf "The upper bound of the elliptic normal basis complexity is %o.\n", UpB;
   return "!!!!!!!!!!!!!!!!! End of computations !!!!!!!!!!!!!!!!!!!!!!!!";
end function;

/*************** Computing the weight of the three remaining rows ****************/
ENB3RowsWeight := function(Fqn, Nelt, x_b, fa, Vr, Vrx)
    
    Fq := BaseField(Fqn); n := Degree(Fqn, Fq); q:=#Fq;
	NB := []; NB[1] := Nelt; 
    for i:= 2 to n do
        NB[i] := Nelt^(q^(i-1));
    end for;	
    //Naive method for the first row CP_0.
    VFqn, Psi := VectorSpace(Fqn, Fq, NB);
    Iota := Psi(x_b);
    CP_0 := Psi(Nelt*Nelt);
    //Using new results for the second, CP_1, and last rows, CP_n1.
    FAv := 0*Vr; FAv[2] := -fa^2;
    Conv1    := ConvolutionProductNaive(Iota, FAv);
    Conv2    := ConvolutionProductNaive(Vrx, FAv); 
    SumConv2 := KCyclicShift(Vr, 1)-Conv2;
    Conv3    := ConvolutionProductNaive(InvConvolutionProduct(Vr, Fq), SumConv2);
    CP_1     := Conv1 + Conv3;
    FAv_n1    := 0*Vr; FAv_n1[1] := -fa^2;
    Conv1_n1  := ConvolutionProductNaive(Iota, FAv_n1);
    Conv2_n1     := ConvolutionProductNaive(Vrx, FAv_n1);
    SumConv2_n1  := KCyclicShift(Vr, n-1)-Conv2_n1;
    Conv3_n1    := ConvolutionProductNaive(InvConvolutionProduct(Vr, Fq), SumConv2_n1);
    CP_n1       := Conv1_n1 + Conv3_n1; 
	
    return CP_0, Iota, CP_1, CP_n1; 
	
end function;
