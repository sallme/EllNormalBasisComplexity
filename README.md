# EllNormalBasisComplexity
This package provide functions that define an elliptic normal basis and give a
bounds on its complexity. The package consists of the file 'ellnbcomplexity.m' 
itself which contains several prelimenary functions that handle operations 
on special vectors define in the follow artcile:
"D. Panario, M. Sall, Q. Wang The complexity of elliptic normal bases" (to appear).
The prelimary functions help also in the computation of elliptic functions of degree
two through which we define an elliptic normal bases. For more details about these
function we refer to the beginning of the file 'ellnbcomplexity.m'.
If an elliptic normal basis exist on the considered extension, the main function 
ENBparamsComputation computes its parameters step by step and returns a format that 
contains all these parameters. These parameters are used by the function 
'ENBcomplexityBounds' that computes the weight of the special vectors and returns
an upper bound and a lower bound on the complexity of the basis.
If there is no elliptic normal basis, the package returns the following statement
"Sorry, we can't define an elliptic normal basis over this finite field extension".
