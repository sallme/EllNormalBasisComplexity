# EllNormalBasisComplexity
This package provides functions that define an elliptic normal basis and give 
bounds on its complexity. By complexity, we mean the number of non-zero elements
of the multiplication table. The package consists of the file 'ellnbcomplexity.m' 
itself which contains several preliminary functions that handle operations 
on special vectors defined in the following article:
"D. Panario, M. Sall, Q. Wang The complexity of elliptic normal bases".
The preliminary functions help also in the computation of elliptic functions of degree
two through which we define an elliptic normal basis. For more details about these
functions, we refer to the beginning of the file 'ellnbcomplexity.m'.
If an elliptic normal basis exists over the considered extension, the main function 
'ENBparamsComputation' computes its parameters step by step and returns a format that 
contains all these parameters. These parameters are used by the function 
'ENBcomplexityBounds' that computes the weight of the special vectors and returns
an upper bound and a lower bound on the complexity of the basis.
If there is no elliptic normal basis, the package returns the following statement
"Sorry, we can't define an elliptic normal basis over this finite field extension".
Finally in this repository, we provide 9 .JPG files that represent the screenshots of some
tests we have made. For each example of elliptic normal basis construction we have 3 .JPG
file. In these three files, we define the process of construction, the parameters that lead to 
the bounds on the complexity of the constructed basis, and the weight of the principal special
vectors. 
