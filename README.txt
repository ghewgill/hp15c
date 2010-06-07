jsmat
=====
Greg Hewgill
http://hewgill.com

jsmat is a direct port of Jama (http://math.nist.gov/javanumerics/jama/) to
Javascript.

The Matrix class provides general purpose matrix arithmetic functions, plus the
following five decompositions:

    - Cholesky Decomposition of symmetric, positive definite matrices
    - LU Decomposition (Gaussian elimination) of rectangular matrices
    - QR Decomposition of rectangular matrices
    - Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices
    - Singular Value Decomposition of rectangular matrices

Example
-------

The following simple example solves a 3x3 linear system Ax=b and computes the 
norm of the residual.

    var array = [[1.,2.,3],[4.,5.,6.],[7.,8.,10.]]; 
    var A = new Matrix(array); 
    var b = Matrix.random(3,1); 
    var x = A.solve(b); 
    var Residual = A.times(x).minus(b); 
    var rnorm = Residual.normInf();

Tests
-----

The test.js module implements a test suite. Run with Rhino as:

    java -cp js.jar org.mozilla.javascript.tools.shell.Main -debug -w test.js

Authors
-------

Jama was originally designed and implemented by Joe Hicklin, Cleve Moler, and
Peter Webb from The MathWorks and Ronald F. Boisvert, Bruce Miller, Roldan
Pozo, and Karin Remington from NIST.

The translation from Java to Javascript was done by Greg Hewgill.
