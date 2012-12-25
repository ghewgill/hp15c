/**
   Jama = Java Matrix class.
<P>
   The Java Matrix Class provides the fundamental operations of numerical
   linear algebra.  Various constructors create Matrices from two dimensional
   arrays of double precision floating point numbers.  Various "gets" and
   "sets" provide access to submatrices and matrix elements.  Several methods 
   implement basic matrix arithmetic, including matrix addition and
   multiplication, matrix norms, and element-by-element array operations.
   Methods for reading and printing matrices are also included.  All the
   operations in this version of the Matrix Class involve real matrices.
   Complex matrices may be handled in a future version.
<P>
   Five fundamental matrix decompositions, which consist of pairs or triples
   of matrices, permutation vectors, and the like, produce results in five
   decomposition classes.  These decompositions are accessed by the Matrix
   class to compute solutions of simultaneous linear equations, determinants,
   inverses and other matrix functions.  The five decompositions are:
<P><UL>
   <LI>Cholesky Decomposition of symmetric, positive definite matrices.
   <LI>LU Decomposition of rectangular matrices.
   <LI>QR Decomposition of rectangular matrices.
   <LI>Singular Value Decomposition of rectangular matrices.
   <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
</UL>
<DL>
<DT><B>Example of use:</B></DT>
<P>
<DD>Solve a linear system A x = b and compute the residual norm, ||b - A x||.
<P><PRE>
      double[][] vals = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}};
      Matrix A = new Matrix(vals);
      Matrix b = Matrix.random(3,1);
      Matrix x = A.solve(b);
      Matrix r = A.times(x).minus(b);
      double rnorm = r.normInf();
</PRE></DD>
</DL>

@author The MathWorks, Inc. and the National Institute of Standards and Technology.
@version 5 August 1998
*/

function ArrayIndexOutOfBoundsException(index) {
    this.name = "ArrayIndexOutOfBoundsException";
    this.message = "index " + index;
    this.getMessage = function() { return message; }
}

function IllegalArgumentException(message) {
    this.name = "IllegalArgumentException";
    this.message = message;
    this.getMessage = function() { return message; }
}

function hypot(a, b) {
    var r;
    if (Math.abs(a) > Math.abs(b)) {
        r = b/a;
        r = Math.abs(a)*Math.sqrt(1+r*r);
    } else if (b != 0) {
        r = a/b;
        r = Math.abs(b)*Math.sqrt(1+r*r);
    } else {
        r = 0.0;
    }
    return r;
}

function LUDecomposition(A) {

/* ------------------------
   Class variables
 * ------------------------ */

    /** Array for internal storage of decomposition.
    @serial internal array storage.
    */
    this.LU = null;

    /** Row and column dimensions, and pivot sign.
    @serial column dimension.
    @serial row dimension.
    @serial pivot sign.
    */
    this.m = 0;
    this.n = 0;
    this.pivsign = 0;

    /** Internal storage of pivot vector.
    @serial pivot vector.
    */
    this.piv = null;

/* ------------------------
   Constructor
 * ------------------------ */

    /** LU Decomposition
    @param  A   Rectangular matrix
    @return     Structure to access L, U and piv.
    */

    // Use a "left-looking", dot-product, Crout/Doolittle algorithm.

    this.LU = A.getArrayCopy();
    this.m = A.getRowDimension();
    this.n = A.getColumnDimension();
    this.piv = new Array(this.m);
    for (var i = 0; i < this.m; i++) {
        this.piv[i] = i;
    }
    this.pivsign = 1;
    var LUrowi = null;
    var LUcolj = new Array(this.m);

    // Outer loop.

    for (var j = 0; j < this.n; j++) {

        // Make a copy of the j-th column to localize references.

        for (var i = 0; i < this.m; i++) {
            LUcolj[i] = this.LU[i][j];
        }

        // Apply previous transformations.

        for (var i = 0; i < this.m; i++) {
            LUrowi = this.LU[i];

            // Most of the time is spent in the following dot product.

            var kmax = Math.min(i,j);
            var s = 0.0;
            for (var k = 0; k < kmax; k++) {
                s += LUrowi[k]*LUcolj[k];
            }

            LUrowi[j] = LUcolj[i] -= s;
        }

        // Find pivot and exchange if necessary.

        var p = j;
        for (var i = j+1; i < this.m; i++) {
            if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
                p = i;
            }
        }
        if (p != j) {
            for (var k = 0; k < this.n; k++) {
                var t = this.LU[p][k]; this.LU[p][k] = this.LU[j][k]; this.LU[j][k] = t;
            }
            var k = this.piv[p]; this.piv[p] = this.piv[j]; this.piv[j] = k;
            this.pivsign = -this.pivsign;
        }

        // Compute multipliers.
     
        if (j < this.m & this.LU[j][j] != 0.0) {
            for (var i = j+1; i < this.m; i++) {
                this.LU[i][j] /= this.LU[j][j];
            }
        }
    }

/* ------------------------
   Public Methods
 * ------------------------ */

    /** Is the matrix nonsingular?
    @return     true if U, and hence A, is nonsingular.
    */

    this.isNonsingular = function() {
        for (var j = 0; j < this.n; j++) {
            if (this.LU[j][j] == 0)
                return false;
        }
        return true;
    }

    /** Return lower triangular factor
    @return     L
    */

    this.getL = function() {
        var X = new Matrix(this.m,this.n);
        var L = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                if (i > j) {
                    L[i][j] = this.LU[i][j];
                } else if (i == j) {
                    L[i][j] = 1.0;
                } else {
                    L[i][j] = 0.0;
                }
            }
        }
        return X;
    }

    /** Return upper triangular factor
    @return     U
    */

    this.getU = function() {
        var X = new Matrix(this.n,this.n);
        var U = X.getArray();
        for (var i = 0; i < this.n; i++) {
            for (var j = 0; j < this.n; j++) {
                if (i <= j) {
                    U[i][j] = this.LU[i][j];
                } else {
                    U[i][j] = 0.0;
                }
            }
        }
        return X;
    }

    /** Return pivot permutation vector
    @return     piv
    */

    this.getPivot = function() {
        var p = new Array(this.m);
        for (var i = 0; i < this.m; i++) {
            p[i] = this.piv[i];
        }
        return p;
    }
//
//   /** Return pivot permutation vector as a one-dimensional double array
//   @return     (double) piv
//   */
//
//   public double[] getDoublePivot () {
//      double[] vals = new double[m];
//      for (int i = 0; i < m; i++) {
//         vals[i] = (double) piv[i];
//      }
//      return vals;
//   }

    /** Determinant
    @return     det(A)
    @exception  IllegalArgumentException  Matrix must be square
    */

    this.det = function() {
        if (this.m != this.n) {
            throw new IllegalArgumentException("Matrix must be square.");
        }
        var d = this.pivsign;
        for (var j = 0; j < this.n; j++) {
            d *= this.LU[j][j];
        }
        return d;
    }

    /** Solve A*X = B
    @param  B   A Matrix with as many rows as A and any number of columns.
    @return     X so that L*U*X = B(piv,:)
    @exception  IllegalArgumentException Matrix row dimensions must agree.
    @exception  RuntimeException  Matrix is singular.
    */

    this.solve = function(B) {
        if (B.getRowDimension() != this.m) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (!this.isNonsingular()) {
            throw new RuntimeException("Matrix is singular.");
        }

        // Copy right hand side with pivoting
        var nx = B.getColumnDimension();
        var Xmat = B.getMatrix(this.piv,0,nx-1);
        var X = Xmat.getArray();

        // Solve L*Y = B(piv,:)
        for (var k = 0; k < this.n; k++) {
            for (var i = k+1; i < this.n; i++) {
                for (var j = 0; j < nx; j++) {
                    X[i][j] -= X[k][j]*this.LU[i][k];
                }
            }
        }
        // Solve U*X = Y;
        for (var k = this.n-1; k >= 0; k--) {
            for (var j = 0; j < nx; j++) {
                X[k][j] /= this.LU[k][k];
            }
            for (var i = 0; i < k; i++) {
                for (var j = 0; j < nx; j++) {
                    X[i][j] -= X[k][j]*this.LU[i][k];
                }
            }
        }
        return Xmat;
    }
}

/** QR Decomposition.
<P>
    For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
    orthogonal matrix Q and an n-by-n upper triangular matrix R so that
    A = Q*R.
<P>
    The QR decompostion always exists, even if the matrix does not have
    full rank, so the constructor will never fail.  The primary use of the
    QR decomposition is in the least squares solution of nonsquare systems
    of simultaneous linear equations.  This will fail if isFullRank()
    returns false.
*/

function QRDecomposition(A) {

/* ------------------------
   Class variables
 * ------------------------ */

    /** Array for internal storage of decomposition.
    @serial internal array storage.
    */
    this.QR = null;

    /** Row and column dimensions.
    @serial column dimension.
    @serial row dimension.
    */
    this.m = 0;
    this.n = 0;

    /** Array for internal storage of diagonal of R.
    @serial diagonal of R.
    */
    this.Rdiag = 0;

/* ------------------------
   Constructor
 * ------------------------ */

    /** QR Decomposition, computed by Householder reflections.
    @param A    Rectangular matrix
    @return     Structure to access R and the Householder vectors and compute Q.
    */

    // Initialize.
    this.QR = A.getArrayCopy();
    this.m = A.getRowDimension();
    this.n = A.getColumnDimension();
    this.Rdiag = new Array(n);

    // Main loop.
    for (var k = 0; k < this.n; k++) {
        // Compute 2-norm of k-th column without under/overflow.
        var nrm = 0;
        for (var i = k; i < this.m; i++) {
            nrm = hypot(nrm,this.QR[i][k]);
        }

        if (nrm != 0.0) {
            // Form k-th Householder vector.
            if (this.QR[k][k] < 0) {
                nrm = -nrm;
            }
            for (var i = k; i < this.m; i++) {
                this.QR[i][k] /= nrm;
            }
            this.QR[k][k] += 1.0;

            // Apply transformation to remaining columns.
            for (var j = k+1; j < this.n; j++) {
                var s = 0.0; 
                for (var i = k; i < this.m; i++) {
                    s += this.QR[i][k]*this.QR[i][j];
                }
                s = -s/this.QR[k][k];
                for (var i = k; i < this.m; i++) {
                    this.QR[i][j] += s*this.QR[i][k];
                }
            }
        }
        this.Rdiag[k] = -nrm;
    }

/* ------------------------
   Public Methods
 * ------------------------ */

//   /** Is the matrix full rank?
//   @return     true if R, and hence A, has full rank.
//   */
//
//   public boolean isFullRank () {
//      for (int j = 0; j < n; j++) {
//         if (Rdiag[j] == 0)
//            return false;
//      }
//      return true;
//   }
//
//   /** Return the Householder vectors
//   @return     Lower trapezoidal matrix whose columns define the reflections
//   */
//
//   public Matrix getH () {
//      Matrix X = new Matrix(m,n);
//      double[][] H = X.getArray();
//      for (int i = 0; i < m; i++) {
//         for (int j = 0; j < n; j++) {
//            if (i >= j) {
//               H[i][j] = QR[i][j];
//            } else {
//               H[i][j] = 0.0;
//            }
//         }
//      }
//      return X;
//   }

    /** Return the upper triangular factor
    @return     R
    */

    this.getR = function() {
        var X = new Matrix(this.n,this.n);
        var R = X.getArray();
        for (var i = 0; i < this.n; i++) {
            for (var j = 0; j < this.n; j++) {
                if (i < j) {
                    R[i][j] = this.QR[i][j];
                } else if (i == j) {
                    R[i][j] = this.Rdiag[i];
                } else {
                    R[i][j] = 0.0;
                }
            }
        }
        return X;
    }

    /** Generate and return the (economy-sized) orthogonal factor
    @return     Q
    */

    this.getQ = function() {
        var X = new Matrix(this.m,this.n);
        var Q = X.getArray();
        for (var k = this.n-1; k >= 0; k--) {
            for (var i = 0; i < this.m; i++) {
                Q[i][k] = 0.0;
            }
            Q[k][k] = 1.0;
            for (var j = k; j < this.n; j++) {
                if (this.QR[k][k] != 0) {
                    var s = 0.0;
                    for (var i = k; i < this.m; i++) {
                        s += this.QR[i][k]*Q[i][j];
                    }
                    s = -s/this.QR[k][k];
                    for (var i = k; i < this.m; i++) {
                        Q[i][j] += s*this.QR[i][k];
                    }
                }
            }
        }
        return X;
    }
//
//   /** Least squares solution of A*X = B
//   @param B    A Matrix with as many rows as A and any number of columns.
//   @return     X that minimizes the two norm of Q*R*X-B.
//   @exception  IllegalArgumentException  Matrix row dimensions must agree.
//   @exception  RuntimeException  Matrix is rank deficient.
//   */
//
//   public Matrix solve (Matrix B) {
//      if (B.getRowDimension() != m) {
//         throw new IllegalArgumentException("Matrix row dimensions must agree.");
//      }
//      if (!this.isFullRank()) {
//         throw new RuntimeException("Matrix is rank deficient.");
//      }
//      
//      // Copy right hand side
//      int nx = B.getColumnDimension();
//      double[][] X = B.getArrayCopy();
//
//      // Compute Y = transpose(Q)*B
//      for (int k = 0; k < n; k++) {
//         for (int j = 0; j < nx; j++) {
//            double s = 0.0; 
//            for (int i = k; i < m; i++) {
//               s += QR[i][k]*X[i][j];
//            }
//            s = -s/QR[k][k];
//            for (int i = k; i < m; i++) {
//               X[i][j] += s*QR[i][k];
//            }
//         }
//      }
//      // Solve R*X = Y;
//      for (int k = n-1; k >= 0; k--) {
//         for (int j = 0; j < nx; j++) {
//            X[k][j] /= Rdiag[k];
//         }
//         for (int i = 0; i < k; i++) {
//            for (int j = 0; j < nx; j++) {
//               X[i][j] -= X[k][j]*QR[i][k];
//            }
//         }
//      }
//      return (new Matrix(X,n,nx).getMatrix(0,n-1,0,nx-1));
//   }
}

    /** Singular Value Decomposition.
    <P>
    For an m-by-n matrix A with m >= n, the singular value decomposition is
    an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
    an n-by-n orthogonal matrix V so that A = U*S*V'.
    <P>
    The singular values, sigma[k] = S[k][k], are ordered so that
    sigma[0] >= sigma[1] >= ... >= sigma[n-1].
    <P>
    The singular value decompostion always exists, so the constructor will
    never fail.  The matrix condition number and the effective numerical
    rank can be computed from this decomposition.
    */

function SingularValueDecomposition(Arg) {

/* ------------------------
   Class variables
 * ------------------------ */

    /** Arrays for internal storage of U and V.
    @serial internal storage of U.
    @serial internal storage of V.
    */
    this.U = null;
    this.V = null;

    /** Array for internal storage of singular values.
    @serial internal storage of singular values.
    */
    this.s = null;

    /** Row and column dimensions.
    @serial row dimension.
    @serial column dimension.
    */
    this.m = 0;
    this.n = 0;

/* ------------------------
   Constructor
 * ------------------------ */

    /** Construct the singular value decomposition
    @param A    Rectangular matrix
    @return     Structure to access U, S and V.
    */

    // Derived from LINPACK code.
    // Initialize.
    var A = Arg.getArrayCopy();
    this.m = Arg.getRowDimension();
    this.n = Arg.getColumnDimension();

    /* Apparently the failing cases are only a proper subset of (m<n), 
        so let's not throw error.  Correct fix to come later?
    if (this.m<this.n) {
        throw new IllegalArgumentException("Jama SVD only works for m >= n"); }
    */
    var nu = Math.min(this.m,this.n);
    this.s = new Array(Math.min(this.m+1,this.n));
    this.U = new Array(this.m);
    for (var i = 0; i < this.m; i++) {
        this.U[i] = new Array(nu);
        for (var j = 0; j < this.n; j++) {
            this.U[i][j] = 0;
        }
    }
    this.V = new Array(this.n);
    for (var i = 0; i < this.n; i++) {
        this.V[i] = new Array(this.n);
        for (var j = 0; j < this.n; j++) {
            this.V[i][j] = 0;
        }
    }
    var e = new Array(this.n);
    var work = new Array(this.m);
    var wantu = true;
    var wantv = true;

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.

    var nct = Math.min(this.m-1,this.n);
    var nrt = Math.max(0,Math.min(this.n-2,this.m));
    for (var k = 0; k < Math.max(nct,nrt); k++) {
        if (k < nct) {

            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            this.s[k] = 0;
            for (var i = k; i < this.m; i++) {
                this.s[k] = hypot(this.s[k],A[i][k]);
            }
            if (this.s[k] != 0.0) {
                if (A[k][k] < 0.0) {
                    this.s[k] = -this.s[k];
                }
                for (var i = k; i < this.m; i++) {
                    A[i][k] /= this.s[k];
                }
                A[k][k] += 1.0;
            }
            this.s[k] = -this.s[k];
        }
        for (var j = k+1; j < this.n; j++) {
            if ((k < nct) & (this.s[k] != 0.0))  {

                // Apply the transformation.

                var t = 0;
                for (var i = k; i < this.m; i++) {
                    t += A[i][k]*A[i][j];
                }
                t = -t/A[k][k];
                for (var i = k; i < this.m; i++) {
                    A[i][j] += t*A[i][k];
                }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e[j] = A[k][j];
        }
        if (wantu & (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (var i = k; i < this.m; i++) {
                this.U[i][k] = A[i][k];
            }
        }
        if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (var i = k+1; i < this.n; i++) {
                e[k] = hypot(e[k],e[i]);
            }
            if (e[k] != 0.0) {
                if (e[k+1] < 0.0) {
                    e[k] = -e[k];
                }
                for (var i = k+1; i < this.n; i++) {
                    e[i] /= e[k];
                }
                e[k+1] += 1.0;
            }
            e[k] = -e[k];
            if ((k+1 < this.m) & (e[k] != 0.0)) {

                // Apply the transformation.

                for (var i = k+1; i < this.m; i++) {
                    work[i] = 0.0;
                }
                for (var j = k+1; j < this.n; j++) {
                    for (var i = k+1; i < this.m; i++) {
                        work[i] += e[j]*A[i][j];
                    }
                }
                for (var j = k+1; j < this.n; j++) {
                    var t = -e[j]/e[k+1];
                    for (var i = k+1; i < this.m; i++) {
                        A[i][j] += t*work[i];
                    }
                }
            }
            if (wantv) {

                // Place the transformation in V for subsequent
                // back multiplication.

                for (var i = k+1; i < this.n; i++) {
                    this.V[i][k] = e[i];
                }
            }
        }
    }

    // Set up the final bidiagonal matrix or order p.

    var p = Math.min(this.n,this.m+1);
    if (nct < this.n) {
        this.s[nct] = A[nct][nct];
    }
    if (this.m < p) {
        this.s[p-1] = 0.0;
    }
    if (nrt+1 < p) {
        e[nrt] = A[nrt][p-1];
    }
    e[p-1] = 0.0;

    // If required, generate U.

    if (wantu) {
        for (var j = nct; j < nu; j++) {
            for (var i = 0; i < this.m; i++) {
                this.U[i][j] = 0.0;
            }
            this.U[j][j] = 1.0;
        }
        for (var k = nct-1; k >= 0; k--) {
            if (this.s[k] != 0.0) {
                for (var j = k+1; j < nu; j++) {
                    var t = 0;
                    for (var i = k; i < this.m; i++) {
                        t += this.U[i][k]*this.U[i][j];
                    }
                    t = -t/this.U[k][k];
                    for (var i = k; i < this.m; i++) {
                        this.U[i][j] += t*this.U[i][k];
                    }
                }
                for (var i = k; i < this.m; i++ ) {
                    this.U[i][k] = -this.U[i][k];
                }
                this.U[k][k] = 1.0 + this.U[k][k];
                for (var i = 0; i < k-1; i++) {
                    this.U[i][k] = 0.0;
                }
            } else {
                for (var i = 0; i < this.m; i++) {
                    this.U[i][k] = 0.0;
                }
                this.U[k][k] = 1.0;
            }
        }
    }

    // If required, generate V.

    if (wantv) {
        for (var k = this.n-1; k >= 0; k--) {
            if ((k < nrt) & (e[k] != 0.0)) {
                for (var j = k+1; j < nu; j++) {
                    var t = 0;
                    for (var i = k+1; i < this.n; i++) {
                        t += this.V[i][k]*this.V[i][j];
                    }
                    t = -t/this.V[k+1][k];
                    for (var i = k+1; i < this.n; i++) {
                        this.V[i][j] += t*this.V[i][k];
                    }
                }
            }
            for (var i = 0; i < this.n; i++) {
                this.V[i][k] = 0.0;
            }
            this.V[k][k] = 1.0;
        }
    }

    // Main iteration loop for the singular values.

    var pp = p-1;
    var iter = 0;
    var eps = Math.pow(2.0,-52.0);
    var tiny = Math.pow(2.0,-966.0);
    while (p > 0) {
        var k,kase;

        // Here is where a test for too many iterations would go.

        // This section of the program inspects for
        // negligible elements in the s and e arrays.  On
        // completion the variables kase and k are set as follows.

        // kase = 1     if s(p) and e[k-1] are negligible and k<p
        // kase = 2     if s(k) is negligible and k<p
        // kase = 3     if e[k-1] is negligible, k<p, and
        //              s(k), ..., s(p) are not negligible (qr step).
        // kase = 4     if e(p-1) is negligible (convergence).

        for (k = p-2; k >= -1; k--) {
            if (k == -1) {
                break;
            }
            if (Math.abs(e[k]) <=
                    tiny + eps*(Math.abs(this.s[k]) + Math.abs(this.s[k+1]))) {
                e[k] = 0.0;
                break;
            }
        }
        if (k == p-2) {
            kase = 4;
        } else {
            var ks;
            for (ks = p-1; ks >= k; ks--) {
                if (ks == k) {
                    break;
                }
                var t = (ks != p ? Math.abs(e[ks]) : 0.) + 
                        (ks != k+1 ? Math.abs(e[ks-1]) : 0.);
                if (Math.abs(this.s[ks]) <= tiny + eps*t)  {
                    this.s[ks] = 0.0;
                    break;
                }
            }
            if (ks == k) {
                kase = 3;
            } else if (ks == p-1) {
                kase = 1;
            } else {
                kase = 2;
                k = ks;
            }
        }
        k++;

        // Perform the task indicated by kase.

        switch (kase) {

            // Deflate negligible s(p).

            case 1: {
                var f = e[p-2];
                e[p-2] = 0.0;
                for (var j = p-2; j >= k; j--) {
                    var t = hypot(this.s[j],f);
                    var cs = this.s[j]/t;
                    var sn = f/t;
                    this.s[j] = t;
                    if (j != k) {
                        f = -sn*e[j-1];
                        e[j-1] = cs*e[j-1];
                    }
                    if (wantv) {
                        for (var i = 0; i < this.n; i++) {
                            t = cs*this.V[i][j] + sn*this.V[i][p-1];
                            this.V[i][p-1] = -sn*this.V[i][j] + cs*this.V[i][p-1];
                            this.V[i][j] = t;
                        }
                    }
                }
            }
            break;

            // Split at negligible s(k).

            case 2: {
                var f = e[k-1];
                e[k-1] = 0.0;
                for (var j = k; j < p; j++) {
                    var t = hypot(this.s[j],f);
                    var cs = this.s[j]/t;
                    var sn = f/t;
                    this.s[j] = t;
                    f = -sn*e[j];
                    e[j] = cs*e[j];
                    if (wantu) {
                        for (var i = 0; i < this.m; i++) {
                            t = cs*this.U[i][j] + sn*this.U[i][k-1];
                            this.U[i][k-1] = -sn*this.U[i][j] + cs*this.U[i][k-1];
                            this.U[i][j] = t;
                        }
                    }
                }
            }
            break;

            // Perform one qr step.

            case 3: {

                // Calculate the shift.

                var scale = Math.max(Math.max(Math.max(Math.max(
                        Math.abs(this.s[p-1]),Math.abs(this.s[p-2])),Math.abs(e[p-2])), 
                        Math.abs(this.s[k])),Math.abs(e[k]));
                var sp = this.s[p-1]/scale;
                var spm1 = this.s[p-2]/scale;
                var epm1 = e[p-2]/scale;
                var sk = this.s[k]/scale;
                var ek = e[k]/scale;
                var b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
                var c = (sp*epm1)*(sp*epm1);
                var shift = 0.0;
                if ((b != 0.0) | (c != 0.0)) {
                    shift = Math.sqrt(b*b + c);
                    if (b < 0.0) {
                        shift = -shift;
                    }
                    shift = c/(b + shift);
                }
                var f = (sk + sp)*(sk - sp) + shift;
                var g = sk*ek;

                // Chase zeros.

                for (var j = k; j < p-1; j++) {
                    var t = hypot(f,g);
                    var cs = f/t;
                    var sn = g/t;
                    if (j != k) {
                        e[j-1] = t;
                    }
                    f = cs*this.s[j] + sn*e[j];
                    e[j] = cs*e[j] - sn*this.s[j];
                    g = sn*this.s[j+1];
                    this.s[j+1] = cs*this.s[j+1];
                    if (wantv) {
                        for (var i = 0; i < this.n; i++) {
                            t = cs*this.V[i][j] + sn*this.V[i][j+1];
                            this.V[i][j+1] = -sn*this.V[i][j] + cs*this.V[i][j+1];
                            this.V[i][j] = t;
                        }
                    }
                    t = hypot(f,g);
                    cs = f/t;
                    sn = g/t;
                    this.s[j] = t;
                    f = cs*e[j] + sn*this.s[j+1];
                    this.s[j+1] = -sn*e[j] + cs*this.s[j+1];
                    g = sn*e[j+1];
                    e[j+1] = cs*e[j+1];
                    if (wantu && (j < this.m-1)) {
                        for (var i = 0; i < this.m; i++) {
                            t = cs*this.U[i][j] + sn*this.U[i][j+1];
                            this.U[i][j+1] = -sn*this.U[i][j] + cs*this.U[i][j+1];
                            this.U[i][j] = t;
                        }
                    }
                }
                e[p-2] = f;
                iter = iter + 1;
            }
            break;

            // Convergence.

            case 4: {

                // Make the singular values positive.

                if (this.s[k] <= 0.0) {
                    this.s[k] = (this.s[k] < 0.0 ? -this.s[k] : 0.0);
                    if (wantv) {
                        for (var i = 0; i <= pp; i++) {
                            this.V[i][k] = -this.V[i][k];
                        }
                    }
                }

                // Order the singular values.

                while (k < pp) {
                    if (this.s[k] >= this.s[k+1]) {
                        break;
                    }
                    var t = this.s[k];
                    this.s[k] = this.s[k+1];
                    this.s[k+1] = t;
                    if (wantv && (k < this.n-1)) {
                        for (var i = 0; i < this.n; i++) {
                            t = this.V[i][k+1]; this.V[i][k+1] = this.V[i][k]; this.V[i][k] = t;
                        }
                    }
                    if (wantu && (k < this.m-1)) {
                        for (var i = 0; i < this.m; i++) {
                            t = this.U[i][k+1]; this.U[i][k+1] = this.U[i][k]; this.U[i][k] = t;
                        }
                    }
                    k++;
                }
                iter = 0;
                p--;
            }
            break;
        }
    }

/* ------------------------
   Public Methods
 * ------------------------ */

    /** Return the left singular vectors
    @return     U
    */

    this.getU = function() {
        return new Matrix(this.U,this.m,Math.min(this.m+1,this.n));
    }

    /** Return the right singular vectors
    @return     V
    */

    this.getV = function() {
        return new Matrix(this.V,this.n,this.n);
    }

    /** Return the one-dimensional array of singular values
    @return     diagonal of S.
    */

    this.getSingularValues = function() {
        return this.s;
    }

    /** Return the diagonal matrix of singular values
    @return     S
    */

    this.getS = function() {
        var X = new Matrix(this.n,this.n);
        var S = X.getArray();
        for (var i = 0; i < this.n; i++) {
            for (var j = 0; j < this.n; j++) {
                S[i][j] = 0.0;
            }
            S[i][i] = this.s[i];
        }
        return X;
    }
//
//   /** Two norm
//   @return     max(S)
//   */
//
//   public double norm2 () {
//      return s[0];
//   }

    /** Two norm condition number
    @return     max(S)/min(S)
    */

    this.cond = function() {
        return this.s[0]/this.s[Math.min(this.m,this.n)-1];
    }

    /** Effective numerical matrix rank
    @return     Number of nonnegligible singular values.
    */

    this.rank = function() {
        var eps = Math.pow(2.0,-52.0);
        var tol = Math.max(this.m,this.n)*this.s[0]*eps;
        var r = 0;
        for (var i = 0; i < this.s.length; i++) {
            if (this.s[i] > tol) {
                r++;
            }
        }
        return r;
    }
}

    /** Cholesky Decomposition.
    <P>
    For a symmetric, positive definite matrix A, the Cholesky decomposition
    is an lower triangular matrix L so that A = L*L'.
    <P>
    If the matrix is not symmetric or positive definite, the constructor
    returns a partial decomposition and sets an internal flag that may
    be queried by the isSPD() method.
    */

function CholeskyDecomposition(Arg) {

/* ------------------------
   Class variables
 * ------------------------ */

    /** Array for internal storage of decomposition.
    @serial internal array storage.
    */
    this.L = null;

    /** Row and column dimension (square matrix).
    @serial matrix dimension.
    */
    this.n = 0;

    /** Symmetric and positive definite flag.
    @serial is symmetric and positive definite flag.
    */
    this.isspd = false;

/* ------------------------
   Constructor
 * ------------------------ */

    /** Cholesky algorithm for symmetric and positive definite matrix.
    @param  A   Square, symmetric matrix.
    @return     Structure to access L and isspd flag.
    */

    // Initialize.
    var A = Arg.getArray();
    this.n = Arg.getRowDimension();
    this.L = new Array(this.n);
    for (var i = 0; i < this.n; i++) {
        this.L[i] = new Array(this.n);
    }
    this.isspd = (Arg.getColumnDimension() == n);
    // Main loop.
    for (var j = 0; j < this.n; j++) {
        var Lrowj = this.L[j];
        var d = 0.0;
        for (var k = 0; k < j; k++) {
            var Lrowk = this.L[k];
            var s = 0.0;
            for (var i = 0; i < k; i++) {
                s += Lrowk[i]*Lrowj[i];
            }
            Lrowj[k] = s = (A[j][k] - s)/this.L[k][k];
            d = d + s*s;
            this.isspd = this.isspd & (A[k][j] == A[j][k]); 
        }
        d = A[j][j] - d;
        this.isspd = this.isspd & (d > 0.0);
        this.L[j][j] = Math.sqrt(Math.max(d,0.0));
        for (var k = j+1; k < this.n; k++) {
            this.L[j][k] = 0.0;
        }
    }

/* ------------------------
   Temporary, experimental code.
 * ------------------------ *\

   \** Right Triangular Cholesky Decomposition.
   <P>
   For a symmetric, positive definite matrix A, the Right Cholesky
   decomposition is an upper triangular matrix R so that A = R'*R.
   This constructor computes R with the Fortran inspired column oriented
   algorithm used in LINPACK and MATLAB.  In Java, we suspect a row oriented,
   lower triangular decomposition is faster.  We have temporarily included
   this constructor here until timing experiments confirm this suspicion.
   *\

   \** Array for internal storage of right triangular decomposition. **\
   private transient double[][] R;

   \** Cholesky algorithm for symmetric and positive definite matrix.
   @param  A           Square, symmetric matrix.
   @param  rightflag   Actual value ignored.
   @return             Structure to access R and isspd flag.
   *\

   public CholeskyDecomposition (Matrix Arg, int rightflag) {
      // Initialize.
      double[][] A = Arg.getArray();
      n = Arg.getColumnDimension();
      R = new double[n][n];
      isspd = (Arg.getColumnDimension() == n);
      // Main loop.
      for (int j = 0; j < n; j++) {
         double d = 0.0;
         for (int k = 0; k < j; k++) {
            double s = A[k][j];
            for (int i = 0; i < k; i++) {
               s = s - R[i][k]*R[i][j];
            }
            R[k][j] = s = s/R[k][k];
            d = d + s*s;
            isspd = isspd & (A[k][j] == A[j][k]); 
         }
         d = A[j][j] - d;
         isspd = isspd & (d > 0.0);
         R[j][j] = Math.sqrt(Math.max(d,0.0));
         for (int k = j+1; k < n; k++) {
            R[k][j] = 0.0;
         }
      }
   }

   \** Return upper triangular factor.
   @return     R
   *\

   public Matrix getR () {
      return new Matrix(R,n,n);
   }

\* ------------------------
   End of temporary code.
 * ------------------------ */

/* ------------------------
   Public Methods
 * ------------------------ */

//   /** Is the matrix symmetric and positive definite?
//   @return     true if A is symmetric and positive definite.
//   */
//
//   public boolean isSPD () {
//      return isspd;
//   }

    /** Return triangular factor.
    @return     L
    */

    this.getL = function() {
        return new Matrix(this.L,this.n,this.n);
    }

    /** Solve A*X = B
    @param  B   A Matrix with as many rows as A and any number of columns.
    @return     X so that L*L'*X = B
    @exception  IllegalArgumentException  Matrix row dimensions must agree.
    @exception  RuntimeException  Matrix is not symmetric positive definite.
    */

    this.solve = function(B) {
        if (B.getRowDimension() != this.n) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (!this.isspd) {
            throw new RuntimeException("Matrix is not symmetric positive definite.");
        }

        // Copy right hand side.
        var X = B.getArrayCopy();
        var nx = B.getColumnDimension();

        // Solve L*Y = B;
        for (var k = 0; k < this.n; k++) {
            for (var j = 0; j < nx; j++) {
                for (var i = 0; i < k ; i++) {
                    X[k][j] -= X[i][j]*this.L[k][i];
                }
                X[k][j] /= this.L[k][k];
            }
        }

        // Solve L'*X = Y;
        for (var k = this.n-1; k >= 0; k--) {
            for (var j = 0; j < nx; j++) {
                for (var i = k+1; i < this.n ; i++) {
                    X[k][j] -= X[i][j]*this.L[i][k];
                }
                X[k][j] /= this.L[k][k];
            }
        }
      
      
        return new Matrix(X,this.n,nx);
    }
}

/** Eigenvalues and eigenvectors of a real matrix. 
<P>
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal.
    I.e. A = V.times(D.times(V.transpose())) and 
    V.times(V.transpose()) equals the identity matrix.
<P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
    columns of V represent the eigenvectors in the sense that A*V = V*D,
    i.e. A.times(V) equals V.times(D).  The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon V.cond().
**/

function EigenvalueDecomposition(Arg) {

/* ------------------------
   Class variables
 * ------------------------ */

    /** Row and column dimension (square matrix).
    @serial matrix dimension.
    */
    this.n = 0;

    /** Symmetry flag.
    @serial internal symmetry flag.
    */
    this.isymmetric = false;

    /** Arrays for internal storage of eigenvalues.
    @serial internal storage of eigenvalues.
    */
    this.d = null;
    this.e = null;

    /** Array for internal storage of eigenvectors.
    @serial internal storage of eigenvectors.
    */
    this.V = null;

    /** Array for internal storage of nonsymmetric Hessenberg form.
    @serial internal storage of nonsymmetric Hessenberg form.
    */
    this.H = null;

    /** Working storage for nonsymmetric algorithm.
    @serial working storage for nonsymmetric algorithm.
    */
    this.ort = null;

/* ------------------------
   Private Methods
 * ------------------------ */

    // Symmetric Householder reduction to tridiagonal form.

    this.tred2 = function() {

        //  This is derived from the Algol procedures tred2 by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.

        for (var j = 0; j < this.n; j++) {
            this.d[j] = this.V[this.n-1][j];
        }

        // Householder reduction to tridiagonal form.
   
        for (var i = this.n-1; i > 0; i--) {
   
            // Scale to avoid under/overflow.
   
            var scale = 0.0;
            var h = 0.0;
            for (var k = 0; k < i; k++) {
                scale = scale + Math.abs(this.d[k]);
            }
            if (scale == 0.0) {
                this.e[i] = this.d[i-1];
                for (var j = 0; j < i; j++) {
                    this.d[j] = this.V[i-1][j];
                    this.V[i][j] = 0.0;
                    this.V[j][i] = 0.0;
                }
            } else {
   
                // Generate Householder vector.
   
                for (var k = 0; k < i; k++) {
                    this.d[k] /= scale;
                    h += this.d[k] * this.d[k];
                }
                var f = this.d[i-1];
                var g = Math.sqrt(h);
                if (f > 0) {
                    g = -g;
                }
                this.e[i] = scale * g;
                h = h - f * g;
                this.d[i-1] = f - g;
                for (var j = 0; j < i; j++) {
                    this.e[j] = 0.0;
                }
   
                // Apply similarity transformation to remaining columns.
   
                for (var j = 0; j < i; j++) {
                    f = this.d[j];
                    this.V[j][i] = f;
                    g = this.e[j] + this.V[j][j] * f;
                    for (var k = j+1; k <= i-1; k++) {
                        g += this.V[k][j] * this.d[k];
                        this.e[k] += this.V[k][j] * f;
                    }
                    this.e[j] = g;
                }
                f = 0.0;
                for (var j = 0; j < i; j++) {
                    this.e[j] /= h;
                    f += this.e[j] * this.d[j];
                }
                var hh = f / (h + h);
                for (var j = 0; j < i; j++) {
                    this.e[j] -= hh * this.d[j];
                }
                for (var j = 0; j < i; j++) {
                    f = this.d[j];
                    g = this.e[j];
                    for (var k = j; k <= i-1; k++) {
                        this.V[k][j] -= (f * this.e[k] + g * this.d[k]);
                    }
                    this.d[j] = this.V[i-1][j];
                    this.V[i][j] = 0.0;
                }
            }
            this.d[i] = h;
        }
   
        // Accumulate transformations.
   
        for (var i = 0; i < this.n-1; i++) {
            this.V[this.n-1][i] = this.V[i][i];
            this.V[i][i] = 1.0;
            var h = this.d[i+1];
            if (h != 0.0) {
                for (var k = 0; k <= i; k++) {
                    this.d[k] = this.V[k][i+1] / h;
                }
                for (var j = 0; j <= i; j++) {
                    var g = 0.0;
                    for (var k = 0; k <= i; k++) {
                        g += this.V[k][i+1] * this.V[k][j];
                    }
                    for (var k = 0; k <= i; k++) {
                        this.V[k][j] -= g * this.d[k];
                    }
                }
            }
            for (var k = 0; k <= i; k++) {
                this.V[k][i+1] = 0.0;
            }
        }
        for (var j = 0; j < this.n; j++) {
            this.d[j] = this.V[this.n-1][j];
            this.V[this.n-1][j] = 0.0;
        }
        this.V[this.n-1][this.n-1] = 1.0;
        this.e[0] = 0.0;
    } 

    // Symmetric tridiagonal QL algorithm.
   
    this.tql2 = function() {

        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
   
        for (var i = 1; i < this.n; i++) {
            this.e[i-1] = this.e[i];
        }
        this.e[n-1] = 0.0;
   
        var f = 0.0;
        var tst1 = 0.0;
        var eps = Math.pow(2.0,-52.0);
        for (var l = 0; l < this.n; l++) {

            // Find small subdiagonal element
   
            tst1 = Math.max(tst1,Math.abs(this.d[l]) + Math.abs(this.e[l]));
            var m = l;
            while (m < this.n) {
                if (Math.abs(this.e[m]) <= eps*tst1) {
                    break;
                }
                m++;
            }
   
            // If m == l, d[l] is an eigenvalue,
            // otherwise, iterate.
   
            if (m > l) {
                var iter = 0;
                do {
                    iter = iter + 1;  // (Could check iteration count here.)
   
                    // Compute implicit shift
   
                    var g = this.d[l];
                    var p = (this.d[l+1] - g) / (2.0 * this.e[l]);
                    var r = hypot(p,1.0);
                    if (p < 0) {
                        r = -r;
                    }
                    this.d[l] = this.e[l] / (p + r);
                    this.d[l+1] = this.e[l] * (p + r);
                    var dl1 = this.d[l+1];
                    var h = g - this.d[l];
                    for (var i = l+2; i < this.n; i++) {
                        this.d[i] -= h;
                    }
                    f = f + h;
   
                    // Implicit QL transformation.
   
                    p = this.d[m];
                    var c = 1.0;
                    var c2 = c;
                    var c3 = c;
                    var el1 = this.e[l+1];
                    var s = 0.0;
                    var s2 = 0.0;
                    for (var i = m-1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * this.e[i];
                        h = c * p;
                        r = hypot(p,this.e[i]);
                        this.e[i+1] = s * r;
                        s = this.e[i] / r;
                        c = p / r;
                        p = c * this.d[i] - s * g;
                        this.d[i+1] = h + s * (c * g + s * this.d[i]);
   
                        // Accumulate transformation.
   
                        for (var k = 0; k < this.n; k++) {
                            h = this.V[k][i+1];
                            this.V[k][i+1] = s * this.V[k][i] + c * h;
                            this.V[k][i] = c * this.V[k][i] - s * h;
                        }
                    }
                    p = -s * s2 * c3 * el1 * this.e[l] / dl1;
                    this.e[l] = s * p;
                    this.d[l] = c * p;
   
                    // Check for convergence.
   
                } while (Math.abs(this.e[l]) > eps*tst1);
            }
            this.d[l] = this.d[l] + f;
            this.e[l] = 0.0;
        }
     
        // Sort eigenvalues and corresponding vectors.
   
        for (var i = 0; i < this.n-1; i++) {
            var k = i;
            var p = this.d[i];
            for (var j = i+1; j < this.n; j++) {
                if (this.d[j] < p) {
                    k = j;
                    p = this.d[j];
                }
            }
            if (k != i) {
                this.d[k] = this.d[i];
                this.d[i] = p;
                for (var j = 0; j < this.n; j++) {
                    p = this.V[j][i];
                    this.V[j][i] = this.V[j][k];
                    this.V[j][k] = p;
                }
            }
        }
    }

    // Nonsymmetric reduction to Hessenberg form.

    this.orthes = function() {
   
        //  This is derived from the Algol procedures orthes and ortran,
        //  by Martin and Wilkinson, Handbook for Auto. Comp.,
        //  Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutines in EISPACK.
   
        var low = 0;
        var high = this.n-1;
   
        for (var m = low+1; m <= high-1; m++) {
   
            // Scale column.
   
            var scale = 0.0;
            for (var i = m; i <= high; i++) {
                scale = scale + Math.abs(this.H[i][m-1]);
            }
            if (scale != 0.0) {
   
                // Compute Householder transformation.
   
                var h = 0.0;
                for (var i = high; i >= m; i--) {
                    this.ort[i] = this.H[i][m-1]/scale;
                    h += this.ort[i] * this.ort[i];
                }
                var g = Math.sqrt(h);
                if (this.ort[m] > 0) {
                    g = -g;
                }
                h = h - this.ort[m] * g;
                this.ort[m] = this.ort[m] - g;
   
                // Apply Householder similarity transformation
                // H = (I-u*u'/h)*H*(I-u*u')/h)
   
                for (var j = m; j < this.n; j++) {
                    var f = 0.0;
                    for (var i = high; i >= m; i--) {
                        f += this.ort[i]*this.H[i][j];
                    }
                    f = f/h;
                    for (var i = m; i <= high; i++) {
                        this.H[i][j] -= f*this.ort[i];
                    }
                }
   
                for (var i = 0; i <= high; i++) {
                    var f = 0.0;
                    for (var j = high; j >= m; j--) {
                        f += this.ort[j]*this.H[i][j];
                    }
                    f = f/h;
                    for (var j = m; j <= high; j++) {
                        this.H[i][j] -= f*this.ort[j];
                    }
                }
                this.ort[m] = scale*this.ort[m];
                this.H[m][m-1] = scale*g;
            }
        }
   
        // Accumulate transformations (Algol's ortran).

        for (var i = 0; i < this.n; i++) {
            for (var j = 0; j < this.n; j++) {
                this.V[i][j] = (i == j ? 1.0 : 0.0);
            }
        }

        for (var m = high-1; m >= low+1; m--) {
            if (this.H[m][m-1] != 0.0) {
                for (var i = m+1; i <= high; i++) {
                    this.ort[i] = this.H[i][m-1];
                }
                for (var j = m; j <= high; j++) {
                    var g = 0.0;
                    for (var i = m; i <= high; i++) {
                        g += this.ort[i] * this.V[i][j];
                    }
                    // Double division avoids possible underflow
                    g = (g / this.ort[m]) / this.H[m][m-1];
                    for (var i = m; i <= high; i++) {
                        this.V[i][j] += g * this.ort[i];
                    }
                }
            }
        }
    }


    // Complex scalar division.

    this.cdiv = function(xr, xi, yr, yi) {
        var r,d;
        if (Math.abs(yr) > Math.abs(yi)) {
            r = yi/yr;
            d = yr + r*yi;
            return [(xr + r*xi)/d, (xi - r*xr)/d];
        } else {
            r = yr/yi;
            d = yi + r*yr;
            return [(r*xr + xi)/d, (r*xi - xr)/d];
        }
    }


    // Nonsymmetric reduction from Hessenberg to real Schur form.

    this.hqr2 = function() {
   
        //  This is derived from the Algol procedure hqr2,
        //  by Martin and Wilkinson, Handbook for Auto. Comp.,
        //  Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
   
        // Initialize
   
        var nn = this.n;
        var n = nn-1;
        var low = 0;
        var high = nn-1;
        var eps = Math.pow(2.0,-52.0);
        var exshift = 0.0;
        var p=0,q=0,r=0,s=0,z=0,t,w,x,y;
   
        // Store roots isolated by balanc and compute matrix norm
   
        var norm = 0.0;
        for (var i = 0; i < nn; i++) {
            if (i < low | i > high) {
                this.d[i] = this.H[i][i];
                this.e[i] = 0.0;
            }
            for (var j = Math.max(i-1,0); j < nn; j++) {
                norm = norm + Math.abs(this.H[i][j]);
            }
        }
   
        // Outer loop over eigenvalue index
   
        var iter = 0;
        while (n >= low) {
   
            // Look for single small sub-diagonal element
   
            var l = n;
            while (l > low) {
                s = Math.abs(this.H[l-1][l-1]) + Math.abs(this.H[l][l]);
                if (s == 0.0) {
                    s = norm;
                }
                if (Math.abs(this.H[l][l-1]) < eps * s) {
                    break;
                }
                l--;
            }
       
            // Check for convergence
            // One root found
   
            if (l == n) {
                this.H[n][n] = this.H[n][n] + exshift;
                this.d[n] = this.H[n][n];
                this.e[n] = 0.0;
                n--;
                iter = 0;
   
            // Two roots found
   
            } else if (l == n-1) {
                w = this.H[n][n-1] * this.H[n-1][n];
                p = (this.H[n-1][n-1] - this.H[n][n]) / 2.0;
                q = p * p + w;
                z = Math.sqrt(Math.abs(q));
                this.H[n][n] = this.H[n][n] + exshift;
                this.H[n-1][n-1] = this.H[n-1][n-1] + exshift;
                x = this.H[n][n];
   
                // Real pair
   
                if (q >= 0) {
                    if (p >= 0) {
                        z = p + z;
                    } else {
                        z = p - z;
                    }
                    this.d[n-1] = x + z;
                    this.d[n] = this.d[n-1];
                    if (z != 0.0) {
                        this.d[n] = x - w / z;
                    }
                    this.e[n-1] = 0.0;
                    this.e[n] = 0.0;
                    x = this.H[n][n-1];
                    s = Math.abs(x) + Math.abs(z);
                    p = x / s;
                    q = z / s;
                    r = Math.sqrt(p * p+q * q);
                    p = p / r;
                    q = q / r;
   
                    // Row modification
   
                    for (var j = n-1; j < nn; j++) {
                        z = this.H[n-1][j];
                        this.H[n-1][j] = q * z + p * this.H[n][j];
                        this.H[n][j] = q * this.H[n][j] - p * z;
                    }
   
                    // Column modification
   
                    for (var i = 0; i <= n; i++) {
                        z = this.H[i][n-1];
                        this.H[i][n-1] = q * z + p * this.H[i][n];
                        this.H[i][n] = q * this.H[i][n] - p * z;
                    }
   
                    // Accumulate transformations
   
                    for (var i = low; i <= high; i++) {
                        z = this.V[i][n-1];
                        this.V[i][n-1] = q * z + p * this.V[i][n];
                        this.V[i][n] = q * this.V[i][n] - p * z;
                    }
   
                // Complex pair
   
                } else {
                    this.d[n-1] = x + p;
                    this.d[n] = x + p;
                    this.e[n-1] = z;
                    this.e[n] = -z;
                }
                n = n - 2;
                iter = 0;
   
            // No convergence yet
   
            } else {
   
                // Form shift
   
                x = this.H[n][n];
                y = 0.0;
                w = 0.0;
                if (l < n) {
                    y = this.H[n-1][n-1];
                    w = this.H[n][n-1] * this.H[n-1][n];
                }
   
                // Wilkinson's original ad hoc shift
   
                if (iter == 10) {
                    exshift += x;
                    for (var i = low; i <= n; i++) {
                        this.H[i][i] -= x;
                    }
                    s = Math.abs(this.H[n][n-1]) + Math.abs(this.H[n-1][n-2]);
                    x = y = 0.75 * s;
                    w = -0.4375 * s * s;
                }

                // MATLAB's new ad hoc shift

                if (iter == 30) {
                    s = (y - x) / 2.0;
                    s = s * s + w;
                    if (s > 0) {
                        s = Math.sqrt(s);
                        if (y < x) {
                            s = -s;
                        }
                        s = x - w / ((y - x) / 2.0 + s);
                        for (var i = low; i <= n; i++) {
                            this.H[i][i] -= s;
                        }
                        exshift += s;
                        x = y = w = 0.964;
                    }
                }
   
                iter = iter + 1;   // (Could check iteration count here.)
   
                // Look for two consecutive small sub-diagonal elements
   
                var m = n-2;
                while (m >= l) {
                    z = this.H[m][m];
                    r = x - z;
                    s = y - z;
                    p = (r * s - w) / this.H[m+1][m] + this.H[m][m+1];
                    q = this.H[m+1][m+1] - z - r - s;
                    r = this.H[m+2][m+1];
                    s = Math.abs(p) + Math.abs(q) + Math.abs(r);
                    p = p / s;
                    q = q / s;
                    r = r / s;
                    if (m == l) {
                        break;
                    }
                    if (Math.abs(this.H[m][m-1]) * (Math.abs(q) + Math.abs(r)) <
                        eps * (Math.abs(p) * (Math.abs(this.H[m-1][m-1]) + Math.abs(z) +
                        Math.abs(this.H[m+1][m+1])))) {
                            break;
                    }
                    m--;
                }
   
                for (var i = m+2; i <= n; i++) {
                    this.H[i][i-2] = 0.0;
                    if (i > m+2) {
                        this.H[i][i-3] = 0.0;
                    }
                }
   
                // Double QR step involving rows l:n and columns m:n
   
                for (var k = m; k <= n-1; k++) {
                    var notlast = (k != n-1);
                    if (k != m) {
                        p = this.H[k][k-1];
                        q = this.H[k+1][k-1];
                        r = (notlast ? this.H[k+2][k-1] : 0.0);
                        x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                        if (x != 0.0) {
                            p = p / x;
                            q = q / x;
                            r = r / x;
                        }
                    }
                    if (x == 0.0) {
                        break;
                    }
                    s = Math.sqrt(p * p + q * q + r * r);
                    if (p < 0) {
                        s = -s;
                    }
                    if (s != 0) {
                        if (k != m) {
                            this.H[k][k-1] = -s * x;
                        } else if (l != m) {
                            this.H[k][k-1] = -this.H[k][k-1];
                        }
                        p = p + s;
                        x = p / s;
                        y = q / s;
                        z = r / s;
                        q = q / p;
                        r = r / p;
   
                        // Row modification
   
                        for (var j = k; j < nn; j++) {
                            p = this.H[k][j] + q * this.H[k+1][j];
                            if (notlast) {
                                p = p + r * this.H[k+2][j];
                                this.H[k+2][j] = this.H[k+2][j] - p * z;
                            }
                            this.H[k][j] = this.H[k][j] - p * x;
                            this.H[k+1][j] = this.H[k+1][j] - p * y;
                        }
   
                        // Column modification
   
                        for (var i = 0; i <= Math.min(n,k+3); i++) {
                            p = x * this.H[i][k] + y * this.H[i][k+1];
                            if (notlast) {
                                p = p + z * this.H[i][k+2];
                                this.H[i][k+2] = this.H[i][k+2] - p * r;
                            }
                            this.H[i][k] = this.H[i][k] - p;
                            this.H[i][k+1] = this.H[i][k+1] - p * q;
                        }
   
                        // Accumulate transformations
   
                        for (var i = low; i <= high; i++) {
                            p = x * this.V[i][k] + y * this.V[i][k+1];
                            if (notlast) {
                                p = p + z * this.V[i][k+2];
                                this.V[i][k+2] = this.V[i][k+2] - p * r;
                            }
                            this.V[i][k] = this.V[i][k] - p;
                            this.V[i][k+1] = this.V[i][k+1] - p * q;
                        }
                    }  // (s != 0)
                }  // k loop
            }  // check convergence
        }  // while (n >= low)
      
        // Backsubstitute to find vectors of upper triangular form

        if (norm == 0.0) {
            return;
        }
   
        for (n = nn-1; n >= 0; n--) {
            p = this.d[n];
            q = this.e[n];
   
            // Real vector
   
            if (q == 0) {
                var l = n;
                this.H[n][n] = 1.0;
                for (var i = n-1; i >= 0; i--) {
                    w = this.H[i][i] - p;
                    r = 0.0;
                    for (var j = l; j <= n; j++) {
                        r = r + this.H[i][j] * this.H[j][n];
                    }
                    if (this.e[i] < 0.0) {
                        z = w;
                        s = r;
                    } else {
                        l = i;
                        if (this.e[i] == 0.0) {
                            if (w != 0.0) {
                                this.H[i][n] = -r / w;
                            } else {
                                this.H[i][n] = -r / (eps * norm);
                            }
   
                        // Solve real equations
   
                        } else {
                            x = this.H[i][i+1];
                            y = this.H[i+1][i];
                            q = (this.d[i] - p) * (this.d[i] - p) + this.e[i] * this.e[i];
                            t = (x * s - z * r) / q;
                            this.H[i][n] = t;
                            if (Math.abs(x) > Math.abs(z)) {
                                this.H[i+1][n] = (-r - w * t) / x;
                            } else {
                                this.H[i+1][n] = (-s - y * t) / z;
                            }
                        }
   
                        // Overflow control
   
                        t = Math.abs(this.H[i][n]);
                        if ((eps * t) * t > 1) {
                            for (var j = i; j <= n; j++) {
                                this.H[j][n] = this.H[j][n] / t;
                            }
                        }
                    }
                }
   
            // Complex vector
   
            } else if (q < 0) {
                var l = n-1;

                // Last vector component imaginary so matrix is triangular
   
                if (Math.abs(this.H[n][n-1]) > Math.abs(this.H[n-1][n])) {
                    this.H[n-1][n-1] = q / this.H[n][n-1];
                    this.H[n-1][n] = -(this.H[n][n] - p) / this.H[n][n-1];
                } else {
                    var zz = this.cdiv(0.0,-this.H[n-1][n],this.H[n-1][n-1]-p,q);
                    this.H[n-1][n-1] = zz[0];
                    this.H[n-1][n] = zz[1];
                }
                this.H[n][n-1] = 0.0;
                this.H[n][n] = 1.0;
                for (var i = n-2; i >= 0; i--) {
                    var ra,sa,vr,vi;
                    ra = 0.0;
                    sa = 0.0;
                    for (var j = l; j <= n; j++) {
                        ra = ra + this.H[i][j] * this.H[j][n-1];
                        sa = sa + this.H[i][j] * this.H[j][n];
                    }
                    w = this.H[i][i] - p;
   
                    if (this.e[i] < 0.0) {
                        z = w;
                        r = ra;
                        s = sa;
                    } else {
                        l = i;
                        if (this.e[i] == 0) {
                            var zz = this.cdiv(-ra,-sa,w,q);
                            this.H[i][n-1] = zz[0];
                            this.H[i][n] = zz[1];
                        } else {
   
                            // Solve complex equations
   
                            x = this.H[i][i+1];
                            y = this.H[i+1][i];
                            vr = (this.d[i] - p) * (this.d[i] - p) + this.e[i] * this.e[i] - q * q;
                            vi = (this.d[i] - p) * 2.0 * q;
                            if (vr == 0.0 & vi == 0.0) {
                                vr = eps * norm * (Math.abs(w) + Math.abs(q) +
                                    Math.abs(x) + Math.abs(y) + Math.abs(z));
                            }
                            var zz = this.cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
                            this.H[i][n-1] = zz[0];
                            this.H[i][n] = zz[1];
                            if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
                                this.H[i+1][n-1] = (-ra - w * this.H[i][n-1] + q * this.H[i][n]) / x;
                                this.H[i+1][n] = (-sa - w * this.H[i][n] - q * this.H[i][n-1]) / x;
                            } else {
                                var zz = this.cdiv(-r-y*this.H[i][n-1],-s-y*this.H[i][n],z,q);
                                this.H[i+1][n-1] = zz[0];
                                this.H[i+1][n] = zz[1];
                            }
                        }
   
                        // Overflow control

                        t = Math.max(Math.abs(this.H[i][n-1]),Math.abs(this.H[i][n]));
                        if ((eps * t) * t > 1) {
                            for (var j = i; j <= n; j++) {
                                this.H[j][n-1] = this.H[j][n-1] / t;
                                this.H[j][n] = this.H[j][n] / t;
                            }
                        }
                    }
                }
            }
        }
   
        // Vectors of isolated roots
   
        for (var i = 0; i < nn; i++) {
            if (i < low | i > high) {
                for (var j = i; j < nn; j++) {
                    this.V[i][j] = this.H[i][j];
                }
            }
        }
   
        // Back transformation to get eigenvectors of original matrix
   
        for (var j = nn-1; j >= low; j--) {
            for (var i = low; i <= high; i++) {
                z = 0.0;
                for (var k = low; k <= Math.min(j,high); k++) {
                    z = z + this.V[i][k] * this.H[k][j];
                }
                this.V[i][j] = z;
            }
        }
    }


/* ------------------------
   Constructor
 * ------------------------ */

    /** Check for symmetry, then construct the eigenvalue decomposition
    @param A    Square matrix
    @return     Structure to access D and V.
    */

    var A = Arg.getArray();
    this.n = Arg.getColumnDimension();
    this.V = new Array(this.n);
    for (var i = 0; i < this.n; i++) {
        this.V[i] = new Array(this.n);
    }
    this.d = new Array(this.n);
    this.e = new Array(this.n);

    this.issymmetric = true;
    for (var j = 0; (j < this.n) & this.issymmetric; j++) {
        for (var i = 0; (i < this.n) & this.issymmetric; i++) {
            this.issymmetric = (A[i][j] == A[j][i]);
        }
    }

    if (this.issymmetric) {
        for (var i = 0; i < this.n; i++) {
            for (var j = 0; j < this.n; j++) {
                this.V[i][j] = A[i][j];
            }
        }

        // Tridiagonalize.
        this.tred2();

        // Diagonalize.
        this.tql2();

    } else {
        this.H = new Array(this.n);
        for (var i = 0; i < this.n; i++) {
            this.H[i] = new Array(this.n);
        }
        this.ort = new Array(n);
     
        for (var j = 0; j < this.n; j++) {
            for (var i = 0; i < this.n; i++) {
                this.H[i][j] = A[i][j];
            }
        }

        // Reduce to Hessenberg form.
        this.orthes();

        // Reduce Hessenberg to real Schur form.
        this.hqr2();
    }

/* ------------------------
   Public Methods
 * ------------------------ */

    /** Return the eigenvector matrix
    @return     V
    */

    this.getV = function() {
        return new Matrix(this.V,this.n,this.n);
    }
//
//   /** Return the real parts of the eigenvalues
//   @return     real(diag(D))
//   */
//
//   public double[] getRealEigenvalues () {
//      return d;
//   }
//
//   /** Return the imaginary parts of the eigenvalues
//   @return     imag(diag(D))
//   */
//
//   public double[] getImagEigenvalues () {
//      return e;
//   }

    /** Return the block diagonal eigenvalue matrix
    @return     D
    */

    this.getD = function() {
        var X = new Matrix(this.n,this.n);
        var D = X.getArray();
        for (var i = 0; i < this.n; i++) {
            for (var j = 0; j < this.n; j++) {
                D[i][j] = 0.0;
            }
            D[i][i] = this.d[i];
            if (this.e[i] > 0) {
                D[i][i+1] = this.e[i];
            } else if (this.e[i] < 0) {
                D[i][i-1] = this.e[i];
            }
        }
        return X;
    }
}


function Matrix() {

/* ------------------------
   Class variables
 * ------------------------ */

    /** Array for internal storage of elements.
    @serial internal array storage.
    */
    this.A = null;

    /** Row and column dimensions.
    @serial row dimension.
    @serial column dimension.
    */
    this.m = 0;
    this.n = 0;

    if (typeof(arguments[0]) === "number" && typeof(arguments[1]) === "number") {

        /** Construct an m-by-n matrix of zeros. 
        @param m    Number of rows.
        @param n    Number of colums.
        */

        var m = arguments[0];
        var n = arguments[1];
        this.m = m;
        this.n = n;
        this.A = Array(m);
        for (var i = 0; i < m; i++) {
            this.A[i] = Array(n);
            for (var j = 0; j < n; j++) {
                this.A[i][j] = 0.;
            }
        }

        /** Construct an m-by-n constant matrix.
        @param m    Number of rows.
        @param n    Number of colums.
        @param s    Fill the matrix with this scalar value.
        */

        if (typeof(arguments[2]) === "number") {
            var s = arguments[2];
            for (var i = 0; i < m; i++) {
                for (var j = 0; j < n; j++) {
                    this.A[i][j] = s;
                }
            }
        }

    } else if (typeof(arguments[0]) === "object" && arguments.length == 1) {

        /** Construct a matrix from a 2-D array.
        @param A    Two-dimensional array of doubles.
        @exception  IllegalArgumentException All rows must have the same length
        @see        #constructWithCopy
        */

        var A = arguments[0];
        this.m = A.length;
        this.n = A[0].length;
        for (var i = 0; i < this.m; i++) {
            if (A[i].length != this.n) {
                throw new IllegalArgumentException("All rows must have the same length.");
            }
        }
        this.A = A;

    } else if (typeof(arguments[0]) === "object" && typeof(arguments[1]) === "number" && typeof(arguments[2]) === "number") {

        /** Construct a matrix quickly without checking arguments.
        @param A    Two-dimensional array of doubles.
        @param m    Number of rows.
        @param n    Number of colums.
        */

        this.A = arguments[0];
        this.m = arguments[1];
        this.n = arguments[2];

    } else if (typeof(arguments[0]) === "object" && typeof(arguments[1]) === "number") {

        /** Construct a matrix from a one-dimensional packed array
        @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
        @param m    Number of rows.
        @exception  IllegalArgumentException Array length must be a multiple of m.
        */

        var vals = arguments[0];
        var m = arguments[1];
        this.m = m;
        this.n = (m != 0 ? Math.floor(vals.length/m) : 0);
        if (m*this.n != vals.length) {
            throw new IllegalArgumentException("Array length must be a multiple of m.");
        }
        this.A = Array(m);
        for (var i = 0; i < m; i++) {
            this.A[i] = Array(this.n);
        }
        for (var i = 0; i < m; i++) {
            for (var j = 0; j < this.n; j++) {
                this.A[i][j] = vals[i+j*m];
            }
        }

    } else {

        throw new IllegalArgumentException("Invalid constructor parameters");

    }

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Construct a matrix from a copy of a 2-D array.
   @param A    Two-dimensional array of doubles.
   @exception  IllegalArgumentException All rows must have the same length
   */

   Matrix.constructWithCopy = function(A) {
      var m = A.length;
      var n = A[0].length;
      var X = new Matrix(m,n);
      var C = X.getArray();
      for (var i = 0; i < m; i++) {
         if (A[i].length != n) {
            throw new IllegalArgumentException
               ("All rows must have the same length.");
         }
         for (var j = 0; j < n; j++) {
            C[i][j] = A[i][j];
         }
      }
      return X;
   }

   /** Make a deep copy of a matrix
   */

   this.copy = function() {
      var X = new Matrix(this.m,this.n);
      var C = X.getArray();
      for (var i = 0; i < this.m; i++) {
         for (var j = 0; j < this.n; j++) {
            C[i][j] = this.A[i][j];
         }
      }
      return X;
   }

//   /** Clone the Matrix object.
//   */
//
//   public Object clone () {
//      return this.copy();
//   }

    /** Access the internal two-dimensional array.
    @return     Pointer to the two-dimensional array of matrix elements.
    */

    this.getArray = function() {
        return this.A;
    }

    /** Copy the internal two-dimensional array.
    @return     Two-dimensional array copy of matrix elements.
    */

    this.getArrayCopy = function() {
        var C = new Array(this.m);
        for (var i = 0; i < this.m; i++) {
            C[i] = new Array(this.n);
            for (var j = 0; j < this.n; j++) {
                C[i][j] = this.A[i][j];
            }
        }
        return C;
    }

    /** Make a one-dimensional column packed copy of the internal array.
    @return     Matrix elements packed in a one-dimensional array by columns.
    */

    this.getColumnPackedCopy = function() {
        var vals = new Array(this.m*this.n);
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                vals[i+j*m] = this.A[i][j];
            }
        }
        return vals;
    }

    /** Make a one-dimensional row packed copy of the internal array.
    @return     Matrix elements packed in a one-dimensional array by rows.
    */

    this.getRowPackedCopy = function() {
        var vals = new Array(this.m*this.n);
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                vals[i*n+j] = this.A[i][j];
            }
        }
        return vals;
    }

    /** Get row dimension.
    @return     m, the number of rows.
    */

    this.getRowDimension = function() {
        return this.m;
    }

    /** Get column dimension.
    @return     n, the number of columns.
    */

    this.getColumnDimension = function() {
        return this.n;
    }

    /** Get a single element.
    @param i    Row index.
    @param j    Column index.
    @return     A(i,j)
    @exception  ArrayIndexOutOfBoundsException
    */

    this.get = function(i, j) {
        if (i < 0 || i >= this.A.length) {
            throw new ArrayIndexOutOfBoundsException(i);
        }
        if (j < 0 || j >= this.A[i].length) {
            throw new ArrayIndexOutOfBoundsException(j);
        }
        return this.A[i][j];
    }

    this.getMatrix = function() {

        /** Get a submatrix.
        @param i0   Initial row index
        @param i1   Final row index
        @param j0   Initial column index
        @param j1   Final column index
        @return     A(i0:i1,j0:j1)
        @exception  ArrayIndexOutOfBoundsException Submatrix indices
        */

        if (typeof(arguments[0]) === "number" && typeof(arguments[1]) === "number" && typeof(arguments[2]) === "number" && typeof(arguments[3]) === "number") {

            var i0 = arguments[0];
            var i1 = arguments[1];
            var j0 = arguments[2];
            var j1 = arguments[3];

            var X = new Matrix(i1-i0+1,j1-j0+1);
            var B = X.getArray();
            for (var i = i0; i <= i1; i++) {
                for (var j = j0; j <= j1; j++) {
                    if (i < 0 || i >= this.m || j < 0 || j >= this.n) {
                        throw new ArrayIndexOutOfBoundsException("Submatrix indices");
                    }
                    B[i-i0][j-j0] = this.A[i][j];
                }
            }
            return X;

        /** Get a submatrix.
        @param r    Array of row indices.
        @param c    Array of column indices.
        @return     A(r(:),c(:))
        @exception  ArrayIndexOutOfBoundsException Submatrix indices
        */

        } else if (typeof(arguments[0]) === "object" && typeof(arguments[1]) === "object") {

            var r = arguments[0];
            var c = arguments[1];

            var X = new Matrix(r.length,c.length);
            var B = X.getArray();
            for (var i = 0; i < r.length; i++) {
                for (var j = 0; j < c.length; j++) {
                    if (r[i] < 0 || r[i] >= this.m || c[j] < 0 || c[j] >= this.n) {
                        throw new ArrayIndexOutOfBoundsException("Submatrix indices");
                    }
                    B[i][j] = this.A[r[i]][c[j]];
                }
            }
            return X;

        /** Get a submatrix.
        @param i0   Initial row index
        @param i1   Final row index
        @param c    Array of column indices.
        @return     A(i0:i1,c(:))
        @exception  ArrayIndexOutOfBoundsException Submatrix indices
        */

        } else if (typeof(arguments[0]) === "number" && typeof(arguments[1]) === "number" && typeof(arguments[2]) === "object") {

            var i0 = arguments[0];
            var i1 = arguments[1];
            var c = arguments[2];

            var X = new Matrix(i1-i0+1,c.length);
            var B = X.getArray();
            for (var i = i0; i <= i1; i++) {
                for (var j = 0; j < c.length; j++) {
                    if (i < 0 || i >= this.m || c[j] < 0 || c[j] >= this.n) {
                        throw new ArrayIndexOutOfBoundsException("Submatrix indices");
                    }
                    B[i-i0][j] = this.A[i][c[j]];
                }
            }
            return X;

        /** Get a submatrix.
        @param r    Array of row indices.
        @param i0   Initial column index
        @param i1   Final column index
        @return     A(r(:),j0:j1)
        @exception  ArrayIndexOutOfBoundsException Submatrix indices
        */

        } else if (typeof(arguments[0]) === "object" && typeof(arguments[1]) === "number" && typeof(arguments[2]) === "number") {

            var r = arguments[0];
            var j0 = arguments[1];
            var j1 = arguments[2];

            var X = new Matrix(r.length,j1-j0+1);
            var B = X.getArray();
            for (var i = 0; i < r.length; i++) {
                for (var j = j0; j <= j1; j++) {
                    if (r[i] < 0 || r[i] >= this.m || j < 0 || j >= this.n) {
                        throw new ArrayIndexOutOfBoundsException("Submatrix indices");
                    }
                    B[i][j-j0] = this.A[r[i]][j];
                }
            }
            return X;

        } else {

            throw IllegalArgumentException("Unexpected constructor parameters");

        }

    }

    /** Set a single element.
    @param i    Row index.
    @param j    Column index.
    @param s    A(i,j).
    @exception  ArrayIndexOutOfBoundsException
    */

    this.set = function(i, j, s) {
        if (i < 0 || i >= this.A.length) {
            throw new ArrayIndexOutOfBoundsException(i);
        }
        if (j < 0 || j >= this.A[i].length) {
            throw new ArrayIndexOutOfBoundsException(j);
        }
        this.A[i][j] = s;
    }

    this.setMatrix = function() {

        /** Set a submatrix.
        @param i0   Initial row index
        @param i1   Final row index
        @param j0   Initial column index
        @param j1   Final column index
        @param X    A(i0:i1,j0:j1)
        @exception  ArrayIndexOutOfBoundsException Submatrix indices
        */

        if (typeof(arguments[0]) === "number" && typeof(arguments[1]) === "number" && typeof(arguments[2]) === "number" && typeof(arguments[3]) === "number" && typeof(arguments[4]) === "object") {

            var i0 = arguments[0];
            var i1 = arguments[1];
            var j0 = arguments[2];
            var j1 = arguments[3];
            var X = arguments[4];

            for (var i = i0; i <= i1; i++) {
                for (var j = j0; j <= j1; j++) {
                    if (i < 0 || i >= this.m || j < 0 || j >= this.n) {
                        throw new ArrayIndexOutOfBoundsException("Submatrix indices");
                    }
                    this.A[i][j] = X.get(i-i0,j-j0);
                }
            }

        /** Set a submatrix.
        @param r    Array of row indices.
        @param c    Array of column indices.
        @param X    A(r(:),c(:))
        @exception  ArrayIndexOutOfBoundsException Submatrix indices
        */

        } else if (typeof(arguments[0]) === "object" && typeof(arguments[1]) === "object" && typeof(arguments[2]) === "object") {

            var r = arguments[0];
            var c = arguments[1];
            var X = arguments[2];

            for (var i = 0; i < r.length; i++) {
                for (var j = 0; j < c.length; j++) {
                    if (r[i] < 0 || r[i] >= this.m || c[j] < 0 || c[j] >= this.n) {
                        throw new ArrayIndexOutOfBoundsException("Submatrix indices");
                    }
                    this.A[r[i]][c[j]] = X.get(i,j);
                }
            }

        /** Set a submatrix.
        @param r    Array of row indices.
        @param j0   Initial column index
        @param j1   Final column index
        @param X    A(r(:),j0:j1)
        @exception  ArrayIndexOutOfBoundsException Submatrix indices
        */

        } else if (typeof(arguments[0]) === "object" && typeof(arguments[1]) === "number" && typeof(arguments[2]) === "number" && typeof(arguments[3]) === "object") {

            var r = arguments[0];
            var j0 = arguments[1];
            var j1 = arguments[2];
            var X = arguments[3];

            for (var i = 0; i < r.length; i++) {
                for (var j = j0; j <= j1; j++) {
                    if (r[i] < 0 || r[i] >= this.m || j < 0 || j >= this.n) {
                        throw new ArrayIndexOutOfBoundsException("Submatrix indices");
                    }
                    this.A[r[i]][j] = X.get(i,j-j0);
                }
            }

        /** Set a submatrix.
        @param i0   Initial row index
        @param i1   Final row index
        @param c    Array of column indices.
        @param X    A(i0:i1,c(:))
        @exception  ArrayIndexOutOfBoundsException Submatrix indices
        */

        } else if (typeof(arguments[0]) === "number" && typeof(arguments[1]) === "number" && typeof(arguments[2]) === "object" && typeof(arguments[2]) === "object") {

            var i0 = arguments[0];
            var i1 = arguments[1];
            var c = arguments[2];
            var X = arguments[3];

            for (var i = i0; i <= i1; i++) {
                for (var j = 0; j < c.length; j++) {
                    if (i < 0 || i >= this.m || c[j] < 0 || c[j] >= this.n) {
                        throw new ArrayIndexOutOfBoundsException("Submatrix indices");
                    }
                    this.A[i][c[j]] = X.get(i-i0,j);
                }
            }

        }
    }

    /** Matrix transpose.
    @return    A'
    */

    this.transpose = function() {
        var X = new Matrix(this.n,this.m);
        var C = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                C[j][i] = this.A[i][j];
            }
        }
        return X;
    }

    /** One norm
    @return    maximum column sum.
    */

    this.norm1 = function() {
        var f = 0;
        for (var j = 0; j < this.n; j++) {
            var s = 0;
            for (var i = 0; i < this.m; i++) {
                s += Math.abs(this.A[i][j]);
            }
            f = Math.max(f,s);
        }
        return f;
    }
//
//   /** Two norm
//   @return    maximum singular value.
//   */
//
//   public double norm2 () {
//      return (new SingularValueDecomposition(this).norm2());
//   }

    /** Infinity norm
    @return    maximum row sum.
    */

    this.normInf = function() {
        var f = 0;
        for (var i = 0; i < this.m; i++) {
            var s = 0;
            for (var j = 0; j < this.n; j++) {
                s += Math.abs(this.A[i][j]);
            }
            f = Math.max(f,s);
        }
        return f;
    }

    /** Frobenius norm
    @return    sqrt of sum of squares of all elements.
    */

    this.normF = function() {
        var f = 0;
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                f = hypot(f,this.A[i][j]);
            }
        }
        return f;
    }

    /**  Unary minus
    @return    -A
    */

    this.uminus = function() {
        var X = new Matrix(this.m,this.n);
        var C = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                C[i][j] = -this.A[i][j];
            }
        }
        return X;
    }

    /** C = A + B
    @param B    another matrix
    @return     A + B
    */

    this.plus = function(B) {
        this.checkMatrixDimensions(B);
        var X = new Matrix(this.m,this.n);
        var C = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                C[i][j] = this.A[i][j] + B.A[i][j];
            }
        }
        return X;
    }

    /** A = A + B
    @param B    another matrix
    @return     A + B
    */

    this.plusEquals = function(B) {
        this.checkMatrixDimensions(B);
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                this.A[i][j] = this.A[i][j] + B.A[i][j];
            }
        }
        return this;
    }

    /** C = A - B
    @param B    another matrix
    @return     A - B
    */

    this.minus = function(B) {
        this.checkMatrixDimensions(B);
        var X = new Matrix(this.m,this.n);
        var C = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                C[i][j] = this.A[i][j] - B.A[i][j];
            }
        }
        return X;
    }

    /** A = A - B
    @param B    another matrix
    @return     A - B
    */

    this.minusEquals = function(B) {
        this.checkMatrixDimensions(B);
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                this.A[i][j] = this.A[i][j] - B.A[i][j];
            }
        }
        return this;
    }

    /** Element-by-element multiplication, C = A.*B
    @param B    another matrix
    @return     A.*B
    */

    this.arrayTimes = function(B) {
        this.checkMatrixDimensions(B);
        var X = new Matrix(this.m,this.n);
        var C = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                C[i][j] = this.A[i][j] * B.A[i][j];
            }
        }
        return X;
    }

    /** Element-by-element multiplication in place, A = A.*B
    @param B    another matrix
    @return     A.*B
    */

    this.arrayTimesEquals = function(B) {
        this.checkMatrixDimensions(B);
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                this.A[i][j] = this.A[i][j] * B.A[i][j];
            }
        }
        return this;
    }

    /** Element-by-element right division, C = A./B
    @param B    another matrix
    @return     A./B
    */

    this.arrayRightDivide = function(B) {
        this.checkMatrixDimensions(B);
        var X = new Matrix(this.m,this.n);
        var C = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                C[i][j] = this.A[i][j] / B.A[i][j];
            }
        }
        return X;
    }

    /** Element-by-element right division in place, A = A./B
    @param B    another matrix
    @return     A./B
    */

    this.arrayRightDivideEquals = function(B) {
        this.checkMatrixDimensions(B);
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                this.A[i][j] = this.A[i][j] / B.A[i][j];
            }
        }
        return this;
    }

    /** Element-by-element left division, C = A.\B
    @param B    another matrix
    @return     A.\B
    */

    this.arrayLeftDivide = function(B) {
        this.checkMatrixDimensions(B);
        var X = new Matrix(this.m,this.n);
        var C = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                C[i][j] = B.A[i][j] / this.A[i][j];
            }
        }
        return X;
    }

    /** Element-by-element left division in place, A = A.\B
    @param B    another matrix
    @return     A.\B
    */

    this.arrayLeftDivideEquals = function(B) {
        this.checkMatrixDimensions(B);
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                this.A[i][j] = B.A[i][j] / this.A[i][j];
            }
        }
        return this;
    }

    /** Multiply a matrix by a scalar, C = s*A
    @param s    scalar
    @return     s*A
    */

    this.timesScalar = function(s) {
        var X = new Matrix(this.m,this.n);
        var C = X.getArray();
        for (var i = 0; i < this.m; i++) {
            for (var j = 0; j < this.n; j++) {
                C[i][j] = s*this.A[i][j];
            }
        }
        return X;
    }
//
//   /** Multiply a matrix by a scalar in place, A = s*A
//   @param s    scalar
//   @return     replace A by s*A
//   */
//
//   public Matrix timesEquals (double s) {
//      for (int i = 0; i < m; i++) {
//         for (int j = 0; j < n; j++) {
//            A[i][j] = s*A[i][j];
//         }
//      }
//      return this;
//   }

    /** Linear algebraic matrix multiplication, A * B
    @param B    another matrix
    @return     Matrix product, A * B
    @exception  IllegalArgumentException Matrix inner dimensions must agree.
    */

    this.times = function(B) {
        if (B.m != this.n) {
            throw new IllegalArgumentException("Matrix inner dimensions must agree.");
        }
        var X = new Matrix(this.m,B.n);
        var C = X.getArray();
        var Bcolj = new Array(n);
        for (var j = 0; j < B.n; j++) {
            for (var k = 0; k < this.n; k++) {
                Bcolj[k] = B.A[k][j];
            }
            for (var i = 0; i < this.m; i++) {
                var Arowi = this.A[i];
                var s = 0;
                for (var k = 0; k < this.n; k++) {
                    s += Arowi[k]*Bcolj[k];
                }
                C[i][j] = s;
            }
        }
        return X;
    }

    /** LU Decomposition
    @return     LUDecomposition
    @see LUDecomposition
    */

    this.lu = function() {
        return new LUDecomposition(this);
    }

    /** QR Decomposition
    @return     QRDecomposition
    @see QRDecomposition
    */

    this.qr = function() {
        return new QRDecomposition(this);
    }

    /** Cholesky Decomposition
    @return     CholeskyDecomposition
    @see CholeskyDecomposition
    */

    this.chol = function() {
        return new CholeskyDecomposition(this);
    }

    /** Singular Value Decomposition
    @return     SingularValueDecomposition
    @see SingularValueDecomposition
    */

    this.svd = function() {
        return new SingularValueDecomposition(this);
    }

    /** Eigenvalue Decomposition
    @return     EigenvalueDecomposition
    @see EigenvalueDecomposition
    */

    this.eig = function() {
        return new EigenvalueDecomposition(this);
    }

    /** Solve A*X = B
    @param B    right hand side
    @return     solution if A is square, least squares solution otherwise
    */

    this.solve = function(B) {
        return (this.m == this.n ? (new LUDecomposition(this)).solve(B) :
                                   (new QRDecomposition(this)).solve(B));
    }
//
//   /** Solve X*A = B, which is also A'*X' = B'
//   @param B    right hand side
//   @return     solution if A is square, least squares solution otherwise.
//   */
//
//   public Matrix solveTranspose (Matrix B) {
//      return transpose().solve(B.transpose());
//   }

    /** Matrix inverse or pseudoinverse
    @return     inverse(A) if A is square, pseudoinverse otherwise.
    */

    this.inverse = function() {
        return this.solve(Matrix.identity(this.m,this.m));
    }

    /** Matrix determinant
    @return     determinant
    */

    this.det = function() {
        return new LUDecomposition(this).det();
    }

    /** Matrix rank
    @return     effective numerical rank, obtained from SVD.
    */

    this.rank = function() {
        return new SingularValueDecomposition(this).rank();
    }

    /** Matrix condition (2 norm)
    @return     ratio of largest to smallest singular value.
    */

    this.cond = function() {
        return new SingularValueDecomposition(this).cond();
    }

    /** Matrix trace.
    @return     sum of the diagonal elements.
    */

    this.trace = function() {
        var t = 0;
        for (var i = 0; i < Math.min(this.m,this.n); i++) {
            t += this.A[i][i];
        }
        return t;
    }

   /** Generate matrix with random elements
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n matrix with uniformly distributed random elements.
   */

   Matrix.random = function(m, n) {
      var A = new Matrix(m,n);
      var X = A.getArray();
      for (var i = 0; i < m; i++) {
         for (var j = 0; j < n; j++) {
            X[i][j] = Math.random();
         }
      }
      return A;
   }

   /** Generate identity matrix
   @param m    Number of rows.
   @param n    Number of colums.
   @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
   */

   Matrix.identity = function(m, n) {
      var A = new Matrix(m,n);
      var X = A.getArray();
      for (var i = 0; i < m; i++) {
         for (var j = 0; j < n; j++) {
            X[i][j] = (i == j ? 1.0 : 0.0);
         }
      }
      return A;
   }
//
//
//   /** Print the matrix to stdout.   Line the elements up in columns
//     * with a Fortran-like 'Fw.d' style format.
//   @param w    Column width.
//   @param d    Number of digits after the decimal.
//   */
//
//   public void print (int w, int d) {
//      print(new PrintWriter(System.out,true),w,d); }
//
//   /** Print the matrix to the output stream.   Line the elements up in
//     * columns with a Fortran-like 'Fw.d' style format.
//   @param output Output stream.
//   @param w      Column width.
//   @param d      Number of digits after the decimal.
//   */
//
//   public void print (PrintWriter output, int w, int d) {
//      DecimalFormat format = new DecimalFormat();
//      format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
//      format.setMinimumIntegerDigits(1);
//      format.setMaximumFractionDigits(d);
//      format.setMinimumFractionDigits(d);
//      format.setGroupingUsed(false);
//      print(output,format,w+2);
//   }
//
//   /** Print the matrix to stdout.  Line the elements up in columns.
//     * Use the format object, and right justify within columns of width
//     * characters.
//     * Note that is the matrix is to be read back in, you probably will want
//     * to use a NumberFormat that is set to US Locale.
//   @param format A  Formatting object for individual elements.
//   @param width     Field width for each column.
//   @see java.text.DecimalFormat#setDecimalFormatSymbols
//   */
//
//   public void print (NumberFormat format, int width) {
//      print(new PrintWriter(System.out,true),format,width); }
//
//   // DecimalFormat is a little disappointing coming from Fortran or C's printf.
//   // Since it doesn't pad on the left, the elements will come out different
//   // widths.  Consequently, we'll pass the desired column width in as an
//   // argument and do the extra padding ourselves.
//
//   /** Print the matrix to the output stream.  Line the elements up in columns.
//     * Use the format object, and right justify within columns of width
//     * characters.
//     * Note that is the matrix is to be read back in, you probably will want
//     * to use a NumberFormat that is set to US Locale.
//   @param output the output stream.
//   @param format A formatting object to format the matrix elements 
//   @param width  Column width.
//   @see java.text.DecimalFormat#setDecimalFormatSymbols
//   */
//
//   public void print (PrintWriter output, NumberFormat format, int width) {
//      output.println();  // start on new line.
//      for (int i = 0; i < m; i++) {
//         for (int j = 0; j < n; j++) {
//            String s = format.format(A[i][j]); // format the number
//            int padding = Math.max(1,width-s.length()); // At _least_ 1 space
//            for (int k = 0; k < padding; k++)
//               output.print(' ');
//            output.print(s);
//         }
//         output.println();
//      }
//      output.println();   // end with blank line.
//   }
//
//   /** Read a matrix from a stream.  The format is the same the print method,
//     * so printed matrices can be read back in (provided they were printed using
//     * US Locale).  Elements are separated by
//     * whitespace, all the elements for each row appear on a single line,
//     * the last row is followed by a blank line.
//   @param input the input stream.
//   */
//
//   public static Matrix read (BufferedReader input) throws java.io.IOException {
//      StreamTokenizer tokenizer= new StreamTokenizer(input);
//
//      // Although StreamTokenizer will parse numbers, it doesn't recognize
//      // scientific notation (E or D); however, Double.valueOf does.
//      // The strategy here is to disable StreamTokenizer's number parsing.
//      // We'll only get whitespace delimited words, EOL's and EOF's.
//      // These words should all be numbers, for Double.valueOf to parse.
//
//      tokenizer.resetSyntax();
//      tokenizer.wordChars(0,255);
//      tokenizer.whitespaceChars(0, ' ');
//      tokenizer.eolIsSignificant(true);
//      java.util.Vector v = new java.util.Vector();
//
//      // Ignore initial empty lines
//      while (tokenizer.nextToken() == StreamTokenizer.TT_EOL);
//      if (tokenizer.ttype == StreamTokenizer.TT_EOF)
//	throw new java.io.IOException("Unexpected EOF on matrix read.");
//      do {
//         v.addElement(Double.valueOf(tokenizer.sval)); // Read & store 1st row.
//      } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
//
//      int n = v.size();  // Now we've got the number of columns!
//      double row[] = new double[n];
//      for (int j=0; j<n; j++)  // extract the elements of the 1st row.
//         row[j]=((Double)v.elementAt(j)).doubleValue();
//      v.removeAllElements();
//      v.addElement(row);  // Start storing rows instead of columns.
//      while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
//         // While non-empty lines
//         v.addElement(row = new double[n]);
//         int j = 0;
//         do {
//            if (j >= n) throw new java.io.IOException
//               ("Row " + v.size() + " is too long.");
//            row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
//         } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
//         if (j < n) throw new java.io.IOException
//            ("Row " + v.size() + " is too short.");
//      }
//      int m = v.size();  // Now we've got the number of rows.
//      double[][] A = new double[m][];
//      v.copyInto(A);  // copy the rows out of the vector
//      return new Matrix(A);
//   }


/* ------------------------
   Private Methods
 * ------------------------ */

    /** Check if size(A) == size(B) **/

    this.checkMatrixDimensions = function(B) {
        if (B.m != this.m || B.n != this.n) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
    }

}
