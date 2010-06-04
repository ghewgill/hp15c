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

//   /** Is the matrix nonsingular?
//   @return     true if U, and hence A, is nonsingular.
//   */
//
//   public boolean isNonsingular () {
//      for (int j = 0; j < n; j++) {
//         if (LU[j][j] == 0)
//            return false;
//      }
//      return true;
//   }
//
//   /** Return lower triangular factor
//   @return     L
//   */
//
//   public Matrix getL () {
//      Matrix X = new Matrix(m,n);
//      double[][] L = X.getArray();
//      for (int i = 0; i < m; i++) {
//         for (int j = 0; j < n; j++) {
//            if (i > j) {
//               L[i][j] = LU[i][j];
//            } else if (i == j) {
//               L[i][j] = 1.0;
//            } else {
//               L[i][j] = 0.0;
//            }
//         }
//      }
//      return X;
//   }
//
//   /** Return upper triangular factor
//   @return     U
//   */
//
//   public Matrix getU () {
//      Matrix X = new Matrix(n,n);
//      double[][] U = X.getArray();
//      for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//            if (i <= j) {
//               U[i][j] = LU[i][j];
//            } else {
//               U[i][j] = 0.0;
//            }
//         }
//      }
//      return X;
//   }
//
//   /** Return pivot permutation vector
//   @return     piv
//   */
//
//   public int[] getPivot () {
//      int[] p = new int[m];
//      for (int i = 0; i < m; i++) {
//         p[i] = piv[i];
//      }
//      return p;
//   }
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
//
//   /** Solve A*X = B
//   @param  B   A Matrix with as many rows as A and any number of columns.
//   @return     X so that L*U*X = B(piv,:)
//   @exception  IllegalArgumentException Matrix row dimensions must agree.
//   @exception  RuntimeException  Matrix is singular.
//   */
//
//   public Matrix solve (Matrix B) {
//      if (B.getRowDimension() != m) {
//         throw new IllegalArgumentException("Matrix row dimensions must agree.");
//      }
//      if (!this.isNonsingular()) {
//         throw new RuntimeException("Matrix is singular.");
//      }
//
//      // Copy right hand side with pivoting
//      int nx = B.getColumnDimension();
//      Matrix Xmat = B.getMatrix(piv,0,nx-1);
//      double[][] X = Xmat.getArray();
//
//      // Solve L*Y = B(piv,:)
//      for (int k = 0; k < n; k++) {
//         for (int i = k+1; i < n; i++) {
//            for (int j = 0; j < nx; j++) {
//               X[i][j] -= X[k][j]*LU[i][k];
//            }
//         }
//      }
//      // Solve U*X = Y;
//      for (int k = n-1; k >= 0; k--) {
//         for (int j = 0; j < nx; j++) {
//            X[k][j] /= LU[k][k];
//         }
//         for (int i = 0; i < k; i++) {
//            for (int j = 0; j < nx; j++) {
//               X[i][j] -= X[k][j]*LU[i][k];
//            }
//         }
//      }
//      return Xmat;
//   }
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

        this.A = A;
        this.m = m;
        this.n = n;

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
      for (var i = 0; i < m; i++) {
         for (var j = 0; j < n; j++) {
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
                    B[i][j] = A[r[i]][c[j]];
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
                    B[i-i0][j] = A[i][c[j]];
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
                    B[i][j-j0] = A[r[i]][j];
                }
            }
            return X;

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
//
//   /** LU Decomposition
//   @return     LUDecomposition
//   @see LUDecomposition
//   */
//
//   public LUDecomposition lu () {
//      return new LUDecomposition(this);
//   }
//
//   /** QR Decomposition
//   @return     QRDecomposition
//   @see QRDecomposition
//   */
//
//   public QRDecomposition qr () {
//      return new QRDecomposition(this);
//   }
//
//   /** Cholesky Decomposition
//   @return     CholeskyDecomposition
//   @see CholeskyDecomposition
//   */
//
//   public CholeskyDecomposition chol () {
//      return new CholeskyDecomposition(this);
//   }
//
//   /** Singular Value Decomposition
//   @return     SingularValueDecomposition
//   @see SingularValueDecomposition
//   */
//
//   public SingularValueDecomposition svd () {
//      return new SingularValueDecomposition(this);
//   }
//
//   /** Eigenvalue Decomposition
//   @return     EigenvalueDecomposition
//   @see EigenvalueDecomposition
//   */
//
//   public EigenvalueDecomposition eig () {
//      return new EigenvalueDecomposition(this);
//   }
//
//   /** Solve A*X = B
//   @param B    right hand side
//   @return     solution if A is square, least squares solution otherwise
//   */
//
//   public Matrix solve (Matrix B) {
//      return (m == n ? (new LUDecomposition(this)).solve(B) :
//                       (new QRDecomposition(this)).solve(B));
//   }
//
//   /** Solve X*A = B, which is also A'*X' = B'
//   @param B    right hand side
//   @return     solution if A is square, least squares solution otherwise.
//   */
//
//   public Matrix solveTranspose (Matrix B) {
//      return transpose().solve(B.transpose());
//   }
//
//   /** Matrix inverse or pseudoinverse
//   @return     inverse(A) if A is square, pseudoinverse otherwise.
//   */
//
//   public Matrix inverse () {
//      return solve(identity(m,m));
//   }

    /** Matrix determinant
    @return     determinant
    */

    this.det = function() {
        return new LUDecomposition(this).det();
    }
//
//   /** Matrix rank
//   @return     effective numerical rank, obtained from SVD.
//   */
//
//   public int rank () {
//      return new SingularValueDecomposition(this).rank();
//   }
//
//   /** Matrix condition (2 norm)
//   @return     ratio of largest to smallest singular value.
//   */
//
//   public double cond () {
//      return new SingularValueDecomposition(this).cond();
//   }

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
