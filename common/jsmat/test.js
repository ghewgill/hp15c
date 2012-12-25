load("matrix.js");

/** TestMatrix tests the functionality of the Jama Matrix class and associated decompositions.
<P>
Run the test from the command line using
<BLOCKQUOTE><PRE><CODE>
 java Jama.test.TestMatrix 
</CODE></PRE></BLOCKQUOTE>
Detailed output is provided indicating the functionality being tested
and whether the functionality is correctly implemented.   Exception handling
is also tested.  
<P>
The test is designed to run to completion and give a summary of any implementation errors
encountered. The final output should be:
<BLOCKQUOTE><PRE><CODE>
      TestMatrix completed.
      Total errors reported: n1
      Total warning reported: n2
</CODE></PRE></BLOCKQUOTE>
If the test does not run to completion, this indicates that there is a 
substantial problem within the implementation that was not anticipated in the test design.  
The stopping point should give an indication of where the problem exists.
**/

function RuntimeException(message) {
    this.name = "RuntimeException";
    this.message = message;
}

/** private utility routines **/

function check() {
    print("you called check()");
}

/** Check magnitude of difference of scalars. **/

function check_scalar(x, y) {
    var eps = Math.pow(2.0,-52.0);
    if (x == 0 & Math.abs(y) < 10*eps) return;
    if (y == 0 & Math.abs(x) < 10*eps) return;
    if (Math.abs(x-y) > 10*eps*Math.max(Math.abs(x),Math.abs(y))) {
        throw new RuntimeException("The difference x-y is too large: x = " + Double.toString(x) + "  y = " + Double.toString(y));
    }
}

/** Check norm of difference of "vectors". **/

function check_vector(x, y) {
    if (x.length == y.length ) {
        for (var i=0;i<x.length;i++) {
            check_scalar(x[i],y[i]);
        } 
    } else {
        throw new RuntimeException("Attempt to compare vectors of different lengths");
    }
}

/** Check norm of difference of arrays. **/

function check_array(x, y) {
    var A = new Matrix(x);
    var B = new Matrix(y);
    check_matrix(A,B);
}

/** Check norm of difference of Matrices. **/

function check_matrix(X, Y) {
    var eps = Math.pow(2.0,-52.0);
    if (X.norm1() == 0. & Y.norm1() < 10*eps) return;
    if (Y.norm1() == 0. & X.norm1() < 10*eps) return;
    if (X.minus(Y).norm1() > 1000*eps*Math.max(X.norm1(),Y.norm1())) {
        throw new RuntimeException("The norm of (X-Y) is too large: " +  Double.toString(X.minus(Y).norm1()));
    }
}

/** Print appropriate messages for successful outcome try **/

function try_success (s, e) {
    print(">    " + s + "success");
    if ( e != "" ) {
        print(">      Message: " + e);
    }
}

/** Print appropriate messages for unsuccessful outcome try **/

function try_failure (count, s, e) {
    print(">    " + s + "*** failure ***\n>      Message: " + e);
    return ++count;
}

/** Print appropriate messages for unsuccessful outcome try **/

function try_warning (count, s, e) {
    print(">    " + s + "*** warning ***\n>      Message: " + e);
    return ++count;
}

/** Print a row vector. **/

function print_row(x, w, d) {
    // Use format Fw.d for all elements.
    print();
    new Matrix(x,1).print(w,d);
    print();
}

new Matrix(0,0); // hack to initialise pseudo-static functions

var A,B,C,Z,O,I,R,S,X,SUB,M,T,SQ,DEF,SOL;
var errorCount=0;
var warningCount=0;
var tmp, s;
var columnwise = [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.];
var rowwise = [1.,4.,7.,10.,2.,5.,8.,11.,3.,6.,9.,12.];
var avals = [[1.,4.,7.,10.],[2.,5.,8.,11.],[3.,6.,9.,12.]];
var rankdef = avals;
var tvals =  [[1.,2.,3.],[4.,5.,6.],[7.,8.,9.],[10.,11.,12.]];
var subavals = [[5.,8.,11.],[6.,9.,12.]];
var rvals = [[1.,4.,7.],[2.,5.,8.,11.],[3.,6.,9.,12.]];
var pvals = [[4.,1.,1.],[1.,2.,3.],[1.,3.,6.]];
var ivals = [[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.]];
var evals = 
   [[0.,1.,0.,0.],[1.,0.,2.e-7,0.],[0.,-2.e-7,0.,1.],[0.,0.,1.,0.]];
var square = [[166.,188.,210.],[188.,214.,240.],[210.,240.,270.]];
var sqSolution = [[13.],[15.]];
var condmat = [[1.,3.],[7.,9.]];
var rows=3,cols=4;
var invalidld=5;/* should trigger bad shape for construction with val */
var raggedr=0; /* (raggedr,raggedc) should be out of bounds in ragged array */
var raggedc=4; 
var validld=3; /* leading dimension of intended test Matrices */
var nonconformld=4; /* leading dimension which is valid, but nonconforming */
var ib=1,ie=2,jb=1,je=3; /* index ranges for sub Matrix */
var rowindexset = [1,2]; 
var badrowindexset = [1,3]; 
var columnindexset = [1,2,3];
var badcolumnindexset = [1,2,4];
var columnsummax = 33.;
var rowsummax = 30.;
var sumofdiagonals = 15;
var sumofsquares = 650;

/** 
Constructors and constructor-like methods:
    double[], int
    double[][]  
    int, int
    int, int, double
    int, int, double[][]
    constructWithCopy(double[][])
    random(int,int)
    identity(int)
**/

print("\nTesting constructors and constructor-like methods...");
try{  
    /** check that exception is thrown in packed constructor with invalid length **/
    A = new Matrix(columnwise,invalidld);
    errorCount = try_failure(errorCount,"Catch invalid length in packed constructor... ",
                "exception not thrown for invalid input");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("Catch invalid length in packed constructor... ",
                    e.getMessage());
    } else {
        throw e;
    }
}
try{ 
    /** check that exception is thrown in default constructor
        if input array is 'ragged' **/
    A = new Matrix(rvals);
    tmp = A.get(raggedr,raggedc);
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("Catch ragged input to default constructor... ",
                        e.getMessage());
    } else if (e.name === "ArrayIndexOutOfBoundsException") {
       errorCount = try_failure(errorCount,"Catch ragged input to constructor... ",
                   "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
    } else {
        throw e;
    }
}
try{ 
    /** check that exception is thrown in constructWithCopy
        if input array is 'ragged' **/
    A = Matrix.constructWithCopy(rvals);
    tmp = A.get(raggedr,raggedc);
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("Catch ragged input to constructWithCopy... ",e.getMessage());
    } else if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"Catch ragged input to constructWithCopy... ","exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
    } else {
        throw e;
    }
}

A = new Matrix(columnwise,validld);
B = new Matrix(avals);
tmp = B.get(0,0);
avals[0][0] = 0.0;
C = B.minus(A);
avals[0][0] = tmp;
B = Matrix.constructWithCopy(avals);
tmp = B.get(0,0);
avals[0][0] = 0.0;
if ( ( tmp - B.get(0,0) ) != 0.0 ) {
    /** check that constructWithCopy behaves properly **/
    errorCount = try_failure(errorCount,"constructWithCopy... ","copy not effected... data visible outside");
} else {
    try_success("constructWithCopy... ","");
}
avals[0][0] = columnwise[0]; 
I = new Matrix(ivals);
try {
    check_matrix(I,Matrix.identity(3,4));
    try_success("identity... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"identity... ","identity Matrix not successfully created");
}   

/**   
Access Methods:
    getColumnDimension()
    getRowDimension()
    getArray()
    getArrayCopy()
    getColumnPackedCopy()
    getRowPackedCopy()
    get(int,int)
    getMatrix(int,int,int,int)
    getMatrix(int,int,int[])
    getMatrix(int[],int,int)
    getMatrix(int[],int[])
    set(int,int,double)
    setMatrix(int,int,int,int,Matrix)
    setMatrix(int,int,int[],Matrix)
    setMatrix(int[],int,int,Matrix)
    setMatrix(int[],int[],Matrix)
**/

print("\nTesting access methods...");

/**
Various get methods:
**/

B = new Matrix(avals);
if (B.getRowDimension() != rows) {
    errorCount = try_failure(errorCount,"getRowDimension... ","");
} else {
    try_success("getRowDimension... ","");
}
if (B.getColumnDimension() != cols) {
    errorCount = try_failure(errorCount,"getColumnDimension... ","");
} else {
    try_success("getColumnDimension... ","");
}
B = new Matrix(avals);
var barray = B.getArray();
if ( barray != avals ) {
    errorCount = try_failure(errorCount,"getArray... ","");
} else {
    try_success("getArray... ","");
}
barray = B.getArrayCopy();
if ( barray == avals ) {
    errorCount = try_failure(errorCount,"getArrayCopy... ","data not (deep) copied");
}
try {
    check_array(barray,avals);
    try_success("getArrayCopy... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"getArrayCopy... ","data not successfully (deep) copied");
}
var bpacked = B.getColumnPackedCopy();
try {
    check_vector(bpacked,columnwise);
    try_success("getColumnPackedCopy... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"getColumnPackedCopy... ","data not successfully (deep) copied by columns");
}
bpacked = B.getRowPackedCopy();
try {
    check_vector(bpacked,rowwise);
    try_success("getRowPackedCopy... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"getRowPackedCopy... ","data not successfully (deep) copied by rows");
}
try {
    tmp = B.get(B.getRowDimension(),B.getColumnDimension()-1);
    errorCount = try_failure(errorCount,"get(int,int)... ","OutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            tmp = B.get(B.getRowDimension()-1,B.getColumnDimension());
            errorCount = try_failure(errorCount,"get(int,int)... ","OutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e1.name === "ArrayIndexOutOfBoundsException") {
                try_success("get(int,int)... OutofBoundsException... ","");
            } else {
                throw e1;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"get(int,int)... ","OutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
}
try {
    if (B.get(B.getRowDimension()-1,B.getColumnDimension()-1) != 
        avals[B.getRowDimension()-1][B.getColumnDimension()-1] ) {
        errorCount = try_failure(errorCount,"get(int,int)... ","Matrix entry (i,j) not successfully retreived");
    } else {
        try_success("get(int,int)... ","");
    }
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"get(int,int)... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}
SUB = new Matrix(subavals);
try {
    M = B.getMatrix(ib,ie+B.getRowDimension()+1,jb,je);
    errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)... ","ArrayIndexOutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            M = B.getMatrix(ib,ie,jb,je+B.getColumnDimension()+1);
            errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)... ","ArrayIndexOutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e1.name === "ArrayIndexOutOfBoundsException") {
                try_success("getMatrix(int,int,int,int)... ArrayIndexOutOfBoundsException... ","");
            } else {
                throw e1;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)... ","ArrayIndexOutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
}
try {
    M = B.getMatrix(ib,ie,jb,je);
    try {
        check_matrix(SUB,M);
        try_success("getMatrix(int,int,int,int)... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)... ","submatrix not successfully retreived");
    }
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"getMatrix(int,int,int,int)... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}

try {
    M = B.getMatrix(ib,ie,badcolumnindexset);
    errorCount = try_failure(errorCount,"getMatrix(int,int,int[])... ","ArrayIndexOutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            M = B.getMatrix(ib,ie+B.getRowDimension()+1,columnindexset);
            errorCount = try_failure(errorCount,"getMatrix(int,int,int[])... ","ArrayIndexOutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e1.name === "ArrayIndexOutOfBoundsException") {
                try_success("getMatrix(int,int,int[])... ArrayIndexOutOfBoundsException... ","");
            } else {
                throw e1;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"getMatrix(int,int,int[])... ","ArrayIndexOutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
} 
try {
    M = B.getMatrix(ib,ie,columnindexset);
    try {
        check_matrix(SUB,M);
        try_success("getMatrix(int,int,int[])... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"getMatrix(int,int,int[])... ","submatrix not successfully retreived");
    }
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"getMatrix(int,int,int[])... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}
try {
    M = B.getMatrix(badrowindexset,jb,je);
    errorCount = try_failure(errorCount,"getMatrix(int[],int,int)... ","ArrayIndexOutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            M = B.getMatrix(rowindexset,jb,je+B.getColumnDimension()+1);
            errorCount = try_failure(errorCount,"getMatrix(int[],int,int)... ","ArrayIndexOutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e1.name === "ArrayIndexOutOfBoundsException") {
                try_success("getMatrix(int[],int,int)... ArrayIndexOutOfBoundsException... ","");
            } else {
                throw e1;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"getMatrix(int[],int,int)... ","ArrayIndexOutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
} 
try {
    M = B.getMatrix(rowindexset,jb,je);
    try {
        check_matrix(SUB,M);
        try_success("getMatrix(int[],int,int)... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"getMatrix(int[],int,int)... ","submatrix not successfully retreived");
    }
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"getMatrix(int[],int,int)... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}
try {
    M = B.getMatrix(badrowindexset,columnindexset);
    errorCount = try_failure(errorCount,"getMatrix(int[],int[])... ","ArrayIndexOutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            M = B.getMatrix(rowindexset,badcolumnindexset);
            errorCount = try_failure(errorCount,"getMatrix(int[],int[])... ","ArrayIndexOutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e1.name === "ArrayIndexOutOfBoundsException") {
                try_success("getMatrix(int[],int[])... ArrayIndexOutOfBoundsException... ","");
            } else {
                throw e1;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"getMatrix(int[],int[])... ","ArrayIndexOutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
} 
try {
    M = B.getMatrix(rowindexset,columnindexset);
    try {
        check_matrix(SUB,M);
        try_success("getMatrix(int[],int[])... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"getMatrix(int[],int[])... ","submatrix not successfully retreived");
    }
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"getMatrix(int[],int[])... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}

/**
Various set methods:
**/

try {
    B.set(B.getRowDimension(),B.getColumnDimension()-1,0.);
    errorCount = try_failure(errorCount,"set(int,int,double)... ","OutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            B.set(B.getRowDimension()-1,B.getColumnDimension(),0.);
            errorCount = try_failure(errorCount,"set(int,int,double)... ","OutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e1.name === "ArrayIndexOutOfBoundsException") {
                try_success("set(int,int,double)... OutofBoundsException... ","");
            } else {
                throw e1;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"set(int,int,double)... ","OutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
}
try {
    B.set(ib,jb,0.);
    tmp = B.get(ib,jb);
    try {
        check_scalar(tmp,0.);
        try_success("set(int,int,double)... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"set(int,int,double)... ","Matrix element not successfully set");
    }
} catch ( e) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"set(int,int,double)... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}
M = new Matrix(2,3,0.);
try {
    B.setMatrix(ib,ie+B.getRowDimension()+1,jb,je,M);
    errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            B.setMatrix(ib,ie,jb,je+B.getColumnDimension()+1,M);
            errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e.name === "ArrayIndexOutOfBoundsException") {
                try_success("setMatrix(int,int,int,int,Matrix)... ArrayIndexOutOfBoundsException... ","");
            } else {
                throw e;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
}
try {
    B.setMatrix(ib,ie,jb,je,M);
    try {
        check_matrix(M.minus(B.getMatrix(ib,ie,jb,je)),M);
        try_success("setMatrix(int,int,int,int,Matrix)... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ","submatrix not successfully set");
    }
    B.setMatrix(ib,ie,jb,je,SUB);
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"setMatrix(int,int,int,int,Matrix)... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}
try {
    B.setMatrix(ib,ie+B.getRowDimension()+1,columnindexset,M);
    errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            B.setMatrix(ib,ie,badcolumnindexset,M);
            errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e1.name === "ArrayIndexOutOfBoundsException") {
                try_success("setMatrix(int,int,int[],Matrix)... ArrayIndexOutOfBoundsException... ","");
            } else {
                throw e1;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
}
try {
    B.setMatrix(ib,ie,columnindexset,M);
    try {
        check_matrix(M.minus(B.getMatrix(ib,ie,columnindexset)),M);
        try_success("setMatrix(int,int,int[],Matrix)... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ","submatrix not successfully set");
    }
    B.setMatrix(ib,ie,jb,je,SUB);
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"setMatrix(int,int,int[],Matrix)... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}
try {
    B.setMatrix(rowindexset,jb,je+B.getColumnDimension()+1,M);
    errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            B.setMatrix(badrowindexset,jb,je,M);
            errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
        } catch ( e ) {
            if (e.name === "ArrayIndexOutOfBoundsException") {
                try_success("setMatrix(int[],int,int,Matrix)... ArrayIndexOutOfBoundsException... ","");
            } else {
                throw e;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
}
try {
    B.setMatrix(rowindexset,jb,je,M);
    try {
        check_matrix(M.minus(B.getMatrix(rowindexset,jb,je)),M);
        try_success("setMatrix(int[],int,int,Matrix)... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ","submatrix not successfully set");
    }
    B.setMatrix(ib,ie,jb,je,SUB);
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"setMatrix(int[],int,int,Matrix)... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}
try {
    B.setMatrix(rowindexset,badcolumnindexset,M);
    errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        try {
            B.setMatrix(badrowindexset,columnindexset,M);
            errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
        } catch ( e1 ) {
            if (e.name === "ArrayIndexOutOfBoundsException") {
                try_success("setMatrix(int[],int[],Matrix)... ArrayIndexOutOfBoundsException... ","");
            } else {
                throw e;
            }
        }
    } else if (e.name === "IllegalArgumentException") {
        errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ","ArrayIndexOutOfBoundsException expected but not thrown");
    } else {
        throw e;
    }
}
try {
    B.setMatrix(rowindexset,columnindexset,M);
    try {
        check_matrix(M.minus(B.getMatrix(rowindexset,columnindexset)),M);
        try_success("setMatrix(int[],int[],Matrix)... ","");
    } catch ( e ) {
        errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ","submatrix not successfully set");
    }
} catch ( e ) {
    if (e.name === "ArrayIndexOutOfBoundsException") {
        errorCount = try_failure(errorCount,"setMatrix(int[],int[],Matrix)... ","Unexpected ArrayIndexOutOfBoundsException");
    } else {
        throw e;
    }
}

/** 
Array-like methods:
    minus
    minusEquals
    plus
    plusEquals
    arrayLeftDivide
    arrayLeftDivideEquals
    arrayRightDivide
    arrayRightDivideEquals
    arrayTimes
    arrayTimesEquals
    uminus
**/

print("\nTesting array-like methods...");
S = new Matrix(columnwise,nonconformld);
R = Matrix.random(A.getRowDimension(),A.getColumnDimension());
A = R;
try {
    S = A.minus(S);
    errorCount = try_failure(errorCount,"minus conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("minus conformance check... ","");
    } else {
        throw e;
    }
}
if (A.minus(R).norm1() != 0.) {
    errorCount = try_failure(errorCount,"minus... ","(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)");
} else {
    try_success("minus... ","");
}
A = R.copy();
A.minusEquals(R);
Z = new Matrix(A.getRowDimension(),A.getColumnDimension());
try {
    A.minusEquals(S);
    errorCount = try_failure(errorCount,"minusEquals conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("minusEquals conformance check... ","");
    } else {
        throw e;
    }
}
if (A.minus(Z).norm1() != 0.) {
    errorCount = try_failure(errorCount,"minusEquals... ","(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)");
} else {
    try_success("minusEquals... ","");
}

A = R.copy();
B = Matrix.random(A.getRowDimension(),A.getColumnDimension());
C = A.minus(B); 
try {
    S = A.plus(S);
    errorCount = try_failure(errorCount,"plus conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("plus conformance check... ","");
    } else {
        throw e;
    }
}
try {
    check_matrix(C.plus(B),A);
    try_success("plus... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"plus... ","(C = A - B, but C + B != A)");
}
C = A.minus(B);
C.plusEquals(B);
try {
    A.plusEquals(S);
    errorCount = try_failure(errorCount,"plusEquals conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("plusEquals conformance check... ","");
    } else {
        throw e;
    }
}
try {
    check_matrix(C,A);
    try_success("plusEquals... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"plusEquals... ","(C = A - B, but C = C + B != A)");
}
A = R.uminus();
try {
    check_matrix(A.plus(R),Z);
    try_success("uminus... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"uminus... ","(-A + A != zeros)");
}
A = R.copy();
O = new Matrix(A.getRowDimension(),A.getColumnDimension(),1.0);
C = A.arrayLeftDivide(R);
try {
    S = A.arrayLeftDivide(S);
    errorCount = try_failure(errorCount,"arrayLeftDivide conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("arrayLeftDivide conformance check... ","");
    } else {
        throw e;
    }
}
try {
    check_matrix(C,O);
    try_success("arrayLeftDivide... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"arrayLeftDivide... ","(M.\\M != ones)");
}
try {
    A.arrayLeftDivideEquals(S);
    errorCount = try_failure(errorCount,"arrayLeftDivideEquals conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("arrayLeftDivideEquals conformance check... ","");
    } else {
        throw e;
    }
}
A.arrayLeftDivideEquals(R);
try {
    check_matrix(A,O);
    try_success("arrayLeftDivideEquals... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"arrayLeftDivideEquals... ","(M.\\M != ones)");
}
A = R.copy();
try {
    A.arrayRightDivide(S);
    errorCount = try_failure(errorCount,"arrayRightDivide conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("arrayRightDivide conformance check... ","");
    } else {
        throw e;
    }
}
C = A.arrayRightDivide(R);
try {
    check_matrix(C,O);
    try_success("arrayRightDivide... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"arrayRightDivide... ","(M./M != ones)");
}
try {
    A.arrayRightDivideEquals(S);
    errorCount = try_failure(errorCount,"arrayRightDivideEquals conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("arrayRightDivideEquals conformance check... ","");
    } else {
        throw e;
    }
}
A.arrayRightDivideEquals(R);
try {
    check_matrix(A,O);
    try_success("arrayRightDivideEquals... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"arrayRightDivideEquals... ","(M./M != ones)");
}
A = R.copy();
B = Matrix.random(A.getRowDimension(),A.getColumnDimension());
try {
    S = A.arrayTimes(S);
    errorCount = try_failure(errorCount,"arrayTimes conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("arrayTimes conformance check... ","");
    } else {
        throw e;
    }
}
C = A.arrayTimes(B);
try {
    check_matrix(C.arrayRightDivideEquals(B),A);
    try_success("arrayTimes... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"arrayTimes... ","(A = R, C = A.*B, but C./B != A)");
}
try {
    A.arrayTimesEquals(S);
    errorCount = try_failure(errorCount,"arrayTimesEquals conformance check... ","nonconformance not raised");
} catch ( e ) {
    if (e.name === "IllegalArgumentException") {
        try_success("arrayTimesEquals conformance check... ","");
    } else {
        throw e;
    }
}
A.arrayTimesEquals(B);
try {
    check_matrix(A.arrayRightDivideEquals(B),R);
    try_success("arrayTimesEquals... ","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"arrayTimesEquals... ","(A = R, A = A.*B, but A./B != R)");
}

///**   
//I/O methods:
//   read
//   print
//   serializable:
//     writeObject
//     readObject
//**/
//print("\nTesting I/O methods...\n");
//   try {
//      DecimalFormat fmt = new DecimalFormat("0.0000E00");
//	    fmt.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
//
//      PrintWriter FILE = new PrintWriter(new FileOutputStream("JamaTestMatrix.out"));
//      A.print(FILE,fmt,10);
//      FILE.close();
//      R = Matrix.read(new BufferedReader(new FileReader("JamaTestMatrix.out")));
//      if (A.minus(R).norm1() < .001 ) {
//         try_success("print()/read()...","");
//      } else {
//         errorCount = try_failure(errorCount,"print()/read()...","Matrix read from file does not match Matrix printed to file");
//      }
//   } catch ( java.io.IOException ioe ) {
//     warningCount = try_warning(warningCount,"print()/read()...","unexpected I/O error, unable to run print/read test;  check write permission in current directory and retry");
//   } catch(Exception e) {
//      try {
//         e.printStackTrace(System.out);
//         warningCount = try_warning(warningCount,"print()/read()...","Formatting error... will try JDK1.1 reformulation...");
//         DecimalFormat fmt = new DecimalFormat("0.0000");
//         PrintWriter FILE = new PrintWriter(new FileOutputStream("JamaTestMatrix.out"));
//         A.print(FILE,fmt,10);
//         FILE.close();
//         R = Matrix.read(new BufferedReader(new FileReader("JamaTestMatrix.out")));
//         if (A.minus(R).norm1() < .001 ) {
//            try_success("print()/read()...","");
//         } else {
//            errorCount = try_failure(errorCount,"print()/read() (2nd attempt) ...","Matrix read from file does not match Matrix printed to file");
//         }
//      } catch ( java.io.IOException ioe ) {
//        warningCount = try_warning(warningCount,"print()/read()...","unexpected I/O error, unable to run print/read test;  check write permission in current directory and retry");
//   }
//}
//
//R = Matrix.random(A.getRowDimension(),A.getColumnDimension());
//String tmpname = "TMPMATRIX.serial";
//try {
//   ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(tmpname));
//   out.writeObject(R);
//   ObjectInputStream sin = new ObjectInputStream(new FileInputStream(tmpname));
//   A = (Matrix) sin.readObject();
// 
//   try {
//      check(A,R);
//      try_success("writeObject(Matrix)/readObject(Matrix)...","");
//   } catch ( java.lang.RuntimeException e ) {
//     errorCount = try_failure(errorCount,"writeObject(Matrix)/readObject(Matrix)...","Matrix not serialized correctly");
//   }
//} catch ( java.io.IOException ioe ) {
//   warningCount = try_warning(warningCount,"writeObject()/readObject()...","unexpected I/O error, unable to run serialization test;  check write permission in current directory and retry");
//} catch(Exception e) {
//   errorCount = try_failure(errorCount,"writeObject(Matrix)/readObject(Matrix)...","unexpected error in serialization test");
//}

/**
LA methods:
    transpose
    times
    cond
    rank
    det
    trace
    norm1
    norm2
    normF
    normInf
    solve
    solveTranspose
    inverse
    chol
    eig
    lu
    qr
    svd 
**/

print("\nTesting linear algebra methods...");
A = new Matrix(columnwise,3);
T = new Matrix(tvals);
T = A.transpose();
try {
    check_matrix(A.transpose(),T);
    try_success("transpose...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"transpose()...","transpose unsuccessful");
}
A.transpose();
try {
    check_scalar(A.norm1(),columnsummax);
    try_success("norm1...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"norm1()...","incorrect norm calculation");
}
try {
    check_scalar(A.normInf(),rowsummax);
    try_success("normInf()...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"normInf()...","incorrect norm calculation");
}
try {
    check_scalar(A.normF(),Math.sqrt(sumofsquares));
    try_success("normF...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"normF()...","incorrect norm calculation");
}
try {
    check_scalar(A.trace(),sumofdiagonals);
    try_success("trace()...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"trace()...","incorrect trace calculation");
}
try {
    check_scalar(A.getMatrix(0,A.getRowDimension()-1,0,A.getRowDimension()-1).det(),0.);
    try_success("det()...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"det()...","incorrect determinant calculation");
}
SQ = new Matrix(square);
try {
    check_matrix(A.times(A.transpose()),SQ);
    try_success("times(Matrix)...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"times(Matrix)...","incorrect Matrix-Matrix product calculation");
}
try {
    check_matrix(A.timesScalar(0.),Z);
    try_success("times(double)...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"times(double)...","incorrect Matrix-scalar product calculation");
}

A = new Matrix(columnwise,4);
var QR = A.qr();
R = QR.getR();
try {
   check_matrix(A,QR.getQ().times(R));
   try_success("QRDecomposition...","");
} catch ( e ) {
   errorCount = try_failure(errorCount,"QRDecomposition...","incorrect QR decomposition calculation");
}
var SVD = A.svd();
try {
    check_matrix(A,SVD.getU().times(SVD.getS().times(SVD.getV().transpose())));
    try_success("SingularValueDecomposition...","");
} catch ( e ) {
    print(e); print(e.name);
    errorCount = try_failure(errorCount,"SingularValueDecomposition...","incorrect singular value decomposition calculation");
}
DEF = new Matrix(rankdef);
try {
    check_scalar(DEF.rank(),Math.min(DEF.getRowDimension(),DEF.getColumnDimension())-1);
    try_success("rank()...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"rank()...","incorrect rank calculation");
}
B = new Matrix(condmat);
SVD = B.svd(); 
var singularvalues = SVD.getSingularValues();
try {
    check_scalar(B.cond(),singularvalues[0]/singularvalues[Math.min(B.getRowDimension(),B.getColumnDimension())-1]);
    try_success("cond()...","");
} catch ( e ) {
    print(e);
    errorCount = try_failure(errorCount,"cond()...","incorrect condition number calculation");
}
var n = A.getColumnDimension();
A = A.getMatrix(0,n-1,0,n-1);
A.set(0,0,0.);
var LU = A.lu();
try {
    check_matrix(A.getMatrix(LU.getPivot(),0,n-1),LU.getL().times(LU.getU()));
    try_success("LUDecomposition...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"LUDecomposition...","incorrect LU decomposition calculation");
    print(e);
}
X = A.inverse();
try {
    check_matrix(A.times(X),Matrix.identity(3,3));
    try_success("inverse()...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"inverse()...","incorrect inverse calculation");
}
O = new Matrix(SUB.getRowDimension(),1,1.0);
SOL = new Matrix(sqSolution);
SQ = SUB.getMatrix(0,SUB.getRowDimension()-1,0,SUB.getRowDimension()-1);
try {
    check_matrix(SQ.solve(SOL),O); 
    try_success("solve()...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"solve()...",e.getMessage());
}
A = new Matrix(pvals);
var Chol = A.chol(); 
var L = Chol.getL();
try {
    check_matrix(A,L.times(L.transpose()));
    try_success("CholeskyDecomposition...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"CholeskyDecomposition...","incorrect Cholesky decomposition calculation");
}
X = Chol.solve(Matrix.identity(3,3));
try {
    check_matrix(A.times(X),Matrix.identity(3,3));
    try_success("CholeskyDecomposition solve()...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"CholeskyDecomposition solve()...","incorrect Choleskydecomposition solve calculation");
}
var Eig = A.eig();
var D = Eig.getD();
var V = Eig.getV();
try {
    check_matrix(A.times(V),V.times(D));
    try_success("EigenvalueDecomposition (symmetric)...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"EigenvalueDecomposition (symmetric)...","incorrect symmetric Eigenvalue decomposition calculation");
}
A = new Matrix(evals);
Eig = A.eig();
D = Eig.getD();
V = Eig.getV();
try {
    check_matrix(A.times(V),V.times(D));
    try_success("EigenvalueDecomposition (nonsymmetric)...","");
} catch ( e ) {
    errorCount = try_failure(errorCount,"EigenvalueDecomposition (nonsymmetric)...","incorrect nonsymmetric Eigenvalue decomposition calculation");
}

print("\nTestMatrix completed.\n");
print("Total errors reported: " + errorCount);
print("Total warnings reported: " + warningCount);
