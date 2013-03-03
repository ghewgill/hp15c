var On = true;
var Display;
var DecimalSwap;
var Stack = [0, 0, 0, 0];
var StackI = [0, 0, 0, 0];
var LastX = 0;
var LastXI = 0;
var Reg = new Array(60);
var DisableKeys;
var Entry;
var DigitEntry = false;
var NewDigitEntry;
var StackLift = false;
var OldStackLift, NewStackLift;
var DelayUpdate = 0;
var DisplayTimeout = 0;
var TemporaryDisplay = false;
var Shift = 0;
var Prefix;
var OldPrefix;
var LcdDisplay;
var DisplayMode = 1; // 1=FIX 2=SCI 3=ENG
var DisplayDigits = 4;
var FullCircle = 360;
var TrigFactor = Math.PI / 180;
var Flags = [false, false, false, false, false, false, false, false, false, false];
var User = false;
var Prgm = false;
var Program = [null];
var PC = 0;
var Running = false;
var RunTimer = null;
var BlinkOn = false;
var Blinker = null;
var ReturnStack = [];
var Result = 0;

var MAX = 9.999999999e99;
var MAX_MAG = 99;

var OpcodeIndex = {'/':10, '*':20, '-':30, '+':40};

function CalcError(n) {
    this.name = "CalcError";
    this.message = "Error " + n;
    this.code = n;
}

function OpcodeInfo(keys, defn, programmable, user) {
    this.keys = keys;
    this.defn = defn;
    this.programmable = programmable !== false;
    this.user = user;
}

var OpSqrt      = new OpcodeInfo([11],      op_sqrt);
var OpA         = new OpcodeInfo([42,11],   op_A);
var OpX2        = new OpcodeInfo([43,11],   op_x2);
var OpEx        = new OpcodeInfo([12],      op_ex);
var OpB         = new OpcodeInfo([42,12],   op_B);
var OpLn        = new OpcodeInfo([43,12],   op_ln);
var Op10x       = new OpcodeInfo([13],      op_10x);
var OpC         = new OpcodeInfo([42,13],   op_C);
var OpLog       = new OpcodeInfo([43,13],   op_log);
var OpYx        = new OpcodeInfo([14],      op_yx);
var OpD         = new OpcodeInfo([42,14],   op_D);
var OpPct       = new OpcodeInfo([43,14],   op_pct);
var Op1x        = new OpcodeInfo([15],      op_1x);
var OpE         = new OpcodeInfo([42,15],   op_E);
var OpDpct      = new OpcodeInfo([43,15],   op_dpct);
var OpChs       = new OpcodeInfo([16],      op_chs);
//var OpMatrix
var OpAbs       = new OpcodeInfo([43,16],   op_abs);
var Op7         = new OpcodeInfo([7],       function() { op_input('7'); });
//var OpFix
var OpDeg       = new OpcodeInfo([43,7],    op_deg);
var Op8         = new OpcodeInfo([8],       function() { op_input('8'); });
//var OpSci
var OpRad       = new OpcodeInfo([43,8],    op_rad);
var Op9         = new OpcodeInfo([9],       function() { op_input('9'); });
//var OpEng
var OpGrd       = new OpcodeInfo([43,9],    op_grd);
var OpDiv       = new OpcodeInfo([10],      op_div);
//var OpSolve
var OpLe        = new OpcodeInfo([43,10],   op_le);
var OpSst       = new OpcodeInfo([21],      op_sst, false);
//var OpLbl
var OpBst       = new OpcodeInfo([43,21],   op_bst, false);
//var OpGto
//var OpHyp
//var OpAhyp
var OpSin       = new OpcodeInfo([23],      op_sin);
//var OpDim
var OpAsin      = new OpcodeInfo([43,23],   op_asin);
var OpCos       = new OpcodeInfo([24],      op_cos);
var OpIndex     = new OpcodeInfo([42,24],   op_index, false);
var OpAcos      = new OpcodeInfo([43,24],   op_acos);
var OpTan       = new OpcodeInfo([25],      op_tan);
var OpI         = new OpcodeInfo([42,25],   op_I);
var OpAtan      = new OpcodeInfo([43,25],   op_atan);
var OpEex       = new OpcodeInfo([26],      op_eex);
//var OpResult
var OpPi        = new OpcodeInfo([43,26],   op_pi);
var Op4         = new OpcodeInfo([4],       function() { op_input('4'); });
//var OpXchg
//var OpSf
var Op5         = new OpcodeInfo([5],       function() { op_input('5'); });
//var OpDse       = new OpcodeInfo([42,5],    op_dse);
//var OpCf
var Op6         = new OpcodeInfo([6],       function() { op_input('6'); });
//var OpIsg       = new OpcodeInfo([42,6],    op_isg);
//var OpFtest
var OpMul       = new OpcodeInfo([20],      op_mul);
//var OpIntegrate
var OpEq        = new OpcodeInfo([43,20],   op_eq);
var OpRs        = new OpcodeInfo([31],      op_rs);
var OpPse       = new OpcodeInfo([42,31],   op_pse);
var OpPr        = new OpcodeInfo([43,31],   op_pr, false);
//var OpGsb
var OpClearStat = new OpcodeInfo([42,32],   op_clear_stat);
var OpRtn       = new OpcodeInfo([43,32],   op_rtn);
var OpRoll      = new OpcodeInfo([33],      op_roll);
var OpClearPrgm = new OpcodeInfo([42,33],   op_clear_prgm, false);
var OpRollup    = new OpcodeInfo([43,33],   op_rollup);
var OpXy        = new OpcodeInfo([34],      op_xy);
var OpClearReg  = new OpcodeInfo([42,34],   op_clear_reg);
var OpRnd       = new OpcodeInfo([43,34],   op_rnd);
var OpBack      = new OpcodeInfo([35],      op_back, false);
var OpClearPrefix=new OpcodeInfo([42,35],   op_clear_prefix, false);
var OpClx       = new OpcodeInfo([43,35],   op_clx);
var OpEnter     = new OpcodeInfo([36],      op_enter);
var OpRand      = new OpcodeInfo([42,36],   op_rand);
var OpLastx     = new OpcodeInfo([43,36],   op_lastx);
var Op1         = new OpcodeInfo([1],       function() { op_input('1'); });
var OpToR       = new OpcodeInfo([42,1],    op_to_r);
var OpToP       = new OpcodeInfo([43,1],    op_to_p);
var Op2         = new OpcodeInfo([2],       function() { op_input('2'); });
var OpToHms     = new OpcodeInfo([42,2],    op_to_hms);
var OpToH       = new OpcodeInfo([43,2],    op_to_h);
var Op3         = new OpcodeInfo([3],       function() { op_input('3'); });
var OpToRad     = new OpcodeInfo([42,3],    op_to_rad);
var OpToDeg     = new OpcodeInfo([43,3],    op_to_deg);
var OpSub       = new OpcodeInfo([30],      op_sub);
var OpReIm      = new OpcodeInfo([42,30],   op_re_im);
//var OpTest
var OpOn        = new OpcodeInfo([41],      op_on, false);
//var OpSto
var OpFrac      = new OpcodeInfo([42,44],   op_frac);
var OpInt       = new OpcodeInfo([43,44],   op_int);
//var OpRcl
var OpUser      = new OpcodeInfo([42,46],   op_user, false);
var OpMem       = new OpcodeInfo([43,46],   op_mem, false);
var Op0         = new OpcodeInfo([0],       function() { op_input('0'); });
var OpFact      = new OpcodeInfo([42,0],    op_fact);
var OpMean      = new OpcodeInfo([43,0],    op_mean);
var OpDot       = new OpcodeInfo([48],      function() { op_input('.'); });
var OpYhat      = new OpcodeInfo([42,48],   op_yhat);
var OpS         = new OpcodeInfo([43,48],   op_s);
var OpSum       = new OpcodeInfo([49],      op_sum);
var OpLr        = new OpcodeInfo([42,49],   op_lr);
var OpSumsub    = new OpcodeInfo([43,49],   op_sumsub);
var OpAdd       = new OpcodeInfo([40],      op_add);
var OpPyx       = new OpcodeInfo([42,40],   op_Pyx);
var OpCyx       = new OpcodeInfo([43,40],   op_Cyx);

function _(x) { return function(k) { return new Opcode(x); }; }

var CharTable = {
    'q': [_(OpSqrt), _(OpA), _(OpX2)],
    'E': [_(OpEx), _(OpB), _(OpLn)],
    ')': [_(Op10x), _(OpC), _(OpLog)],
    '^': [_(OpYx), _(OpD), _(OpPct)],
    '\\':[_(Op1x), _(OpE), _(OpDpct)],
    '_': [_(OpChs), decode_matrix, _(OpAbs)],
    '7': [_(Op7), decode_fix, _(OpDeg)],
    '8': [_(Op8), decode_sci, _(OpRad)],
    '9': [_(Op9), decode_eng, _(OpGrd)],
    '/': [_(OpDiv), decode_solve, _(OpLe)],
    'T': [_(OpSst), decode_lbl, _(OpBst)],
    'G': [decode_gto, decode_hyp, decode_ahyp],
    's': [_(OpSin), decode_dim, _(OpAsin)],
    'c': [_(OpCos), _(OpIndex), _(OpAcos)],
    't': [_(OpTan), _(OpI), _(OpAtan)],
    'e': [_(OpEex), decode_result, _(OpPi)],
    '4': [_(Op4), decode_xchg, decode_sf],
    '5': [_(Op5), decode_dse, decode_cf],
    '6': [_(Op6), decode_isg, decode_ftest],
    '*': [_(OpMul), decode_integrate, _(OpEq)],
    'P': [_(OpRs), _(OpPse), _(OpPr)],
    'U': [decode_gsb, _(OpClearStat), _(OpRtn)],
    'r': [_(OpRoll), _(OpClearPrgm), _(OpRollup)],
    'x': [_(OpXy), _(OpClearReg), _(OpRnd)],
    '\b': [_(OpBack), _(OpClearPrefix), _(OpClx)],
    '\r': [_(OpEnter), _(OpRand), _(OpLastx)],
    '\n': [_(OpEnter), _(OpRand), _(OpLastx)],
    '1': [_(Op1), _(OpToR), _(OpToP)],
    '2': [_(Op2), _(OpToHms), _(OpToH)],
    '3': [_(Op3), _(OpToRad), _(OpToDeg)],
    '-': [_(OpSub), _(OpReIm), decode_test],
    '\x1b': [_(OpOn), _(OpOn), _(OpOn)],
    'f': [decode_f, decode_f, decode_f],
    'g': [decode_g, decode_g, decode_g],
    'S': [decode_sto, _(OpFrac), _(OpInt)],
    'R': [decode_rcl, _(OpUser), _(OpMem)],
    '0': [_(Op0), _(OpFact), _(OpMean)],
    '.': [_(OpDot), _(OpYhat), _(OpS)],
    ';': [_(OpSum), _(OpLr), _(OpSumsub)],
    '+': [_(OpAdd), _(OpPyx), _(OpCyx)],
    ' ': _(OpEnter),
    '!': _(OpFact),
    '@': _(OpX2),
    '%': _(OpPct),
    'A': _(OpA),
    'B': _(OpB),
    'C': _(OpC),
    'D': _(OpD),
    'L': _(OpLastx),
    'a': _(OpAbs),
    'i': _(OpInt),
    'I': _(OpI),
    'l': _(OpLn),
    'p': _(OpPi),
    '\x12': _(OpRand)
};

var KeyTable = [
    ['q', 'E', ')', '^', '\\','_', '7', '8', '9', '/'],
    ['T', 'G', 's', 'c', 't', 'e', '4', '5', '6', '*'],
    ['P', 'U', 'r', 'x', '\b','\r','1', '2', '3', '-'],
    ['\x1b', 'f', 'g', 'S', 'R', '\r','0', '.', ';', '+']
];
var ExtraKeyTable = [
    [3, 6, -1, '!'],
    [0, 0,  1, '@'],
    [0, 3,  1, '%'],
    [0, 0, -1, 'A'],
    [0, 1, -1, 'B'],
    [0, 2, -1, 'C'],
    [0, 3, -1, 'D'],
    [3, 5,  1, 'L'],
    [0, 5,  1, 'a'],
    [3, 3,  1, 'i'],
    [1, 4, -1, 'I'],
    [0, 1,  1, 'l'],
    [1, 5,  1, 'p']
];

function sign(x) {
    if (x === 0) {
        return 0;
    }
    return (x > 0) ? 1 : -1;
}

function sin_drg_mode(x) {
    if (FullCircle === 360 || FullCircle === 400) {
        // check for exact values in deg/grad mode to avoid roundoff errors
        var t = Math.abs(x % FullCircle);
        if (t === 0 || t === FullCircle/2) {
            return 0;
        }
    }
    return Math.sin(x * TrigFactor);
}

function cos_drg_mode(x) {
    if (FullCircle === 360 || FullCircle === 400) {
        // check for exact values in deg/grad mode to avoid roundoff errors
        var t = Math.abs(x % FullCircle);
        if (t === FullCircle/4 || t === FullCircle*3/4) {
            return 0;
        }
    }
    return Math.cos(x * TrigFactor);
}

function tan_drg_mode(x) {
    if (FullCircle === 360 || FullCircle === 400) {
        // check for exact values in deg/grad mode to avoid roundoff errors
        var t = Math.abs(x % FullCircle);
        if (t === 0 || t === FullCircle/2) {
            return 0;
        }
        if (t === FullCircle/4 || t === FullCircle*3/4) {
            Flags[9] = true;
            return MAX;
        }
    }
    return Math.tan(x * TrigFactor);
}

function sinh(x) {
    return (Math.exp(x) - Math.exp(-x)) / 2;
}

function cosh(x) {
    return (Math.exp(x) + Math.exp(-x)) / 2;
}

function tanh(x) {
    return (Math.exp(2*x) - 1) / (Math.exp(2*x) + 1);
}

function Complex(re, im) {
    this.re = re;
    this.im = im;

    this.acos = function() {
        var s1 = new Complex(1 - this.re, -this.im).sqrt();
        var s2 = new Complex(1 + this.re, this.im).sqrt();
        var r1 = 2 * Math.atan2(s1.re, s2.re);
        var i1 = (s2.re*s1.im) - (s2.im*s1.re);
        i1 = sign(i1) * Math.log(Math.abs(i1) + Math.sqrt(i1*i1 + 1));
        return new Complex(r1, i1);
    };

    this.acosh = function() {
        var s1 = new Complex(this.re - 1, this.im).sqrt();
        var s2 = new Complex(this.re + 1, this.im).sqrt();
        var r1 = (s1.re*s2.re) + (s1.im*s2.im);
        r1 = sign(r1) * Math.log(Math.abs(r1) + Math.sqrt(r1*r1 + 1));
        var i1 = 2 * Math.atan2(s1.im, s2.re);
        return new Complex(r1, i1);
    };

    this.abs = function() {
        return Math.sqrt(this.re*this.re + this.im*this.im);
    };

    this.add = function(that) {
        return new Complex(this.re + that.re, this.im + that.im);
    };

    this.arg = function() {
        return Math.atan2(this.im, this.re);
    };

    this.asin = function() {
        var s1 = new Complex(1 + this.re, this.im).sqrt();
        var s2 = new Complex(1 - this.re, -this.im).sqrt();
        var r1 = (s1.re*s2.im) - (s2.re*s1.im);
        r1 = sign(r1) * Math.log(Math.abs(r1) + Math.sqrt(r1*r1 + 1));
        var i1 = Math.atan2(this.re, (s1.re*s2.re) - (s1.im*s2.im));
        return new Complex(i1, -r1);
    };

    this.asinh = function() {
        var s1 = new Complex(1 + this.im, -this.re).sqrt();
        var s2 = new Complex(1 - this.im, this.re).sqrt();
        var r1 = (s1.re * s2.im) - (s2.re * s1.im);
        r1 = sign(r1) * Math.log(Math.abs(r1) + Math.sqrt(r1 * r1 + 1));
        var i1 = Math.atan2(this.im, (s1.re * s2.re) - (s1.im * s2.im));
        return new Complex(r1, i1);
    };

    this.atan = function() {
        if (this.re === 0 && Math.abs(this.im) === 1) {
            Flags[9] = true;
            return new Complex(0, sign(this.im) * MAX);
        }
        var rsign = 1;
        if (this.re === 0 && Math.abs(this.im) > 1) {
            rsign = -1;
        }
        var u = Complex.i.add(this).div(Complex.i.sub(this));
        var w = u.log();
        return new Complex(rsign * -w.im/2, w.re/2);
    };

    this.atanh = function() {
        if (this.im === 0 && Math.abs(this.re) === 1) {
            Flags[9] = true;
            return sign(this.re) * MAX;
        }
        var u = Complex.one.add(this).div(Complex.one.sub(this)).log();
        return new Complex(u.re/2, u.im/2);
    };

    this.cos = function() {
        return new Complex(Math.cos(this.re)*cosh(this.im), -Math.sin(this.re)*sinh(this.im));
    };

    this.cosh = function() {
        return new Complex(cosh(this.re)*Math.cos(this.im), sinh(this.re)*Math.sin(this.im));
    };

    this.div = function(that) {
        var d = that.re*that.re + that.im*that.im;
        return new Complex((this.re*that.re + this.im*that.im) / d, (this.im*that.re - this.re*that.im) / d);
    };

    this.exp = function() {
        var r = Math.exp(this.re);
        return new Complex(r * Math.cos(this.im), r * Math.sin(this.im));
    };

    this.exp10 = function() {
        var r = Math.pow(10, this.re);
        var t = this.im * Math.LN10;
        return new Complex(r * Math.cos(t), r * Math.sin(t));
    };

    this.inv = function() {
        var d = this.re*this.re + this.im*this.im;
        return new Complex(this.re/d, -this.im/d);
    };

    this.log = function() {
        return new Complex(Math.log(this.re*this.re + this.im*this.im)/2, this.arg());
    };

    this.log10 = function() {
        var u = this.log();
        return new Complex(u.re / Math.LN10, u.im / Math.LN10);
    };

    this.mul = function(that) {
        return new Complex(this.re * that.re - this.im * that.im, this.re * that.im + this.im * that.re);
    };

    this.pow = function(that) {
        var p = this.arg();
        var a = this.abs();
        var r = Math.pow(a, that.re) * Math.exp(-that.im * p);
        var t = that.re * p + that.im * Math.log(a);
        return new Complex(r * Math.cos(t), r * Math.sin(t));
    };

    this.sin = function() {
        return new Complex(Math.sin(this.re)*cosh(this.im), Math.cos(this.re)*sinh(this.im));
    };

    this.sinh = function() {
        return new Complex(sinh(this.re)*Math.cos(this.im), cosh(this.re)*Math.sin(this.im));
    };

    this.sub = function(that) {
        return new Complex(this.re - that.re, this.im - that.im);
    };

    this.sqrt = function() {
        var a = this.abs();
        return new Complex(Math.sqrt((this.re + a) / 2), (this.im < 0 ? -1 : 1) * Math.sqrt((-this.re + a) / 2));
    };

    this.square = function() {
        return new Complex(this.re*this.re - this.im*this.im, 2*this.re*this.im);
    };

    this.tan = function() {
        var u = new Complex(Math.tan(this.re), tanh(this.im));
        return u.div(new Complex(1, -u.re*u.im));
    };

    this.tanh = function() {
        var u = new Complex(tanh(this.re), Math.tan(this.im));
        return u.div(new Complex(1, u.re*u.im));
    };

    this.toString = function() {
        return "(" + this.re + "," + this.im + ")";
    };
}

Complex.one = new Complex(1, 0);
Complex.i = new Complex(0, 1);

var A = 0;
var B = 1;
var C = 2;
var D = 3;
var E = 4;

function Mat() {
    if (typeof(arguments[0]) === "number" && typeof(arguments[1]) === "number") {
        this.rows = arguments[0];
        this.cols = arguments[1];
        this.m = new Matrix(this.rows, this.cols, arguments[2]);
    } else if (typeof(arguments[0]) === "object") {
        this.m = arguments[0];
        this.rows = this.m.getRowDimension();
        this.cols = this.m.getColumnDimension();
    }

    this.complex2 = function() {
        var r = new Matrix(this.rows, this.cols*2);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.cols; j++) {
                r.set(i, j, this.m.get(i, j));
                if (i < this.rows/2) {
                    r.set(this.rows/2+i, this.cols+j, this.m.get(i, j));
                } else {
                    r.set(i-this.rows/2, this.cols+j, -this.m.get(i, j));
                }
            }
        }
        return new Mat(r);
    };

    this.complex3 = function() {
        return new Mat(this.m.getMatrix(0, this.rows-1, 0, this.cols/2-1));
    };

    this.copy = function() {
        return new Mat(this.m.copy());
    };

    this.det = function() {
        return this.m.det();
    };

    this.get = function(row, col) {
        return this.m.get(row-1, col-1);
    };

    this.inverse = function() {
        return new Mat(this.m.inverse());
    };

    this.minus = function(B) {
        return new Mat(this.m.minus(B.m));
    };

    this.norm = function() {
        return this.m.normInf();
    };

    this.normF = function() {
        return this.m.normF();
    };

    this.partition = function() {
        var r = new Matrix(this.rows*2, this.cols/2);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.cols; j += 2) {
                r.set(i, j/2, this.m.get(i, j));
                r.set(this.rows+i, j/2, this.m.get(i, j+1));
            }
        }
        return new Mat(r);
    };

    this.plus = function(B) {
        return new Mat(this.m.plus(B.m));
    };

    this.residual = function(Y, X) {
        return new Mat(this.m.minus(Y.m.times(X.m)));
    };

    this.set = function(row, col, value) {
        this.m.set(row-1, col-1, value);
    };

    this.times = function(B) {
        return new Mat(this.m.times(B.m));
    };

    this.timesScalar = function(s) {
        return new Mat(this.m.timesScalar(s));
    };

    this.transpose = function() {
        return new Mat(this.m.transpose());
    };

    this.toString = function() {
        return "<Mat " + this.rows + "," + this.cols + ">";
    };

    this.unpartition = function() {
        var r = new Matrix(this.rows/2, this.cols*2);
        for (var i = 0; i < this.rows/2; i++) {
            for (var j = 0; j < this.cols; j++) {
                r.set(i, j*2, this.m.get(i, j));
                r.set(i, j*2+1, this.m.get(this.rows/2+i, j));
            }
        }
        return new Mat(r);
    };
}

var g_Matrix = [new Mat(0, 0),
                new Mat(0, 0),
                new Mat(0, 0),
                new Mat(0, 0),
                new Mat(0, 0)];

function Descriptor(label) {
    this.label = label;

    this.toString = function() {
        return "<Descriptor " + this.label + " (" + g_Matrix[this.label].rows + "," + g_Matrix[this.label].cols + ")>";
    };
}

function Opcode(info, fn) {
    this.info = info;
    this.fn = fn;

    this.exec = function() {
        OldStackLift = StackLift;
        NewDigitEntry = false;
        NewStackLift = true;
        try {
            if (this.fn === undefined) {
                this.info.defn();
            } else {
                this.fn();
            }
        } finally {
            DigitEntry = NewDigitEntry;
            StackLift = NewStackLift;
            if (DigitEntry) {
                if (Entry.length > 0 && Entry.charAt(Entry.length-1) === 'e') {
                    Entry += "0";
                }
                Stack[0] = Number(Entry);
            }
        }
    };
}

function trunc(x) {
    if (x < 0) {
        return -Math.floor(-x);
    } else {
        return Math.floor(x);
    }
}

function log10int(x) {
    var mag = 0;
    x = Math.abs(x);
    if (x >= 1) {
        while (x >= 10) {
            mag++;
            x /= 10;
        }
    } else if (x > 0) {
        while (x < 1) {
            mag--;
            x *= 10;
        }
    }
    return mag;
}

function update_lcd(s) {
    LcdDisplay = s;
    Display.clear_digits();
    if (!On) {
        return;
    }
    if (Flags[9] && !BlinkOn) {
        return;
    }
    var eex = null;
    var d = 0;
    for (var i = 0; i < s.length; i++) {
        var c = s.charAt(i);
        if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'E') || c === 'o' || c === 'r' || c === 'u') {
            if (eex !== null) {
                eex = eex * 10 + c.charCodeAt(0) - "0".charCodeAt(0);
                var t = sprintf("%02d", eex);
                Display.set_digit(8, t.charAt(0));
                Display.set_digit(9, t.charAt(1));
            } else if (d < 10) {
                Display.set_digit(d, c);
                d++;
            }
        } else if (c === '.') {
            Display.set_decimal(d-1);
        } else if (c === ',') {
            Display.set_comma(d-1);
        } else if (c === '-') {
            if (eex !== null) {
                Display.set_digit(7, "-");
            } else if (d > 0) {
                Display.set_digit(d, "-");
                d++;
            } else {
                Display.set_neg();
            }
        } else if (c === 'e') {
            eex = 0;
            d = 9;
            Display.clear_digit(7);
            Display.set_digit(8, "0");
            Display.set_digit(9, "0");
        } else if (c === ' ') {
            d++;
        }
    }
}

function insert_commas(s) {
    var sign = "";
    if (s.charAt(0) === "-") {
        sign = s.charAt(0);
        s = s.substr(1);
    }
    var d = s.indexOf(".");
    if (d < 0) {
        d = s.indexOf("e");
        if (d < 0) {
            d = s.length;
        }
    }
    while (true) {
        d -= 3;
        if (d <= 0) {
            break;
        }
        s = s.substr(0, d) + "," + s.substr(d);
    }
    if (DecimalSwap) {
        for (var i = 0; i < s.length; i++) {
            if (s.charAt(i) === ".") {
                s = s.substr(0, i) + "," + s.substr(i+1);
            } else if (s.charAt(i) === ",") {
                s = s.substr(0, i) + "." + s.substr(i+1);
            }
        }
    }
    return sign + s;
}

function format_fix(n) {
    var x = Math.round(n * Math.pow(10, DisplayDigits));
    // TODO: var s
    s = x.toString();
    while (s.length < DisplayDigits+1) {
        s = '0' + s;
    }
    s = s.substr(0, s.length-DisplayDigits) + '.' + s.substr(s.length-DisplayDigits);
    s = insert_commas(s);
    return s;
}

function format_sci(n, mag) {
    var x = Math.round(n * Math.pow(10, DisplayDigits - mag));
    // Not sure what case the following loop addresses
    while (log10int(x) > DisplayDigits) {
        if (mag >= MAX_MAG) {
            x = Math.floor(n * Math.pow(10, DisplayDigits - mag));
            break;
        }
        mag++;
        x /= 10;
    }
    s = x.toString();
    while (s.length < DisplayDigits+1) {
        s = '0' + s;
    }
    s = s.substr(0, 1) + '.' + s.substr(1);
    s += "e" + mag;
    return s;
}

function format_eng(n, mag) {
    var x = Math.round(n * Math.pow(10, DisplayDigits - mag));
    s = x.toString();
    while (s.length < DisplayDigits+1) {
        s = '0' + s;
    }
    var ilen = 1;
    while (mag % 3) {
        ilen++;
        if (s.length < ilen) {
            s += '0';
        }
        mag--;
    }
    s = s.substr(0, ilen) + '.' + s.substr(ilen);
    s += "e" + mag;
    return s;
}

function update_display_num(n) {
    if (n instanceof Descriptor) {
        update_lcd(sprintf("%c    %2d %2d",
            "A".charCodeAt(0) + n.label,
            g_Matrix[n.label].rows,
            g_Matrix[n.label].cols));
    } else {
        if (n > MAX) {
            n = MAX;
            Flags[9] = true;
        } else if (n < -MAX) {
            n = -MAX;
            Flags[9] = true;
        }
        var s = n.toString();
        var mag = log10int(n);
        var sign = "";
        if (n < 0) {
            n = -n;
            sign = "-";
        }
        var dm = DisplayMode;
        if (dm === 1 && (mag >= 10 || mag < -DisplayDigits)) {
            dm = 2;
        }
        switch (dm) {
        case 1:
            s = format_fix(n);
            break;
        case 2:
            s = format_sci(n, mag);
            break;
        case 3:
            s = format_eng(n, mag);
            break;
        }
        update_lcd(sign + s);
    }
}

function update_display() {
    if (Prgm) {
        var s = sprintf("%03d-", PC);
        if (PC > 0) {
            if (Program[PC].info.user) {
                s = s.substr(0, 3) + "u";
            }
            var keys = Program[PC].info.keys;
            switch (keys.length) {
                case 1:
                    s += sprintf("    %2d", keys[0]);
                    break;
                case 2:
                    if (Math.floor(keys[1]) === keys[1]) {
                        s += sprintf(" %2d %2d", keys[0], keys[1]);
                    } else {
                        s += sprintf(" %2d  .%d", keys[0], keys[1]*10);
                    }
                    break;
                case 3:
                    if (Math.floor(keys[2]) === keys[2]) {
                        s += sprintf("%2d,%2d,%2d", keys[0], keys[1], keys[2]);
                    } else {
                        s += sprintf("%2d,%2d,.%d", keys[0], keys[1], keys[2]*10);
                    }
                    break;
            }
        }
        update_lcd(s);
    } else {
        if (DigitEntry) {
            if (Entry !== "") {
                update_lcd(insert_commas(Entry));
            } else {
                update_display_num(0);
            }
        } else {
            update_display_num(Stack[0]);
        }
    }
    Display.set_complex(Flags[8]);
    Display.set_prgm(Prgm);
    if (Flags[9]) {
        if (Blinker === null) {
            Blinker = setInterval(function() {
                BlinkOn = !BlinkOn;
                update_display();
            }, 300);
        }
    } else {
        if (Blinker !== null) {
            clearInterval(Blinker);
            Blinker = null;
        }
    }
}

function push(x, forcelift) {
    if (forcelift || StackLift) {
        Stack[3] = Stack[2]; StackI[3] = StackI[2];
        Stack[2] = Stack[1]; StackI[2] = StackI[1];
        Stack[1] = Stack[0]; StackI[1] = StackI[0];
    }
    Stack[0] = x;
    StackI[0] = 0;
    StackLift = true;
}

function fill(x) {
    for (var i in Stack) {
        Stack[i] = x;
    }
}

function unop(f) {
    LastX = Stack[0];
    LastXI = StackI[0];
    var r = f(Stack[0]);
    if (isNaN(r) || r === Infinity || r === -Infinity) {
        throw new CalcError(0);
    }
    Stack[0] = r;
}

function binop(f) {
    LastX = Stack[0];
    LastXI = StackI[0];
    var r = f(Stack[1], Stack[0]);
    if (isNaN(r) || r === Infinity || r === -Infinity) {
        throw new CalcError(0);
    }
    Stack[0] = r;
    Stack[1] = Stack[2]; StackI[1] = StackI[2];
    Stack[2] = Stack[3]; StackI[2] = StackI[3];
}

function unopc(f) {
    LastX = Stack[0];
    LastXI = StackI[0];
    var r = f(new Complex(Stack[0], StackI[0]));
    if (r instanceof Complex) {
        Stack[0] = r.re;
        StackI[0] = r.im;
    } else {
        Stack[0] = r;
        StackI[0] = 0;
    }
}

function binopc(f) {
    LastX = Stack[0];
    LastXI = StackI[0];
    var r = f(new Complex(Stack[1], StackI[1]), new Complex(Stack[0], StackI[0]));
    Stack[0] = r.re;
    StackI[0] = r.im;
    Stack[1] = Stack[2]; StackI[1] = StackI[2];
    Stack[2] = Stack[3]; StackI[2] = StackI[3];
}

function unopm(f) {
    LastX = Stack[0];
    LastXI = StackI[0];
    var r = f(Stack[0]);
    g_Matrix[Result] = r;
    Stack[0] = new Descriptor(Result);
}

function binopm(f) {
    LastX = Stack[0];
    LastXI = StackI[0];
    var r = f(Stack[1], Stack[0]);
    g_Matrix[Result] = r;
    Stack[0] = new Descriptor(Result);
    Stack[1] = Stack[2]; StackI[1] = StackI[2];
    Stack[2] = Stack[3]; StackI[2] = StackI[3];
}

function op_sqrt() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.sqrt();
        });
    } else {
        unop(Math.sqrt);
    }
}

function op_A() {
    op_gsb(11);
}

function op_x2() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.square();
        });
    } else {
        unop(function(x) { return x * x; });
    }
}

function op_ex() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.exp();
        });
    } else {
        unop(Math.exp);
    }
}

function op_B() {
    op_gsb(12);
}

function op_ln() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.log();
        });
    } else {
        unop(Math.log);
    }
}

function op_10x() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.exp10();
        });
    } else {
        unop(function(x) { return Math.pow(10, x); });
    }
}

function op_C() {
    op_gsb(13);
}

function op_log() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.log10();
        });
    } else {
        unop(function(x) { return Math.log(x) / Math.LN10; });
    }
}

function op_yx() {
    if (Flags[8]) {
        binopc(function(y, x) {
            return y.pow(x);
        });
    } else {
        binop(Math.pow);
    }
}

function op_D() {
    op_gsb(14);
}

function op_pct() {
    LastX = Stack[0];
    LastXI = StackI[0];
    Stack[0] = Stack[1] * Stack[0] / 100;
}

function op_1x() {
    if (Stack[0] instanceof Descriptor) {
        unopm(function(x) {
            return g_Matrix[x.label].inverse();
        });
    } else if (Flags[8]) {
        unopc(function(x) {
            return x.inv();
        });
    } else {
        unop(function(x) { return 1 / x; });
    }
}

function op_E() {
    op_gsb(15);
}

function op_dpct() {
    binop(function(y, x) {
        return (x - y) / y * 100;
    });
}

function op_chs() {
    if (DigitEntry) {
        var i = Entry.indexOf('e');
        if (i >= 0) {
            if (i+1 < Entry.length && Entry.charAt(i+1) === '-') {
                Entry = Entry.substr(0, i+1) + Entry.substr(i+2);
            } else {
                Entry = Entry.substr(0, i+1) + "-" + Entry.substr(i+1);
            }
        } else if (Entry.charAt(0) === '-') {
            Entry = Entry.substr(1);
        } else {
            Entry = '-' + Entry;
        }
        NewDigitEntry = true;
    } else if (Stack[0] instanceof Descriptor) {
        var x = Stack[0].label;
        g_Matrix[x] = g_Matrix[x].timesScalar(-1);
    } else {
        Stack[0] = -Stack[0];
    }
}

function op_matrix_clear() {
    g_Matrix = [new Mat(0, 0),
                new Mat(0, 0),
                new Mat(0, 0),
                new Mat(0, 0),
                new Mat(0, 0)];
}

function op_matrix_home() {
    Reg[0] = 1;
    Reg[1] = 1;
}

function op_matrix_complex2() {
    var m = Stack[0].label;
    g_Matrix[m] = g_Matrix[m].complex2();
}

function op_matrix_complex3() {
    var m = Stack[0].label;
    g_Matrix[m] = g_Matrix[m].complex3();
}

function op_matrix_transpose() {
    var m = Stack[0].label;
    g_Matrix[m] = g_Matrix[m].transpose();
}

function op_matrix_transmul() {
    binopm(function(y, x) {
        return g_Matrix[y.label].transpose().times(g_Matrix[x.label]);
    });
}

function op_matrix_residual() {
    binopm(function(y, x) {
        var Y = g_Matrix[y.label];
        var X = g_Matrix[x.label];
        var R = g_Matrix[Result];
        return R.residual(Y, X);
    });
}

function op_matrix_norm() {
    unop(function(x) {
        return g_Matrix[x.label].norm();
    });
}

function op_matrix_normf() {
    unop(function(x) {
        return g_Matrix[x.label].normF();
    });
}

function op_matrix_det() {
    unop(function(x) {
        return g_Matrix[x.label].det();
    });
}

function op_abs() {
    if (Stack[0] instanceof Descriptor) {
        throw new CalcError(1);
    }
    if (Flags[8]) {
        unopc(function(x) {
            return new Complex(x.abs(), 0);
        });
    } else {
        unop(Math.abs);
    }
}

function op_fix(n) {
    DisplayMode = 1;
    DisplayDigits = n;
    StackLift = OldStackLift;
}

function op_fix_index() {
    op_fix(trunc(Reg.I));
}

function op_deg() {
    FullCircle = 360;
    TrigFactor = Math.PI / 180;
    Display.set_trigmode(null);
    StackLift = OldStackLift;
}

function op_sci(n) {
    DisplayMode = 2;
    DisplayDigits = n;
    StackLift = OldStackLift;
}

function op_sci_index() {
    op_sci(trunc(Reg.I));
}

function op_rad() {
    FullCircle = Math.PI * 2; // for consistency, but we will not use this value
    TrigFactor = 1;
    Display.set_trigmode("RAD");
    StackLift = OldStackLift;
}

function op_eng(n) {
    DisplayMode = 3;
    DisplayDigits = n;
    StackLift = OldStackLift;
}

function op_eng_index() {
    op_eng(trunc(Reg.I));
}

function op_grd() {
    FullCircle = 400;
    TrigFactor = Math.PI / 200;
    Display.set_trigmode("GRAD");
    StackLift = OldStackLift;
}

function op_div() {
    if (Stack[0] instanceof Descriptor && Stack[1] instanceof Descriptor) {
        binopm(function(y, x) {
            return g_Matrix[x.label].inverse().times(g_Matrix[y.label]);
        });
    } else if (Stack[0] instanceof Descriptor) {
        binopm(function(y, x) {
            return g_Matrix[x.label].inverse().timesScalar(y);
        });
    } else if (Stack[1] instanceof Descriptor) {
        binopm(function(y, x) {
            return g_Matrix[y.label].timesScalar(1/x);
        });
    } else if (Flags[8]) {
        binopc(function(y, x) {
            return y.div(x);
        });
    } else {
        binop(function(y, x) { return y / x; });
    }
}

function op_solve(n) {
    var call = function(n) {
        var r = ReturnStack.length;
        op_gsb(n);
        while (ReturnStack.length > r) {
            step();
        }
    };
    // This is http://mathworld.wolfram.com/SecantMethod.html
    var eps = 1e-9;
    var maxiter = 100;
    var x0, x1, x2, y0, y1, y2;
    x0 = Stack[1];
    x1 = Stack[0];
    if (x1 === x0) {
        x1 = x0 * 1.01;
    }
    while (true) {
        fill(x0);
        call(n);
        y0 = Stack[0];
        fill(x1);
        call(n);
        y1 = Stack[0];
        x2 = x1 - y1 * ((x1 - x0) / (y1 - y0));
        if (isNaN(x2) || x2 === Infinity || x2 === -Infinity) {
            if (Running) {
                PC++;
                return;
            } else {
                throw new CalcError(8);
            }
        }
        fill(x2);
        call(n);
        y2 = Stack[0];
        //alert("x0=" + x0 + " y0=" + y0 + "\n" +
        //      "x1=" + x1 + " y1=" + y1 + "\n" +
        //      "x2=" + x2 + " y2=" + y2 + "\n");
        x0 = x1;
        x1 = x2;
        if (Math.abs(x1 - x0) < eps) {
            break;
        }
        if (--maxiter <= 0) {
            if (Running) {
                PC++;
                return;
            } else {
                throw new CalcError(8);
            }
        }
    }
    push(y2);
    push(x1);
    push(x2);
}

function op_le() {
    if (Stack[0] <= Stack[1]) {
        // execute next opcode
    } else {
        PC++;
    }
}

function op_sst() {
    if (Prgm) {
        PC++;
        if (PC >= Program.length) {
            PC = 0;
        }
    } else {
        step();
    }
    StackLift = OldStackLift;
}

function op_bst() {
    if (PC > 0) {
        PC--;
    }
    StackLift = OldStackLift;
}

function op_gto_immediate(n) {
    PC = n;
}

function op_gto_label(n) {
    var p = PC + 1;
    while (true) {
        if (p >= Program.length) {
            p = 0;
        } else {
            //print(p + ": " + Program[p].info.keys);
            if (Program[p].info.keys.toString() === [42,21,n].toString()) {
                break;
            }
        }
        if (p === PC) {
            throw new CalcError(4);
        }
        p++;
    }
    PC = p;
}

function op_gto_index() {
    var i = trunc(Reg.I);
    if (i < 0) {
        PC = -i;
    } else if (i < 10) {
        op_gto_label(i);
    } else if (i < 20) {
        op_gto_label((i - 10) / 10);
    } else if (i < 24) {
        op_gto_label(i - 9);
    } else {
        alert("error 4 ");
    }
}

function op_sinh() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.sinh();
        });
    } else {
        unop(sinh);
    }
}

function op_cosh() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.cosh();
        });
    } else {
        unop(cosh);
    }
}

function op_tanh() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.tanh();
        });
    } else {
        unop(tanh);
    }
}

function op_asinh() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.asinh();
        });
    } else {
        unop(function(x) {
            return sign(x) * Math.log(Math.abs(x) + Math.sqrt(x*x + 1));
        });
    }
}

function op_acosh() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.acosh();
        });
    } else {
        unop(function(x) {
            return Math.log(x + Math.sqrt(x*x - 1));
        });
    }
}

function op_atanh() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.atanh();
        });
    } else {
        unop(function(x) {
            return Math.log((1 + x) / (1 - x)) / 2;
        });
    }
}

function op_sin() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.sin();
        });
    } else {
        unop(sin_drg_mode);
    }
}

function op_dim(m) {
    var oldmat = g_Matrix[m].m;
    var r = Stack[1];
    var c = Stack[0];
    var i;
    if (oldmat.getRowDimension() > 0 && r > 0 && c > 0) {
        var a = oldmat.getArray();
        if (c < a[0].length) {
            for (i = 0; i < a.length; i++) {
                a[i] = a[i].slice(0, c);
            }
        } else if (c > a[0].length) {
            for (i = 0; i < a.length; i++) {
                for (j = a[i].length; j < c; j++) {
                    a[i][j] = 0;
                }
            }
        }
        if (r < a.length) {
            a = a.slice(0, r);
        } else {
            while (r > a.length) {
                var b = new Array(a[0].length);
                for (i = 0; i < c; i++) {
                    b[i] = 0;
                }
                a[a.length] = b;
            }
        }
        g_Matrix[m] = new Mat(new Matrix(a));
    } else {
        g_Matrix[m] = new Mat(Stack[1], Stack[0]);
    }
}

function op_dim_i() {
    // noop
}

function op_asin() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.asin();
        });
    } else {
        unop(function(x) {
            return Math.asin(x) / TrigFactor;
        });
    }
}

function op_cos() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.cos();
        });
    } else {
        unop(cos_drg_mode);
    }
}

function op_index() {
    if (Flags[8]) {
        update_display_num(StackI[0]);
        DelayUpdate = 1000;
    }
    StackLift = OldStackLift;
}

function op_acos() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.acos();
        });
    } else {
        unop(function(x) {
            return Math.acos(x) / TrigFactor;
        });
    }
}

function op_tan() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.tan();
        });
    } else {
        unop(tan_drg_mode);
    }
}

function op_I() {
    Flags[8] = true;
    StackI[0] = Stack[0];
    Stack[0] = Stack[1];
    Stack[1] = Stack[2]; StackI[1] = StackI[2];
    Stack[2] = Stack[3]; StackI[2] = StackI[3];
}

function op_atan() {
    if (Flags[8]) {
        unopc(function(x) {
            return x.atan();
        });
    } else {
        unop(function(x) {
            return Math.atan(x) / TrigFactor;
        });
    }
}

function op_eex() {
    op_input('e');
}

function op_result(m) {
    Result = m;
}

function op_pi() {
    push(Math.PI);
}

function op_xchg(r) {
    var t = Reg[r];
    Reg[r] = Stack[0];
    Stack[0] = t;
}

function op_xchg_index() {
    op_xchg(Math.floor(Math.abs(Reg.I)));
}

function op_sf(f) {
    Flags[f] = true;
}

function op_sf_index() {
    op_sf(Math.floor(Math.abs(Reg.I)));
}

function op_dse(r) {
    var n = trunc(Reg[r]);
    var f = Reg[r] - n;
    var s = sprintf("%.5f", f);
    var x = +s.substr(2, 3);
    var y = +s.substr(5, 2);
    if (y === 0) {
        y = 1;
    }
    n -= y;
    Reg[r] = n + f;
    if (n <= x) {
        PC++;
    }
}

function op_dse_index() {
    op_dse(Math.floor(Math.abs(Reg.I)));
}

function op_cf(f) {
    Flags[f] = false;
}

function op_cf_index() {
    op_cf(Math.floor(Math.abs(Reg.I)));
}

function op_isg(r) {
    var n = trunc(Reg[r]);
    var f = Reg[r] - n;
    var s = sprintf("%.5f", f);
    var x = +s.substr(2, 3);
    var y = +s.substr(5, 2);
    if (y === 0) {
        y = 1;
    }
    n += y;
    Reg[r] = n + f;
    if (n > x) {
        PC++;
    }
}

function op_isg_index() {
    op_isg(Math.floor(Math.abs(Reg.I)));
}

function op_ftest(f) {
    if (!Flags[f]) {
        PC++;
    }
}

function op_ftest_index() {
    op_ftest(Math.floor(Math.abs(Reg.I)));
}

function op_mul() {
    if (Stack[0] instanceof Descriptor && Stack[1] instanceof Descriptor) {
        binopm(function(y, x) {
            return g_Matrix[y.label].times(g_Matrix[x.label]);
        });
    } else if (Stack[0] instanceof Descriptor) {
        binopm(function(y, x) {
            return g_Matrix[x.label].timesScalar(Stack[1]);
        });
    } else if (Stack[1] instanceof Descriptor) {
        binopm(function(y, x) {
            return g_Matrix[y.label].timesScalar(Stack[0]);
        });
    } else if (Flags[8]) {
        binopc(function(y, x) {
            return y.mul(x);
        });
    } else {
        binop(function(y, x) { return y * x; });
    }
}

function op_integrate(n) {
    var call = function(n) {
        var r = ReturnStack.length;
        op_gsb(n);
        while (ReturnStack.length > r) {
            step();
        }
    };
    // This is http://mathworld.wolfram.com/SimpsonsRule.html
    var eps = 1e-9;
    var x0 = Stack[1];
    var x1 = Stack[0];
    var d = 0;
    var prev = 0;
    var steps = 32;
    var r;
    while (true) {
        r = 0;
        for (var j = 0; j < steps; j += 2) {
            fill(x0 + (x1-x0)*j/steps);
            call(n);
            r += Stack[0];
            fill(x0 + (x1-x0)*(j+1)/steps);
            call(n);
            r += 4 * Stack[0];
            fill(x0 + (x1-x0)*(j+2)/steps);
            call(n);
            r += Stack[0];
        }
        r *= ((x1-x0)/steps) / 3;
        d = Math.abs(r - prev);
        if (d < eps) {
            break;
        }
        prev = r;
        steps *= 2;
    }
    Stack[3] = x0;
    Stack[2] = x1;
    Stack[1] = d;
    Stack[0] = r;
}

function op_eq() {
    if (Stack[0] === 0) {
        // execute next opcode
    } else {
        PC++;
    }
}

function op_rs() {
    Running = !Running;
    StackLift = OldStackLift;
}

function op_pse() {
    //alert("Unimplemented: PSE");
    // TODO set RunningPause?
    update_display();
    //if (!confirm("pause. keep going?")) {
    //    Running = false;
    //}
    StackLift = OldStackLift;
}

function op_pr() {
    Prgm = !Prgm;
    StackLift = OldStackLift;
}

function op_gsb(n) {
    if (Running) {
        ReturnStack.push(PC);
    } else {
        ReturnStack.push(0);
    }
    op_gto_label(n);
    Running = true;
}

function op_gsb_index() {
    if (Running) {
        ReturnStack.push(PC);
    } else {
        ReturnStack.push(0);
    }
    op_gto_index();
    Running = true;
}

function op_clear_stat() {
    var i;
    for (i in Stack) {
        Stack[i] = 0;
    }
    for (i = 2; i <= 7; i++) {
        Reg[i] = 0;
    }
    StackLift = OldStackLift;
}

function op_rtn() {
    if (ReturnStack.length > 0) {
        PC = ReturnStack.pop();
    } else {
        PC = 0;
    }
    if (PC === 0) {
        Running = false;
    }
}

function op_roll() {
    var t = Stack[0];
    var ti = StackI[0];
    Stack[0] = Stack[1]; StackI[0] = StackI[1];
    Stack[1] = Stack[2]; StackI[1] = StackI[2];
    Stack[2] = Stack[3]; StackI[2] = StackI[3];
    Stack[3] = t;
    StackI[3] = ti;
}

function op_clear_prgm() {
    if (Prgm) {
        Program = [null];
    }
    PC = 0;
}

function op_rollup() {
    var t = Stack[3];
    var ti = StackI[3];
    Stack[3] = Stack[2]; StackI[3] = StackI[2];
    Stack[2] = Stack[1]; StackI[2] = StackI[1];
    Stack[1] = Stack[0]; StackI[1] = StackI[0];
    Stack[0] = t;
    StackI[0] = ti;
}

function op_xy() {
    var t = Stack[0];
    var ti = StackI[0];
    Stack[0] = Stack[1]; StackI[0] = StackI[1];
    Stack[1] = t;
    StackI[1] = ti;
}

function op_clear_reg() {
    for (var i in Reg) {
        Reg[i] = 0;
    }
    StackLift = OldStackLift;
}

function op_rnd() {
    unop(function(x) {
        var factor = Math.pow(10, DisplayDigits + log10int(x));
        return Math.round(x * factor) / factor;
    });
}

function op_back() {
    if (Prgm) {
        if (PC > 0) {
            Program.splice(PC, 1);
            PC--;
        }
        return;
    } else if (Flags[9]) {
        Flags[9] = false;
        if (Stack[0] > MAX) {
            Stack[0] = MAX;
        } else if (Stack[0] < -MAX) {
            Stack[0] = -MAX;
        }
        return;
    } else if (DigitEntry && Entry.length > 0) {
        Entry = Entry.substr(0, Entry.length-1);
    } else {
        op_clx();
        Entry = "";
    }
    NewDigitEntry = true;
}

function op_clear_prefix() {
    Prefix = null;
    var x = Math.abs(Stack[0]);
    var s = sprintf("%.9e", x).replace(".", "").substr(0, 10);
    while (s.length < 10) {
        s += '0';
    }
    update_lcd(s);
    DelayUpdate = 1000;
    StackLift = OldStackLift;
}

function op_clx() {
    Stack[0] = 0;
    NewStackLift = false;
}

function op_enter() {
    push(Stack[0], true);
    // push() only pushes the real part,
    // so copy the imaginary part too
    StackI[0] = StackI[1];
    NewStackLift = false;
}

function op_rand() {
    push(Math.random());
}

function op_lastx() {
    push(LastX);
    StackI[0] = LastXI;
}

function op_to_r() {
    LastX = Stack[0];
    LastXI = StackI[0];
    var t, r;
    if (Flags[8]) {
        t = StackI[0];
        r = Stack[0];
        StackI[0] = r * sin_drg_mode(t);
        Stack[0] = r * cos_drg_mode(t);
    } else {
        t = Stack[1];
        r = Stack[0];
        Stack[1] = r * sin_drg_mode(t);
        Stack[0] = r * cos_drg_mode(t);
    }
}

function op_to_p() {
    LastX = Stack[0];
    LastXI = StackI[0];
    var x, y;
    if (Flags[8]) {
        y = StackI[0];
        x = Stack[0];
        StackI[0] = Math.atan2(y, x) / TrigFactor;
        Stack[0] = Math.sqrt(x*x + y*y);
    } else {
        y = Stack[1];
        x = Stack[0];
        Stack[1] = Math.atan2(y, x) / TrigFactor;
        Stack[0] = Math.sqrt(x*x + y*y);
    }
}

function op_to_hms() {
    unop(function(x) {
        var r = Math.floor(x);
        x -= r;
        x *= 60;
        r += Math.floor(x) / 100;
        x -= Math.floor(x);
        x *= 60;
        r += x / 10000;
        return r;
    });
}

function op_to_h() {
    unop(function(x) {
        var r = Math.floor(x);
        x -= Math.floor(x);
        x *= 100;
        r += Math.floor(x) / 60;
        x -= Math.floor(x);
        x *= 100;
        r += x / 3600;
        return r;
    });
}

function op_to_rad() {
    unop(function(x) { return x * Math.PI / 180; });
}

function op_to_deg() {
    unop(function(x) { return x * 180 / Math.PI; });
}

function op_sub() {
    if (Stack[0] instanceof Descriptor && Stack[1] instanceof Descriptor) {
        binopm(function(y, x) {
            return g_Matrix[y.label].minus(g_Matrix[x.label]);
        });
    } else if (Stack[0] instanceof Descriptor) {
        binopm(function(y, x) {
            var m = g_Matrix[x.label];
            return new Mat(m.rows, m.cols, Stack[1]).minus(m);
        });
    } else if (Stack[1] instanceof Descriptor) {
        binopm(function(y, x) {
            var m = g_Matrix[y.label];
            return m.minus(new Mat(m.rows, m.cols, Stack[0]));
        });
    } else if (Flags[8]) {
        binopc(function(y, x) {
            return y.sub(x);
        });
    } else {
        binop(function(y, x) { return y - x; });
    }
}

function op_re_im() {
    Flags[8] = true;
    var t = StackI[0];
    StackI[0] = Stack[0];
    Stack[0] = t;
}

function op_test(t) {
    var b;
    switch (t) {
        case 0: b = Stack[0] !== 0; break;
        case 1: b = Stack[0]  > 0; break;
        case 2: b = Stack[0]  < 0; break;
        case 3: b = Stack[0] >= 0; break;
        case 4: b = Stack[0] <= 0; break;
        case 5: b = Stack[0] == Stack[1]; break;
        case 6: b = Stack[0] != Stack[1]; break;
        case 7: b = Stack[0]  > Stack[1]; break;
        case 8: b = Stack[0]  < Stack[1]; break;
        case 9: b = Stack[0] >= Stack[1]; break;
    }
    if (!b) {
        PC++;
    }
}

function op_on() {
    On = !On;
}

function op_sto_reg(n) {
    Reg[n] = Stack[0];
}

function op_sto_op_reg(op, n) {
    switch (op) {
        case '+': Reg[n] += Stack[0]; break;
        case '-': Reg[n] -= Stack[0]; break;
        case '*': Reg[n] *= Stack[0]; break;
        case '/': Reg[n] /= Stack[0]; break;
    }
}

function op_sto_index(user) {
    if (Reg.I instanceof Descriptor) {
        op_sto_matrix(Reg.I.label, user);
    } else {
        op_sto_reg(Math.floor(Math.abs(Reg.I)));
    }
}

function op_sto_op_index(op) {
    op_sto_op_reg(op, Math.floor(Math.abs(Reg.I)));
}

function op_sto_matrix(m, user) {
    g_Matrix[m].set(Reg[0], Reg[1], Stack[0]);
    update_lcd(sprintf("%c %2d,%d", "A".charCodeAt(0) + m, Reg[0], Reg[1]));
    DelayUpdate = 1000;
    if (user) {
        if (Reg[1] < g_Matrix[m].cols) {
            Reg[1]++;
        } else if (Reg[0] < g_Matrix[m].rows) {
            Reg[0]++;
            Reg[1] = 1;
        } else {
            Reg[0] = 1;
            Reg[1] = 1;
            if (Running) {
                PC++;
            }
        }
    }
}

function op_sto_matrix_imm(m) {
    var x = Stack[0];
    var y = Stack[1];
    Stack[0] = Stack[2]; StackI[0] = StackI[2];
    Stack[1] = Stack[3]; StackI[1] = StackI[3];
    Stack[2] = Stack[3]; StackI[2] = StackI[3];
    g_Matrix[m].set(y, x, Stack[0]);
}

function op_sto_matrix_all(m) {
    if (Stack[0] instanceof Descriptor) {
        g_Matrix[m] = g_Matrix[Stack[0].label].copy();
    } else {
        for (var j = 1; j <= g_Matrix[m].rows; j++) {
            for (var i = 1; i <= g_Matrix[m].cols; i++) {
                g_Matrix[m].set(j, i, Stack[0]);
            }
        }
    }
}

function op_sto_result() {
    if (Stack[0] instanceof Descriptor) {
        Result = Stack[0].label;
    } else {
        // TODO: error
    }
}

function op_frac() {
    unop(function(x) {
        return x - trunc(x);
    });
}

function op_int() {
    unop(trunc);
}

function op_rcl_reg(r) {
    if (StackLift) {
        push(Reg[r]);
    } else {
        Stack[0] = Reg[r];
        StackI[0] = 0;
        StackLift = true;
    }
}

function op_rcl_op_reg(op, r) {
    switch (op) {
        case '+': Stack[0] += Reg[r]; break;
        case '-': Stack[0] -= Reg[r]; break;
        case '*': Stack[0] *= Reg[r]; break;
        case '/': Stack[0] /= Reg[r]; break;
    }
}

function op_rcl_index(user) {
    if (Reg.I instanceof Descriptor) {
        op_rcl_matrix(Reg.I.label, user);
    } else {
        op_rcl_reg(Math.floor(Math.abs(Reg.I)));
    }
}

function op_rcl_op_index(op) {
    op_rcl_op_reg(op, Math.floor(Math.abs(Reg.I)));
}

function op_rcl_descriptor(m) {
    push(new Descriptor(m));
}

function op_rcl_dim(m) {
    push(g_Matrix[m].rows);
    push(g_Matrix[m].cols);
}

function op_rcl_matrix(m, user) {
    if (StackLift) {
        push(g_Matrix[m].get(Reg[0], Reg[1]));
    } else {
        Stack[0] = g_Matrix[m].get(Reg[0], Reg[1]);
        StackLift = true;
    }
    update_lcd(sprintf("%c %2d,%d", "A".charCodeAt(0) + m, Reg[0], Reg[1]));
    DelayUpdate = 1000;
    if (user) {
        if (Reg[1] < g_Matrix[m].cols) {
            Reg[1]++;
        } else if (Reg[0] < g_Matrix[m].rows) {
            Reg[0]++;
            Reg[1] = 1;
        } else {
            Reg[0] = 1;
            Reg[1] = 1;
            if (Running) {
                PC++;
            }
        }
    }
}

function op_rcl_matrix_imm(m) {
    binop(function(y, x) {
        return g_Matrix[m].get(y, x);
    });
}

function op_rcl_result() {
    push(new Descriptor(Result));
}

function op_user() {
    User = !User;
    Display.set_user(User);
    StackLift = OldStackLift;
}

function op_mem() {
    alert("Unimplemented: MEM");
    StackLift = OldStackLift;
}

function op_fact() {
    unop(function(x) {
        if (x >= 0 && x === Math.floor(x)) {
            if (x > 69) {
                Flags[9] = true;
                return MAX;
            }
            var r = 1;
            while (x > 1) {
                r *= x;
                x -= 1;
            }
            return r;
        } else {
            x += 1;
            var gamma = function(z) {
                // Spouge approximation
                // https://en.wikipedia.org/wiki/Spouge's_approximation
                var kc = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
                var kf = 1.0;
                kc[0] = Math.sqrt(2.0 * Math.PI);
                for (var k = 1; k < 12; k++) {
                    kc[k] = Math.exp(12.0 - k) * Math.pow(12.0 - k, k - 0.5) / kf;
                    kf *= -k;
                }
                var acc = kc[0];
                for (k = 1; k < 12; k++) {
                    acc += kc[k] / (z + k);
                }
                acc *= Math.exp(-(z + 12)) * Math.pow(z + 12, z + 0.5);
                return acc / z;
            };
            if (x <= 0 && x === Math.floor(x)) {
                Flags[9] = true;
                return -MAX;
            }
            if (x < -70.06400563) {
                return 0;
            }
            if (x > 70.95757445) {
                Flags[9] = true;
                return MAX;
            }
            if (x >= -10) {
                return gamma(x);
            } else {
                return Math.PI / (gamma(1-x) * Math.sin(Math.PI * x));
            }
        }
    });
}

function op_mean() {
    push(Reg[5] / Reg[2]);
    push(Reg[3] / Reg[2]);
}

function op_yhat() {
    LastX = Stack[0];
    LastXI = StackI[0];
    var M = Reg[2] * Reg[4] - Reg[3] * Reg[3];
    var N = Reg[2] * Reg[6] - Reg[5] * Reg[5];
    var P = Reg[2] * Reg[7] - Reg[3] * Reg[5];
    push(P / Math.sqrt(M * N));
    push((M * Reg[5] + P * (Reg[2] * LastX - Reg[3])) / (Reg[2] * M));
}

function op_s() {
    var M = Reg[2] * Reg[4] - Reg[3] * Reg[3];
    var N = Reg[2] * Reg[6] - Reg[5] * Reg[5];
    push(Math.sqrt(N / (Reg[2] * (Reg[2] - 1))));
    push(Math.sqrt(M / (Reg[2] * (Reg[2] - 1))));
}

function op_sum() {
    LastX = Stack[0];
    LastXI = StackI[0];
    Reg[2] += 1;
    Reg[3] += Stack[0];
    Reg[4] += Stack[0] * Stack[0];
    Reg[5] += Stack[1];
    Reg[6] += Stack[1] * Stack[1];
    Reg[7] += Stack[0] * Stack[1];
    Stack[0] = Reg[2];
    NewStackLift = false;
}

function op_lr() {
    var M = Reg[2] * Reg[4] - Reg[3] * Reg[3];
    var N = Reg[2] * Reg[6] - Reg[5] * Reg[5];
    var P = Reg[2] * Reg[7] - Reg[3] * Reg[5];
    push(P / M);
    push((M * Reg[5] - P * Reg[3]) / (Reg[2] * M));
}

function op_sumsub() {
    LastX = Stack[0];
    LastXI = StackI[0];
    Reg[2] -= 1;
    Reg[3] -= Stack[0];
    Reg[4] -= Stack[0] * Stack[0];
    Reg[5] -= Stack[1];
    Reg[6] -= Stack[1] * Stack[1];
    Reg[7] -= Stack[0] * Stack[1];
    Stack[0] = Reg[2];
    NewStackLift = false;
}

function op_add() {
    if (Stack[0] instanceof Descriptor && Stack[1] instanceof Descriptor) {
        binopm(function(y, x) {
            return g_Matrix[y.label].plus(g_Matrix[x.label]);
        });
    } else if (Stack[0] instanceof Descriptor) {
        binopm(function(y, x) {
            var m = g_Matrix[x.label];
            return m.plus(new Mat(m.rows, m.cols, Stack[1]));
        });
    } else if (Stack[1] instanceof Descriptor) {
        binopm(function(y, x) {
            var m = g_Matrix[y.label];
            return m.plus(new Mat(m.rows, m.cols, Stack[0]));
        });
    } else if (Flags[8]) {
        binopc(function(y, x) {
            return y.add(x);
        });
    } else {
        binop(function(y, x) { return y + x; });
    }
}

function op_Pyx() {
    if (Stack[0] instanceof Descriptor) {
        var m = Stack[0].label;
        g_Matrix[m] = g_Matrix[m].partition();
    } else {
        binop(function(y, x) {
            var r = 1;
            var t = y - x;
            while (y > t) {
                r *= y;
                y--;
            }
            return r;
        });
    }
}

function op_Cyx() {
    if (Stack[0] instanceof Descriptor) {
        var m = Stack[0].label;
        g_Matrix[m] = g_Matrix[m].unpartition();
    } else {
        binop(function(y, x) {
            var r = 1;
            var t = y - x;
            while (y > t) {
                r *= y;
                y--;
            }
            while (x > 1) {
                r /= x;
                x--;
            }
            return r;
        });
    }
}

function op_input(c) {
    if (!DigitEntry) {
        if (StackLift) {
            push("");
        }
        Entry = "";
    }
    if (Entry.length === 0) {
        switch (c) {
            case 'e':
                Entry = "1";
                break;
            case '.':
                Entry = "0";
                break;
        }
    }
    Entry = Entry + c;
    NewDigitEntry = true;
}

function decode_matrix(k) {
    Prefix = function(k) {
        switch (k) {
            case '0': return new Opcode(new OpcodeInfo([42,16,0], op_matrix_clear));
            case '1': return new Opcode(new OpcodeInfo([42,16,1], op_matrix_home));
            case '2': return new Opcode(new OpcodeInfo([42,16,2], op_matrix_complex2));
            case '3': return new Opcode(new OpcodeInfo([42,16,3], op_matrix_complex3));
            case '4': return new Opcode(new OpcodeInfo([42,16,4], op_matrix_transpose));
            case '5': return new Opcode(new OpcodeInfo([42,16,5], op_matrix_transmul));
            case '6': return new Opcode(new OpcodeInfo([42,16,6], op_matrix_residual));
            case '7': return new Opcode(new OpcodeInfo([42,16,7], op_matrix_norm));
            case '8': return new Opcode(new OpcodeInfo([42,16,8], op_matrix_normf));
            case '9': return new Opcode(new OpcodeInfo([42,16,9], op_matrix_det));
        }
    };
    return null;
}

function decode_fix(k) {
    Prefix = function(k) {
        if (k >= '0' && k <= '9') {
            var i = Number(k);
            return new Opcode(new OpcodeInfo([42,7,i]), function() { op_fix(i); });
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([42,7,25]), op_fix_index);
        }
    };
    return null;
}

function decode_sci(k) {
    Prefix = function(k) {
        if (k >= '0' && k <= '9') {
            var i = Number(k);
            return new Opcode(new OpcodeInfo([42,8,i]), function() { op_sci(i); });
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([42,8,25]), op_sci_index);
        }
    };
    return null;
}

function decode_eng(k) {
    Prefix = function(k) {
        if (k >= '0' && k <= '9') {
            var i = Number(k);
            return new Opcode(new OpcodeInfo([42,9,i]), function() { op_eng(i); });
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([42,9,25]), op_eng_index);
        }
    };
    return null;
}

function decode_solve(k) {
    var f = 1;
    Prefix = function(k) {
        var i;
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            i = Number(k) / f;
            return new Opcode(new OpcodeInfo([42,10,i]), function() { op_solve(i); });
        } else {
            i = "qE)^\\".indexOf(k);
            if (i >= 0) {
                return new Opcode(new OpcodeInfo([42,10,11+i]), function() { op_solve(11+i); });
            }
        }
    };
    return null;
}

function decode_lbl(k) {
    var f = 1;
    Prefix = function(k) {
        var i;
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            i = Number(k) / f;
            return new Opcode(new OpcodeInfo([42,21,i]), function() {});
        } else {
            i = "qE)^\\".indexOf(k);
            if (i >= 0) {
                return new Opcode(new OpcodeInfo([42,21,11+i]), function() {});
            }
        }
    };
    return null;
}

function decode_gto() {
    var f = 1;
    var immediate = false;
    var n = 0;
    var i = 0;
    Prefix = function(k) {
        var x;
        if (k === '_' && n === 0) {
            immediate = true;
            Prefix = OldPrefix;
            return null;
        } else if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            x = Number(k);
            if (immediate) {
                n = n * 10 + x;
                i++;
                if (i === 3) {
                    return new Opcode(new OpcodeInfo([22], null, false), function() { op_gto_immediate(n); });
                }
                Prefix = OldPrefix;
                return null;
            } else {
                x /= f;
                return new Opcode(new OpcodeInfo([22,x]), function() { op_gto_label(x); });
            }
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([22,25]), function() { op_gto_index(); });
        } else {
            x = "qE)^\\".indexOf(k);
            if (x >= 0) {
                return new Opcode(new OpcodeInfo([22,11+x]), function() { op_gto_label(11+x); });
            }
        }
    };
    return null;
}

function decode_hyp(k) {
    Prefix = function(k) {
        switch (k) {
            case 's': return new Opcode(new OpcodeInfo([42,22,23], op_sinh));
            case 'c': return new Opcode(new OpcodeInfo([42,22,24], op_cosh));
            case 't': return new Opcode(new OpcodeInfo([42,22,25], op_tanh));
        }
    };
    return null;
}

function decode_ahyp(k) {
    Prefix = function(k) {
        switch (k) {
            case 's': return new Opcode(new OpcodeInfo([43,22,23], op_asinh));
            case 'c': return new Opcode(new OpcodeInfo([43,22,24], op_acosh));
            case 't': return new Opcode(new OpcodeInfo([43,22,25], op_atanh));
        }
    };
    return null;
}

function decode_dim(k) {
    Prefix = function(k) {
        if (k == 'c') {
            return new Opcode(new OpcodeInfo([42,23,24]), op_dim_i);
        } else {
            var i = "qE)^\\".indexOf(k);
            if (i >= 0) {
                return new Opcode(new OpcodeInfo([42,23,11+i]), function() { op_dim(i); });
            }
        }
    };
    return null;
}

function decode_result(k) {
    Prefix = function(k) {
        var i = "qE)^\\".indexOf(k);
        if (i >= 0) {
            return new Opcode(new OpcodeInfo([42,26,11+i]), function() { op_result(i); });
        }
    };
    return null;
}

function decode_xchg(k) {
    var f = 0;
    Prefix = function(k) {
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            var i = Number(k) + f;
            return new Opcode(new OpcodeInfo([42,4,i]), function() { op_xchg(i); });
        } else if (k === 'c') {
            return new Opcode(new OpcodeInfo([42,4,24]), op_xchg_index);
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([42,4,25]), function() { op_xchg('I'); });
        }
    };
    return null;
}

function decode_sf(k) {
    Prefix = function(k) {
        if (k >= '0' && k <= '9') {
            var i = Number(k);
            return new Opcode(new OpcodeInfo([43,4,i]), function() { op_sf(i); });
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([43,4,25]), op_sf_index);
        }
    };
    return null;
}

function decode_dse(k) {
    var f = 0;
    Prefix = function(k) {
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            var i = Number(k) + f;
            return new Opcode(new OpcodeInfo([42,5,i]), function() { op_dse(i); });
        } else if (k === 'c') {
            return new Opcode(new OpcodeInfo([42,5,24]), op_dse_index);
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([42,5,25]), function() { op_dse('I'); });
        }
    };
    return null;
}

function decode_cf(k) {
    Prefix = function(k) {
        if (k >= '0' && k <= '9') {
            var i = Number(k);
            return new Opcode(new OpcodeInfo([43,5,i]), function() { op_cf(i); });
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([43,5,25]), op_cf_index);
        }
    };
    return null;
}

function decode_isg(k) {
    var f = 0;
    Prefix = function(k) {
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            var i = Number(k) + f;
            return new Opcode(new OpcodeInfo([42,6,i]), function() { op_isg(i); });
        } else if (k === 'c') {
            return new Opcode(new OpcodeInfo([42,6,24]), op_isg_index);
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([42,6,25]), function() { op_isg('I'); });
        }
    };
    return null;
}

function decode_ftest(k) {
    Prefix = function(k) {
        if (k >= '0' && k <= '9') {
            var i = Number(k);
            return new Opcode(new OpcodeInfo([43,6,i]), function() { op_ftest(i); });
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([43,6,25]), op_ftest_index);
        }
    };
    return null;
}

function decode_integrate(k) {
    var f = 1;
    Prefix = function(k) {
        var i;
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            i = Number(k) / f;
            return new Opcode(new OpcodeInfo([42,20,i]), function() { op_integrate(i); });
        } else {
            i = "qE)^\\".indexOf(k);
            if (i >= 0) {
                return new Opcode(new OpcodeInfo([42,20,11+i]), function() { op_integrate(11+i); });
            }
        }
    };
    return null;
}

function decode_gsb() {
    var f = 1;
    Prefix = function(k) {
        var i;
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            i = Number(k) / f;
            return new Opcode(new OpcodeInfo([32,i]), function() { op_gsb(i); });
        } else if (k === 't') {
            return new Opcode(new OpcodeInfo([32,25]), function() { op_gsb_index(); });
        } else {
            i = "qE)^\\".indexOf(k);
            if (i >= 0) {
                return new Opcode(new OpcodeInfo([32,11+i]), function() { op_gsb(11+i); });
            }
        }
    };
    return null;
}

function decode_test(k) {
    Prefix = function(k) {
        if (k >= '0' && k <= '9') {
            var i = Number(k);
            return new Opcode(new OpcodeInfo([43,30,i]), function() { op_test(i); });
        }
    };
    return null;
}

function decode_f() {
    Shift = 1;
    Display.clear_shift();
    Display.set_shift("f");
    return null;
}

function decode_g() {
    Shift = 2;
    Display.clear_shift();
    Display.set_shift("g");
    return null;
}

function decode_sto(k) {
    var f = 0;
    var op = null;
    var g = false;
    Prefix = function(k) {
        var i, u;
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            i = Number(k) + f;
            if (op !== null) {
                return new Opcode(new OpcodeInfo([44,OpcodeIndex[op],i]), function() { op_sto_op_reg(op, i); });
            } else {
                return new Opcode(new OpcodeInfo([44,i]), function() { op_sto_reg(i); });
            }
        } else if (k === '+' || k === '-' || k === '*' || k === '/') {
            op = k;
            Prefix = OldPrefix;
            return null;
        } else if (k === '_') {
            Prefix = function(k) {
                i = "qE)^\\".indexOf(k);
                if (i >= 0) {
                    return new Opcode(new OpcodeInfo([44,16,11+i]), function() { op_sto_matrix_all(i); });
                }
            };
            return null;
        } else if (k === 'c') {
            if (op !== null) {
                return new Opcode(new OpcodeInfo([44,OpcodeIndex[op],24]), function() { op_sto_op_index(op); });
            } else {
                u = User;
                return new Opcode(new OpcodeInfo([44,24]), function() { op_sto_index(u); });
            }
        } else if (k === 't') {
            if (op !== null) {
                return new Opcode(new OpcodeInfo([44,OpcodeIndex[op],25]), function() { op_sto_op_reg(op, 'I'); });
            } else {
                return new Opcode(new OpcodeInfo([44,25]), function() { op_sto_reg('I'); });
            }
        } else if (k === 'e') {
            return new Opcode(new OpcodeInfo([44,26]), op_sto_result);
        } else if (k === 'g') {
            g = true;
            Prefix = OldPrefix;
            return null;
        } else {
            i = "qE)^\\".indexOf(k);
            if (i >= 0) {
                if (g) {
                    return new Opcode(new OpcodeInfo([44,43,11+i]), function() { op_sto_matrix_imm(i); });
                } else {
                    u = User; // capture current value
                    return new Opcode(new OpcodeInfo([44,11+i], null, true, User), function() { op_sto_matrix(i, u); });
                }
            }
        }
    };
    return null;
}

function decode_rcl(k) {
    var f = 0;
    var op = null;
    var g = false;
    Prefix = function(k) {
        var i, u;
        if (k === '.') {
            f = 10;
            Prefix = OldPrefix;
            return null;
        } else if (k >= '0' && k <= '9') {
            i = Number(k) + f;
            if (op !== null) {
                return new Opcode(new OpcodeInfo([45,OpcodeIndex[op],i]), function() { op_rcl_op_reg(op, i); });
            } else {
                return new Opcode(new OpcodeInfo([45,i]), function() { op_rcl_reg(i); });
            }
        } else if (k === '+' || k === '-' || k === '*' || k === '/') {
            op = k;
            Prefix = OldPrefix;
            return null;
        } else if (k === '_') {
            Prefix = function(k) {
                i = "qE)^\\".indexOf(k);
                if (i >= 0) {
                    return new Opcode(new OpcodeInfo([45,16,11+i]), function() { op_rcl_descriptor(i); });
                }
            };
            return null;
        } else if (k === 's') {
            Prefix = function(k) {
                i = "qE)^\\".indexOf(k);
                if (i >= 0) {
                    return new Opcode(new OpcodeInfo([45,23,11+i]), function() { op_rcl_dim(i); });
                }
            };
            return null;
        } else if (k === 'c') {
            if (op !== null) {
                return new Opcode(new OpcodeInfo([45,OpcodeIndex[op],24]), function() { op_rcl_op_index(op); });
            } else {
                u = User;
                return new Opcode(new OpcodeInfo([45,24]), function() { op_rcl_index(u); });
            }
        } else if (k === 't') {
            if (op !== null) {
                return new Opcode(new OpcodeInfo([45,OpcodeIndex[op],25]), function() { op_rcl_op_reg(op, 'I'); });
            } else {
                return new Opcode(new OpcodeInfo([45,25]), function() { op_rcl_reg('I'); });
            }
        } else if (k === 'e') {
            return new Opcode(new OpcodeInfo([45,26]), op_rcl_result);
        } else if (k === 'g') {
            g = true;
            Prefix = OldPrefix;
            return null;
        } else {
            i = "qE)^\\".indexOf(k);
            if (i >= 0) {
                if (g) {
                    return new Opcode(new OpcodeInfo([45,43,11+i]), function() { op_rcl_matrix_imm(i); });
                } else {
                    u = User; // capture current value
                    return new Opcode(new OpcodeInfo([45,11+i], null, true, User), function() { op_rcl_matrix(i, u); });
                }
            }
        }
    };
    return null;
}

function decode(k) {
    var d = null;
    var s = Shift;
    Shift = -1;
    if (Prefix != null) {
        d = Prefix;
    } else if (typeof(CharTable[k]) === "object") {
        if (User && "qE)^\\".indexOf(k) >= 0) {
            switch (s) {
                case 0: s = 1; break;
                case 1: s = 0; break;
            }
        }
        d = CharTable[k][s];
    } else {
        d = CharTable[k];
        if (d === undefined) {
            return null;
        }
    }
    OldPrefix = Prefix;
    Prefix = null;
    var r = d(k);
    if (Shift === -1) {
        Shift = 0;
        Display.clear_shift();
    }
    return r;
}

function step() {
    if (PC === 0) {
        PC = 1;
    }
    if (PC < Program.length) {
        //console.log("PC", PC, Program[PC].info.defn);
        var p = PC;
        PC++;
        try {
            Program[p].exec();
        } catch (e) {
            Running = false;
            if (e.name === "CalcError") {
                update_lcd("Error " + e.code);
                DelayUpdate = -1;
            } else {
                throw e;
            }
        }
    } else {
        op_rtn();
    }
}

function run() {
    RunTimer = null;
    if (!Running) {
        alert("run() called when not Running");
        return;
    }
    step();
    if (Running) {
        RunTimer = setTimeout(run, 0);
    } else {
        update_display();
    }
}

function delay_update_timeout() {
    if (!TemporaryDisplay) {
        update_display();
    }
    DisplayTimeout = 0;
}

function key_up() {
    if (TemporaryDisplay) {
        TemporaryDisplay = false;
        if (DisplayTimeout === 0) {
            delay_update_timeout();
        }
    }
}

function key_down(k, override) {
    if (!On && k != '\x1b') {
        return;
    }
    if (DisableKeys && !override) {
        return;
    }
    var op = decode(k);
    if (op === undefined) {
        //console.log("undefined decode: "+k);
        return;
    }
    if (Running) {
        clearTimeout(RunTimer);
        Running = false;
        return;
    }
    if (op !== null) {
        try {
            if (Prgm && op.info.programmable) {
                PC++;
                Program.splice(PC, 0, op);
            } else {
                op.exec();
                if (Running) {
                    RunTimer = setTimeout(run, 0);
                }
            }
        } catch (e) {
            if (e.name === "CalcError") {
                update_lcd("Error " + e.code);
                DelayUpdate = -1;
            } else {
                throw e;
            }
        }
    }
    if (DelayUpdate === 0) {
        update_display();
    } else {
        if (DelayUpdate > 0) {
            TemporaryDisplay = true;
            DisplayTimeout = setTimeout(delay_update_timeout, DelayUpdate);
        }
        DelayUpdate = 0;
    }
}

function key(k, override) {
    key_down(k, override);
    key_up();
}

function paste(s) {
    if (!Prgm) {
        new Opcode(null, function() { push(s); }).exec();
        update_display();
    }
}

function init() {
    for (var i = 0; i < Reg.length; i++) {
        Reg[i] = 0;
    }
    update_display();
}
