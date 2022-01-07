
#load "xerbla.jai";
#load "i_amax.jai";

#load "Level1/_axpy.jai";
#load "Level1/_copy.jai";
#load "Level1/_dotc.jai";
#load "Level1/_dotu.jai";
#load "Level1/_scal.jai";
#load "Level1/_rot.jai";
#load "Level1/_swap.jai";
#load "Level1/_asum.jai";
#load "Level1/_abs.jai";
#load "Level1/_rotm.jai";
#load "Level1/_rotmg.jai";
#load "Level1/_nrm2.jai";

#load "Level2/_gbmv.jai";
#load "Level2/_gemv.jai";
#load "Level2/_geru.jai";
#load "Level2/_hbmv.jai";
#load "Level2/_hemv.jai";
#load "Level2/_her.jai";
#load "Level2/_her2.jai";
#load "Level2/_hpmv.jai";
#load "Level2/_hpr.jai";
#load "Level2/_hpr2.jai";
#load "Level2/_tbmv.jai";
#load "Level2/_tbsv.jai";
#load "Level2/_tpmv.jai";
#load "Level2/_tpsv.jai";
#load "Level2/_trmv.jai";
#load "Level2/_trsv.jai";

#load "Level3/_gemm.jai";
#load "Level3/_hemm.jai";
#load "Level3/_her2k.jai";
#load "Level3/_herk.jai";
#load "Level3/_trmm.jai";
#load "Level3/_trsm.jai";

#load "Testing/cblat1.jai";


// Regex to replace Fortran comple syntax (a,b) with C32.{a,b}: 
//                 \(([^\)]*)\)   ->   C32.{$1}


// ? tests if val equals any of lower_case_chars
char_eq :: inline ($$val: s64, $lower_case_chars: ..s64) -> bool {
    #if lower_case_chars.count == 1 {
        return val == lower_case_chars[0] || val == lower_case_chars[0]-32;
    } else {
        for_unroll(0, lower_case_chars.count-1, #code {
            if lower_case_chars[it] == val || lower_case_chars[it] == val+32 then return true;
        });
        return false;
    }
}

err :: (msg: string, args: ..Any, $exit_program: bool=false) {
    #if exit_program {
        defer assert(false);
    }
    s := sprint("ERROR: %", msg);
    defer free(s);
    println(s, ..args, colo=.FG_RED);
}

data_array_2d :: ($T: Type, $d1: int, $d2: int, vals: ..$U, $row_major: bool=true) -> [d1][d2]T {
    res : [d1][d2]T;
    i1 := 0;
    i2 := 0;
    for vals {
        res[i1][i2] = ncast(T,it);

        #if row_major {
            i2 += 1;
            if i2 == d2 {
                i2 = 0;
                i1 += 1;
            }
        } else {
            i1 += 1;
            if i1 == d1 {
                i1 = 0;
                i2 += 1;
            }
        }
    }
    return res;
}
data_array_3d :: ($T: Type, $d1: int, $d2: int, $d3: int, vals: ..$U, $row_major: bool=true) -> [d1][d2][d3]T {
    res : [d1][d2][d3]T;
    i1 := 0;
    i2 := 0;
    i3 := 0;
    for vals {
        res[i1][i2][i3] = ncast(T,it);
        
        #if row_major {
            i3 += 1;
            if i3 == d3 {
                i3 = 0;
                i2 += 1;
            }
            if i2 == d2 {
                i2 = 0;
                i1 += 1;
            }
        } else {
            i1 += 1;
            if i1 == d1 {
                i1 = 0;
                i2 += 1;
            }
            if i2 == d2 {
                i2 = 0;
                i3 += 1;
            }
        }
    }
    return res;
}

// Anderson E. (2017)
// Algorithm 978: Safe Scaling in the Level 1 BLAS
// ACM Trans Math Softw 44:1--28
// https://doi.org/10.1145/3061665
EPS_32    : float32 : 0.5960464478E-07;
ULP_32    : float32 : 0.1192092896E-06;
SAFMIN_32 : float32 : 0.1175494351E-37;
SAFMAX_32 : float32 : 0.8507059173E+38;
SMLNUM_32 : float32 : 0.9860761315E-31;
BIGNUM_32 : float32 : 0.1014120480E+32;
RTMIN_32  : float32 : 0.3140184864E-15;
RTMAX_32  : float32 : 0.3184525782E+16;
// Blue's scaling constants
TSML_32 : float32 : 0.1084202172E-18;
TBIG_32 : float32 : 0.4503599627E+16;
SSML_32 : float32 : 0.3777893186E+23;
SBIG_32 : float32 : 0.1323488980E-22;
// 64-bit
EPS_64    : float64 : 0.11102230246251565404E-015;
ULP_64    : float64 : 0.22204460492503130808E-015;
SAFMIN_64 : float64 : 0.22250738585072013831E-307;
SAFMAX_64 : float64 : 0.44942328371557897693E+308;
SMLNUM_64 : float64 : 0.10020841800044863890E-291;
BIGNUM_64 : float64 : 0.99792015476735990583E+292;
RTMIN_64  : float64 : 0.10010415475915504622E-145;
RTMAX_64  : float64 : 0.99895953610111751404E+146;
// Blue's scaling constants
TSML_64 : float64 : 0.14916681462400413487E-153;
TBIG_64 : float64 : 0.19979190722022350281E+147;
SSML_64 : float64 : 0.44989137945431963828E+162;
SBIG_64 : float64 : 0.11113793747425387417E-161;



test_lapack :: () {
    println_push("LAPACK", color = .FG_WHITE);

    {
        println_push("char_eq", color = .FG_GREEN);
        val := #char "a";
        VAL := #char "A";
        println("a = %", val);
        println("A = %", VAL);
        println("\"a\" == \"a\" == \"A\": %, %", char_eq(val, #char"a"), char_eq(VAL, #char"a"));
    }

    // test_lapack_axpy();
    // test_lapack_copy();
    // test_lapack_dotc();
    // test_lapack_dotu();
    // test_lapack_scal();
    // test_lapack_rot();
    // test_lapack_swap();
    // test_lapack_asum();
}