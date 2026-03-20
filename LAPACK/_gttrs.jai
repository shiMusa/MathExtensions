// ? LAPACK *gttrs()
// ? http://www.netlib.org/lapack/explore-html/
//
// ? original comment sgttrs/cgttrs/dgttrs/zgttrs
// *GTTRS solves one of the systems of equations
//    A*X = B  or  A**T*X = B,
// with a tridiagonal matrix A using the LU factorization computed
// by *GTTRF.
//
sgttrs :: #bake_arguments _gttrs(RT=float32, MT=float32);
dgttrs :: #bake_arguments _gttrs(RT=float64, MT=float64);
cgttrs :: #bake_arguments _gttrs(RT=float32, MT=C32);
zgttrs :: #bake_arguments _gttrs(RT=float64, MT=C64);

#scope_file
_gttrs_forward_substitution :: (
    $$itrans: int,
    $$n: int,
    dl: $DL/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$nrhs: int,
    $MT: Type
) #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    // This helper isn't currently used; kept for future factoring.
    pass;
}

_gttrs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGTTRS"
    } else if MT == float64 {
        "DGTTRS"
    } else if MT == C32 {
        "CGTTRS"
    } else {
        "ZGTTRS"
    };

    info.* = 0;

    notran : bool = trans == #char "N" || trans == #char "n";
    if !notran &&
       !(trans == #char "T" || trans == #char "t" || trans == #char "C" || trans == #char "c") {
        info.* = -1;
        xerbla(routine_name, -info.*);
        return;
    } else if n < 0 {
        info.* = -2;
        xerbla(routine_name, -info.*);
        return;
    } else if nrhs < 0 {
        info.* = -3;
        xerbla(routine_name, -info.*);
        return;
    } else if ldb < max(n, 1) {
        info.* = -10;
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 || nrhs == 0 {
        return;
    }

    // Decode TRANS
    iTRANS : int = 0;
    if MT == float32 || MT == float64 {
        if notran {
            iTRANS = 0;
        } else {
            iTRANS = 1; // For real matrices, 'C' behaves like transpose
        }
    } else {
        if notran {
            iTRANS = 0;
        } else if trans == #char "T" || trans == #char "t" {
            iTRANS = 1;
        } else {
            iTRANS = 2; // conjugate-transpose
        }
    }

    // Solve using the LAPACK helper (*GTTS2).
    #if MT == float32 {
        sgtts2(iTRANS, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
    } else if MT == float64 {
        dgtts2(iTRANS, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
    } else if MT == C32 {
        cgtts2(iTRANS, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
    } else {
        zgtts2(iTRANS, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
    }
}

