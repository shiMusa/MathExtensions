// ? LAPACK *gbtrs()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbtrs/cgbtrs/dgbtrs/zgbtrs
// *GBTRS solves a system of linear equations with a general band
// matrix using the LU factorization computed by *GBTRF.

sgbtrs :: #bake_arguments _gbtrs(MT=float32);
dgbtrs :: #bake_arguments _gbtrs(MT=float64);
cgbtrs :: #bake_arguments _gbtrs(MT=C32);
zgbtrs :: #bake_arguments _gbtrs(MT=C64);

#scope_file
_gbtrs_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    #if CHECKS {
        assert(ld >= 1);
        assert(row >= 1);
        assert(col >= 1);
        assert(count >= 0);
    }

    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld + (row - 1);
    view.count = count;
    return view;
}

_gbtrs_tbsv :: (
    $$uplo: u8,
    $$trans: u8,
    $$diag: u8,
    $$n: int,
    $$k: int,
    a: $A/ArrayType2D,
    $$lda: int,
    x: *$X/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT;
} {
    #if MT == float32 {
        stbsv(uplo, trans, diag, n, k, a, lda, x, 1);
    } else if MT == float64 {
        dtbsv(uplo, trans, diag, n, k, a, lda, x, 1);
    } else if MT == C32 {
        ctbsv(uplo, trans, diag, n, k, a, lda, x, 1);
    } else {
        ztbsv(uplo, trans, diag, n, k, a, lda, x, 1);
    }
}

_gbtrs_swap_rows :: (
    $$nrhs: int,
    b: *$B/ArrayType2D,
    $$row1: int,
    $$row2: int
) {
    if row1 == row2 || nrhs <= 0 then return;

    for rhs:1..nrhs {
        temp := fget(b.*, row1, rhs);
        fset(b, row1, rhs, fget(b.*, row2, rhs));
        fset(b, row2, rhs, temp);
    }
}

_gbtrs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGBTRS"
    } else if MT == float64 {
        "DGBTRS"
    } else if MT == C32 {
        "CGBTRS"
    } else {
        "ZGBTRS"
    };

    kd, l, lm, rhs : int = ---;
    lnoti, notran : bool = ---;
    temp : MT = ---;

    info.* = 0;
    notran = lsame(trans, #char "N");
    if !notran && !lsame(trans, #char "T") && !lsame(trans, #char "C") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if kl < 0 {
        info.* = -3;
    } else if ku < 0 {
        info.* = -4;
    } else if nrhs < 0 {
        info.* = -5;
    } else if ldab < 2*kl + ku + 1 {
        info.* = -7;
    } else if ldb < max(1, n) {
        info.* = -10;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 || nrhs == 0 then return;

    kd = ku + kl + 1;
    lnoti = kl > 0;

    if notran {
        if lnoti {
            for j:1..n-1 {
                lm = min(kl, n-j);
                l = fget(ipiv, j);
                if l != j {
                    _gbtrs_swap_rows(nrhs, b, l, j);
                }

                if lm > 0 {
                    for rhs:1..nrhs {
                        temp = fget(b.*, j, rhs);
                        for i:1..lm {
                            fset(b, j + i, rhs, fget(b.*, j + i, rhs) - fget(ab, kd + i, j)*temp);
                        }
                    }
                }
            }
        }

        for rhs:1..nrhs {
            x_view := _gbtrs_matrix_column(b.*, ldb, 1, rhs, n);
            _gbtrs_tbsv(#char "U", #char "N", #char "N", n, kl + ku, ab, ldab, *x_view, MT);
        }
    } else if lsame(trans, #char "T") {
        for rhs:1..nrhs {
            x_view := _gbtrs_matrix_column(b.*, ldb, 1, rhs, n);
            _gbtrs_tbsv(#char "U", #char "T", #char "N", n, kl + ku, ab, ldab, *x_view, MT);
        }

        if lnoti {
            for <j:1..n-1 {
                lm = min(kl, n-j);
                if lm > 0 {
                    for rhs:1..nrhs {
                        temp = fget(b.*, j, rhs);
                        for i:1..lm {
                            temp -= fget(b.*, j + i, rhs) * fget(ab, kd + i, j);
                        }
                        fset(b, j, rhs, temp);
                    }
                }

                l = fget(ipiv, j);
                if l != j {
                    _gbtrs_swap_rows(nrhs, b, l, j);
                }
            }
        }
    } else {
        for rhs:1..nrhs {
            x_view := _gbtrs_matrix_column(b.*, ldb, 1, rhs, n);
            _gbtrs_tbsv(#char "U", #char "C", #char "N", n, kl + ku, ab, ldab, *x_view, MT);
        }

        if lnoti {
            for <j:1..n-1 {
                lm = min(kl, n-j);
                if lm > 0 {
                    for rhs:1..nrhs {
                        temp = fget(b.*, j, rhs);
                        for i:1..lm {
                            temp -= fget(b.*, j + i, rhs) * conjugate(fget(ab, kd + i, j));
                        }
                        fset(b, j, rhs, temp);
                    }
                }

                l = fget(ipiv, j);
                if l != j {
                    _gbtrs_swap_rows(nrhs, b, l, j);
                }
            }
        }
    }
}
