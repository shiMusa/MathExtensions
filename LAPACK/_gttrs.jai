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

    // Solve each RHS column independently (matches DGTTS2/CGTTS2 semantics).
    if iTRANS == 0 {
        // Solve A * X = B.
        for j:1..nrhs {
            // Solve L * x = b.
            for i:1..n-1 {
                if fget(ipiv, i) == i {
                    fset(b, i + 1, j, fget(b.*, i + 1, j) - fget(dl, i)*fget(b.*, i, j));
                } else {
                    // Pivot implies IPIV(i) == i+1
                    temp_i := fget(b.*, i, j);
                    temp_ip := fget(b.*, i + 1, j);
                    fset(b, i, j, temp_ip);
                    fset(b, i + 1, j, temp_i - fget(dl, i)*temp_ip);
                }
            }

            // Solve U * x = b.
            fset(b, n, j, fget(b.*, n, j) / fget(d, n));
            if n > 1 {
                fset(b, n - 1, j,
                    (fget(b.*, n - 1, j) - fget(du, n - 1)*fget(b.*, n, j)) / fget(d, n - 1));
            }
            if n > 2 {
                for <i:n - 2..1 {
                    fset(b, i, j,
                        (fget(b.*, i, j)
                         - fget(du, i)*fget(b.*, i + 1, j)
                         - fget(du2, i)*fget(b.*, i + 2, j)) / fget(d, i));
                }
            }
        }
    } else if iTRANS == 1 {
        // Solve A^T * X = B.
        for j:1..nrhs {
            // Solve U^T * x = b.
            fset(b, 1, j, fget(b.*, 1, j) / fget(d, 1));
            if n > 1 {
                fset(b, 2, j,
                    (fget(b.*, 2, j) - fget(du, 1)*fget(b.*, 1, j)) / fget(d, 2));
            }
            if n > 2 {
                for i:3..n {
                    fset(b, i, j,
                        (fget(b.*, i, j)
                         - fget(du, i - 1)*fget(b.*, i - 1, j)
                         - fget(du2, i - 2)*fget(b.*, i - 2, j)) / fget(d, i));
                }
            }

            // Solve L^T * x = b.
            if n > 1 {
                for <i:n - 1..1 {
                    if fget(ipiv, i) == i {
                        fset(b, i, j, fget(b.*, i, j) - fget(dl, i)*fget(b.*, i + 1, j));
                    } else {
                        // Pivot implies IPIV(i) == i+1
                        temp_ip1 := fget(b.*, i + 1, j);
                        fset(b, i + 1, j, fget(b.*, i, j) - fget(dl, i)*temp_ip1);
                        fset(b, i, j, temp_ip1);
                    }
                }
            }
        }
    } else {
        // Solve A^H * X = B (complex only, ITRANS == 2).
        for j:1..nrhs {
            // Solve U^H * x = b.
            fset(b, 1, j, fget(b.*, 1, j) / conjugate(fget(d, 1)));
            if n > 1 {
                fset(b, 2, j,
                    (fget(b.*, 2, j) - conjugate(fget(du, 1))*fget(b.*, 1, j)) / conjugate(fget(d, 2)));
            }
            if n > 2 {
                for i:3..n {
                    fset(b, i, j,
                        (fget(b.*, i, j)
                         - conjugate(fget(du, i - 1))*fget(b.*, i - 1, j)
                         - conjugate(fget(du2, i - 2))*fget(b.*, i - 2, j)) / conjugate(fget(d, i)));
                }
            }

            // Solve L^H * x = b.
            if n > 1 {
                for <i:n - 1..1 {
                    if fget(ipiv, i) == i {
                        fset(b, i, j, fget(b.*, i, j) - conjugate(fget(dl, i))*fget(b.*, i + 1, j));
                    } else {
                        temp_ip1 := fget(b.*, i + 1, j);
                        fset(b, i + 1, j, fget(b.*, i, j) - conjugate(fget(dl, i))*temp_ip1);
                        fset(b, i, j, temp_ip1);
                    }
                }
            }
        }
    }
}

