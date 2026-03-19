// ? LAPACK *gtsv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgtsv/cgtsv/dgtsv/zgtsv
// *GTSV computes the solution to a tridiagonal system by Gaussian
// elimination with partial pivoting.

sgtsv :: (
    $$n: int,
    $$nrhs: int,
    dl: *$DL/ArrayType1D,
    d: *$D/ArrayType1D,
    du: *$DU/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gtsv_real(n, nrhs, dl, d, du, b, ldb, info, RT=float32);
}

dgtsv :: (
    $$n: int,
    $$nrhs: int,
    dl: *$DL/ArrayType1D,
    d: *$D/ArrayType1D,
    du: *$DU/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gtsv_real(n, nrhs, dl, d, du, b, ldb, info, RT=float64);
}

cgtsv :: (
    $$n: int,
    $$nrhs: int,
    dl: *$DL/ArrayType1D,
    d: *$D/ArrayType1D,
    du: *$DU/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gtsv_complex(n, nrhs, dl, d, du, b, ldb, info, RT=float32, MT=C32);
}

zgtsv :: (
    $$n: int,
    $$nrhs: int,
    dl: *$DL/ArrayType1D,
    d: *$D/ArrayType1D,
    du: *$DU/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gtsv_complex(n, nrhs, dl, d, du, b, ldb, info, RT=float64, MT=C64);
}

#scope_file
_gtsv_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_gtsv_real :: (
    $$n: int,
    $$nrhs: int,
    dl: *$DL/ArrayType1D,
    d: *$D/ArrayType1D,
    du: *$DU/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           DL.data_type == RT &&
           D.data_type == RT &&
           DU.data_type == RT &&
           B.data_type == RT;
} {
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { "SGTSV " } else { "DGTSV " };

    fact, temp : RT = ---;

    info.* = 0;
    if n < 0 {
        info.* = -1;
    } else if nrhs < 0 {
        info.* = -2;
    } else if ldb < max(1, n) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    if nrhs == 1 {
        for i:1..n - 2 {
            if abs(fget(d.*, i)) >= abs(fget(dl.*, i)) {
                if fget(d.*, i) != zero_r {
                    fact = fget(dl.*, i) / fget(d.*, i);
                    fset(d, i + 1, fget(d.*, i + 1) - fact*fget(du.*, i));
                    fset(b, i + 1, 1, fget(b.*, i + 1, 1) - fact*fget(b.*, i, 1));
                } else {
                    info.* = i;
                    return;
                }
                fset(dl, i, zero_r);
            } else {
                fact = fget(d.*, i) / fget(dl.*, i);
                fset(d, i, fget(dl.*, i));
                temp = fget(d.*, i + 1);
                fset(d, i + 1, fget(du.*, i) - fact*temp);
                fset(dl, i, fget(du.*, i + 1));
                fset(du, i + 1, -fact*fget(dl.*, i));
                fset(du, i, temp);
                temp = fget(b.*, i, 1);
                fset(b, i, 1, fget(b.*, i + 1, 1));
                fset(b, i + 1, 1, temp - fact*fget(b.*, i + 1, 1));
            }
        }

        if n > 1 {
            i := n - 1;
            if abs(fget(d.*, i)) >= abs(fget(dl.*, i)) {
                if fget(d.*, i) != zero_r {
                    fact = fget(dl.*, i) / fget(d.*, i);
                    fset(d, i + 1, fget(d.*, i + 1) - fact*fget(du.*, i));
                    fset(b, i + 1, 1, fget(b.*, i + 1, 1) - fact*fget(b.*, i, 1));
                } else {
                    info.* = i;
                    return;
                }
            } else {
                fact = fget(d.*, i) / fget(dl.*, i);
                fset(d, i, fget(dl.*, i));
                temp = fget(d.*, i + 1);
                fset(d, i + 1, fget(du.*, i) - fact*temp);
                fset(du, i, temp);
                temp = fget(b.*, i, 1);
                fset(b, i, 1, fget(b.*, i + 1, 1));
                fset(b, i + 1, 1, temp - fact*fget(b.*, i + 1, 1));
            }
        }

        if fget(d.*, n) == zero_r {
            info.* = n;
            return;
        }
    } else {
        for i:1..n - 2 {
            if abs(fget(d.*, i)) >= abs(fget(dl.*, i)) {
                if fget(d.*, i) != zero_r {
                    fact = fget(dl.*, i) / fget(d.*, i);
                    fset(d, i + 1, fget(d.*, i + 1) - fact*fget(du.*, i));
                    for j:1..nrhs {
                        fset(b, i + 1, j, fget(b.*, i + 1, j) - fact*fget(b.*, i, j));
                    }
                } else {
                    info.* = i;
                    return;
                }
                fset(dl, i, zero_r);
            } else {
                fact = fget(d.*, i) / fget(dl.*, i);
                fset(d, i, fget(dl.*, i));
                temp = fget(d.*, i + 1);
                fset(d, i + 1, fget(du.*, i) - fact*temp);
                fset(dl, i, fget(du.*, i + 1));
                fset(du, i + 1, -fact*fget(dl.*, i));
                fset(du, i, temp);
                for j:1..nrhs {
                    temp = fget(b.*, i, j);
                    fset(b, i, j, fget(b.*, i + 1, j));
                    fset(b, i + 1, j, temp - fact*fget(b.*, i + 1, j));
                }
            }
        }

        if n > 1 {
            i := n - 1;
            if abs(fget(d.*, i)) >= abs(fget(dl.*, i)) {
                if fget(d.*, i) != zero_r {
                    fact = fget(dl.*, i) / fget(d.*, i);
                    fset(d, i + 1, fget(d.*, i + 1) - fact*fget(du.*, i));
                    for j:1..nrhs {
                        fset(b, i + 1, j, fget(b.*, i + 1, j) - fact*fget(b.*, i, j));
                    }
                } else {
                    info.* = i;
                    return;
                }
            } else {
                fact = fget(d.*, i) / fget(dl.*, i);
                fset(d, i, fget(dl.*, i));
                temp = fget(d.*, i + 1);
                fset(d, i + 1, fget(du.*, i) - fact*temp);
                fset(du, i, temp);
                for j:1..nrhs {
                    temp = fget(b.*, i, j);
                    fset(b, i, j, fget(b.*, i + 1, j));
                    fset(b, i + 1, j, temp - fact*fget(b.*, i + 1, j));
                }
            }
        }

        if fget(d.*, n) == zero_r {
            info.* = n;
            return;
        }
    }

    for j:1..nrhs {
        fset(b, n, j, fget(b.*, n, j) / fget(d.*, n));
        if n > 1 {
            fset(b, n - 1, j, (fget(b.*, n - 1, j) - fget(du.*, n - 1)*fget(b.*, n, j)) / fget(d.*, n - 1));
        }
        for <k:n - 2..1 {
            fset(b, k, j, (fget(b.*, k, j) - fget(du.*, k)*fget(b.*, k + 1, j) - fget(dl.*, k)*fget(b.*, k + 2, j)) / fget(d.*, k));
        }
    }
}

_gtsv_complex :: (
    $$n: int,
    $$nrhs: int,
    dl: *$DL/ArrayType1D,
    d: *$D/ArrayType1D,
    du: *$DU/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           B.data_type == MT;
} {
    zero_t :: #run cpx(zero(RT), zero(RT));
    routine_name :: #run if MT == C32 { "CGTSV " } else { "ZGTSV " };

    mult, temp : MT = ---;

    info.* = 0;
    if n < 0 {
        info.* = -1;
    } else if nrhs < 0 {
        info.* = -2;
    } else if ldb < max(1, n) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    for k:1..n - 1 {
        if fget(dl.*, k) == zero_t {
            if fget(d.*, k) == zero_t {
                info.* = k;
                return;
            }
        } else if _gtsv_abs1(fget(d.*, k), RT, MT) >= _gtsv_abs1(fget(dl.*, k), RT, MT) {
            mult = fget(dl.*, k) / fget(d.*, k);
            fset(d, k + 1, fget(d.*, k + 1) - mult*fget(du.*, k));
            for j:1..nrhs {
                fset(b, k + 1, j, fget(b.*, k + 1, j) - mult*fget(b.*, k, j));
            }
            if k < n - 1 {
                fset(dl, k, zero_t);
            }
        } else {
            mult = fget(d.*, k) / fget(dl.*, k);
            fset(d, k, fget(dl.*, k));
            temp = fget(d.*, k + 1);
            fset(d, k + 1, fget(du.*, k) - mult*temp);
            if k < n - 1 {
                fset(dl, k, fget(du.*, k + 1));
                fset(du, k + 1, -mult*fget(dl.*, k));
            }
            fset(du, k, temp);
            for j:1..nrhs {
                temp = fget(b.*, k, j);
                fset(b, k, j, fget(b.*, k + 1, j));
                fset(b, k + 1, j, temp - mult*fget(b.*, k + 1, j));
            }
        }
    }

    if fget(d.*, n) == zero_t {
        info.* = n;
        return;
    }

    for j:1..nrhs {
        fset(b, n, j, fget(b.*, n, j) / fget(d.*, n));
        if n > 1 {
            fset(b, n - 1, j, (fget(b.*, n - 1, j) - fget(du.*, n - 1)*fget(b.*, n, j)) / fget(d.*, n - 1));
        }
        for <k:n - 2..1 {
            fset(b, k, j, (fget(b.*, k, j) - fget(du.*, k)*fget(b.*, k + 1, j) - fget(dl.*, k)*fget(b.*, k + 2, j)) / fget(d.*, k));
        }
    }
}
