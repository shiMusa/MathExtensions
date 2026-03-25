// ? LAPACK *laqp2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqp2/dlaqp2/claqp2/zlaqp2
// *LAQP2 computes a QR factorization with column pivoting of the matrix
// block.

slaqp2 :: #bake_arguments _laqp2(RT=float32, MT=float32);
dlaqp2 :: #bake_arguments _laqp2(RT=float64, MT=float64);
claqp2 :: #bake_arguments _laqp2(RT=float32, MT=C32);
zlaqp2 :: #bake_arguments _laqp2(RT=float64, MT=C64);

#scope_file
_laqp2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_laqp2_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = max(0, count);
    return view;
}

_laqp2_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _laqp2_raw_slice(a.data + (col - 1)*lda + (row - 1), max(0, count));
}

_laqp2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _laqp2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqp2_iamax :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx);
    } else {
        return idamax(n, x, incx);
    }
}

_laqp2_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "E");
    } else {
        return dlamch(#char "E");
    }
}

_laqp2_nrm2 :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return snrm2(n, x, incx);
    } else if MT == float64 {
        return dnrm2(n, x, incx);
    } else if MT == C32 {
        return scnrm2(n, x, incx);
    } else {
        return dznrm2(n, x, incx);
    }
}

_laqp2_swap :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sswap(n, x, incx, y, incy);
    } else if MT == float64 {
        dswap(n, x, incx, y, incy);
    } else if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_laqp2_larfg :: (
    $$n: int,
    alpha: *$MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    tau: *$TAU,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           TAU == MT;
} {
    #if MT == float32 {
        slarfg(n, alpha, x, incx, tau);
    } else if MT == float64 {
        dlarfg(n, alpha, x, incx, tau);
    } else if MT == C32 {
        clarfg(n, alpha, x, incx, tau);
    } else {
        zlarfg(n, alpha, x, incx, tau);
    }
}

_laqp2_larf1f :: (
    $$side: u8,
    $$m: int,
    $$n: int,
    v: $V/ArrayType1D,
    $$incv: int,
    $$tau: $MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == float64 {
        dlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == C32 {
        clarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else {
        zlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    }
}

_laqp2_apply_tau :: inline (tau: $MT, $MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == float32 || MT == float64 {
        return tau;
    } else {
        return conj(tau);
    }
}

_laqp2 :: (
    $$m: int,
    $$n: int,
    $$offset: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    jpvt: *$JPVT/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D,
    vn2: *$VN2/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           JPVT.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == MT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);

    i, itemp, j, mn, offpi, pvt : int = ---;
    temp, temp2, tol3z : RT = ---;
    tau_i : MT = ---;

    mn = min(m - offset, n);
    tol3z = sqrt(_laqp2_lamch_epsilon(RT));

    for i:1..mn {
        offpi = offset + i;

        pvt = (i - 1) + _laqp2_iamax(n - i + 1, _laqp2_subvector(vn1.*, i, n - i + 1), 1, RT);
        if pvt != i {
            _laqp2_swap(m, _laqp2_matrix_column(a.*, lda, 1, pvt, m), 1, _laqp2_matrix_column(a.*, lda, 1, i, m), 1, MT);
            itemp = fget(jpvt.*, pvt);
            fset(jpvt, pvt, fget(jpvt.*, i));
            fset(jpvt, i, itemp);
            fset(vn1, pvt, fget(vn1.*, i));
            fset(vn2, pvt, fget(vn2.*, i));
        }

        if offpi < m {
            _laqp2_larfg(m - offpi + 1, a.data + (i - 1)*lda + (offpi - 1), *_laqp2_matrix_column(a.*, lda, offpi + 1, i, m - offpi), 1, tau.data + (i - 1), RT, MT);
        } else {
            _laqp2_larfg(1, a.data + (i - 1)*lda + (m - 1), *_laqp2_matrix_column(a.*, lda, m, i, 0), 1, tau.data + (i - 1), RT, MT);
        }

        if i < n {
            tau_i = _laqp2_apply_tau(fget(tau.*, i), MT);
            _laqp2_larf1f(#char "L", m - offpi + 1, n - i, _laqp2_matrix_column(a.*, lda, offpi, i, m - offpi + 1), 1,
                          tau_i, *_laqp2_matrix_storage(a.*, lda, offpi, i + 1, m - offpi + 1, n - i), lda, work, RT, MT);
        }

        for j:i + 1..n {
            if fget(vn1.*, j) != zero_r {
                temp = one_r - (abs(fget(a.*, offpi, j)) / fget(vn1.*, j))*(abs(fget(a.*, offpi, j)) / fget(vn1.*, j));
                temp = max(temp, zero_r);
                temp2 = temp*(fget(vn1.*, j) / fget(vn2.*, j))*(fget(vn1.*, j) / fget(vn2.*, j));
                if temp2 <= tol3z {
                    if offpi < m {
                        fset(vn1, j, _laqp2_nrm2(m - offpi, _laqp2_matrix_column(a.*, lda, offpi + 1, j, m - offpi), 1, RT, MT));
                        fset(vn2, j, fget(vn1.*, j));
                    } else {
                        fset(vn1, j, zero_r);
                        fset(vn2, j, zero_r);
                    }
                } else {
                    fset(vn1, j, fget(vn1.*, j)*sqrt(temp));
                }
            }
        }
    }
}
