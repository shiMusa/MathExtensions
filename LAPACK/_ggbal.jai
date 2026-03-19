// ? LAPACK *ggbal()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggbal/cggbal/dggbal/zggbal
// *GGBAL balances a pair of general matrices by permutation and/or scaling.

sggbal :: #bake_arguments _ggbal(RT=float32, MT=float32);
dggbal :: #bake_arguments _ggbal(RT=float64, MT=float64);
cggbal :: #bake_arguments _ggbal(RT=float32, MT=C32);
zggbal :: #bake_arguments _ggbal(RT=float64, MT=C64);

#scope_file
_ggbal_matrix_col_prefix :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda + (row - 1);
    view.count = count;
    return view;
}

_ggbal_matrix_row_suffix :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda + (row - 1);
    view.count = if count <= 0 then 0 else 1 + (count - 1)*lda;
    return view;
}

_ggbal_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_ggbal_swap_matrix_vec :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
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

_ggbal_scal_matrix_vec :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, incx);
    } else if MT == float64 {
        dscal(n, alpha, x, incx);
    } else if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_ggbal_scal_real_vec :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        sscal(n, alpha, x, incx);
    } else {
        dscal(n, alpha, x, incx);
    }
}

_ggbal_axpy_real :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT;
} {
    #if RT == float32 {
        saxpy(n, alpha, x, incx, y, incy);
    } else {
        daxpy(n, alpha, x, incx, y, incy);
    }
}

_ggbal_dot_real :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $RT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT;
} {
    #if RT == float32 {
        return sdot(n, x, incx, y, incy);
    } else {
        return ddot(n, x, incx, y, incy);
    }
}

_ggbal_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return isamax(n, x, incx) + 1;
    } else if MT == float64 {
        return idamax(n, x, incx) + 1;
    } else if MT == C32 {
        return icamax(n, x, incx) + 1;
    } else {
        return izamax(n, x, incx) + 1;
    }
}

_ggbal_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_ggbal_is_nonzero :: inline (z: $MT, $RT: Type, $MT: Type) -> bool #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return z != zero(MT);
    } else {
        return real(z) != zero(RT) || imag(z) != zero(RT);
    }
}

_ggbal_abs1_for_logs :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_ggbal_abs_value :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    return abs(z);
}

_ggbal_pow_sclfac :: inline ($$exponent: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    base :: #run ncast(RT, 10.0);
    result: RT = one(RT);

    if exponent > 0 {
        for _:1..exponent {
            result *= base;
        }
    } else if exponent < 0 {
        for _:1..-exponent {
            result /= base;
        }
    }

    return result;
}

_ggbal_swap_cols_prefix :: (
    $$len: int,
    a: *$A/ArrayType2D,
    b: *$B/ArrayType2D,
    $$lda: int,
    $$ldb: int,
    $$col1: int,
    $$col2: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    if col1 == col2 then return;

    a_col1 := _ggbal_matrix_col_prefix(a.*, lda, 1, col1, len);
    a_col2 := _ggbal_matrix_col_prefix(a.*, lda, 1, col2, len);
    b_col1 := _ggbal_matrix_col_prefix(b.*, ldb, 1, col1, len);
    b_col2 := _ggbal_matrix_col_prefix(b.*, ldb, 1, col2, len);

    _ggbal_swap_matrix_vec(len, *a_col1, 1, *a_col2, 1, MT);
    _ggbal_swap_matrix_vec(len, *b_col1, 1, *b_col2, 1, MT);
}

_ggbal_swap_rows_suffix :: (
    $$n: int,
    $$k: int,
    a: *$A/ArrayType2D,
    b: *$B/ArrayType2D,
    $$lda: int,
    $$ldb: int,
    $$row1: int,
    $$row2: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    if row1 == row2 then return;

    a_row1 := _ggbal_matrix_row_suffix(a.*, lda, row1, k, n - k + 1);
    a_row2 := _ggbal_matrix_row_suffix(a.*, lda, row2, k, n - k + 1);
    b_row1 := _ggbal_matrix_row_suffix(b.*, ldb, row1, k, n - k + 1);
    b_row2 := _ggbal_matrix_row_suffix(b.*, ldb, row2, k, n - k + 1);

    _ggbal_swap_matrix_vec(n - k + 1, *a_row1, lda, *a_row2, lda, MT);
    _ggbal_swap_matrix_vec(n - k + 1, *b_row1, ldb, *b_row2, ldb, MT);
}

_ggbal_scale_row_suffix :: (
    $$n: int,
    $$k: int,
    g: $RT,
    a: *$A/ArrayType2D,
    b: *$B/ArrayType2D,
    $$lda: int,
    $$ldb: int,
    $$row: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    a_row := _ggbal_matrix_row_suffix(a.*, lda, row, k, n - k + 1);
    b_row := _ggbal_matrix_row_suffix(b.*, ldb, row, k, n - k + 1);
    _ggbal_scal_matrix_vec(n - k + 1, g, *a_row, lda, RT, MT);
    _ggbal_scal_matrix_vec(n - k + 1, g, *b_row, ldb, RT, MT);
}

_ggbal_scale_col_prefix :: (
    $$l: int,
    f: $RT,
    a: *$A/ArrayType2D,
    b: *$B/ArrayType2D,
    $$lda: int,
    $$ldb: int,
    $$col: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    a_col := _ggbal_matrix_col_prefix(a.*, lda, 1, col, l);
    b_col := _ggbal_matrix_col_prefix(b.*, ldb, 1, col, l);
    _ggbal_scal_matrix_vec(l, f, *a_col, 1, RT, MT);
    _ggbal_scal_matrix_vec(l, f, *b_col, 1, RT, MT);
}

_ggbal :: (
    $$job: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    ilo: *int,
    ihi: *int,
    lscale: *$LSCALE/ArrayType1D,
    rscale: *$RSCALE/ArrayType1D,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           B.data_type == MT &&
           LSCALE.data_type == RT &&
           RSCALE.data_type == RT &&
           WORK.data_type == RT;
} {
    zero_r  :: #run zero(RT);
    half_r  :: #run ncast(RT, 0.5);
    one_r   :: #run one(RT);
    three_r :: #run ncast(RT, 3.0);
    sclfac  :: #run ncast(RT, 10.0);
    routine_name :: #run if MT == float32 {
        "SGGBAL"
    } else if MT == float64 {
        "DGGBAL"
    } else if MT == C32 {
        "CGGBAL"
    } else {
        "ZGGBAL"
    };

    alpha, basl, beta, cab, cmax, coef, coef2, coef5, cor, ew, ewc, gamma, pgamma, rab, sfmax, sfmin, sum, t, ta, tb, tc : RT = ---;
    found, no_nonzero : bool = ---;
    i, icab, iflow, ip1, ir, irab, it, j, jc, jp1, k, kount, l, lcab, lm1, lrab, lsfmax, lsfmin, m, nr, nrp2, row_pick, col_pick : int = ---;

    w1 := _ggbal_subvector(work.*, 1, n);
    w2 := _ggbal_subvector(work.*, n + 1, n);
    w3 := _ggbal_subvector(work.*, 2*n + 1, n);
    w4 := _ggbal_subvector(work.*, 3*n + 1, n);
    w5 := _ggbal_subvector(work.*, 4*n + 1, n);
    w6 := _ggbal_subvector(work.*, 5*n + 1, n);

    info.* = 0;
    if !lsame(job, #char "N") &&
       !lsame(job, #char "P") &&
       !lsame(job, #char "S") &&
       !lsame(job, #char "B") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if ldb < max(1, n) {
        info.* = -6;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        ilo.* = 1;
        ihi.* = n;
        return;
    }

    if n == 1 {
        ilo.* = 1;
        ihi.* = n;
        fset(lscale, 1, one_r);
        fset(rscale, 1, one_r);
        return;
    }

    if lsame(job, #char "N") {
        ilo.* = 1;
        ihi.* = n;
        for i:1..n {
            fset(lscale, i, one_r);
            fset(rscale, i, one_r);
        }
        return;
    }

    k = 1;
    l = n;

    if !lsame(job, #char "S") {
        while true {
            if l == 1 {
                fset(rscale, 1, one_r);
                fset(lscale, 1, one_r);
                break;
            }

            lm1 = l - 1;
            found = false;
            row_pick = 0;
            col_pick = 0;

            for <i:1..l {
                no_nonzero = true;
                jp1 = 1;
                for j:1..lm1 {
                    jp1 = j + 1;
                    if _ggbal_is_nonzero(fget(a.*, i, j), RT, MT) ||
                       _ggbal_is_nonzero(fget(b.*, i, j), RT, MT) {
                        no_nonzero = false;
                        break;
                    }
                }

                if no_nonzero {
                    row_pick = i;
                    col_pick = l;
                    m = l;
                    iflow = 1;
                    found = true;
                    break;
                }

                no_nonzero = true;
                for j:jp1..l {
                    if _ggbal_is_nonzero(fget(a.*, i, j), RT, MT) ||
                       _ggbal_is_nonzero(fget(b.*, i, j), RT, MT) {
                        no_nonzero = false;
                        break;
                    }
                }

                if no_nonzero {
                    row_pick = i;
                    col_pick = jp1 - 1;
                    m = l;
                    iflow = 1;
                    found = true;
                    break;
                }
            }

            if found {
                fset(lscale, m, ncast(RT, row_pick));
                if row_pick != m {
                    _ggbal_swap_rows_suffix(n, k, a, b, lda, ldb, row_pick, m, MT);
                }

                fset(rscale, m, ncast(RT, col_pick));
                if col_pick != m {
                    _ggbal_swap_cols_prefix(l, a, b, lda, ldb, col_pick, m, MT);
                }

                if iflow == 1 {
                    l = lm1;
                    continue;
                }
            }

            found = false;
            lm1 = l - 1;
            for j:k..l {
                no_nonzero = true;
                ip1 = k;
                for i:k..lm1 {
                    ip1 = i + 1;
                    if _ggbal_is_nonzero(fget(a.*, i, j), RT, MT) ||
                       _ggbal_is_nonzero(fget(b.*, i, j), RT, MT) {
                        no_nonzero = false;
                        break;
                    }
                }

                if no_nonzero {
                    row_pick = l;
                    col_pick = j;
                    m = k;
                    iflow = 2;
                    found = true;
                    break;
                }

                no_nonzero = true;
                for i:ip1..l {
                    if _ggbal_is_nonzero(fget(a.*, i, j), RT, MT) ||
                       _ggbal_is_nonzero(fget(b.*, i, j), RT, MT) {
                        no_nonzero = false;
                        break;
                    }
                }

                if no_nonzero {
                    row_pick = ip1 - 1;
                    col_pick = j;
                    m = k;
                    iflow = 2;
                    found = true;
                    break;
                }
            }

            if !found {
                break;
            }

            fset(lscale, m, ncast(RT, row_pick));
            if row_pick != m {
                _ggbal_swap_rows_suffix(n, k, a, b, lda, ldb, row_pick, m, MT);
            }

            fset(rscale, m, ncast(RT, col_pick));
            if col_pick != m {
                _ggbal_swap_cols_prefix(l, a, b, lda, ldb, col_pick, m, MT);
            }

            k += 1;
        }
    }

    ilo.* = k;
    ihi.* = l;

    if lsame(job, #char "P") {
        for i:ilo.*..ihi.* {
            fset(lscale, i, one_r);
            fset(rscale, i, one_r);
        }
        return;
    }

    if ilo.* == ihi.* then return;

    nr = ihi.* - ilo.* + 1;
    for i:ilo.*..ihi.* {
        fset(rscale, i, zero_r);
        fset(lscale, i, zero_r);
        fset(w1, i, zero_r);
        fset(w2, i, zero_r);
        fset(w3, i, zero_r);
        fset(w4, i, zero_r);
        fset(w5, i, zero_r);
        fset(w6, i, zero_r);
    }

    basl = log10(sclfac);
    for i:ilo.*..ihi.* {
        for j:ilo.*..ihi.* {
            if _ggbal_is_nonzero(fget(a.*, i, j), RT, MT) {
                ta = log10(_ggbal_abs1_for_logs(fget(a.*, i, j), RT, MT)) / basl;
            } else {
                ta = zero_r;
            }

            if _ggbal_is_nonzero(fget(b.*, i, j), RT, MT) {
                tb = log10(_ggbal_abs1_for_logs(fget(b.*, i, j), RT, MT)) / basl;
            } else {
                tb = zero_r;
            }

            fset(w5, i, fget(w5, i) - ta - tb);
            fset(w6, j, fget(w6, j) - ta - tb);
        }
    }

    coef = one_r / ncast(RT, 2*nr);
    coef2 = coef*coef;
    coef5 = half_r*coef2;
    nrp2 = nr + 2;
    beta = zero_r;
    it = 1;

    while true {
        w1_active := _ggbal_subvector(w1, ilo.*, nr);
        w2_active := _ggbal_subvector(w2, ilo.*, nr);
        w3_active := _ggbal_subvector(w3, ilo.*, nr);
        w4_active := _ggbal_subvector(w4, ilo.*, nr);
        w5_active := _ggbal_subvector(w5, ilo.*, nr);
        w6_active := _ggbal_subvector(w6, ilo.*, nr);

        gamma =
            _ggbal_dot_real(nr, *w5_active, 1, *w5_active, 1, RT) +
            _ggbal_dot_real(nr, *w6_active, 1, *w6_active, 1, RT);

        ew = zero_r;
        ewc = zero_r;
        for i:ilo.*..ihi.* {
            ew += fget(w5, i);
            ewc += fget(w6, i);
        }

        gamma = coef*gamma - coef2*(ew*ew + ewc*ewc) - coef5*(ew - ewc)*(ew - ewc);
        if gamma == zero_r then break;

        if it != 1 {
            beta = gamma / pgamma;
        }
        t  = coef5*(ewc - three_r*ew);
        tc = coef5*(ew - three_r*ewc);

        _ggbal_scal_real_vec(nr, beta, *w1_active, 1, RT);
        _ggbal_scal_real_vec(nr, beta, *w2_active, 1, RT);
        _ggbal_axpy_real(nr, coef, *w5_active, 1, *w2_active, 1, RT);
        _ggbal_axpy_real(nr, coef, *w6_active, 1, *w1_active, 1, RT);

        for i:ilo.*..ihi.* {
            fset(w1, i, fget(w1, i) + tc);
            fset(w2, i, fget(w2, i) + t);
        }

        for i:ilo.*..ihi.* {
            kount = 0;
            sum = zero_r;
            for j:ilo.*..ihi.* {
                if _ggbal_is_nonzero(fget(a.*, i, j), RT, MT) {
                    kount += 1;
                    sum += fget(w1, j);
                }
                if _ggbal_is_nonzero(fget(b.*, i, j), RT, MT) {
                    kount += 1;
                    sum += fget(w1, j);
                }
            }
            fset(w3, i, ncast(RT, kount)*fget(w2, i) + sum);
        }

        for j:ilo.*..ihi.* {
            kount = 0;
            sum = zero_r;
            for i:ilo.*..ihi.* {
                if _ggbal_is_nonzero(fget(a.*, i, j), RT, MT) {
                    kount += 1;
                    sum += fget(w2, i);
                }
                if _ggbal_is_nonzero(fget(b.*, i, j), RT, MT) {
                    kount += 1;
                    sum += fget(w2, i);
                }
            }
            fset(w4, j, ncast(RT, kount)*fget(w1, j) + sum);
        }

        sum =
            _ggbal_dot_real(nr, *w2_active, 1, *w3_active, 1, RT) +
            _ggbal_dot_real(nr, *w1_active, 1, *w4_active, 1, RT);
        alpha = gamma / sum;

        cmax = zero_r;
        for i:ilo.*..ihi.* {
            cor = alpha*fget(w2, i);
            if abs(cor) > cmax {
                cmax = abs(cor);
            }
            fset(lscale, i, fget(lscale.*, i) + cor);

            cor = alpha*fget(w1, i);
            if abs(cor) > cmax {
                cmax = abs(cor);
            }
            fset(rscale, i, fget(rscale.*, i) + cor);
        }

        if cmax < half_r then break;

        _ggbal_axpy_real(nr, -alpha, *w3_active, 1, *w5_active, 1, RT);
        _ggbal_axpy_real(nr, -alpha, *w4_active, 1, *w6_active, 1, RT);

        pgamma = gamma;
        it += 1;
        if it > nrp2 then break;
    }

    sfmin = _ggbal_lamch_safe_minimum(RT);
    sfmax = one_r / sfmin;
    lsfmin = cast(int) (log10(sfmin)/basl + one_r);
    lsfmax = cast(int) (log10(sfmax)/basl);

    for i:ilo.*..ihi.* {
        irab = _ggbal_iamax_1based(n - ilo.* + 1, _ggbal_matrix_row_suffix(a.*, lda, i, ilo.*, n - ilo.* + 1), lda, RT, MT);
        rab = _ggbal_abs_value(fget(a, i, irab + ilo.* - 1), RT, MT);
        irab = _ggbal_iamax_1based(n - ilo.* + 1, _ggbal_matrix_row_suffix(b.*, ldb, i, ilo.*, n - ilo.* + 1), ldb, RT, MT);
        rab = max(rab, _ggbal_abs_value(fget(b, i, irab + ilo.* - 1), RT, MT));
        lrab = cast(int) (log10(rab + sfmin)/basl + one_r);

        ir = cast(int) (fget(lscale.*, i) + sign(half_r, fget(lscale.*, i)));
        ir = max(ir, lsfmin);
        ir = min(ir, lsfmax);
        ir = min(ir, lsfmax - lrab);
        fset(lscale, i, _ggbal_pow_sclfac(ir, RT));

        icab = _ggbal_iamax_1based(ihi.*, _ggbal_matrix_col_prefix(a.*, lda, 1, i, ihi.*), 1, RT, MT);
        cab = _ggbal_abs_value(fget(a, icab, i), RT, MT);
        icab = _ggbal_iamax_1based(ihi.*, _ggbal_matrix_col_prefix(b.*, ldb, 1, i, ihi.*), 1, RT, MT);
        cab = max(cab, _ggbal_abs_value(fget(b, icab, i), RT, MT));
        lcab = cast(int) (log10(cab + sfmin)/basl + one_r);

        jc = cast(int) (fget(rscale.*, i) + sign(half_r, fget(rscale.*, i)));
        jc = max(jc, lsfmin);
        jc = min(jc, lsfmax);
        jc = min(jc, lsfmax - lcab);
        fset(rscale, i, _ggbal_pow_sclfac(jc, RT));
    }

    for i:ilo.*..ihi.* {
        _ggbal_scale_row_suffix(n, ilo.*, fget(lscale.*, i), a, b, lda, ldb, i, RT, MT);
    }

    for j:ilo.*..ihi.* {
        _ggbal_scale_col_prefix(ihi.*, fget(rscale.*, j), a, b, lda, ldb, j, RT, MT);
    }
}
