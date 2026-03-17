// ? LAPACK *gebal()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgebal/cgebal/dgebal/zgebal
// *GEBAL balances a general matrix by permutation and/or scaling.

sgebal :: #bake_arguments _gebal(RT=float32, MT=float32);
dgebal :: #bake_arguments _gebal(RT=float64, MT=float64);
cgebal :: #bake_arguments _gebal(RT=float32, MT=C32);
zgebal :: #bake_arguments _gebal(RT=float64, MT=C64);

#scope_file
_gebal_matrix_col_prefix :: inline (
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

_gebal_matrix_row_suffix :: inline (
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

_gebal_swap :: (
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

_gebal_scal :: (
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

_gebal_swap_cols_prefix :: (
    $$len: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$col1: int,
    $$col2: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT;
} {
    if col1 == col2 then return;

    col1_view := _gebal_matrix_col_prefix(a.*, lda, 1, col1, len);
    col2_view := _gebal_matrix_col_prefix(a.*, lda, 1, col2, len);
    _gebal_swap(len, *col1_view, 1, *col2_view, 1, MT);
}

_gebal_swap_rows_suffix :: (
    $$n: int,
    $$k: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row1: int,
    $$row2: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT;
} {
    if row1 == row2 then return;

    row1_view := _gebal_matrix_row_suffix(a.*, lda, row1, k, n-k+1);
    row2_view := _gebal_matrix_row_suffix(a.*, lda, row2, k, n-k+1);
    _gebal_swap(n-k+1, *row1_view, lda, *row2_view, lda, MT);
}

_gebal_scale_row_suffix :: (
    $$n: int,
    $$k: int,
    g: $RT,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    row_view := _gebal_matrix_row_suffix(a.*, lda, row, k, n-k+1);
    _gebal_scal(n-k+1, g, *row_view, lda, RT, MT);
}

_gebal_scale_col_prefix :: (
    $$l: int,
    f: $RT,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$col: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    col_view := _gebal_matrix_col_prefix(a.*, lda, 1, col, l);
    _gebal_scal(l, f, *col_view, 1, RT, MT);
}

_gebal_is_nonzero :: inline (z: $MT, $RT: Type, $MT: Type) -> bool #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return z != zero(MT);
    } else {
        return real(z) != zero(RT) || imag(z) != zero(RT);
    }
}

_gebal_abs_value :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    return abs(z);
}

_gebal_iamax_1based :: (
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

_gebal_nrm2 :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) -> RT #modify {
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

_gebal_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gebal_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_gebal_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_gebal :: (
    $$job: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ilo: *int,
    ihi: *int,
    scale: *$SCALE/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           SCALE.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    sclfac :: #run ncast(RT, 2.0);
    factor :: #run ncast(RT, 0.95);
    routine_name :: #run if MT == float32 {
        "SGEBAL"
    } else if MT == float64 {
        "DGEBAL"
    } else if MT == C32 {
        "CGEBAL"
    } else {
        "ZGEBAL"
    };

    noconv, canswap : bool = ---;
    i, ica, ira, j, k, l : int = ---;
    c, ca, f, g, r, ra, s, sfmax1, sfmax2, sfmin1, sfmin2 : RT = ---;

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
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        ilo.* = 1;
        ihi.* = 0;
        return;
    }

    if lsame(job, #char "N") {
        for i:1..n {
            fset(scale, i, one_r);
        }
        ilo.* = 1;
        ihi.* = n;
        return;
    }

    k = 1;
    l = n;

    if !lsame(job, #char "S") {
        noconv = true;
        while noconv {
            noconv = false;
            for <i:1..l {
                canswap = true;
                for j:1..l {
                    if i != j && _gebal_is_nonzero(fget(a.*, i, j), RT, MT) {
                        canswap = false;
                        break;
                    }
                }

                if canswap {
                    fset(scale, l, ncast(RT, i));
                    if i != l {
                        _gebal_swap_cols_prefix(l, a, lda, i, l, MT);
                        _gebal_swap_rows_suffix(n, k, a, lda, i, l, MT);
                    }
                    noconv = true;

                    if l == 1 {
                        ilo.* = 1;
                        ihi.* = 1;
                        return;
                    }

                    l -= 1;
                }
            }
        }

        noconv = true;
        while noconv {
            noconv = false;
            for j:k..l {
                canswap = true;
                for i:k..l {
                    if i != j && _gebal_is_nonzero(fget(a.*, i, j), RT, MT) {
                        canswap = false;
                        break;
                    }
                }

                if canswap {
                    fset(scale, k, ncast(RT, j));
                    if j != k {
                        _gebal_swap_cols_prefix(l, a, lda, j, k, MT);
                        _gebal_swap_rows_suffix(n, k, a, lda, j, k, MT);
                    }
                    noconv = true;
                    k += 1;
                }
            }
        }
    }

    for i:k..l {
        fset(scale, i, one_r);
    }

    if lsame(job, #char "P") {
        ilo.* = k;
        ihi.* = l;
        return;
    }

    sfmin1 = _gebal_lamch_safe_minimum(RT) / _gebal_lamch_precision(RT);
    sfmax1 = one_r / sfmin1;
    sfmin2 = sfmin1 * sclfac;
    sfmax2 = one_r / sfmin2;

    noconv = true;
    while noconv {
        noconv = false;

        for i:k..l {
            col_view := _gebal_matrix_col_prefix(a.*, lda, k, i, l-k+1);
            row_view := _gebal_matrix_row_suffix(a.*, lda, i, k, l-k+1);
            c = _gebal_nrm2(l-k+1, col_view, 1, RT, MT);
            r = _gebal_nrm2(l-k+1, row_view, lda, RT, MT);

            ica = _gebal_iamax_1based(l, _gebal_matrix_col_prefix(a.*, lda, 1, i, l), 1, RT, MT);
            ca = _gebal_abs_value(fget(a.*, ica, i), RT, MT);
            ira = _gebal_iamax_1based(n-k+1, _gebal_matrix_row_suffix(a.*, lda, i, k, n-k+1), lda, RT, MT);
            ra = _gebal_abs_value(fget(a.*, i, ira+k-1), RT, MT);

            if c == zero_r || r == zero_r then continue;

            if _gebal_isnan(c + ca + r + ra, RT) {
                info.* = -3;
                xerbla(routine_name, -info.*);
                return;
            }

            g = r / sclfac;
            f = one_r;
            s = c + r;

            while c < g && max(f, c, ca) < sfmax2 && min(r, g, ra) > sfmin2 {
                f *= sclfac;
                c *= sclfac;
                ca *= sclfac;
                r /= sclfac;
                g /= sclfac;
                ra /= sclfac;
            }

            g = c / sclfac;
            while g >= r && max(r, ra) < sfmax2 && min(f, c, g, ca) > sfmin2 {
                f /= sclfac;
                c /= sclfac;
                g /= sclfac;
                ca /= sclfac;
                r *= sclfac;
                ra *= sclfac;
            }

            if (c + r) >= factor * s then continue;
            if f < one_r && fget(scale.*, i) < one_r {
                if f * fget(scale.*, i) <= sfmin1 then continue;
            }
            if f > one_r && fget(scale.*, i) > one_r {
                if fget(scale.*, i) >= sfmax1 / f then continue;
            }

            g = one_r / f;
            fset(scale, i, fget(scale.*, i) * f);
            noconv = true;

            _gebal_scale_row_suffix(n, k, g, a, lda, i, RT, MT);
            _gebal_scale_col_prefix(l, f, a, lda, i, RT, MT);
        }
    }

    ilo.* = k;
    ihi.* = l;
}
