// ? LAPACK *gbtrf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbtrf/cgbtrf/dgbtrf/zgbtrf
// *GBTRF computes an LU factorization of a general band matrix
// using partial pivoting with row interchanges.

sgbtrf :: #bake_arguments _gbtrf(MT=float32);
dgbtrf :: #bake_arguments _gbtrf(MT=float64);
cgbtrf :: #bake_arguments _gbtrf(MT=C32);
zgbtrf :: #bake_arguments _gbtrf(MT=C64);

#scope_file
_gbtrf_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gbtrf_matrix_segment :: inline (a: $A/ArrayType2D, $$ld: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
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

_gbtrf_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$ld: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int,
    $$lda_view: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda_view + m;
    return _gbtrf_raw_slice(a.data + (col - 1)*ld + (row - 1), count);
}

_gbtrf_vector_storage :: inline (
    a: $A/ArrayType2D,
    $$ld: int,
    $$row: int,
    $$col: int,
    $$n: int,
    $$inc: int
) -> [] A.data_type {
    count := ifx n <= 0 then 0 else 1 + (n - 1)*abs(inc);
    return _gbtrf_raw_slice(a.data + (col - 1)*ld + (row - 1), count);
}

_gbtrf_ipiv_slice :: inline (ipiv: $IPIV/ArrayType1D, $$offset: int, $$count: int) -> [] IPIV.data_type {
    view: [] IPIV.data_type = ipiv;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gbtrf_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $MT: Type
) -> int #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return isamax(n, x, 1) + 1;
    } else if MT == float64 {
        return idamax(n, x, 1) + 1;
    } else if MT == C32 {
        return icamax(n, x, 1) + 1;
    } else {
        return izamax(n, x, 1) + 1;
    }
}

_gbtrf_swap :: (
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

_gbtrf_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
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
        scopy(n, x, incx, y, incy);
    } else if MT == float64 {
        dcopy(n, x, incx, y, incy);
    } else if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_gbtrf_scal :: (
    $$n: int,
    alpha: $MT,
    x: *$X/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, 1);
    } else if MT == float64 {
        dscal(n, alpha, x, 1);
    } else if MT == C32 {
        cscal(n, alpha, x, 1);
    } else {
        zscal(n, alpha, x, 1);
    }
}

_gbtrf_laswp :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$k1: int,
    $$k2: int,
    ipiv: $IPIV/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    #if MT == float32 {
        slaswp(n, a, lda, k1, k2, ipiv, incx);
    } else if MT == float64 {
        dlaswp(n, a, lda, k1, k2, ipiv, incx);
    } else if MT == C32 {
        claswp(n, a, lda, k1, k2, ipiv, incx);
    } else {
        zlaswp(n, a, lda, k1, k2, ipiv, incx);
    }
}

_gbtrf_trsm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == float64 {
        dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_gbtrf_gemm :: (
    $$transa: u8,
    $$transb: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == float32 {
        sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == float64 {
        dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_gbtrf_gbtf2 :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT &&
           IPIV.data_type == int;
} {
    #if MT == float32 {
        sgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
    } else if MT == float64 {
        dgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
    } else if MT == C32 {
        cgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
    } else {
        zgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
    }
}

_gbtrf_rank1_update :: (
    $$km: int,
    $$jj: int,
    $$jm: int,
    ab: *$AB/ArrayType2D,
    $$kv: int,
    alpha: $MT
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT;
} {
    zero_t :: #run zero(MT);

    if km <= 0 || jm <= jj then return;

    temp : MT = ---;
    offset, row, col : int = ---;

    for col:jj+1..jm {
        offset = col - (jj + 1);
        temp = alpha * fget(ab.*, kv - offset, col);
        if temp != zero_t {
            for i:1..km {
                row = kv + i - offset;
                fset(ab, row, col, fget(ab.*, row, col) + fget(ab.*, kv + 1 + i, jj) * temp);
            }
        }
    }
}

_gbtrf :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT &&
           IPIV.data_type == int;
} {
    one_t  :: #run one(MT);
    zero_t :: #run zero(MT);
    minus_one_t :: #run -one(MT);
    nbmax :: 64;
    ldwork :: nbmax + 1;
    routine_name :: #run if MT == float32 {
        "SGBTRF"
    } else if MT == float64 {
        "DGBTRF"
    } else if MT == C32 {
        "CGBTRF"
    } else {
        "ZGBTRF"
    };

    work13_flat: [ldwork*nbmax] MT;
    work31_flat: [ldwork*nbmax] MT;
    work13 := array2d_view(work13_flat, ldwork, row_major=false);
    work31 := array2d_view(work31_flat, ldwork, row_major=false);

    i, i2, i3, ii, ip, j2, j3, jb, jj, jm, jp, ju, k2, km, kv, nb, nw : int = ---;
    j : int = ---;
    minmn :: min(m, n);
    temp : MT = ---;

    kv = ku + kl;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if kl < 0 {
        info.* = -3;
    } else if ku < 0 {
        info.* = -4;
    } else if ldab < kl + kv + 1 {
        info.* = -6;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m == 0 || n == 0 then return;

    nb = ilaenv(1, routine_name, " "s, m, n, kl, ku);
    nb = min(nb, nbmax);

    if nb <= 1 || nb > kl {
        _gbtrf_gbtf2(m, n, kl, ku, ab, ldab, ipiv, info, MT);
        return;
    }

    for j:1..nb {
        for i:1..j-1 {
            fset(work13, i, j, zero_t);
        }
    }

    for j:1..nb {
        for i:j+1..nb {
            fset(work31, i, j, zero_t);
        }
    }

    for j:ku+2..min(kv, n) {
        for i:kv-j+2..kl {
            fset(ab, i, j, zero_t);
        }
    }

    ju = 1;
    j = 1;
    while j <= minmn {
        jb = min(nb, minmn - j + 1);
        i2 = min(kl - jb, m - j - jb + 1);
        i3 = min(jb, m - j - kl + 1);

        for jj:j..j+jb-1 {
            if jj + kv <= n {
                for i:1..kl {
                    fset(ab, i, jj + kv, zero_t);
                }
            }

            km = min(kl, m - jj);
            column_view := _gbtrf_matrix_segment(ab.*, ldab, kv + 1, jj, km + 1);
            jp = _gbtrf_iamax_1based(km + 1, column_view, MT);
            fset(ipiv, jj, jp + jj - j);

            if fget(ab.*, kv + jp, jj) != zero_t {
                ju = max(ju, min(jj + ku + jp - 1, n));
                if jp != 1 {
                    if jp + jj - 1 < j + kl {
                        x_view := _gbtrf_vector_storage(ab.*, ldab, kv + 1 + jj - j, j, jb, ldab - 1);
                        y_view := _gbtrf_vector_storage(ab.*, ldab, kv + jp + jj - j, j, jb, ldab - 1);
                        _gbtrf_swap(jb, *x_view, ldab - 1, *y_view, ldab - 1, MT);
                    } else {
                        if jj > j {
                            x_view := _gbtrf_vector_storage(ab.*, ldab, kv + 1 + jj - j, j, jj - j, ldab - 1);
                            y_view := _gbtrf_vector_storage(work31, ldwork, jp + jj - j - kl, 1, jj - j, ldwork);
                            _gbtrf_swap(jj - j, *x_view, ldab - 1, *y_view, ldwork, MT);
                        }
                        if j + jb - jj > 0 {
                            x_view := _gbtrf_vector_storage(ab.*, ldab, kv + 1, jj, j + jb - jj, ldab - 1);
                            y_view := _gbtrf_vector_storage(ab.*, ldab, kv + jp, jj, j + jb - jj, ldab - 1);
                            _gbtrf_swap(j + jb - jj, *x_view, ldab - 1, *y_view, ldab - 1, MT);
                        }
                    }
                }

                multipliers_view := _gbtrf_matrix_segment(ab.*, ldab, kv + 2, jj, km);
                _gbtrf_scal(km, one_t / fget(ab.*, kv + 1, jj), *multipliers_view, MT);

                jm = min(ju, j + jb - 1);
                if jm > jj {
                    _gbtrf_rank1_update(km, jj, jm, ab, kv, minus_one_t);
                }
            } else if info.* == 0 {
                info.* = jj;
            }

            nw = min(jj - j + 1, i3);
            if nw > 0 {
                src_view := _gbtrf_matrix_segment(ab.*, ldab, kv + kl + 1 - jj + j, jj, nw);
                dst_view := _gbtrf_matrix_segment(work31, ldwork, 1, jj - j + 1, nw);
                _gbtrf_copy(nw, src_view, 1, *dst_view, 1, MT);
            }
        }

        if j + jb <= n {
            j2 = min(ju - j + 1, kv) - jb;
            j3 = max(0, ju - j - kv + 1);

            if j2 > 0 {
                ipiv_view := _gbtrf_ipiv_slice(ipiv.*, j, jb);
                a12_raw := _gbtrf_matrix_storage(ab.*, ldab, kv + 1 - jb, j + jb, jb, j2, ldab - 1);
                a12_view := array2d_view(a12_raw, ldab - 1, row_major=false);
                _gbtrf_laswp(j2, *a12_view, ldab - 1, 1, jb, ipiv_view, 1, MT);
            }

            for i:j..j+jb-1 {
                fset(ipiv, i, fget(ipiv.*, i) + j - 1);
            }

            k2 = j - 1 + jb + j2;
            for i:1..j3 {
                jj = k2 + i;
                for ii:j+i-1..j+jb-1 {
                    ip = fget(ipiv.*, ii);
                    if ip != ii {
                        temp = fget(ab.*, kv + 1 + ii - jj, jj);
                        fset(ab, kv + 1 + ii - jj, jj, fget(ab.*, kv + 1 + ip - jj, jj));
                        fset(ab, kv + 1 + ip - jj, jj, temp);
                    }
                }
            }

            if j2 > 0 {
                l11_raw := _gbtrf_matrix_storage(ab.*, ldab, kv + 1, j, jb, jb, ldab - 1);
                l11_view := array2d_view(l11_raw, ldab - 1, row_major=false);
                a12_raw := _gbtrf_matrix_storage(ab.*, ldab, kv + 1 - jb, j + jb, jb, j2, ldab - 1);
                a12_view := array2d_view(a12_raw, ldab - 1, row_major=false);
                _gbtrf_trsm(#char "L", #char "L", #char "N", #char "U", jb, j2, one_t, l11_view, ldab - 1, *a12_view, ldab - 1, MT);

                if i2 > 0 {
                    a21_raw := _gbtrf_matrix_storage(ab.*, ldab, kv + 1 + jb, j, i2, jb, ldab - 1);
                    a21_view := array2d_view(a21_raw, ldab - 1, row_major=false);
                    a22_raw := _gbtrf_matrix_storage(ab.*, ldab, kv + 1, j + jb, i2, j2, ldab - 1);
                    a22_view := array2d_view(a22_raw, ldab - 1, row_major=false);
                    _gbtrf_gemm(#char "N", #char "N", i2, j2, jb, minus_one_t, a21_view, ldab - 1, a12_view, ldab - 1, one_t, *a22_view, ldab - 1, MT);
                }

                if i3 > 0 {
                    a32_raw := _gbtrf_matrix_storage(ab.*, ldab, kv + kl + 1 - jb, j + jb, i3, j2, ldab - 1);
                    a32_view := array2d_view(a32_raw, ldab - 1, row_major=false);
                    work31_block_raw := _gbtrf_matrix_storage(work31, ldwork, 1, 1, i3, jb, ldwork);
                    work31_block_view := array2d_view(work31_block_raw, ldwork, row_major=false);
                    _gbtrf_gemm(#char "N", #char "N", i3, j2, jb, minus_one_t, work31_block_view, ldwork, a12_view, ldab - 1, one_t, *a32_view, ldab - 1, MT);
                }
            }

            if j3 > 0 {
                for jj:1..j3 {
                    for ii:jj..jb {
                        fset(work13, ii, jj, fget(ab.*, ii - jj + 1, jj + j + kv - 1));
                    }
                }

                l11_raw := _gbtrf_matrix_storage(ab.*, ldab, kv + 1, j, jb, jb, ldab - 1);
                l11_view := array2d_view(l11_raw, ldab - 1, row_major=false);
                work13_raw := _gbtrf_matrix_storage(work13, ldwork, 1, 1, jb, j3, ldwork);
                work13_view := array2d_view(work13_raw, ldwork, row_major=false);
                _gbtrf_trsm(#char "L", #char "L", #char "N", #char "U", jb, j3, one_t, l11_view, ldab - 1, *work13_view, ldwork, MT);

                if i2 > 0 {
                    a21_raw := _gbtrf_matrix_storage(ab.*, ldab, kv + 1 + jb, j, i2, jb, ldab - 1);
                    a21_view := array2d_view(a21_raw, ldab - 1, row_major=false);
                    a23_raw := _gbtrf_matrix_storage(ab.*, ldab, 1 + jb, j + kv, i2, j3, ldab - 1);
                    a23_view := array2d_view(a23_raw, ldab - 1, row_major=false);
                    _gbtrf_gemm(#char "N", #char "N", i2, j3, jb, minus_one_t, a21_view, ldab - 1, work13_view, ldwork, one_t, *a23_view, ldab - 1, MT);
                }

                if i3 > 0 {
                    work31_block_raw := _gbtrf_matrix_storage(work31, ldwork, 1, 1, i3, jb, ldwork);
                    work31_block_view := array2d_view(work31_block_raw, ldwork, row_major=false);
                    a33_raw := _gbtrf_matrix_storage(ab.*, ldab, 1 + kl, j + kv, i3, j3, ldab - 1);
                    a33_view := array2d_view(a33_raw, ldab - 1, row_major=false);
                    _gbtrf_gemm(#char "N", #char "N", i3, j3, jb, minus_one_t, work31_block_view, ldwork, work13_view, ldwork, one_t, *a33_view, ldab - 1, MT);
                }

                for jj:1..j3 {
                    for ii:jj..jb {
                        fset(ab, ii - jj + 1, jj + j + kv - 1, fget(work13, ii, jj));
                    }
                }
            }
        } else {
            for i:j..j+jb-1 {
                fset(ipiv, i, fget(ipiv.*, i) + j - 1);
            }
        }

        for <jj:j..j+jb-1 {
            jp = fget(ipiv.*, jj) - jj + 1;
            if jp != 1 {
                if jp + jj - 1 < j + kl {
                    if jj > j {
                        x_view := _gbtrf_vector_storage(ab.*, ldab, kv + 1 + jj - j, j, jj - j, ldab - 1);
                        y_view := _gbtrf_vector_storage(ab.*, ldab, kv + jp + jj - j, j, jj - j, ldab - 1);
                        _gbtrf_swap(jj - j, *x_view, ldab - 1, *y_view, ldab - 1, MT);
                    }
                } else {
                    if jj > j {
                        x_view := _gbtrf_vector_storage(ab.*, ldab, kv + 1 + jj - j, j, jj - j, ldab - 1);
                        y_view := _gbtrf_vector_storage(work31, ldwork, jp + jj - j - kl, 1, jj - j, ldwork);
                        _gbtrf_swap(jj - j, *x_view, ldab - 1, *y_view, ldwork, MT);
                    }
                }
            }

            nw = min(i3, jj - j + 1);
            if nw > 0 {
                src_view := _gbtrf_matrix_segment(work31, ldwork, 1, jj - j + 1, nw);
                dst_view := _gbtrf_matrix_segment(ab.*, ldab, kv + kl + 1 - jj + j, jj, nw);
                _gbtrf_copy(nw, src_view, 1, *dst_view, 1, MT);
            }
        }

        j += nb;
    }
}
