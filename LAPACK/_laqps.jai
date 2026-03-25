// ? LAPACK *laqps()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqps/dlaqps/claqps/zlaqps
// *LAQPS computes a step of QR factorization with column pivoting using
// BLAS level 3.

slaqps :: #bake_arguments _laqps(RT=float32, MT=float32);
dlaqps :: #bake_arguments _laqps(RT=float64, MT=float64);
claqps :: #bake_arguments _laqps(RT=float32, MT=C32);
zlaqps :: #bake_arguments _laqps(RT=float64, MT=C64);

#scope_file
_laqps_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_laqps_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = max(0, count);
    return view;
}

_laqps_matrix_column :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    return _laqps_raw_slice(a.data + (col - 1)*lda + (row - 1), max(0, count));
}

_laqps_matrix_row :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    storage_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _laqps_raw_slice(a.data + (col - 1)*lda + (row - 1), storage_count);
}

_laqps_matrix_storage :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _laqps_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqps_iamax :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx);
    } else {
        return idamax(n, x, incx);
    }
}

_laqps_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "E");
    } else {
        return dlamch(#char "E");
    }
}

_laqps_nrm2 :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
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

_laqps_swap :: ($$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $MT: Type) #modify {
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

_laqps_larfg :: ($$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *$TAU, $RT: Type, $MT: Type) #modify {
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

_laqps_gemv :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    $$incx: int,
    beta: $MT,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == float64 {
        dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_laqps_gemm :: (
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
    beta: $MT,
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

_laqps_zero_mt :: inline ($RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return zero(MT);
    } else {
        return cpx(zero(RT), zero(RT));
    }
}

_laqps_one_mt :: inline ($RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return one(MT);
    } else {
        return cpx(one(RT), zero(RT));
    }
}

_laqps_conj_if_needed :: inline (z: $MT, $MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == float32 || MT == float64 {
        return z;
    } else {
        return conj(z);
    }
}

_laqps :: (
    $$m: int,
    $$n: int,
    $$offset: int,
    $$nb: int,
    kb: *int,
    a: *$A/ArrayType2D,
    $$lda: int,
    jpvt: *$JPVT/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D,
    vn2: *$VN2/ArrayType1D,
    auxv: *$AUXV/ArrayType1D,
    f: *$F/ArrayType2D,
    $$ldf: int,
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
           AUXV.data_type == MT &&
           F.data_type == MT;
} {
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    zero_mt :: #run _laqps_zero_mt(RT, MT);
    one_mt  :: #run _laqps_one_mt(RT, MT);

    itemp, j, k, lastrk, lsticc, pvt, rk : int = ---;
    temp, temp2, tol3z : RT = ---;
    akk : MT = ---;

    lastrk = min(m, n + offset);
    lsticc = 0;
    k = 0;
    tol3z = sqrt(_laqps_lamch_epsilon(RT));

    while k < nb && lsticc == 0 {
        k = k + 1;
        rk = offset + k;

        pvt = (k - 1) + _laqps_iamax(n - k + 1, _laqps_subvector(vn1.*, k, n - k + 1), 1, RT);
        if pvt != k {
            _laqps_swap(m, _laqps_matrix_column(a.*, lda, 1, pvt, m), 1, _laqps_matrix_column(a.*, lda, 1, k, m), 1, MT);
            _laqps_swap(k - 1, _laqps_matrix_row(f.*, ldf, pvt, 1, k - 1), ldf, _laqps_matrix_row(f.*, ldf, k, 1, k - 1), ldf, MT);
            itemp = fget(jpvt.*, pvt);
            fset(jpvt, pvt, fget(jpvt.*, k));
            fset(jpvt, k, itemp);
            fset(vn1, pvt, fget(vn1.*, k));
            fset(vn2, pvt, fget(vn2.*, k));
        }

        if k > 1 {
            #if MT == C32 || MT == C64 {
                for j:1..k - 1 {
                    fset(f, k, j, conj(fget(f.*, k, j)));
                }
            }
            _laqps_gemv(#char "N", m - rk + 1, k - 1, -one_mt, _laqps_matrix_storage(a.*, lda, rk, 1, m - rk + 1, k - 1), lda,
                        _laqps_matrix_row(f.*, ldf, k, 1, k - 1), ldf, one_mt, *_laqps_matrix_column(a.*, lda, rk, k, m - rk + 1), 1, MT);
            #if MT == C32 || MT == C64 {
                for j:1..k - 1 {
                    fset(f, k, j, conj(fget(f.*, k, j)));
                }
            }
        }

        if rk < m {
            _laqps_larfg(m - rk + 1, a.data + (k - 1)*lda + (rk - 1), *_laqps_matrix_column(a.*, lda, rk + 1, k, m - rk), 1, tau.data + (k - 1), RT, MT);
        } else {
            _laqps_larfg(1, a.data + (k - 1)*lda + (rk - 1), *_laqps_matrix_column(a.*, lda, rk, k, 0), 1, tau.data + (k - 1), RT, MT);
        }

        akk = fget(a.*, rk, k);
        fset(a, rk, k, one_mt);

        if k < n {
            #if MT == float32 || MT == float64 {
                _laqps_gemv(#char "T", m - rk + 1, n - k, fget(tau.*, k), _laqps_matrix_storage(a.*, lda, rk, k + 1, m - rk + 1, n - k), lda,
                            _laqps_matrix_column(a.*, lda, rk, k, m - rk + 1), 1, zero_mt, *_laqps_matrix_column(f.*, ldf, k + 1, k, n - k), 1, MT);
            } else {
                _laqps_gemv(#char "C", m - rk + 1, n - k, fget(tau.*, k), _laqps_matrix_storage(a.*, lda, rk, k + 1, m - rk + 1, n - k), lda,
                            _laqps_matrix_column(a.*, lda, rk, k, m - rk + 1), 1, zero_mt, *_laqps_matrix_column(f.*, ldf, k + 1, k, n - k), 1, MT);
            }
        }

        for j:1..k {
            fset(f, j, k, zero_mt);
        }

        if k > 1 {
            #if MT == float32 || MT == float64 {
                _laqps_gemv(#char "T", m - rk + 1, k - 1, -fget(tau.*, k), _laqps_matrix_storage(a.*, lda, rk, 1, m - rk + 1, k - 1), lda,
                            _laqps_matrix_column(a.*, lda, rk, k, m - rk + 1), 1, zero_mt, *auxv, 1, MT);
            } else {
                _laqps_gemv(#char "C", m - rk + 1, k - 1, -fget(tau.*, k), _laqps_matrix_storage(a.*, lda, rk, 1, m - rk + 1, k - 1), lda,
                            _laqps_matrix_column(a.*, lda, rk, k, m - rk + 1), 1, zero_mt, *auxv, 1, MT);
            }
            _laqps_gemv(#char "N", n, k - 1, one_mt, _laqps_matrix_storage(f.*, ldf, 1, 1, n, k - 1), ldf,
                        auxv.*, 1, one_mt, *_laqps_matrix_column(f.*, ldf, 1, k, n), 1, MT);
        }

        if k < n {
            #if MT == float32 || MT == float64 {
                _laqps_gemv(#char "N", n - k, k, -one_mt, _laqps_matrix_storage(f.*, ldf, k + 1, 1, n - k, k), ldf,
                            _laqps_matrix_row(a.*, lda, rk, 1, k), lda, one_mt, *_laqps_matrix_row(a.*, lda, rk, k + 1, n - k), lda, MT);
            } else {
                _laqps_gemm(#char "N", #char "C", 1, n - k, k, -one_mt, _laqps_matrix_storage(a.*, lda, rk, 1, 1, k), lda,
                            _laqps_matrix_storage(f.*, ldf, k + 1, 1, n - k, k), ldf, one_mt,
                            *_laqps_matrix_storage(a.*, lda, rk, k + 1, 1, n - k), lda, MT);
            }
        }

        if rk < lastrk {
            for j:k + 1..n {
                if fget(vn1.*, j) != zero_r {
                    temp = abs(fget(a.*, rk, j)) / fget(vn1.*, j);
                    temp = max(zero_r, (one_r + temp)*(one_r - temp));
                    temp2 = temp*(fget(vn1.*, j) / fget(vn2.*, j))*(fget(vn1.*, j) / fget(vn2.*, j));
                    if temp2 <= tol3z {
                        fset(vn2, j, ncast(RT, lsticc));
                        lsticc = j;
                    } else {
                        fset(vn1, j, fget(vn1.*, j)*sqrt(temp));
                    }
                }
            }
        }

        fset(a, rk, k, akk);
    }

    kb.* = k;
    rk = offset + kb.*;

    if kb.* < min(n, m - offset) {
        _laqps_gemm(#char "N", #char "C", m - rk, n - kb.*, kb.*, -one_mt,
                    _laqps_matrix_storage(a.*, lda, rk + 1, 1, m - rk, kb.*), lda,
                    _laqps_matrix_storage(f.*, ldf, kb.* + 1, 1, n - kb.*, kb.*), ldf, one_mt,
                    *_laqps_matrix_storage(a.*, lda, rk + 1, kb.* + 1, m - rk, n - kb.*), lda, MT);
    }

    while lsticc > 0 {
        itemp = ncast(int, fget(vn2.*, lsticc));
        fset(vn1, lsticc, _laqps_nrm2(m - rk, _laqps_matrix_column(a.*, lda, rk + 1, lsticc, m - rk), 1, RT, MT));
        fset(vn2, lsticc, fget(vn1.*, lsticc));
        lsticc = itemp;
    }
}
