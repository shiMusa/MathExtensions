// ? LAPACK *hpgst()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpgst/zhpgst
// *HPGST reduces a complex Hermitian-definite generalized eigenproblem
// to standard form using packed storage.

chpgst :: #bake_arguments _hpgst(RT=float32, MT=C32);
zhpgst :: #bake_arguments _hpgst(RT=float64, MT=C64);

#scope_file
_hpgst_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hpgst_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hpgst_raw_slice(x.data + offset - 1, count);
}

_hpgst_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    return _hpgst_raw_slice(a.data + (col - 1)*ld, count);
}

_hpgst_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hpgst_upper_index :: inline ($$i: int, $$j: int) -> int {
    return i + (j - 1)*j/2;
}

_hpgst_lower_index :: inline ($$n: int, $$i: int, $$j: int) -> int {
    return i + (j - 1)*(2*n - j)/2;
}

_hpgst_get :: inline ($$uplo: u8, $$n: int, ap: $AP/ArrayType1D, $$i: int, $$j: int) -> AP.data_type {
    if lsame(uplo, #char "U") {
        return fget(ap, _hpgst_upper_index(i, j));
    } else {
        return fget(ap, _hpgst_lower_index(n, i, j));
    }
}

_hpgst_set :: inline ($$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, $$i: int, $$j: int, value: AP.data_type) {
    if lsame(uplo, #char "U") {
        fset(ap, _hpgst_upper_index(i, j), value);
    } else {
        fset(ap, _hpgst_lower_index(n, i, j), value);
    }
}

_hpgst_axpy :: (
    $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { caxpy(n, alpha, x, incx, y, incy); } else { zaxpy(n, alpha, x, incx, y, incy); }
}

_hpgst_hpmv :: (
    $$uplo: u8, $$n: int, alpha: $MT, ap: $AP/ArrayType1D, x: $X/ArrayType1D, $$incx: int,
    beta: MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); } else { zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); }
}

_hpgst_hpr2 :: (
    $$uplo: u8, $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int,
    ap: *$AP/ArrayType1D, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT && AP.data_type == MT;
} {
    #if MT == C32 { chpr2(uplo, n, alpha, x, incx, y, incy, ap); } else { zhpr2(uplo, n, alpha, x, incx, y, incy, ap); }
}

_hpgst_scal_real :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hpgst_tpmv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, ap: $AP/ArrayType1D, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { ctpmv(uplo, trans, diag, n, ap, x, incx); } else { ztpmv(uplo, trans, diag, n, ap, x, incx); }
}

_hpgst_tpsv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, ap: $AP/ArrayType1D, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { ctpsv(uplo, trans, diag, n, ap, x, incx); } else { ztpsv(uplo, trans, diag, n, ap, x, incx); }
}

_hpgst_dotc :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $MT: Type
) -> MT #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { return cdotc(n, x, incx, y, incy); } else { return zdotc(n, x, incx, y, incy); }
}

_hpgst :: (
    $$itype: int, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, bp: $BP/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && BP.data_type == MT;
} {
    one_r  :: #run one(RT);
    half_r :: #run cast(RT) 0.5;
    cone_t :: #run _hpgst_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHPGST" } else { "ZHPGST" };

    upper : bool = ---;
    j, j1, j1j1, jj, k, k1, k1k1, kk : int = ---;
    ajj, akk, bjj, bkk : RT = ---;
    ct : MT = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if itype == 1 {
        if upper {
            jj = 0;
            for j:1..n {
                j1 = jj + 1;
                jj = jj + j;

                _hpgst_set(uplo, n, ap, j, j, _hpgst_complex_from_real(real(fget(ap.*, jj)), RT, MT));
                bjj = real(fget(bp, jj));

                colj := _hpgst_subvector(ap.*, j1, j);
                _hpgst_tpsv(uplo, #char "C", #char "N", j, bp, *colj, 1, MT);

                if j > 1 {
                    lead_ap := _hpgst_subvector(ap.*, 1, (j - 1)*j/2);
                    bp_col := _hpgst_subvector(bp, j1, j - 1);
                    colj_head := _hpgst_subvector(ap.*, j1, j - 1);
                    _hpgst_hpmv(uplo, j - 1, -cone_t, lead_ap, bp_col, 1, cone_t, *colj_head, 1, MT);
                    _hpgst_scal_real(j - 1, one_r / bjj, *colj_head, 1, RT, MT);
                    diag_val := (fget(ap.*, jj) - _hpgst_dotc(j - 1, colj_head, 1, bp_col, 1, MT)) / _hpgst_complex_from_real(bjj, RT, MT);
                    fset(ap, jj, diag_val);
                } else {
                    fset(ap, jj, fget(ap.*, jj) / _hpgst_complex_from_real(bjj, RT, MT));
                }
            }
        } else {
            kk = 1;
            for k:1..n {
                k1k1 = kk + n - k + 1;

                akk = real(fget(ap.*, kk));
                bkk = real(fget(bp, kk));
                akk = akk / (bkk*bkk);
                fset(ap, kk, _hpgst_complex_from_real(akk, RT, MT));

                if k < n {
                    ntrail := n - k;
                    acol := _hpgst_subvector(ap.*, kk + 1, ntrail);
                    bcol := _hpgst_subvector(bp, kk + 1, ntrail);
                    tail_ap := _hpgst_subvector(ap.*, k1k1, ntrail*(ntrail + 1)/2);

                    _hpgst_scal_real(ntrail, one_r / bkk, *acol, 1, RT, MT);
                    ct = _hpgst_complex_from_real(-half_r*akk, RT, MT);
                    _hpgst_axpy(ntrail, ct, bcol, 1, *acol, 1, MT);
                    _hpgst_hpr2(uplo, ntrail, -cone_t, acol, 1, bcol, 1, *tail_ap, MT);
                    _hpgst_axpy(ntrail, ct, bcol, 1, *acol, 1, MT);

                    tail_bp := _hpgst_subvector(bp, k1k1, ntrail*(ntrail + 1)/2);
                    _hpgst_tpsv(uplo, #char "N", #char "N", ntrail, tail_bp, *acol, 1, MT);
                }
                kk = k1k1;
            }
        }
    } else {
        if upper {
            kk = 0;
            for k:1..n {
                k1 = kk + 1;
                kk = kk + k;

                akk = real(fget(ap.*, kk));
                bkk = real(fget(bp, kk));

                if k > 1 {
                    kprev := k - 1;
                    lead_bp := _hpgst_subvector(bp, 1, kprev*k/2);
                    lead_ap := _hpgst_subvector(ap.*, 1, kprev*k/2);
                    acol := _hpgst_subvector(ap.*, k1, kprev);
                    bcol := _hpgst_subvector(bp, k1, kprev);

                    _hpgst_tpmv(uplo, #char "N", #char "N", kprev, lead_bp, *acol, 1, MT);
                    ct = _hpgst_complex_from_real(half_r*akk, RT, MT);
                    _hpgst_axpy(kprev, ct, bcol, 1, *acol, 1, MT);
                    _hpgst_hpr2(uplo, kprev, cone_t, acol, 1, bcol, 1, *lead_ap, MT);
                    _hpgst_axpy(kprev, ct, bcol, 1, *acol, 1, MT);
                    _hpgst_scal_real(kprev, bkk, *acol, 1, RT, MT);
                }

                fset(ap, kk, _hpgst_complex_from_real(akk*bkk*bkk, RT, MT));
            }
        } else {
            jj = 1;
            for j:1..n {
                j1j1 = jj + n - j + 1;

                ajj = real(fget(ap.*, jj));
                bjj = real(fget(bp, jj));
                if j < n {
                    ntrail := n - j;
                    acol := _hpgst_subvector(ap.*, jj + 1, ntrail);
                    bcol := _hpgst_subvector(bp, jj + 1, ntrail);
                    tail_ap := _hpgst_subvector(ap.*, j1j1, ntrail*(ntrail + 1)/2);
                    tail_bp := _hpgst_subvector(bp, jj, (ntrail + 1)*(ntrail + 2)/2);

                    diag_val := _hpgst_complex_from_real(ajj*bjj, RT, MT) + _hpgst_dotc(ntrail, acol, 1, bcol, 1, MT);
                    fset(ap, jj, diag_val);
                    _hpgst_scal_real(ntrail, bjj, *acol, 1, RT, MT);
                    _hpgst_hpmv(uplo, ntrail, cone_t, tail_ap, bcol, 1, cone_t, *acol, 1, MT);

                    colj := _hpgst_subvector(ap.*, jj, ntrail + 1);
                    _hpgst_tpmv(uplo, #char "C", #char "N", ntrail + 1, tail_bp, *colj, 1, MT);
                } else {
                    fset(ap, jj, _hpgst_complex_from_real(ajj*bjj*bjj, RT, MT));
                }

                jj = j1j1;
            }
        }
    }
}
