// ? LAPACK *hetrd_he2hb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrd_he2hb/zhetrd_he2hb
// *HETRD_HE2HB reduces a Hermitian matrix to Hermitian band form.

chetrd_he2hb :: #bake_arguments _hetrd_he2hb(RT=float32, MT=C32);
zhetrd_he2hb :: #bake_arguments _hetrd_he2hb(RT=float64, MT=C64);

#scope_file
_hetrd_he2hb_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrd_he2hb_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetrd_he2hb_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrd_he2hb_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetrd_he2hb_raw_slice(x.data + (offset - 1), count);
}

_hetrd_he2hb_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrd_he2hb_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrd_he2hb_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int,
    a: $A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 { clacpy(uplo, m, n, a, lda, b, ldb); } else { zlacpy(uplo, m, n, a, lda, b, ldb); }
}

_hetrd_he2hb_laset :: (
    $$uplo: u8, $$m: int, $$n: int,
    alpha: $MT, beta: $MT,
    a: *$A/ArrayType2D, $$lda: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 { claset(uplo, m, n, alpha, beta, a, lda); } else { zlaset(uplo, m, n, alpha, beta, a, lda); }
}

_hetrd_he2hb_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { ccopy(n, x, incx, y, incy); } else { zcopy(n, x, incx, y, incy); }
}

_hetrd_he2hb_gelqf :: (
    $$m: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 { cgelqf(m, n, a, lda, tau, work, lwork, info); } else { zgelqf(m, n, a, lda, tau, work, lwork, info); }
}

_hetrd_he2hb_geqrf :: (
    $$m: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 { cgeqrf(m, n, a, lda, tau, work, lwork, info); } else { zgeqrf(m, n, a, lda, tau, work, lwork, info); }
}

_hetrd_he2hb_larft :: (
    $$direct: u8, $$storev: u8, $$n: int, $$k: int,
    v: $V/ArrayType2D, $$ldv: int,
    tau: $TAU/ArrayType1D,
    t: *$T/ArrayType2D, $$ldt: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && V.data_type == MT && TAU.data_type == MT && T.data_type == MT;
} {
    #if MT == C32 { clarft(direct, storev, n, k, v, ldv, tau, t, ldt); } else { zlarft(direct, storev, n, k, v, ldv, tau, t, ldt); }
}

_hetrd_he2hb_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    b: $B/ArrayType2D, $$ldb: int,
    beta: $MT, c: *$C/ArrayType2D, $$ldc: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 { cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); } else { zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
}

_hetrd_he2hb_hemm :: (
    $$side: u8, $$uplo: u8, $$m: int, $$n: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    b: $B/ArrayType2D, $$ldb: int,
    beta: $MT, c: *$C/ArrayType2D, $$ldc: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 { chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); } else { zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); }
}

_hetrd_he2hb_her2k :: (
    $$uplo: u8, $$trans: u8, $$n: int, $$k: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    b: $B/ArrayType2D, $$ldb: int,
    beta: $RT, c: *$C/ArrayType2D, $$ldc: int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 { cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); } else { zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
}

_hetrd_he2hb_copy_upper_band :: (
    $$j: int, $$lk: int,
    a: $A/ArrayType2D, $$lda: int,
    ab: *$AB/ArrayType2D, $$ldab: int
) {
    for t:0..lk - 1 {
        fset(ab, ldab - t, j, fget(a, j, j + t));
    }
}

_hetrd_he2hb_copy_lower_band :: (
    $$j: int, $$lk: int,
    a: $A/ArrayType2D,
    ab: *$AB/ArrayType2D
) {
    for t:0..lk - 1 {
        fset(ab, 1 + t, j, fget(a, j + t, j));
    }
}

_hetrd_he2hb :: (
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AB.data_type == MT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    rone :: #run one(RT);
    zero_t :: #run zero(MT);
    one_t  :: #run _hetrd_he2hb_complex_from_real(rone, RT, MT);
    half_t :: #run _hetrd_he2hb_complex_from_real(ncast(RT) 0.5, RT, MT);
    routine_name :: #run if MT == C32 { "CHETRD_HE2HB" } else { "ZHETRD_HE2HB" };

    i, j, iinfo, lwmin, pn, pk, lk, ldt, ldw, lds2, lds1, ls2, ls1, lw, lt, tpos, wpos, s2pos, s1pos : int = ---;
    lquery, upper : bool = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;
    if n <= kd + 1 {
        lwmin = 1;
    } else {
        lwmin = ilaenv2stage(4, ifx MT == C32 then "CHETRD_HE2HB" else "ZHETRD_HE2HB", "", n, kd, -1, -1);
    }

    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if kd < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldab < max(1, kd + 1) {
        info.* = -7;
    } else if lwork < lwmin && !lquery {
        info.* = -10;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        fset(work, 1, _hetrd_he2hb_complex_from_real(_hetrd_he2hb_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
        return;
    }

    if n <= kd + 1 {
        if upper {
            for i:1..n {
                lk = min(kd + 1, i);
                for j:0..lk - 1 {
                    fset(ab, kd + 1 - (lk - 1) + j, i, fget(a.*, i - lk + 1 + j, i));
                }
            }
        } else {
            for i:1..n {
                lk = min(kd + 1, n - i + 1);
                for j:0..lk - 1 {
                    fset(ab, 1 + j, i, fget(a.*, i + j, i));
                }
            }
        }
        fset(work, 1, one_t);
        return;
    }

    ldt = kd;
    lds1 = kd;
    lt = ldt*kd;
    lw = n*kd;
    ls1 = lds1*kd;
    ls2 = lwmin - lt - lw - ls1;
    tpos = 1;
    wpos = tpos + lt;
    s1pos = wpos + lw;
    s2pos = s1pos + ls1;
    if upper {
        ldw = kd;
        lds2 = kd;
    } else {
        ldw = n;
        lds2 = n;
    }

    t_raw := _hetrd_he2hb_raw_slice(work.data + (tpos - 1), ldt*kd);
    t_view := array2d_view(t_raw, ldt, row_major=false);
    _hetrd_he2hb_laset(#char "A", ldt, kd, zero_t, zero_t, *t_view, ldt, MT);

    if upper {
        i = 1;
        while i <= n - kd {
            pn = n - i - kd + 1;
            pk = min(n - i - kd + 1, kd);

            a_block_raw := _hetrd_he2hb_matrix_storage(a.*, lda, i, i + kd, kd, pn);
            a_block := array2d_view(a_block_raw, lda, row_major=false);
            tau_view := _hetrd_he2hb_vector_slice(tau, i, pk);
            work_s2 := _hetrd_he2hb_vector_slice(work, s2pos, ls2);
            _hetrd_he2hb_gelqf(kd, pn, *a_block, lda, *tau_view, *work_s2, ls2, info, MT);

            for j:i..i + pk - 1 {
                lk = min(kd, n - j) + 1;
                _hetrd_he2hb_copy_upper_band(j, lk, a.*, lda, ab, ldab);
            }

            _hetrd_he2hb_laset(#char "L", pk, pk, zero_t, one_t, *a_block, lda, MT);
            t_panel_raw := _hetrd_he2hb_raw_slice(work.data + (tpos - 1), ldt*pk);
            t_panel := array2d_view(t_panel_raw, ldt, row_major=false);
            _hetrd_he2hb_larft(#char "F", #char "R", pn, pk, a_block, lda, tau_view, *t_panel, ldt, MT);

            s2_raw := _hetrd_he2hb_raw_slice(work.data + (s2pos - 1), lds2*pn);
            s2_view := array2d_view(s2_raw, lds2, row_major=false);
            _hetrd_he2hb_gemm(#char "C", #char "N", pk, pn, pk, one_t, t_panel, ldt, a_block, lda, zero_t, *s2_view, lds2, MT);

            trailing_raw := _hetrd_he2hb_matrix_storage(a.*, lda, i + kd, i + kd, pn, pn);
            trailing := array2d_view(trailing_raw, lda, row_major=false);
            w_raw := _hetrd_he2hb_raw_slice(work.data + (wpos - 1), ldw*pn);
            w_view := array2d_view(w_raw, ldw, row_major=false);
            _hetrd_he2hb_hemm(#char "R", uplo, pk, pn, one_t, trailing, lda, s2_view, lds2, zero_t, *w_view, ldw, MT);

            s1_raw := _hetrd_he2hb_raw_slice(work.data + (s1pos - 1), lds1*pk);
            s1_view := array2d_view(s1_raw, lds1, row_major=false);
            _hetrd_he2hb_gemm(#char "N", #char "C", pk, pk, pn, one_t, w_view, ldw, s2_view, lds2, zero_t, *s1_view, lds1, MT);
            _hetrd_he2hb_gemm(#char "N", #char "N", pk, pn, pk, -half_t, s1_view, lds1, a_block, lda, one_t, *w_view, ldw, MT);
            _hetrd_he2hb_her2k(uplo, #char "C", pn, pk, -one_t, a_block, lda, w_view, ldw, rone, *trailing, lda, RT, MT);

            i += kd;
        }

        for j:n - kd + 1..n {
            lk = min(kd, n - j) + 1;
            _hetrd_he2hb_copy_upper_band(j, lk, a.*, lda, ab, ldab);
        }
    } else {
        i = 1;
        while i <= n - kd {
            pn = n - i - kd + 1;
            pk = min(n - i - kd + 1, kd);

            a_block_raw := _hetrd_he2hb_matrix_storage(a.*, lda, i + kd, i, pn, kd);
            a_block := array2d_view(a_block_raw, lda, row_major=false);
            tau_view := _hetrd_he2hb_vector_slice(tau, i, pk);
            work_s2 := _hetrd_he2hb_vector_slice(work, s2pos, ls2);
            _hetrd_he2hb_geqrf(pn, kd, *a_block, lda, *tau_view, *work_s2, ls2, info, MT);

            for j:i..i + pk - 1 {
                lk = min(kd, n - j) + 1;
                _hetrd_he2hb_copy_lower_band(j, lk, a.*, ab);
            }

            _hetrd_he2hb_laset(#char "U", pk, pk, zero_t, one_t, *a_block, lda, MT);
            t_panel_raw := _hetrd_he2hb_raw_slice(work.data + (tpos - 1), ldt*pk);
            t_panel := array2d_view(t_panel_raw, ldt, row_major=false);
            _hetrd_he2hb_larft(#char "F", #char "C", pn, pk, a_block, lda, tau_view, *t_panel, ldt, MT);

            s2_raw := _hetrd_he2hb_raw_slice(work.data + (s2pos - 1), lds2*pk);
            s2_view := array2d_view(s2_raw, lds2, row_major=false);
            _hetrd_he2hb_gemm(#char "N", #char "N", pn, pk, pk, one_t, a_block, lda, t_panel, ldt, zero_t, *s2_view, lds2, MT);

            trailing_raw := _hetrd_he2hb_matrix_storage(a.*, lda, i + kd, i + kd, pn, pn);
            trailing := array2d_view(trailing_raw, lda, row_major=false);
            w_raw := _hetrd_he2hb_raw_slice(work.data + (wpos - 1), ldw*pk);
            w_view := array2d_view(w_raw, ldw, row_major=false);
            _hetrd_he2hb_hemm(#char "L", uplo, pn, pk, one_t, trailing, lda, s2_view, lds2, zero_t, *w_view, ldw, MT);

            s1_raw := _hetrd_he2hb_raw_slice(work.data + (s1pos - 1), lds1*pk);
            s1_view := array2d_view(s1_raw, lds1, row_major=false);
            _hetrd_he2hb_gemm(#char "C", #char "N", pk, pk, pn, one_t, s2_view, lds2, w_view, ldw, zero_t, *s1_view, lds1, MT);
            _hetrd_he2hb_gemm(#char "N", #char "N", pn, pk, pk, -half_t, a_block, lda, s1_view, lds1, one_t, *w_view, ldw, MT);
            _hetrd_he2hb_her2k(uplo, #char "N", pn, pk, -one_t, a_block, lda, w_view, ldw, rone, *trailing, lda, RT, MT);

            i += kd;
        }

        for j:n - kd + 1..n {
            lk = min(kd, n - j) + 1;
            _hetrd_he2hb_copy_lower_band(j, lk, a.*, ab);
        }
    }

    fset(work, 1, _hetrd_he2hb_complex_from_real(_hetrd_he2hb_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
}
