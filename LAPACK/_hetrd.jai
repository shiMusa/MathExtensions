// ? LAPACK *hetrd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrd/zhetrd
// *HETRD reduces a Hermitian matrix to real symmetric tridiagonal form.

chetrd :: #bake_arguments _hetrd(RT=float32, MT=C32);
zhetrd :: #bake_arguments _hetrd(RT=float64, MT=C64);

#scope_file
_hetrd_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrd_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetrd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrd_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetrd_raw_slice(x.data + (offset - 1), count);
}

_hetrd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrd_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrd_latrd :: (
    $$uplo: u8,
    $$n: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    w: *$W/ArrayType2D,
    $$ldw: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           E.data_type == RT &&
           TAU.data_type == MT &&
           W.data_type == MT;
} {
    #if MT == C32 {
        clatrd(uplo, n, nb, a, lda, e, tau, w, ldw);
    } else {
        zlatrd(uplo, n, nb, a, lda, e, tau, w, ldw);
    }
}

_hetrd_her2k :: (
    $$uplo: u8,
    $$trans: u8,
    $$n: int,
    $$k: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: $RT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == C32 {
        cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_hetrd :: (
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
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
           D.data_type == RT &&
           E.data_type == RT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    one_r :: #run one(RT);
    cone  :: #run _hetrd_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHETRD" } else { "ZHETRD" };

    i, iinfo, iws, j, kk, ldwork, lwkopt, nb, nbmin, nx : int = ---;
    lquery, upper : bool = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if lwork < 1 && !lquery {
        info.* = -9;
    }

    if info.* == 0 {
        nb = ilaenv(1, "CHETRD", cast(*u8) *uplo, n, -1, -1, -1);
        lwkopt = max(1, n*nb);
        fset(work, 1, _hetrd_complex_from_real(_hetrd_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        fset(work, 1, cone);
        return;
    }

    nx = n;
    iws = 1;
    ldwork = n;
    if nb > 1 && nb < n {
        nx = max(nb, ilaenv(3, "CHETRD", cast(*u8) *uplo, n, -1, -1, -1));
        if nx < n {
            ldwork = n;
            iws = ldwork*nb;
            if lwork < iws {
                nb = max(lwork / ldwork, 1);
                nbmin = ilaenv(2, "CHETRD", cast(*u8) *uplo, n, -1, -1, -1);
                if nb < nbmin {
                    nx = n;
                }
            }
        } else {
            nx = n;
        }
    } else {
        nb = 1;
    }

    if upper {
        kk = n - ((n - nx + nb - 1) / nb)*nb;
        for <i:kk + 1..n - nb + 1 {
            w_raw := _hetrd_raw_slice(work.data, ldwork*nb);
            w_view := array2d_view(w_raw, ldwork, row_major=false);
            _hetrd_latrd(uplo, i + nb - 1, nb, a, lda, e, tau, *w_view, ldwork, RT, MT);

            left_raw := _hetrd_matrix_storage(a.*, lda, 1, 1, i - 1, i - 1);
            left_view := array2d_view(left_raw, lda, row_major=false);
            v_raw := _hetrd_matrix_storage(a.*, lda, 1, i, i - 1, nb);
            v_view := array2d_view(v_raw, lda, row_major=false);
            _hetrd_her2k(uplo, #char "N", i - 1, nb, -cone, v_view, lda, w_view, ldwork, one_r, *left_view, lda, RT, MT);

            for j:i..i + nb - 1 {
                fset(a, j - 1, j, _hetrd_complex_from_real(fget(e.*, j - 1), RT, MT));
                fset(d, j, real(fget(a.*, j, j)));
            }
        }

        a_tail_raw := _hetrd_matrix_storage(a.*, lda, 1, 1, kk, kk);
        a_tail := array2d_view(a_tail_raw, lda, row_major=false);
        d_tail := _hetrd_vector_slice(d, 1, kk);
        e_tail := _hetrd_vector_slice(e, 1, max(0, kk - 1));
        tau_tail := _hetrd_vector_slice(tau, 1, max(0, kk - 1));
        #if MT == C32 {
            chetd2(uplo, kk, *a_tail, lda, *d_tail, *e_tail, *tau_tail, *iinfo);
        } else {
            zhetd2(uplo, kk, *a_tail, lda, *d_tail, *e_tail, *tau_tail, *iinfo);
        }
    } else {
        i = 1;
        while i <= n - nx {
            w_raw := _hetrd_raw_slice(work.data, ldwork*nb);
            w_view := array2d_view(w_raw, ldwork, row_major=false);
            a_sub_raw := _hetrd_matrix_storage(a.*, lda, i, i, n - i + 1, n - i + 1);
            a_sub := array2d_view(a_sub_raw, lda, row_major=false);
            e_sub := _hetrd_vector_slice(e, i, nb);
            tau_sub := _hetrd_vector_slice(tau, i, nb);
            _hetrd_latrd(uplo, n - i + 1, nb, *a_sub, lda, *e_sub, *tau_sub, *w_view, ldwork, RT, MT);

            trail_raw := _hetrd_matrix_storage(a.*, lda, i + nb, i + nb, n - i - nb + 1, n - i - nb + 1);
            trail_view := array2d_view(trail_raw, lda, row_major=false);
            v_raw := _hetrd_matrix_storage(a.*, lda, i + nb, i, n - i - nb + 1, nb);
            v_view := array2d_view(v_raw, lda, row_major=false);
            w_off_raw := _hetrd_raw_slice(work.data + nb, ldwork*nb - nb);
            w_off := array2d_view(w_off_raw, ldwork, row_major=false);
            _hetrd_her2k(uplo, #char "N", n - i - nb + 1, nb, -cone, v_view, lda, w_off, ldwork, one_r, *trail_view, lda, RT, MT);

            for j:i..i + nb - 1 {
                fset(a, j + 1, j, _hetrd_complex_from_real(fget(e.*, j), RT, MT));
                fset(d, j, real(fget(a.*, j, j)));
            }
            i += nb;
        }

        a_tail_raw := _hetrd_matrix_storage(a.*, lda, i, i, n - i + 1, n - i + 1);
        a_tail := array2d_view(a_tail_raw, lda, row_major=false);
        d_tail := _hetrd_vector_slice(d, i, n - i + 1);
        e_tail := _hetrd_vector_slice(e, i, max(0, n - i));
        tau_tail := _hetrd_vector_slice(tau, i, max(0, n - i));
        #if MT == C32 {
            chetd2(uplo, n - i + 1, *a_tail, lda, *d_tail, *e_tail, *tau_tail, *iinfo);
        } else {
            zhetd2(uplo, n - i + 1, *a_tail, lda, *d_tail, *e_tail, *tau_tail, *iinfo);
        }
    }

    fset(work, 1, _hetrd_complex_from_real(_hetrd_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
