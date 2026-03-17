// ? LAPACK *gebrd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgebrd/cgebrd/dgebrd/zgebrd
// *GEBRD reduces a general matrix to bidiagonal form.

sgebrd :: #bake_arguments _gebrd(RT=float32, MT=float32);
dgebrd :: #bake_arguments _gebrd(RT=float64, MT=float64);
cgebrd :: #bake_arguments _gebrd(RT=float32, MT=C32);
zgebrd :: #bake_arguments _gebrd(RT=float64, MT=C64);

#scope_file
_gebrd_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gebrd_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gebrd_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gebrd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gebrd_store_real :: inline (
    a: *$A/ArrayType2D,
    $$row: int,
    $$col: int,
    value: $RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    fset(a, row, col, ncast(MT, value));
}

_gebrd_set_work1 :: inline (work: *$WORK/ArrayType1D, $$value: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           WORK.data_type == MT;
} {
    fset(work, 1, ncast(MT, ncast(RT, value)));
}

_gebrd_gemm :: (
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

_gebrd_labrd :: (
    $$m: int,
    $$n: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    tauq: *$TAUQ/ArrayType1D,
    taup: *$TAUP/ArrayType1D,
    x: *$X/ArrayType2D,
    $$ldx: int,
    y: *$Y/ArrayType2D,
    $$ldy: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           TAUQ.data_type == MT &&
           TAUP.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        slabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
    } else if MT == float64 {
        dlabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
    } else if MT == C32 {
        clabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
    } else {
        zlabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
    }
}

_gebrd_gebd2 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    tauq: *$TAUQ/ArrayType1D,
    taup: *$TAUP/ArrayType1D,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           TAUQ.data_type == MT &&
           TAUP.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgebd2(m, n, a, lda, d, e, tauq, taup, work, info);
    } else if MT == float64 {
        dgebd2(m, n, a, lda, d, e, tauq, taup, work, info);
    } else if MT == C32 {
        cgebd2(m, n, a, lda, d, e, tauq, taup, work, info);
    } else {
        zgebd2(m, n, a, lda, d, e, tauq, taup, work, info);
    }
}

_gebrd :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    tauq: *$TAUQ/ArrayType1D,
    taup: *$TAUP/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           TAUQ.data_type == MT &&
           TAUP.data_type == MT &&
           WORK.data_type == MT;
} {
    one_t :: #run one(MT);
    minus_one_t :: #run -one(MT);
    routine_name :: #run if MT == float32 {
        "SGEBRD"
    } else if MT == float64 {
        "DGEBRD"
    } else if MT == C32 {
        "CGEBRD"
    } else {
        "ZGEBRD"
    };
    trans_y :: #run if MT == float32 || MT == float64 {
        #char "T"
    } else {
        #char "C"
    };

    lquery : bool = ---;
    i, iinfo, j, ldwrkx, ldwrky, lwkmin, lwkopt, minmn, nb, nbmin, nx, ws : int = ---;

    info.* = 0;
    minmn = min(m, n);
    if minmn == 0 {
        lwkmin = 1;
        lwkopt = 1;
    } else {
        lwkmin = max(m, n);
        nb = max(1, ilaenv(1, routine_name, " "s, m, n, -1, -1));
        lwkopt = (m + n) * nb;
    }
    _gebrd_set_work1(work, lwkopt, RT, MT);

    lquery = lwork == -1;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    } else if lwork < lwkmin && !lquery {
        info.* = -10;
    }

    if info.* < 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if minmn == 0 {
        _gebrd_set_work1(work, 1, RT, MT);
        return;
    }

    ws = max(m, n);
    ldwrkx = m;
    ldwrky = n;

    if nb > 1 && nb < minmn {
        nx = max(nb, ilaenv(3, routine_name, " "s, m, n, -1, -1));
        if nx < minmn {
            ws = lwkopt;
            if lwork < ws {
                nbmin = ilaenv(2, routine_name, " "s, m, n, -1, -1);
                if lwork >= (m + n) * nbmin {
                    nb = lwork / (m + n);
                } else {
                    nb = 1;
                    nx = minmn;
                }
            }
        }
    } else {
        nx = minmn;
    }

    i = 1;
    while i <= minmn - nx {
        a_panel_raw := _gebrd_matrix_storage(a.*, lda, i, i, m - i + 1, n - i + 1);
        a_panel_view := array2d_view(a_panel_raw, lda, row_major=false);

        d_view := _gebrd_subvector(d.*, i, nb);
        e_view := _gebrd_subvector(e.*, i, nb);
        tauq_view := _gebrd_subvector(tauq.*, i, nb);
        taup_view := _gebrd_subvector(taup.*, i, nb);

        x_raw := _gebrd_subvector(work.*, 1, ldwrkx * nb);
        y_raw := _gebrd_subvector(work.*, ldwrkx * nb + 1, ldwrky * nb);
        x_view := array2d_view(x_raw, ldwrkx, row_major=false);
        y_view := array2d_view(y_raw, ldwrky, row_major=false);

        _gebrd_labrd(m - i + 1, n - i + 1, nb, *a_panel_view, lda, *d_view, *e_view,
                     *tauq_view, *taup_view, *x_view, ldwrkx, *y_view, ldwrky, RT, MT);

        a21_raw := _gebrd_matrix_storage(a.*, lda, i + nb, i, m - i - nb + 1, nb);
        y2_raw  := _gebrd_matrix_storage(y_view, ldwrky, nb + 1, 1, n - i - nb + 1, nb);
        a22_raw := _gebrd_matrix_storage(a.*, lda, i + nb, i + nb, m - i - nb + 1, n - i - nb + 1);
        a21_view := array2d_view(a21_raw, lda, row_major=false);
        y2_view  := array2d_view(y2_raw, ldwrky, row_major=false);
        a22_view := array2d_view(a22_raw, lda, row_major=false);

        _gebrd_gemm(#char "N", trans_y, m - i - nb + 1, n - i - nb + 1, nb, minus_one_t,
                    a21_view, lda, y2_view, ldwrky, one_t, *a22_view, lda, MT);

        x2_raw  := _gebrd_matrix_storage(x_view, ldwrkx, nb + 1, 1, m - i - nb + 1, nb);
        a12_raw := _gebrd_matrix_storage(a.*, lda, i, i + nb, nb, n - i - nb + 1);
        x2_view  := array2d_view(x2_raw, ldwrkx, row_major=false);
        a12_view := array2d_view(a12_raw, lda, row_major=false);

        _gebrd_gemm(#char "N", #char "N", m - i - nb + 1, n - i - nb + 1, nb, minus_one_t,
                    x2_view, ldwrkx, a12_view, lda, one_t, *a22_view, lda, MT);

        if m >= n {
            for j:i..i + nb - 1 {
                _gebrd_store_real(a, j, j, fget(d.*, j), RT, MT);
                _gebrd_store_real(a, j, j + 1, fget(e.*, j), RT, MT);
            }
        } else {
            for j:i..i + nb - 1 {
                _gebrd_store_real(a, j, j, fget(d.*, j), RT, MT);
                _gebrd_store_real(a, j + 1, j, fget(e.*, j), RT, MT);
            }
        }

        i += nb;
    }

    a_rem_raw := _gebrd_matrix_storage(a.*, lda, i, i, m - i + 1, n - i + 1);
    a_rem_view := array2d_view(a_rem_raw, lda, row_major=false);
    d_rem := _gebrd_subvector(d.*, i, minmn - i + 1);
    e_rem := _gebrd_subvector(e.*, i, max(0, minmn - i));
    tauq_rem := _gebrd_subvector(tauq.*, i, minmn - i + 1);
    taup_rem := _gebrd_subvector(taup.*, i, minmn - i + 1);
    _gebrd_gebd2(m - i + 1, n - i + 1, *a_rem_view, lda, *d_rem, *e_rem, *tauq_rem,
                 *taup_rem, work, *iinfo, RT, MT);

    _gebrd_set_work1(work, ws, RT, MT);
}
