// ? LAPACK *gelqt()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgelqt/cgelqt/dgelqt/zgelqt
// *GELQT computes a blocked compact-WY LQ factorization.

sgelqt :: #bake_arguments _gelqt(RT=float32, MT=float32);
dgelqt :: #bake_arguments _gelqt(RT=float64, MT=float64);
cgelqt :: #bake_arguments _gelqt(RT=float32, MT=C32);
zgelqt :: #bake_arguments _gelqt(RT=float64, MT=C64);

#scope_file
_gelqt_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gelqt_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gelqt_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gelqt_larfb :: (
    $$side: u8,
    $$trans: u8,
    $$direct: u8,
    $$storev: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    v: $V/ArrayType2D,
    $$ldv: int,
    t: $T/ArrayType2D,
    $$ldt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType2D,
    $$ldwork: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else if MT == float64 {
        dlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else if MT == C32 {
        clarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else {
        zlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    }
}

_gelqt_gelqt3 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           T.data_type == MT;
} {
    #if MT == float32 {
        sgelqt3(m, n, a, lda, t, ldt, info);
    } else if MT == float64 {
        dgelqt3(m, n, a, lda, t, ldt, info);
    } else if MT == C32 {
        cgelqt3(m, n, a, lda, t, ldt, info);
    } else {
        zgelqt3(m, n, a, lda, t, ldt, info);
    }
}

_gelqt :: (
    $$m: int,
    $$n: int,
    $$mb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGELQT"
    } else if MT == float64 {
        "DGELQT"
    } else if MT == C32 {
        "CGELQT"
    } else {
        "ZGELQT"
    };

    i, ib, iinfo, k : int = ---;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if mb < 1 || (mb > min(m, n) && min(m, n) > 0) {
        info.* = -3;
    } else if lda < max(1, m) {
        info.* = -5;
    } else if ldt < mb {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    k = min(m, n);
    if k == 0 {
        return;
    }

    i = 1;
    while i <= k {
        ib = min(k - i + 1, mb);

        a_block_raw := _gelqt_matrix_storage(a.*, lda, i, i, ib, n - i + 1);
        a_block := array2d_view(a_block_raw, lda, row_major=false);
        t_block_raw := _gelqt_matrix_storage(t.*, ldt, 1, i, ib, ib);
        t_block := array2d_view(t_block_raw, ldt, row_major=false);
        _gelqt_gelqt3(ib, n - i + 1, *a_block, lda, *t_block, ldt, *iinfo, RT, MT);

        if i + ib <= m {
            c_raw := _gelqt_matrix_storage(a.*, lda, i + ib, i, m - i - ib + 1, n - i + 1);
            c_view := array2d_view(c_raw, lda, row_major=false);
            work_raw := _gelqt_raw_slice(work.data, (n - i + 1)*(m - i - ib + 1));
            work_view := array2d_view(work_raw, m - i - ib + 1, row_major=false);
            _gelqt_larfb(#char "R", #char "N", #char "F", #char "R", m - i - ib + 1, n - i + 1, ib,
                         a_block.*, lda, t_block.*, ldt, *c_view, lda, *work_view, m - i - ib + 1, MT);
        }

        i += mb;
    }
}
