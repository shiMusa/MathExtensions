// ? LAPACK *geqrt()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeqrt/cgeqrt/dgeqrt/zgeqrt
// *GEQRT computes a blocked compact-WY QR factorization.

sgeqrt :: #bake_arguments _geqrt(RT=float32, MT=float32);
dgeqrt :: #bake_arguments _geqrt(RT=float64, MT=float64);
cgeqrt :: #bake_arguments _geqrt(RT=float32, MT=C32);
zgeqrt :: #bake_arguments _geqrt(RT=float64, MT=C64);

#scope_file
_geqrt_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geqrt_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _geqrt_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geqrt_larfb :: (
    $$side: u8, $$trans: u8, $$direct: u8, $$storev: u8,
    $$m: int, $$n: int, $$k: int,
    v: $V/ArrayType2D, $$ldv: int, t: $T/ArrayType2D, $$ldt: int,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType2D, $$ldwork: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT && T.data_type == MT && C.data_type == MT && WORK.data_type == MT;
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

_geqrt_geqrt2 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, t: *$T/ArrayType2D, $$ldt: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && T.data_type == MT;
} {
    #if MT == float32 {
        sgeqrt2(m, n, a, lda, t, ldt, info);
    } else if MT == float64 {
        dgeqrt2(m, n, a, lda, t, ldt, info);
    } else if MT == C32 {
        cgeqrt2(m, n, a, lda, t, ldt, info);
    } else {
        zgeqrt2(m, n, a, lda, t, ldt, info);
    }
}

_geqrt_geqrt3 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, t: *$T/ArrayType2D, $$ldt: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && T.data_type == MT;
} {
    #if MT == float32 {
        sgeqrt3(m, n, a, lda, t, ldt, info);
    } else if MT == float64 {
        dgeqrt3(m, n, a, lda, t, ldt, info);
    } else if MT == C32 {
        cgeqrt3(m, n, a, lda, t, ldt, info);
    } else {
        zgeqrt3(m, n, a, lda, t, ldt, info);
    }
}

_geqrt :: (
    $$m: int, $$n: int, $$nb: int, a: *$A/ArrayType2D, $$lda: int,
    t: *$T/ArrayType2D, $$ldt: int, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && T.data_type == MT && WORK.data_type == MT;
} {
    trans_t :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
    routine_name :: #run if MT == float32 {
        "SGEQRT"
    } else if MT == float64 {
        "DGEQRT"
    } else if MT == C32 {
        "CGEQRT"
    } else {
        "ZGEQRT"
    };

    i, ib, iinfo, k : int = ---;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nb < 1 || (nb > min(m, n) && min(m, n) > 0) {
        info.* = -3;
    } else if lda < max(1, m) {
        info.* = -5;
    } else if ldt < nb {
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
        ib = min(k - i + 1, nb);

        a_block_raw := _geqrt_matrix_storage(a.*, lda, i, i, m - i + 1, ib);
        a_block := array2d_view(a_block_raw, lda, row_major=false);
        t_block_raw := _geqrt_matrix_storage(t.*, ldt, 1, i, ib, ib);
        t_block := array2d_view(t_block_raw, ldt, row_major=false);
        _geqrt_geqrt3(m - i + 1, ib, *a_block, lda, *t_block, ldt, *iinfo, RT, MT);

        if i + ib <= n {
            c_raw := _geqrt_matrix_storage(a.*, lda, i, i + ib, m - i + 1, n - i - ib + 1);
            c_view := array2d_view(c_raw, lda, row_major=false);
            work_raw := _geqrt_raw_slice(work.data, (n - i - ib + 1)*ib);
            work_view := array2d_view(work_raw, n - i - ib + 1, row_major=false);
            _geqrt_larfb(#char "L", trans_t, #char "F", #char "C", m - i + 1, n - i - ib + 1, ib,
                         a_block.*, lda, t_block.*, ldt, *c_view, lda, *work_view, n - i - ib + 1, MT);
        }

        i += nb;
    }
}
