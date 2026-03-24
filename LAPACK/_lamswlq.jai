// ? LAPACK *lamswlq()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slamswlq/dlamswlq/clamswlq/zlamswlq
// *LAMSWLQ applies Q from the short-wide LQ factorization.

slamswlq :: #bake_arguments _lamswlq(RT=float32, MT=float32);
dlamswlq :: #bake_arguments _lamswlq(RT=float64, MT=float64);
clamswlq :: #bake_arguments _lamswlq(RT=float32, MT=C32);
zlamswlq :: #bake_arguments _lamswlq(RT=float64, MT=C64);

#scope_file
_lamswlq_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lamswlq_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _lamswlq_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lamswlq_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_lamswlq_roundup_lwork :: ($RT: Type, x: RT) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_lamswlq_gemlqt :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    $$mb: int,
    a: $A/ArrayType2D,
    $$lda: int,
    t: $T/ArrayType2D,
    $$ldt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgemlqt(side, trans, m, n, k, mb, a, lda, t, ldt, c, ldc, work, info);
    } else if MT == float64 {
        dgemlqt(side, trans, m, n, k, mb, a, lda, t, ldt, c, ldc, work, info);
    } else if MT == C32 {
        cgemlqt(side, trans, m, n, k, mb, a, lda, t, ldt, c, ldc, work, info);
    } else {
        zgemlqt(side, trans, m, n, k, mb, a, lda, t, ldt, c, ldc, work, info);
    }
}

_lamswlq_tpmlqt :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$l: int,
    $$ib: int,
    $$mb: int,
    a: $A/ArrayType2D,
    $$lda: int,
    t: $T/ArrayType2D,
    $$ldt: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           V.data_type == MT &&
           B.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        stpmlqt(side, trans, m, n, ib, l, mb, a, lda, t, ldt, v, ldv, b, ldb, work, info);
    } else if MT == float64 {
        dtpmlqt(side, trans, m, n, ib, l, mb, a, lda, t, ldt, v, ldv, b, ldb, work, info);
    } else if MT == C32 {
        ctpmlqt(side, trans, m, n, ib, l, mb, a, lda, t, ldt, v, ldv, b, ldb, work, info);
    } else {
        ztpmlqt(side, trans, m, n, ib, l, mb, a, lda, t, ldt, v, ldv, b, ldb, work, info);
    }
}

_lamswlq :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    $$mb: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SLAMSWLQ"
    } else if MT == float64 {
        "DLAMSWLQ"
    } else if MT == C32 {
        "CLAMSWLQ"
    } else {
        "ZLAMSWLQ"
    };

    ctr, i, ii, kk, lw, lwmin, minmnk : int = ---;
    left, lquery, notran, right, tran : bool = ---;

    lquery = lwork == -1;
    #if MT == float32 || MT == float64 {
        tran = lsame(trans, #char "T");
    } else {
        tran = lsame(trans, #char "C");
    }
    notran = lsame(trans, #char "N");
    left = lsame(side, #char "L");
    right = lsame(side, #char "R");

    if left {
        lw = n * mb;
    } else {
        lw = m * mb;
    }

    minmnk = min(m, n, k);
    if minmnk == 0 {
        lwmin = 1;
    } else {
        lwmin = max(1, lw);
    }

    info.* = 0;
    if !left && !right {
        info.* = -1;
    } else if !tran && !notran {
        info.* = -2;
    } else if m < k {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if k < 0 {
        info.* = -5;
    } else if k < mb || mb < 1 {
        info.* = -6;
    } else if lda < max(1, k) {
        info.* = -9;
    } else if ldt < max(1, mb) {
        info.* = -11;
    } else if ldc < max(1, m) {
        info.* = -13;
    } else if lwork < lwmin && !lquery {
        info.* = -15;
    }

    if info.* == 0 {
        fset(work, 1, _lamswlq_complex_from_real(_lamswlq_roundup_lwork(RT, ncast(RT, lwmin)), RT, MT));
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if minmnk == 0 {
        return;
    }

    if (nb <= k) || (nb >= max(m, n, k)) {
        _lamswlq_gemlqt(side, trans, m, n, k, mb, a.*, lda, t.*, ldt, c, ldc, work, info, MT);
        return;
    }

    if left && tran {
        kk = mod(m - k, nb - k);
        ctr = (m - k) / (nb - k);

        if kk > 0 {
            ii = m - kk + 1;
            a_blk := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, ii, k, kk), lda, row_major=false);
            t_blk := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, ctr*k + 1, mb, k), ldt, row_major=false);
            top_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, k, n), ldc, row_major=false);
            bottom_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, ii, 1, kk, n), ldc, row_major=false);
            _lamswlq_tpmlqt(#char "L", trans, kk, n, 0, k, mb, a_blk, lda, t_blk, ldt, *top_blk, ldc, *bottom_blk, ldc, work, info, MT);
        } else {
            ii = m + 1;
        }

        i = ii - (nb - k);
        while i >= nb + 1 {
            ctr -= 1;
            a_blk := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, i, k, nb - k), lda, row_major=false);
            t_blk := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, ctr*k + 1, mb, k), ldt, row_major=false);
            top_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, k, n), ldc, row_major=false);
            bottom_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, i, 1, nb - k, n), ldc, row_major=false);
            _lamswlq_tpmlqt(#char "L", trans, nb - k, n, 0, k, mb, a_blk, lda, t_blk, ldt, *top_blk, ldc, *bottom_blk, ldc, work, info, MT);
            i -= nb - k;
        }

        first_c := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, nb, n), ldc, row_major=false);
        first_a := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, 1, k, nb), lda, row_major=false);
        first_t := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, 1, mb, k), ldt, row_major=false);
        _lamswlq_gemlqt(#char "L", trans, nb, n, k, mb, first_a, lda, first_t, ldt, *first_c, ldc, work, info, MT);
    } else if left && notran {
        kk = mod(m - k, nb - k);
        ii = m - kk + 1;
        ctr = 1;

        first_c := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, nb, n), ldc, row_major=false);
        first_a := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, 1, k, nb), lda, row_major=false);
        first_t := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, 1, mb, k), ldt, row_major=false);
        _lamswlq_gemlqt(#char "L", trans, nb, n, k, mb, first_a, lda, first_t, ldt, *first_c, ldc, work, info, MT);

        i = nb + 1;
        while i <= ii - nb + k {
            a_blk := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, i, k, nb - k), lda, row_major=false);
            t_blk := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, ctr*k + 1, mb, k), ldt, row_major=false);
            top_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, k, n), ldc, row_major=false);
            bottom_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, i, 1, nb - k, n), ldc, row_major=false);
            _lamswlq_tpmlqt(#char "L", trans, nb - k, n, 0, k, mb, a_blk, lda, t_blk, ldt, *top_blk, ldc, *bottom_blk, ldc, work, info, MT);
            ctr += 1;
            i += nb - k;
        }
        if ii <= m {
            a_blk := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, ii, k, kk), lda, row_major=false);
            t_blk := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, ctr*k + 1, mb, k), ldt, row_major=false);
            top_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, k, n), ldc, row_major=false);
            bottom_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, ii, 1, kk, n), ldc, row_major=false);
            _lamswlq_tpmlqt(#char "L", trans, kk, n, 0, k, mb, a_blk, lda, t_blk, ldt, *top_blk, ldc, *bottom_blk, ldc, work, info, MT);
        }
    } else if right && notran {
        kk = mod(n - k, nb - k);
        ctr = (n - k) / (nb - k);
        if kk > 0 {
            ii = n - kk + 1;
            a_blk := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, ii, k, kk), lda, row_major=false);
            t_blk := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, ctr*k + 1, mb, k), ldt, row_major=false);
            left_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, m, k), ldc, row_major=false);
            right_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, ii, m, kk), ldc, row_major=false);
            _lamswlq_tpmlqt(#char "R", trans, m, kk, 0, k, mb, a_blk, lda, t_blk, ldt, *left_blk, ldc, *right_blk, ldc, work, info, MT);
        } else {
            ii = n + 1;
        }

        i = ii - (nb - k);
        while i >= nb + 1 {
            ctr -= 1;
            a_blk := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, i, k, nb - k), lda, row_major=false);
            t_blk := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, ctr*k + 1, mb, k), ldt, row_major=false);
            left_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, m, k), ldc, row_major=false);
            right_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, i, m, nb - k), ldc, row_major=false);
            _lamswlq_tpmlqt(#char "R", trans, m, nb - k, 0, k, mb, a_blk, lda, t_blk, ldt, *left_blk, ldc, *right_blk, ldc, work, info, MT);
            i -= nb - k;
        }

        first_c := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, m, nb), ldc, row_major=false);
        first_a := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, 1, k, nb), lda, row_major=false);
        first_t := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, 1, mb, k), ldt, row_major=false);
        _lamswlq_gemlqt(#char "R", trans, m, nb, k, mb, first_a, lda, first_t, ldt, *first_c, ldc, work, info, MT);
    } else if right && tran {
        kk = mod(n - k, nb - k);
        ii = n - kk + 1;
        ctr = 1;

        first_c := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, m, nb), ldc, row_major=false);
        first_a := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, 1, k, nb), lda, row_major=false);
        first_t := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, 1, mb, k), ldt, row_major=false);
        _lamswlq_gemlqt(#char "R", trans, m, nb, k, mb, first_a, lda, first_t, ldt, *first_c, ldc, work, info, MT);

        i = nb + 1;
        while i <= ii - nb + k {
            a_blk := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, i, k, nb - k), lda, row_major=false);
            t_blk := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, ctr*k + 1, mb, k), ldt, row_major=false);
            left_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, m, k), ldc, row_major=false);
            right_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, i, m, nb - k), ldc, row_major=false);
            _lamswlq_tpmlqt(#char "R", trans, m, nb - k, 0, k, mb, a_blk, lda, t_blk, ldt, *left_blk, ldc, *right_blk, ldc, work, info, MT);
            ctr += 1;
            i += nb - k;
        }
        if ii <= n {
            a_blk := array2d_view(_lamswlq_matrix_storage(a.*, lda, 1, ii, k, kk), lda, row_major=false);
            t_blk := array2d_view(_lamswlq_matrix_storage(t.*, ldt, 1, ctr*k + 1, mb, k), ldt, row_major=false);
            left_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, 1, m, k), ldc, row_major=false);
            right_blk := array2d_view(_lamswlq_matrix_storage(c.*, ldc, 1, ii, m, kk), ldc, row_major=false);
            _lamswlq_tpmlqt(#char "R", trans, m, kk, 0, k, mb, a_blk, lda, t_blk, ldt, *left_blk, ldc, *right_blk, ldc, work, info, MT);
        }
    }

    fset(work, 1, _lamswlq_complex_from_real(_lamswlq_roundup_lwork(RT, ncast(RT, lwmin)), RT, MT));
}
