// ? LAPACK *gemlqt()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgemlqt/cgemlqt/dgemlqt/zgemlqt
// *GEMLQT applies an LQ compact-WY block reflector to a matrix.

sgemlqt :: #bake_arguments _gemlqt(RT=float32, MT=float32);
dgemlqt :: #bake_arguments _gemlqt(RT=float64, MT=float64);
cgemlqt :: #bake_arguments _gemlqt(RT=float32, MT=C32);
zgemlqt :: #bake_arguments _gemlqt(RT=float64, MT=C64);

#scope_file
_gemlqt_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gemlqt_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gemlqt_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gemlqt_larfb :: (
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

_gemlqt :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    $$mb: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGEMLQT"
    } else if MT == float64 {
        "DGEMLQT"
    } else if MT == C32 {
        "CGEMLQT"
    } else {
        "ZGEMLQT"
    };

    left, notran, right, tran : bool = ---;
    i, ib, kf, ldwork, q : int = ---;

    info.* = 0;
    left   = lsame(side,  #char "L");
    right  = lsame(side,  #char "R");
    #if MT == float32 || MT == float64 {
        tran   = lsame(trans, #char "T");
    } else {
        tran   = lsame(trans, #char "C");
    }
    notran = lsame(trans, #char "N");

    if left {
        ldwork = max(1, n);
        q = m;
    } else if right {
        ldwork = max(1, m);
        q = n;
    } else {
        ldwork = 1;
        q = 0;
    }

    if !left && !right {
        info.* = -1;
    } else if !tran && !notran {
        info.* = -2;
    } else if m < 0 {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if k < 0 || k > q {
        info.* = -5;
    } else if mb < 1 || (mb > k && k > 0) {
        info.* = -6;
    } else if ldv < max(1, k) {
        info.* = -8;
    } else if ldt < mb {
        info.* = -10;
    } else if ldc < max(1, m) {
        info.* = -12;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m == 0 || n == 0 || k == 0 {
        return;
    }

    if left && notran {
        for i:1..k..mb {
            ib = min(mb, k - i + 1);
            v_blk := array2d_view(_gemlqt_matrix_storage(v.*, ldv, i, i, ib, q - i + 1), ldv, row_major=false);
            t_blk := array2d_view(_gemlqt_matrix_storage(t.*, ldt, 1, i, mb, ib), ldt, row_major=false);
            c_blk := array2d_view(_gemlqt_matrix_storage(c.*, ldc, i, 1, m - i + 1, n), ldc, row_major=false);
            work_blk := array2d_view(_gemlqt_raw_slice(work.data, ldwork*ib), ldwork, row_major=false);
            #if MT == float32 || MT == float64 {
                _gemlqt_larfb(#char "L", #char "T", #char "F", #char "R", m - i + 1, n, ib, v_blk.*, ldv, t_blk.*, ldt, *c_blk, ldc, *work_blk, ldwork, MT);
            } else {
                _gemlqt_larfb(#char "L", #char "C", #char "F", #char "R", m - i + 1, n, ib, v_blk.*, ldv, t_blk.*, ldt, *c_blk, ldc, *work_blk, ldwork, MT);
            }
        }
    } else if right && tran {
        for i:1..k..mb {
            ib = min(mb, k - i + 1);
            v_blk := array2d_view(_gemlqt_matrix_storage(v.*, ldv, i, i, ib, q - i + 1), ldv, row_major=false);
            t_blk := array2d_view(_gemlqt_matrix_storage(t.*, ldt, 1, i, mb, ib), ldt, row_major=false);
            c_blk := array2d_view(_gemlqt_matrix_storage(c.*, ldc, 1, i, m, n - i + 1), ldc, row_major=false);
            work_blk := array2d_view(_gemlqt_raw_slice(work.data, ldwork*ib), ldwork, row_major=false);
            _gemlqt_larfb(#char "R", #char "N", #char "F", #char "R", m, n - i + 1, ib, v_blk.*, ldv, t_blk.*, ldt, *c_blk, ldc, *work_blk, ldwork, MT);
        }
    } else if left && tran {
        kf = ((k - 1)/mb)*mb + 1;
        i = kf;
        while i >= 1 {
            ib = min(mb, k - i + 1);
            v_blk := array2d_view(_gemlqt_matrix_storage(v.*, ldv, i, i, ib, q - i + 1), ldv, row_major=false);
            t_blk := array2d_view(_gemlqt_matrix_storage(t.*, ldt, 1, i, mb, ib), ldt, row_major=false);
            c_blk := array2d_view(_gemlqt_matrix_storage(c.*, ldc, i, 1, m - i + 1, n), ldc, row_major=false);
            work_blk := array2d_view(_gemlqt_raw_slice(work.data, ldwork*ib), ldwork, row_major=false);
            _gemlqt_larfb(#char "L", #char "N", #char "F", #char "R", m - i + 1, n, ib, v_blk.*, ldv, t_blk.*, ldt, *c_blk, ldc, *work_blk, ldwork, MT);
            i -= mb;
        }
    } else if right && notran {
        kf = ((k - 1)/mb)*mb + 1;
        i = kf;
        while i >= 1 {
            ib = min(mb, k - i + 1);
            v_blk := array2d_view(_gemlqt_matrix_storage(v.*, ldv, i, i, ib, q - i + 1), ldv, row_major=false);
            t_blk := array2d_view(_gemlqt_matrix_storage(t.*, ldt, 1, i, mb, ib), ldt, row_major=false);
            c_blk := array2d_view(_gemlqt_matrix_storage(c.*, ldc, 1, i, m, n - i + 1), ldc, row_major=false);
            work_blk := array2d_view(_gemlqt_raw_slice(work.data, ldwork*ib), ldwork, row_major=false);
            #if MT == float32 || MT == float64 {
                _gemlqt_larfb(#char "R", #char "T", #char "F", #char "R", m, n - i + 1, ib, v_blk.*, ldv, t_blk.*, ldt, *c_blk, ldc, *work_blk, ldwork, MT);
            } else {
                _gemlqt_larfb(#char "R", #char "C", #char "F", #char "R", m, n - i + 1, ib, v_blk.*, ldv, t_blk.*, ldt, *c_blk, ldc, *work_blk, ldwork, MT);
            }
            i -= mb;
        }
    }
}
