// ? LAPACK *gemqr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgemqr/cgemqr/dgemqr/zgemqr
// *GEMQR applies Q from the modern GEQR interface.

sgemqr :: #bake_arguments _gemqr(RT=float32, MT=float32);
dgemqr :: #bake_arguments _gemqr(RT=float64, MT=float64);
cgemqr :: #bake_arguments _gemqr(RT=float32, MT=C32);
zgemqr :: #bake_arguments _gemqr(RT=float64, MT=C64);

#scope_file
_gemqr_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gemqr_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gemqr_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gemqr_roundup_lwork :: ($RT: Type, x: RT) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_gemqr_t_meta_int :: inline (t: $T/ArrayType1D, $$index: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           T.data_type == MT;
} {
    #if MT == float32 || MT == float64 {
        return cast(int) fget(t, index);
    } else {
        return cast(int) real(fget(t, index));
    }
}

_gemqr_lamtsqr :: (
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
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dlamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
    } else if MT == C32 {
        clamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
    } else {
        zlamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
    }
}

_gemqr :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType1D,
    $$tsize: int,
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
        "SGEMQR"
    } else if MT == float64 {
        "DGEMQR"
    } else if MT == C32 {
        "CGEMQR"
    } else {
        "ZGEMQR"
    };

    left, lquery, notran, right, tran : bool = ---;
    lw, lwmin, mb, mn, minmnk, nb, nblcks : int = ---;

    lquery = lwork == -1;
    #if MT == float32 || MT == float64 {
        tran   = lsame(trans, #char "T");
    } else {
        tran   = lsame(trans, #char "C");
    }
    notran = lsame(trans, #char "N");
    left   = lsame(side,  #char "L");
    right  = lsame(side,  #char "R");

    if tsize >= 3 {
        mb = _gemqr_t_meta_int(t.*, 2, RT, MT);
        nb = _gemqr_t_meta_int(t.*, 3, RT, MT);
    } else {
        mb = 1;
        nb = 1;
    }

    if left {
        lw = n*nb;
        mn = m;
    } else {
        lw = mb*nb;
        mn = n;
    }

    minmnk = min(m, n, k);
    if minmnk == 0 {
        lwmin = 1;
    } else {
        lwmin = max(1, lw);
    }

    if (mb > k) && (mn > k) {
        if mod(mn - k, mb - k) == 0 {
            nblcks = (mn - k) / (mb - k);
        } else {
            nblcks = (mn - k) / (mb - k) + 1;
        }
    } else {
        nblcks = 1;
    }
    _ = nblcks;

    info.* = 0;
    if !left && !right {
        info.* = -1;
    } else if !tran && !notran {
        info.* = -2;
    } else if m < 0 {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if k < 0 || k > mn {
        info.* = -5;
    } else if lda < max(1, mn) {
        info.* = -7;
    } else if tsize < 5 {
        info.* = -9;
    } else if ldc < max(1, m) {
        info.* = -11;
    } else if lwork < lwmin && !lquery {
        info.* = -13;
    }

    if info.* == 0 {
        fset(work, 1, _gemqr_complex_from_real(_gemqr_roundup_lwork(RT, ncast(RT, lwmin)), RT, MT));
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

    t_data := _gemqr_subvector(t.*, 6, max(0, tsize - 5));
    t_view := array2d_view(t_data, nb, row_major=false);

    if (left && m <= k) || (right && n <= k) || (mb <= k) || (mb >= max(m, n, k)) {
        _gemqrt(side, trans, m, n, k, nb, a, lda, *t_view, nb, c, ldc, work, info, RT, MT);
    } else {
        _gemqr_lamtsqr(side, trans, m, n, k, mb, nb, a, lda, *t_view, nb, c, ldc, work, lwork, info, MT);
    }

    fset(work, 1, _gemqr_complex_from_real(_gemqr_roundup_lwork(RT, ncast(RT, lwmin)), RT, MT));
}
