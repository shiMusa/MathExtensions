// ? LAPACK *ggsvp3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggsvp3/cggsvp3/dggsvp3/zggsvp3
// *GGSVP3 computes orthogonal or unitary matrices U, V and Q for GSVD preprocessing.

sggsvp3 :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$p: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $$tola: float32, $$tolb: float32,
    k: *int, l: *int,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    iwork: *$IWORK/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int
) {
    _ggsvp3_real(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, lwork, info, RT=float32);
}

dggsvp3 :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$p: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $$tola: float64, $$tolb: float64,
    k: *int, l: *int,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    iwork: *$IWORK/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int
) {
    _ggsvp3_real(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, lwork, info, RT=float64);
}

cggsvp3 :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$p: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $$tola: float32, $$tolb: float32,
    k: *int, l: *int,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    iwork: *$IWORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int
) {
    _ggsvp3_complex(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info, RT=float32, MT=C32);
}

zggsvp3 :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$p: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $$tola: float64, $$tolb: float64,
    k: *int, l: *int,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    iwork: *$IWORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int
) {
    _ggsvp3_complex(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info, RT=float64, MT=C64);
}

#scope_file
_ggsvp3_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_ggsvp3_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _ggsvp3_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_ggsvp3_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && (MT == C32 || MT == C64);
} {
    return cpx(x, zero(RT));
}

_ggsvp3_query_to_int :: inline (x: $T, $RT: Type, $T: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (T == RT || T == complex_with_inner_type(RT));
} {
    #if T == float32 || T == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_ggsvp3_roundup_lwork :: inline (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_ggsvp3_laset :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: MT,
    a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        slaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == float64 {
        dlaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_ggsvp3_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT;
} {
    #if MT == float32 {
        slacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == float64 {
        dlacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_ggsvp3_lapmt :: (
    $$forwrd: bool, $$m: int, $$n: int, x: *$X/ArrayType2D, $$ldx: int,
    k: *$K/ArrayType1D, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && K.data_type == int;
} {
    #if MT == float32 {
        slapmt(forwrd, m, n, x, ldx, k);
    } else if MT == float64 {
        dlapmt(forwrd, m, n, x, ldx, k);
    } else if MT == C32 {
        clapmt(forwrd, m, n, x, ldx, k);
    } else {
        zlapmt(forwrd, m, n, x, ldx, k);
    }
}

_ggsvp3_geqp3_real :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && JPVT.data_type == int &&
           TAU.data_type == RT && WORK.data_type == RT;
} {
    #if RT == float32 {
        sgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info);
    } else {
        dgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info);
    }
}

_ggsvp3_geqp3_complex :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, rwork: *$RWORK/ArrayType1D, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT && JPVT.data_type == int &&
           RWORK.data_type == RT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
    } else {
        zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
    }
}

_ggsvp3_geqr2 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sgeqr2(m, n, a, lda, tau, work, info);
    } else if MT == float64 {
        dgeqr2(m, n, a, lda, tau, work, info);
    } else if MT == C32 {
        cgeqr2(m, n, a, lda, tau, work, info);
    } else {
        zgeqr2(m, n, a, lda, tau, work, info);
    }
}

_ggsvp3_gerq2 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sgerq2(m, n, a, lda, tau, work, info);
    } else if MT == float64 {
        dgerq2(m, n, a, lda, tau, work, info);
    } else if MT == C32 {
        cgerq2(m, n, a, lda, tau, work, info);
    } else {
        zgerq2(m, n, a, lda, tau, work, info);
    }
}

_ggsvp3_g2r :: (
    $$m: int, $$n: int, $$k: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sorg2r(m, n, k, a, lda, tau, work, info);
    } else if MT == float64 {
        dorg2r(m, n, k, a, lda, tau, work, info);
    } else if MT == C32 {
        cung2r(m, n, k, a, lda, tau, work, info);
    } else {
        zung2r(m, n, k, a, lda, tau, work, info);
    }
}

_ggsvp3_m2r :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$k: int,
    a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT &&
           C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    } else if MT == float64 {
        dorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    } else if MT == C32 {
        cunm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    } else {
        zunm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
}

_ggsvp3_mr2 :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$k: int,
    a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT &&
           C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sormr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    } else if MT == float64 {
        dormr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    } else if MT == C32 {
        cunmr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    } else {
        zunmr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
    }
}

_ggsvp3_real :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$p: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $$tola: $RT, $$tolb: RT,
    k: *int, l: *int,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    iwork: *$IWORK/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT &&
           U.data_type == RT && V.data_type == RT && Q.data_type == RT &&
           IWORK.data_type == int && TAU.data_type == RT && WORK.data_type == RT;
} {
    zero_t :: #run zero(RT);
    one_t  :: #run one(RT);
    routine_name :: #run if RT == float32 then "SGGSVP3" else "DGGSVP3";
    trans_t :: #run #char "T";

    wantu, wantv, wantq, forwrd, lquery : bool = ---;
    i, j, lwkopt : int = ---;

    info.* = 0;
    wantu = lsame(jobu, #char "U");
    wantv = lsame(jobv, #char "V");
    wantq = lsame(jobq, #char "Q");
    forwrd = true;
    lquery = lwork == -1;
    lwkopt = 1;

    if !(wantu || lsame(jobu, #char "N")) {
        info.* = -1;
    } else if !(wantv || lsame(jobv, #char "N")) {
        info.* = -2;
    } else if !(wantq || lsame(jobq, #char "N")) {
        info.* = -3;
    } else if m < 0 {
        info.* = -4;
    } else if p < 0 {
        info.* = -5;
    } else if n < 0 {
        info.* = -6;
    } else if lda < max(1, m) {
        info.* = -8;
    } else if ldb < max(1, p) {
        info.* = -10;
    } else if ldu < 1 || (wantu && ldu < m) {
        info.* = -16;
    } else if ldv < 1 || (wantv && ldv < p) {
        info.* = -18;
    } else if ldq < 1 || (wantq && ldq < n) {
        info.* = -20;
    } else if lwork < 1 && !lquery {
        info.* = -24;
    }

    if info.* == 0 {
        _ggsvp3_geqp3_real(p, n, b, ldb, iwork, tau, work, -1, info, RT);
        lwkopt = _ggsvp3_query_to_int(fget(work.*, 1), RT, RT);
        if wantv {
            lwkopt = max(lwkopt, p);
        }
        lwkopt = max(lwkopt, min(n, p));
        lwkopt = max(lwkopt, m);
        if wantq {
            lwkopt = max(lwkopt, n);
        }

        _ggsvp3_geqp3_real(m, n, a, lda, iwork, tau, work, -1, info, RT);
        lwkopt = max(lwkopt, _ggsvp3_query_to_int(fget(work.*, 1), RT, RT));
        lwkopt = max(1, lwkopt);
        fset(work, 1, _ggsvp3_roundup_lwork(ncast(RT, lwkopt), RT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lquery return;

    for i:1..n {
        fset(iwork, i, 0);
    }
    _ggsvp3_geqp3_real(p, n, b, ldb, iwork, tau, work, lwork, info, RT);

    _ggsvp3_lapmt(forwrd, m, n, a, lda, iwork, RT);

    l.* = 0;
    for i:1..min(p, n) {
        if abs(fget(b.*, i, i)) > tolb {
            l.* += 1;
        }
    }

    if wantv {
        _ggsvp3_laset(#char "F", p, p, zero_t, zero_t, v, ldv, RT);
        if p > 1 {
            b_sub_raw := _ggsvp3_matrix_storage(b.*, ldb, 2, 1, p - 1, n);
            b_sub := array2d_view(b_sub_raw, ldb, row_major=false);
            v_sub_raw := _ggsvp3_matrix_storage(v.*, ldv, 2, 1, p - 1, n);
            v_sub := array2d_view(v_sub_raw, ldv, row_major=false);
            _ggsvp3_lacpy(#char "L", p - 1, n, *b_sub, ldb, *v_sub, ldv, RT);
        }
        _ggsvp3_g2r(p, p, min(p, n), v, ldv, tau, work, info, RT, RT);
    }

    for j:1..l.* - 1 {
        for i:j + 1..l.* {
            fset(b, i, j, zero_t);
        }
    }
    if p > l.* {
        b_tail_raw := _ggsvp3_matrix_storage(b.*, ldb, l.* + 1, 1, p - l.*, n);
        b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
        _ggsvp3_laset(#char "F", p - l.*, n, zero_t, zero_t, *b_tail, ldb, RT);
    }

    if wantq {
        _ggsvp3_laset(#char "F", n, n, zero_t, one_t, q, ldq, RT);
        _ggsvp3_lapmt(forwrd, n, n, q, ldq, iwork, RT);
    }

    if p >= l.* && n != l.* {
        _ggsvp3_gerq2(l.*, n, b, ldb, tau, work, info, RT);
        _ggsvp3_mr2(#char "R", trans_t, m, n, l.*, b, ldb, tau, a, lda, work, info, RT, RT);
        if wantq {
            _ggsvp3_mr2(#char "R", trans_t, n, n, l.*, b, ldb, tau, q, ldq, work, info, RT, RT);
        }

        _ggsvp3_laset(#char "F", l.*, n - l.*, zero_t, zero_t, b, ldb, RT);
        for j:n - l.* + 1..n {
            for i:j - n + l.* + 1..l.* {
                fset(b, i, j, zero_t);
            }
        }
    }

    for i:1..n - l.* {
        fset(iwork, i, 0);
    }
    _ggsvp3_geqp3_real(m, n - l.*, a, lda, iwork, tau, work, lwork, info, RT);

    k.* = 0;
    for i:1..min(m, n - l.*) {
        if abs(fget(a.*, i, i)) > tola {
            k.* += 1;
        }
    }

    a12_raw := _ggsvp3_matrix_storage(a.*, lda, 1, n - l.* + 1, m, l.*);
    a12 := array2d_view(a12_raw, lda, row_major=false);
    _ggsvp3_m2r(#char "L", trans_t, m, l.*, min(m, n - l.*), a, lda, tau, *a12, lda, work, info, RT, RT);

    if wantu {
        _ggsvp3_laset(#char "F", m, m, zero_t, zero_t, u, ldu, RT);
        if m > 1 {
            a_sub_raw := _ggsvp3_matrix_storage(a.*, lda, 2, 1, m - 1, n - l.*);
            a_sub := array2d_view(a_sub_raw, lda, row_major=false);
            u_sub_raw := _ggsvp3_matrix_storage(u.*, ldu, 2, 1, m - 1, n - l.*);
            u_sub := array2d_view(u_sub_raw, ldu, row_major=false);
            _ggsvp3_lacpy(#char "L", m - 1, n - l.*, *a_sub, lda, *u_sub, ldu, RT);
        }
        _ggsvp3_g2r(m, m, min(m, n - l.*), u, ldu, tau, work, info, RT, RT);
    }

    if wantq && n - l.* > 0 {
        q_left_raw := _ggsvp3_matrix_storage(q.*, ldq, 1, 1, n, n - l.*);
        q_left := array2d_view(q_left_raw, ldq, row_major=false);
        _ggsvp3_lapmt(forwrd, n, n - l.*, *q_left, ldq, iwork, RT);
    }

    for j:1..k.* - 1 {
        for i:j + 1..k.* {
            fset(a, i, j, zero_t);
        }
    }
    if m > k.* {
        a_tail_raw := _ggsvp3_matrix_storage(a.*, lda, k.* + 1, 1, m - k.*, n - l.*);
        a_tail := array2d_view(a_tail_raw, lda, row_major=false);
        _ggsvp3_laset(#char "F", m - k.*, n - l.*, zero_t, zero_t, *a_tail, lda, RT);
    }

    if n - l.* > k.* {
        a11_raw := _ggsvp3_matrix_storage(a.*, lda, 1, 1, k.*, n - l.*);
        a11 := array2d_view(a11_raw, lda, row_major=false);
        _ggsvp3_gerq2(k.*, n - l.*, *a11, lda, tau, work, info, RT);

        if wantq {
            q_left_raw := _ggsvp3_matrix_storage(q.*, ldq, 1, 1, n, n - l.*);
            q_left := array2d_view(q_left_raw, ldq, row_major=false);
            _ggsvp3_mr2(#char "R", trans_t, n, n - l.*, k.*, *a11, lda, tau, *q_left, ldq, work, info, RT, RT);
        }

        _ggsvp3_laset(#char "F", k.*, n - l.* - k.*, zero_t, zero_t, *a11, lda, RT);
        for j:n - l.* - k.* + 1..n - l.* {
            for i:j - n + l.* + k.* + 1..k.* {
                fset(a, i, j, zero_t);
            }
        }
    }

    if m > k.* {
        a23_raw := _ggsvp3_matrix_storage(a.*, lda, k.* + 1, n - l.* + 1, m - k.*, l.*);
        a23 := array2d_view(a23_raw, lda, row_major=false);
        _ggsvp3_geqr2(m - k.*, l.*, *a23, lda, tau, work, info, RT);

        if wantu {
            u_right_raw := _ggsvp3_matrix_storage(u.*, ldu, 1, k.* + 1, m, m - k.*);
            u_right := array2d_view(u_right_raw, ldu, row_major=false);
            _ggsvp3_m2r(#char "R", #char "N", m, m - k.*, min(m - k.*, l.*), *a23, lda, tau, *u_right, ldu, work, info, RT, RT);
        }

        for j:n - l.* + 1..n {
            for i:j - n + k.* + l.* + 1..m {
                fset(a, i, j, zero_t);
            }
        }
    }

    fset(work, 1, _ggsvp3_roundup_lwork(ncast(RT, lwkopt), RT));
}

_ggsvp3_complex :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$p: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $$tola: $RT, $$tolb: RT,
    k: *int, l: *int,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    iwork: *$IWORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT &&
           U.data_type == MT && V.data_type == MT && Q.data_type == MT &&
           IWORK.data_type == int && RWORK.data_type == RT &&
           TAU.data_type == MT && WORK.data_type == MT;
} {
    zero_t :: #run _ggsvp3_complex_from_real(zero(RT), RT, MT);
    one_t  :: #run _ggsvp3_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 then "CGGSVP3" else "ZGGSVP3";
    trans_t :: #run #char "C";

    wantu, wantv, wantq, forwrd, lquery : bool = ---;
    i, j, lwkopt : int = ---;

    info.* = 0;
    wantu = lsame(jobu, #char "U");
    wantv = lsame(jobv, #char "V");
    wantq = lsame(jobq, #char "Q");
    forwrd = true;
    lquery = lwork == -1;
    lwkopt = 1;

    if !(wantu || lsame(jobu, #char "N")) {
        info.* = -1;
    } else if !(wantv || lsame(jobv, #char "N")) {
        info.* = -2;
    } else if !(wantq || lsame(jobq, #char "N")) {
        info.* = -3;
    } else if m < 0 {
        info.* = -4;
    } else if p < 0 {
        info.* = -5;
    } else if n < 0 {
        info.* = -6;
    } else if lda < max(1, m) {
        info.* = -8;
    } else if ldb < max(1, p) {
        info.* = -10;
    } else if ldu < 1 || (wantu && ldu < m) {
        info.* = -16;
    } else if ldv < 1 || (wantv && ldv < p) {
        info.* = -18;
    } else if ldq < 1 || (wantq && ldq < n) {
        info.* = -20;
    } else if lwork < 1 && !lquery {
        info.* = -24;
    }

    if info.* == 0 {
        _ggsvp3_geqp3_complex(p, n, b, ldb, iwork, rwork, tau, work, -1, info, RT, MT);
        lwkopt = _ggsvp3_query_to_int(fget(work.*, 1), RT, MT);
        if wantv {
            lwkopt = max(lwkopt, p);
        }
        lwkopt = max(lwkopt, min(n, p));
        lwkopt = max(lwkopt, m);
        if wantq {
            lwkopt = max(lwkopt, n);
        }

        _ggsvp3_geqp3_complex(m, n, a, lda, iwork, rwork, tau, work, -1, info, RT, MT);
        lwkopt = max(lwkopt, _ggsvp3_query_to_int(fget(work.*, 1), RT, MT));
        lwkopt = max(1, lwkopt);
        fset(work, 1, _ggsvp3_complex_from_real(ncast(RT, lwkopt), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lquery return;

    for i:1..n {
        fset(iwork, i, 0);
    }
    _ggsvp3_geqp3_complex(p, n, b, ldb, iwork, rwork, tau, work, lwork, info, RT, MT);

    _ggsvp3_lapmt(forwrd, m, n, a, lda, iwork, MT);

    l.* = 0;
    for i:1..min(p, n) {
        if abs(fget(b.*, i, i)) > tolb {
            l.* += 1;
        }
    }

    if wantv {
        _ggsvp3_laset(#char "F", p, p, zero_t, zero_t, v, ldv, MT);
        if p > 1 {
            b_sub_raw := _ggsvp3_matrix_storage(b.*, ldb, 2, 1, p - 1, n);
            b_sub := array2d_view(b_sub_raw, ldb, row_major=false);
            v_sub_raw := _ggsvp3_matrix_storage(v.*, ldv, 2, 1, p - 1, n);
            v_sub := array2d_view(v_sub_raw, ldv, row_major=false);
            _ggsvp3_lacpy(#char "L", p - 1, n, *b_sub, ldb, *v_sub, ldv, MT);
        }
        _ggsvp3_g2r(p, p, min(p, n), v, ldv, tau, work, info, RT, MT);
    }

    for j:1..l.* - 1 {
        for i:j + 1..l.* {
            fset(b, i, j, zero_t);
        }
    }
    if p > l.* {
        b_tail_raw := _ggsvp3_matrix_storage(b.*, ldb, l.* + 1, 1, p - l.*, n);
        b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
        _ggsvp3_laset(#char "F", p - l.*, n, zero_t, zero_t, *b_tail, ldb, MT);
    }

    if wantq {
        _ggsvp3_laset(#char "F", n, n, zero_t, one_t, q, ldq, MT);
        _ggsvp3_lapmt(forwrd, n, n, q, ldq, iwork, MT);
    }

    if p >= l.* && n != l.* {
        _ggsvp3_gerq2(l.*, n, b, ldb, tau, work, info, MT);
        _ggsvp3_mr2(#char "R", trans_t, m, n, l.*, b, ldb, tau, a, lda, work, info, RT, MT);
        if wantq {
            _ggsvp3_mr2(#char "R", trans_t, n, n, l.*, b, ldb, tau, q, ldq, work, info, RT, MT);
        }

        _ggsvp3_laset(#char "F", l.*, n - l.*, zero_t, zero_t, b, ldb, MT);
        for j:n - l.* + 1..n {
            for i:j - n + l.* + 1..l.* {
                fset(b, i, j, zero_t);
            }
        }
    }

    for i:1..n - l.* {
        fset(iwork, i, 0);
    }
    _ggsvp3_geqp3_complex(m, n - l.*, a, lda, iwork, rwork, tau, work, lwork, info, RT, MT);

    k.* = 0;
    for i:1..min(m, n - l.*) {
        if abs(fget(a.*, i, i)) > tola {
            k.* += 1;
        }
    }

    a12_raw := _ggsvp3_matrix_storage(a.*, lda, 1, n - l.* + 1, m, l.*);
    a12 := array2d_view(a12_raw, lda, row_major=false);
    _ggsvp3_m2r(#char "L", trans_t, m, l.*, min(m, n - l.*), a, lda, tau, *a12, lda, work, info, RT, MT);

    if wantu {
        _ggsvp3_laset(#char "F", m, m, zero_t, zero_t, u, ldu, MT);
        if m > 1 {
            a_sub_raw := _ggsvp3_matrix_storage(a.*, lda, 2, 1, m - 1, n - l.*);
            a_sub := array2d_view(a_sub_raw, lda, row_major=false);
            u_sub_raw := _ggsvp3_matrix_storage(u.*, ldu, 2, 1, m - 1, n - l.*);
            u_sub := array2d_view(u_sub_raw, ldu, row_major=false);
            _ggsvp3_lacpy(#char "L", m - 1, n - l.*, *a_sub, lda, *u_sub, ldu, MT);
        }
        _ggsvp3_g2r(m, m, min(m, n - l.*), u, ldu, tau, work, info, RT, MT);
    }

    if wantq && n - l.* > 0 {
        q_left_raw := _ggsvp3_matrix_storage(q.*, ldq, 1, 1, n, n - l.*);
        q_left := array2d_view(q_left_raw, ldq, row_major=false);
        _ggsvp3_lapmt(forwrd, n, n - l.*, *q_left, ldq, iwork, MT);
    }

    for j:1..k.* - 1 {
        for i:j + 1..k.* {
            fset(a, i, j, zero_t);
        }
    }
    if m > k.* {
        a_tail_raw := _ggsvp3_matrix_storage(a.*, lda, k.* + 1, 1, m - k.*, n - l.*);
        a_tail := array2d_view(a_tail_raw, lda, row_major=false);
        _ggsvp3_laset(#char "F", m - k.*, n - l.*, zero_t, zero_t, *a_tail, lda, MT);
    }

    if n - l.* > k.* {
        a11_raw := _ggsvp3_matrix_storage(a.*, lda, 1, 1, k.*, n - l.*);
        a11 := array2d_view(a11_raw, lda, row_major=false);
        _ggsvp3_gerq2(k.*, n - l.*, *a11, lda, tau, work, info, MT);

        if wantq {
            q_left_raw := _ggsvp3_matrix_storage(q.*, ldq, 1, 1, n, n - l.*);
            q_left := array2d_view(q_left_raw, ldq, row_major=false);
            _ggsvp3_mr2(#char "R", trans_t, n, n - l.*, k.*, *a11, lda, tau, *q_left, ldq, work, info, RT, MT);
        }

        _ggsvp3_laset(#char "F", k.*, n - l.* - k.*, zero_t, zero_t, *a11, lda, MT);
        for j:n - l.* - k.* + 1..n - l.* {
            for i:j - n + l.* + k.* + 1..k.* {
                fset(a, i, j, zero_t);
            }
        }
    }

    if m > k.* {
        a23_raw := _ggsvp3_matrix_storage(a.*, lda, k.* + 1, n - l.* + 1, m - k.*, l.*);
        a23 := array2d_view(a23_raw, lda, row_major=false);
        _ggsvp3_geqr2(m - k.*, l.*, *a23, lda, tau, work, info, MT);

        if wantu {
            u_right_raw := _ggsvp3_matrix_storage(u.*, ldu, 1, k.* + 1, m, m - k.*);
            u_right := array2d_view(u_right_raw, ldu, row_major=false);
            _ggsvp3_m2r(#char "R", #char "N", m, m - k.*, min(m - k.*, l.*), *a23, lda, tau, *u_right, ldu, work, info, RT, MT);
        }

        for j:n - l.* + 1..n {
            for i:j - n + k.* + l.* + 1..m {
                fset(a, i, j, zero_t);
            }
        }
    }

    fset(work, 1, _ggsvp3_complex_from_real(ncast(RT, lwkopt), RT, MT));
}
