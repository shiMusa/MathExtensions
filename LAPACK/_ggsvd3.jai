// ? LAPACK *ggsvd3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggsvd3/cggsvd3/dggsvd3/zggsvd3
// *GGSVD3 computes the generalized singular value decomposition.

sggsvd3 :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$n: int, $$p: int,
    k: *int, l: *int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, info: *int
) {
    _ggsvd3_real(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info, RT=float32);
}

dggsvd3 :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$n: int, $$p: int,
    k: *int, l: *int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, info: *int
) {
    _ggsvd3_real(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info, RT=float64);
}

cggsvd3 :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$n: int, $$p: int,
    k: *int, l: *int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, info: *int
) {
    _ggsvd3_complex(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, rwork, iwork, info, RT=float32, MT=C32);
}

zggsvd3 :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$n: int, $$p: int,
    k: *int, l: *int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, info: *int
) {
    _ggsvd3_complex(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, rwork, iwork, info, RT=float64, MT=C64);
}

#scope_file
_ggsvd3_roundup_lwork :: inline (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_ggsvd3_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && (MT == C32 || MT == C64);
} {
    return cpx(x, zero(RT));
}

_ggsvd3_query_to_int :: inline (x: $T, $RT: Type, $T: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (T == RT || T == complex_with_inner_type(RT));
} {
    #if T == float32 || T == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_ggsvd3_lange :: (
    $$norm: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int,
    work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && WORK.data_type == RT;
} {
    #if MT == float32 {
        return slange(norm, m, n, a, lda, work);
    } else if MT == float64 {
        return dlange(norm, m, n, a, lda, work);
    } else if MT == C32 {
        return clange(norm, m, n, a, lda, work);
    } else {
        return zlange(norm, m, n, a, lda, work);
    }
}

_ggsvd3_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach);
    } else {
        return dlamch(cmach);
    }
}

_ggsvd3_copy_real :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 {
        scopy(n, x, incx, y, incy);
    } else {
        dcopy(n, x, incx, y, incy);
    }
}

_ggsvd3_ggsvp3_real :: (
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
    info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT &&
           U.data_type == RT && V.data_type == RT && Q.data_type == RT &&
           IWORK.data_type == int && TAU.data_type == RT && WORK.data_type == RT;
} {
    #if RT == float32 {
        sggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, lwork, info);
    } else {
        dggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, lwork, info);
    }
}

_ggsvd3_ggsvp3_complex :: (
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
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT &&
           U.data_type == MT && V.data_type == MT && Q.data_type == MT &&
           IWORK.data_type == int && RWORK.data_type == RT &&
           TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info);
    } else {
        zggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info);
    }
}

_ggsvd3_tgsja :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$p: int, $$n: int, $$k: int, $$l: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $$tola: $RT, $$tolb: RT,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    work: *$WORK/ArrayType1D, ncycle: *int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && B.data_type == MT &&
           ALPHA.data_type == RT && BETA.data_type == RT &&
           U.data_type == MT && V.data_type == MT && Q.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        stgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
    } else if MT == float64 {
        dtgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
    } else if MT == C32 {
        ctgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
    } else {
        ztgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
    }
}

_ggsvd3_real :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$n: int, $$p: int,
    k: *int, l: *int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT &&
           ALPHA.data_type == RT && BETA.data_type == RT &&
           U.data_type == RT && V.data_type == RT && Q.data_type == RT &&
           WORK.data_type == RT && IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 then "SGGSVD3" else "DGGSVD3";

    wantu, wantv, wantq, lquery : bool = ---;
    anorm, bnorm, smax, temp, tola, tolb, ulp, unfl : RT = ---;
    i, ibnd, isub, j, lwkopt, ncycle : int = ---;

    info.* = 0;
    wantu = lsame(jobu, #char "U");
    wantv = lsame(jobv, #char "V");
    wantq = lsame(jobq, #char "Q");
    lquery = lwork == -1;
    lwkopt = 1;
    tola = zero(RT);
    tolb = zero(RT);

    if !(wantu || lsame(jobu, #char "N")) {
        info.* = -1;
    } else if !(wantv || lsame(jobv, #char "N")) {
        info.* = -2;
    } else if !(wantq || lsame(jobq, #char "N")) {
        info.* = -3;
    } else if m < 0 {
        info.* = -4;
    } else if n < 0 {
        info.* = -5;
    } else if p < 0 {
        info.* = -6;
    } else if lda < max(1, m) {
        info.* = -10;
    } else if ldb < max(1, p) {
        info.* = -12;
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
        _ggsvd3_ggsvp3_real(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, work, work, -1, info, RT);
        lwkopt = n + _ggsvd3_query_to_int(fget(work.*, 1), RT, RT);
        lwkopt = max(2*n, lwkopt);
        lwkopt = max(1, lwkopt);
        fset(work, 1, _ggsvd3_roundup_lwork(ncast(RT, lwkopt), RT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lquery return;

    anorm = _ggsvd3_lange(#char "1", m, n, a.*, lda, work, RT, RT);
    bnorm = _ggsvd3_lange(#char "1", p, n, b.*, ldb, work, RT, RT);

    ulp = _ggsvd3_lamch("Precision"s, RT);
    unfl = _ggsvd3_lamch("Safe Minimum"s, RT);
    tola = ncast(RT, max(m, n))*max(anorm, unfl)*ulp;
    tolb = ncast(RT, max(p, n))*max(bnorm, unfl)*ulp;

    work_tail: [] RT;
    work_tail.data = work.data + n;
    work_tail.count = max(0, lwork - n);
    _ggsvd3_ggsvp3_real(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, work, *work_tail, lwork - n, info, RT);

    _ggsvd3_tgsja(jobu, jobv, jobq, m, p, n, k.*, l.*, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, *ncycle, info, RT, RT);

    _ggsvd3_copy_real(n, alpha.*, 1, work, 1, RT);
    ibnd = min(l.*, m - k.*);
    for i:1..ibnd {
        isub = i;
        smax = fget(work.*, k.* + i);
        for j:i + 1..ibnd {
            temp = fget(work.*, k.* + j);
            if temp > smax {
                isub = j;
                smax = temp;
            }
        }

        if isub != i {
            fset(work, k.* + isub, fget(work.*, k.* + i));
            fset(work, k.* + i, smax);
            fset(iwork, k.* + i, k.* + isub);
        } else {
            fset(iwork, k.* + i, k.* + i);
        }
    }

    fset(work, 1, _ggsvd3_roundup_lwork(ncast(RT, lwkopt), RT));
}

_ggsvd3_complex :: (
    $$jobu: u8, $$jobv: u8, $$jobq: u8,
    $$m: int, $$n: int, $$p: int,
    k: *int, l: *int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int,
    v: *$V/ArrayType2D, $$ldv: int,
    q: *$Q/ArrayType2D, $$ldq: int,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT &&
           ALPHA.data_type == RT && BETA.data_type == RT &&
           U.data_type == MT && V.data_type == MT && Q.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 then "CGGSVD3" else "ZGGSVD3";

    wantu, wantv, wantq, lquery : bool = ---;
    anorm, bnorm, smax, temp, tola, tolb, ulp, unfl : RT = ---;
    i, ibnd, isub, j, lwkopt, ncycle : int = ---;

    info.* = 0;
    wantu = lsame(jobu, #char "U");
    wantv = lsame(jobv, #char "V");
    wantq = lsame(jobq, #char "Q");
    lquery = lwork == -1;
    lwkopt = 1;
    tola = zero(RT);
    tolb = zero(RT);

    if !(wantu || lsame(jobu, #char "N")) {
        info.* = -1;
    } else if !(wantv || lsame(jobv, #char "N")) {
        info.* = -2;
    } else if !(wantq || lsame(jobq, #char "N")) {
        info.* = -3;
    } else if m < 0 {
        info.* = -4;
    } else if n < 0 {
        info.* = -5;
    } else if p < 0 {
        info.* = -6;
    } else if lda < max(1, m) {
        info.* = -10;
    } else if ldb < max(1, p) {
        info.* = -12;
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
        _ggsvd3_ggsvp3_complex(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, work, work, -1, info, RT, MT);
        lwkopt = n + _ggsvd3_query_to_int(fget(work.*, 1), RT, MT);
        lwkopt = max(2*n, lwkopt);
        lwkopt = max(1, lwkopt);
        fset(work, 1, _ggsvd3_complex_from_real(ncast(RT, lwkopt), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lquery return;

    anorm = _ggsvd3_lange(#char "1", m, n, a.*, lda, rwork, RT, MT);
    bnorm = _ggsvd3_lange(#char "1", p, n, b.*, ldb, rwork, RT, MT);

    ulp = _ggsvd3_lamch("Precision"s, RT);
    unfl = _ggsvd3_lamch("Safe Minimum"s, RT);
    tola = ncast(RT, max(m, n))*max(anorm, unfl)*ulp;
    tolb = ncast(RT, max(p, n))*max(bnorm, unfl)*ulp;

    work_tail: [] MT;
    work_tail.data = work.data + n;
    work_tail.count = max(0, lwork - n);
    _ggsvd3_ggsvp3_complex(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, work, *work_tail, lwork - n, info, RT, MT);

    _ggsvd3_tgsja(jobu, jobv, jobq, m, p, n, k.*, l.*, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, *ncycle, info, RT, MT);

    _ggsvd3_copy_real(n, alpha.*, 1, rwork, 1, RT);
    ibnd = min(l.*, m - k.*);
    for i:1..ibnd {
        isub = i;
        smax = fget(rwork.*, k.* + i);
        for j:i + 1..ibnd {
            temp = fget(rwork.*, k.* + j);
            if temp > smax {
                isub = j;
                smax = temp;
            }
        }

        if isub != i {
            fset(rwork, k.* + isub, fget(rwork.*, k.* + i));
            fset(rwork, k.* + i, smax);
            fset(iwork, k.* + i, k.* + isub);
        } else {
            fset(iwork, k.* + i, k.* + i);
        }
    }

    fset(work, 1, _ggsvd3_complex_from_real(ncast(RT, lwkopt), RT, MT));
}
