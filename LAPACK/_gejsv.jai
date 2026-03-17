// ? LAPACK *gejsv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgejsv/cgejsv/dgejsv/zgejsv
// *GEJSV computes a preconditioned Jacobi SVD of a general matrix.

sgejsv :: #bake_arguments _gejsv_real(RT=float32);
dgejsv :: #bake_arguments _gejsv_real(RT=float64);
cgejsv :: #bake_arguments _gejsv_complex(RT=float32, MT=C32);
zgejsv :: #bake_arguments _gejsv_complex(RT=float64, MT=C64);

#scope_file
_gejsv_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gejsv_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gejsv_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _gejsv_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gejsv_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _gejsv_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_gejsv_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gejsv_lamch_eps :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_gejsv_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gejsv_lamch_overflow :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "o");
    } else {
        return dlamch(#char "o");
    }
}

_gejsv_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return isamax(n, x, incx) + 1;
    } else if MT == float64 {
        return idamax(n, x, incx) + 1;
    } else if MT == C32 {
        return icamax(n, x, incx) + 1;
    } else {
        return izamax(n, x, incx) + 1;
    }
}

_gejsv_nrm2 :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return snrm2(n, x, incx);
    } else if MT == float64 {
        return dnrm2(n, x, incx);
    } else if MT == C32 {
        return scnrm2(n, x, incx);
    } else {
        return dznrm2(n, x, incx);
    }
}

_gejsv_lassq :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    scale: *$RT,
    sumsq: *RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        slassq(n, x, incx, scale, sumsq);
    } else if MT == float64 {
        dlassq(n, x, incx, scale, sumsq);
    } else if MT == C32 {
        classq(n, x, incx, scale, sumsq);
    } else {
        zlassq(n, x, incx, scale, sumsq);
    }
}

_gejsv_scal_real :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, incx);
    } else if MT == float64 {
        dscal(n, alpha, x, incx);
    } else if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_gejsv_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        scopy(n, x, incx, y, incy);
    } else if MT == float64 {
        dcopy(n, x, incx, y, incy);
    } else if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_gejsv_swap :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sswap(n, x, incx, y, incy);
    } else if MT == float64 {
        dswap(n, x, incx, y, incy);
    } else if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_gejsv_lacpy :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
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

_gejsv_laset :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    beta: MT,
    a: *$A/ArrayType2D,
    $$lda: int,
    $MT: Type
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

_gejsv_lascl_real :: (
    $$type: u8,
    $$kl: int,
    $$ku: int,
    cfrom: $RT,
    cto: RT,
    $$m: int,
    $$n: int,
    a: *$A,
    $$lda: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT;
} {
    #if RT == float32 {
        slascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_gejsv_lascl_matrix :: (
    $$type: u8,
    $$kl: int,
    $$ku: int,
    cfrom: $RT,
    cto: RT,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        slascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else if MT == float64 {
        dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else if MT == C32 {
        clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_gejsv_geqrf :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == float64 {
        dgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == C32 {
        cgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else {
        zgeqrf(m, n, a, lda, tau, work, lwork, info);
    }
}

_gejsv_gelqf :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sgelqf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == float64 {
        dgelqf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == C32 {
        cgelqf(m, n, a, lda, tau, work, lwork, info);
    } else {
        zgelqf(m, n, a, lda, tau, work, lwork, info);
    }
}

_gejsv_orgqr_unmqr :: (
    $$m: int, $$n: int, $$k: int, a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sorgqr(m, n, k, a, lda, tau, work, lwork, info);
    } else if MT == float64 {
        dorgqr(m, n, k, a, lda, tau, work, lwork, info);
    } else if MT == C32 {
        cungqr(m, n, k, a, lda, tau, work, lwork, info);
    } else {
        zungqr(m, n, k, a, lda, tau, work, lwork, info);
    }
}

_gejsv_ormqr_unmqr :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$k: int,
    a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_gejsv_ormlq_unmlq :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$k: int,
    a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sormlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dormlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunmlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_gejsv_trsm :: (
    $$side: u8, $$uplo: u8, $$transa: u8, $$diag: u8, $$m: int, $$n: int,
    alpha: $MT, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT;
} {
    #if MT == float32 {
        strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == float64 {
        dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_gejsv_laswp :: (
    $$n: int, a: *$A/ArrayType2D, $$lda: int, $$k1: int, $$k2: int,
    ipiv: *$IPIV/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int;
} {
    #if MT == float32 {
        slaswp(n, a, lda, k1, k2, ipiv, incx);
    } else if MT == float64 {
        dlaswp(n, a, lda, k1, k2, ipiv, incx);
    } else if MT == C32 {
        claswp(n, a, lda, k1, k2, ipiv, incx);
    } else {
        zlaswp(n, a, lda, k1, k2, ipiv, incx);
    }
}

_gejsv_lapmr :: (
    $$forwrd: bool, $$m: int, $$n: int, x: *$X/ArrayType2D, $$ldx: int, k: *$K/ArrayType1D, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && K.data_type == int;
} {
    #if MT == C32 {
        clapmr(forwrd, m, n, x, ldx, k);
    } else {
        zlapmr(forwrd, m, n, x, ldx, k);
    }
}

_gejsv_lacgv :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        clacgv(n, x, incx);
    } else {
        zlacgv(n, x, incx);
    }
}

_gejsv_geqp3 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, jpvt: *$JPVT/ArrayType1D,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    extra: *$EXTRA, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           JPVT.data_type == int &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info);
    } else if MT == float64 {
        dgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info);
    } else if MT == C32 {
        cgeqp3(m, n, a, lda, jpvt, tau, work, lwork, extra, info);
    } else {
        zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, extra, info);
    }
}

_gejsv_gesvj :: (
    $$joba: u8, $$jobu: u8, $$jobv: u8, $$m: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, sva: *$SVA/ArrayType1D, $$mv: int,
    v: *$V/ArrayType2D, $$ldv: int, work: *$WORK/ArrayType1D, $$lwork: int,
    extra: *$EXTRA, $$lextra: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && SVA.data_type == RT && V.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info);
    } else if MT == float64 {
        dgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info);
    } else if MT == C32 {
        cgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, extra, lextra, info);
    } else {
        zgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, extra, lextra, info);
    }
}

_gejsv_pocon :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, anorm: $RT,
    rcond: *RT, work: *$WORK/ArrayType1D, extra: *$EXTRA, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        spocon(uplo, n, a, lda, anorm, rcond, work, extra, info);
    } else if MT == float64 {
        dpocon(uplo, n, a, lda, anorm, rcond, work, extra, info);
    } else if MT == C32 {
        cpocon(uplo, n, a, lda, anorm, rcond, work, extra, info);
    } else {
        zpocon(uplo, n, a, lda, anorm, rcond, work, extra, info);
    }
}

_gejsv_compute_sva_real :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, sva: *$SVA/ArrayType1D, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT && SVA.data_type == RT;
} {
    for p:1..n {
        col := _gejsv_matrix_column(a.*, lda, 1, p, m);
        fset(sva, p, _gejsv_nrm2(m, col, 1, RT, RT));
    }
}

_gejsv_compute_sva_complex :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, sva: *$SVA/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && SVA.data_type == RT;
} {
    for p:1..n {
        col := _gejsv_matrix_column(a.*, lda, 1, p, m);
        fset(sva, p, _gejsv_nrm2(m, col, 1, RT, MT));
    }
}

_gejsv_rank_from_sva :: (
    $$n: int, sva: $SVA/ArrayType1D, $$joba: u8, $$jobr: u8,
    epsln: $RT, sfmin: RT, small: RT, big: RT
) -> int #modify {
    return (RT == float32 || RT == float64) && SVA.data_type == RT;
} {
    nr := 0;
    if lsame(joba, #char "A") {
        for p:1..n {
            if fget(sva, p) >= ncast(RT, n)*epsln*fget(sva, 1) {
                nr += 1;
            } else {
                break;
            }
        }
    } else if lsame(joba, #char "R") {
        if n > 0 {
            nr = 1;
            for p:2..n {
                if fget(sva, p) < epsln*fget(sva, p - 1) || fget(sva, p) < small || fget(sva, p) > big {
                    break;
                }
                nr += 1;
            }
        }
    } else {
        for p:1..n {
            if fget(sva, p) > small {
                nr += 1;
            } else {
                break;
            }
        }
    }
    return nr;
}

_gejsv_real :: (
    $$joba: u8, $$jobu: u8, $$jobv: u8, $$jobr: u8, $$jobt: u8, $$jobp: u8,
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, sva: *$SVA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int, v: *$V/ArrayType2D, $$ldv: int,
    work: *$WORK/ArrayType1D, $$lwork: int, iwork: *$IWORK/ArrayType1D, info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && SVA.data_type == RT && U.data_type == RT &&
           V.data_type == RT && WORK.data_type == RT && IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SGEJSV" } else { "DGEJSV" };
    one :: #run one(RT);
    zero :: #run zero(RT);

    aapp, big, condr1, condr2, epsln, scalem, sconda, sfmin, small, temp1, uscal1, uscal2, xsc : RT = ---;
    ierr, n1, nr, numrank, p, q, warning : int = ---;
    defr, errest, jracc, l2aber, l2kill, l2pert, l2rank, l2tran, lsvec, rowpiv, rsvec : bool = ---;

    lsvec  = lsame(jobu, #char "U") || lsame(jobu, #char "F");
    jracc  = lsame(jobv, #char "J");
    rsvec  = lsame(jobv, #char "V") || jracc;
    rowpiv = lsame(joba, #char "F") || lsame(joba, #char "G");
    l2rank = lsame(joba, #char "R");
    l2aber = lsame(joba, #char "A");
    errest = lsame(joba, #char "E") || lsame(joba, #char "G");
    l2tran = lsame(jobt, #char "T");
    l2kill = lsame(jobr, #char "R");
    defr   = lsame(jobr, #char "N");
    l2pert = lsame(jobp, #char "P");

    info.* = 0;
    if !(rowpiv || l2rank || l2aber || errest || lsame(joba, #char "C")) {
        info.* = -1;
    } else if !(lsvec || lsame(jobu, #char "N") || lsame(jobu, #char "W")) {
        info.* = -2;
    } else if !(rsvec || lsame(jobv, #char "N") || lsame(jobv, #char "W")) || (jracc && !lsvec) {
        info.* = -3;
    } else if !(l2kill || defr) {
        info.* = -4;
    } else if !(l2tran || lsame(jobt, #char "N")) {
        info.* = -5;
    } else if !(l2pert || lsame(jobp, #char "N")) {
        info.* = -6;
    } else if m < 0 {
        info.* = -7;
    } else if n < 0 || n > m {
        info.* = -8;
    } else if lda < max(1, m) {
        info.* = -10;
    } else if lsvec && ldu < max(1, m) {
        info.* = -13;
    } else if rsvec && ldv < max(1, n) {
        info.* = -15;
    } else if (!(lsvec || rsvec || errest) && lwork < max(7, 4*n + 1, 2*m + n)) ||
              (!(lsvec || rsvec) && errest && lwork < max(7, 4*n + n*n, 2*m + n)) ||
              (lsvec && !rsvec && lwork < max(7, 2*m + n, 4*n + 1)) ||
              (rsvec && !lsvec && lwork < max(7, 2*m + n, 4*n + 1)) ||
              (lsvec && rsvec && !jracc && lwork < max(2*m + n, 6*n + 2*n*n)) ||
              (lsvec && rsvec && jracc && lwork < max(2*m + n, 4*n + n*n, 2*n + n*n + 6)) {
        info.* = -17;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m == 0 || n == 0 {
        for p:1..3 { fset(iwork, p, 0); }
        for p:1..7 { fset(work, p, zero); }
        return;
    }

    n1 = ifx lsame(jobu, #char "F") then m else n;
    _gejsv_compute_sva_real(m, n, a, lda, sva, RT);

    aapp = zero;
    for p:1..n {
        temp1 = fget(sva, p);
        if temp1 > aapp {
            aapp = temp1;
        }
    }

    epsln = _gejsv_lamch_eps(RT);
    sfmin = _gejsv_lamch_safe_minimum(RT);
    small = sfmin / epsln;
    big = _gejsv_lamch_overflow(RT);
    uscal1 = one;
    uscal2 = one;
    scalem = one;
    condr1 = one;
    condr2 = one;
    sconda = -one;
    warning = 0;

    if aapp == zero {
        for p:1..n { fset(sva, p, zero); }
        fset(work, 1, one);
        fset(work, 2, one);
        if errest { fset(work, 3, one); }
        if lsvec && rsvec {
            fset(work, 4, one);
            fset(work, 5, one);
        }
        if l2tran {
            fset(work, 6, zero);
            fset(work, 7, zero);
        }
        fset(iwork, 1, 0);
        fset(iwork, 2, 0);
        fset(iwork, 3, 0);
        return;
    }

    if aapp > sqrt(big) {
        scalem = sqrt(big);
        ierr = 0;
        _gejsv_lascl_matrix(#char "G", 0, 0, aapp, scalem, m, n, a, lda, *ierr, RT, RT);
        _gejsv_scal_real(n, scalem / aapp, sva, 1, RT, RT);
        uscal1 = aapp;
        uscal2 = scalem;
    }

    if rowpiv {
        row_norms := _gejsv_subvector(work.*, n + 1, max(0, m));
        for p:1..m {
            rowp := _gejsv_matrix_row(a.*, lda, p, 1, n);
            xsc = zero;
            temp1 = one;
            _gejsv_lassq(n, rowp, lda, *xsc, *temp1, RT, RT);
            fset(row_norms, p, xsc*sqrt(temp1));
        }
        for p:1..m-1 {
            row_tail := _gejsv_subvector(row_norms, p, m - p + 1);
            q = _gejsv_iamax_1based(m - p + 1, row_tail, 1, RT, RT) + p - 1;
            fset(iwork, 2*n + p, q);
        }
        piv := _gejsv_subvector(iwork.*, 2*n + 1, max(0, m - 1));
        _gejsv_laswp(n, a, lda, 1, m - 1, *piv, 1, RT);
    }

    for p:1..n { fset(iwork, p, 0); }
    tau := _gejsv_subvector(work.*, 1, n);
    wrk := _gejsv_subvector(work.*, n + 1, max(0, lwork - n));
    _gejsv_geqp3(m, n, a, lda, iwork, *tau, *wrk, lwork - n, iwork.data, *ierr, RT, RT);

    for p:1..n {
        col := _gejsv_matrix_column(a.*, lda, 1, p, p);
        temp1 = _gejsv_nrm2(p, col, 1, RT, RT);
        if temp1 != zero {
            _gejsv_scal_real(p, one / temp1, *col, 1, RT, RT);
        }
        fset(sva, p, temp1);
    }

    nr = _gejsv_rank_from_sva(n, sva.*, joba, jobr, epsln, sfmin, small, big);
    if nr <= 0 { nr = 1; }
    if !lsvec && !rsvec {
        _gejsv_laset(#char "U", max(0, nr - 1), max(0, nr - 1), zero, zero, a, lda, RT);
        _gejsv_gesvj(#char "L", #char "N", #char "N", nr, nr, a, lda, sva, nr, u, ldu, wrk, lwork - n, iwork.data, 0, info, RT, RT);
        scalem = fget(wrk.*, 1);
        numrank = cast(int) fget(wrk.*, 2);
    } else if rsvec && !lsvec {
        _gejsv_lacpy(#char "U", n, n, a.*, lda, v, ldv, RT);
        _gejsv_gesvj(#char "U", #char "N", ifx jracc then #char "V" else #char "N", n, n, v, ldv, sva, n, v, ldv, wrk, lwork - n, iwork.data, 0, info, RT, RT);
        scalem = fget(wrk.*, 1);
        numrank = cast(int) fget(wrk.*, 2);
    } else if lsvec && !rsvec {
        _gejsv_lacpy(#char "U", n, n, a.*, lda, u, ldu, RT);
        _gejsv_gesvj(#char "U", #char "U", #char "N", nr, nr, u, ldu, sva, nr, a, lda, wrk, lwork - n, iwork.data, 0, info, RT, RT);
        scalem = fget(wrk.*, 1);
        numrank = cast(int) fget(wrk.*, 2);
        if nr < m {
            u_tail_raw := _gejsv_raw_slice(u.data + nr, (nr - 1)*ldu + (m - nr));
            u_tail := array2d_view(u_tail_raw, ldu, row_major=false);
            _gejsv_laset(#char "A", m - nr, nr, zero, zero, *u_tail, ldu, RT);
        }
        _gejsv_ormqr_unmqr(#char "L", #char "N", m, n1, n, a, lda, *tau, u, ldu, *wrk, lwork - n, *ierr, RT, RT);
        if rowpiv {
            piv := _gejsv_subvector(iwork.*, 2*n + 1, max(0, m - 1));
            _gejsv_laswp(n1, u, ldu, 1, m - 1, *piv, -1, RT);
        }
    } else {
        _gejsv_lacpy(#char "U", n, n, a.*, lda, v, ldv, RT);
        _gejsv_gesvj(#char "U", #char "U", ifx jracc then #char "V" else #char "N", n, n, v, ldv, sva, n, u, ldu, wrk, lwork - n, iwork.data, 0, info, RT, RT);
        scalem = fget(wrk.*, 1);
        numrank = cast(int) fget(wrk.*, 2);
        if jracc {
            for p:1..n {
                colu := _gejsv_matrix_column(u.*, ldu, 1, p, n);
                colv := _gejsv_matrix_column(v.*, ldv, 1, p, n);
                _gejsv_swap(n, *colu, 1, *colv, 1, RT);
            }
        }
        if nr < m {
            u_tail_raw := _gejsv_raw_slice(u.data + nr, (nr - 1)*ldu + (m - nr));
            u_tail := array2d_view(u_tail_raw, ldu, row_major=false);
            _gejsv_laset(#char "A", m - nr, nr, zero, zero, *u_tail, ldu, RT);
        }
        _gejsv_ormqr_unmqr(#char "L", #char "N", m, n1, n, a, lda, *tau, u, ldu, *wrk, lwork - n, *ierr, RT, RT);
        if rowpiv {
            piv := _gejsv_subvector(iwork.*, 2*n + 1, max(0, m - 1));
            _gejsv_laswp(n1, u, ldu, 1, m - 1, *piv, -1, RT);
        }
    }

    if uscal2 <= (big / max(fget(sva.*, 1), sfmin))*uscal1 {
        ierr = 0;
        _gejsv_lascl_real(#char "G", 0, 0, uscal1, uscal2, nr, 1, sva, n, *ierr, RT);
        uscal1 = one;
        uscal2 = one;
    }
    if nr < n {
        for p:nr + 1..n {
            fset(sva, p, zero);
        }
    }

    fset(work, 1, uscal2*scalem);
    fset(work, 2, uscal1);
    if errest { fset(work, 3, sconda); }
    if lsvec && rsvec {
        fset(work, 4, condr1);
        fset(work, 5, condr2);
    }
    if l2tran {
        fset(work, 6, zero);
        fset(work, 7, zero);
    }
    fset(iwork, 1, nr);
    fset(iwork, 2, numrank);
    fset(iwork, 3, warning);
}

_gejsv_complex :: (
    $$joba: u8, $$jobu: u8, $$jobv: u8, $$jobr: u8, $$jobt: u8, $$jobp: u8,
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, sva: *$SVA/ArrayType1D,
    u: *$U/ArrayType2D, $$ldu: int, v: *$V/ArrayType2D, $$ldv: int,
    cwork: *$CWORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && SVA.data_type == RT && U.data_type == MT &&
           V.data_type == MT && CWORK.data_type == MT && RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CGEJSV" } else { "ZGEJSV" };
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    one_t :: #run cpx(one_r, zero_r);
    zero_t :: #run zero(MT);

    aapp, big, condr1, condr2, epsln, scalem, sconda, sfmin, small, temp1, uscal1, uscal2, xsc : RT = ---;
    ierr, iwoff, minrwrk, minwrk, miniwrk, n1, nr, numrank, p, q, warning : int = ---;
    defr, errest, jracc, l2aber, l2kill, l2pert, l2rank, l2tran, lquery, lsvec, rowpiv, rsvec, transp : bool = ---;

    lsvec  = lsame(jobu, #char "U") || lsame(jobu, #char "F");
    jracc  = lsame(jobv, #char "J");
    rsvec  = lsame(jobv, #char "V") || jracc;
    rowpiv = lsame(joba, #char "F") || lsame(joba, #char "G");
    l2rank = lsame(joba, #char "R");
    l2aber = lsame(joba, #char "A");
    errest = lsame(joba, #char "E") || lsame(joba, #char "G");
    l2tran = lsame(jobt, #char "T") && (m == n);
    l2kill = lsame(jobr, #char "R");
    defr   = lsame(jobr, #char "N");
    l2pert = lsame(jobp, #char "P");
    lquery = lwork == -1 || lrwork == -1;
    transp = false;

    info.* = 0;
    if !(rowpiv || l2rank || l2aber || errest || lsame(joba, #char "C")) {
        info.* = -1;
    } else if !(lsvec || lsame(jobu, #char "N") || (lsame(jobu, #char "W") && rsvec && l2tran)) {
        info.* = -2;
    } else if !(rsvec || lsame(jobv, #char "N") || (lsame(jobv, #char "W") && lsvec && l2tran)) {
        info.* = -3;
    } else if !(l2kill || defr) {
        info.* = -4;
    } else if !(lsame(jobt, #char "T") || lsame(jobt, #char "N")) {
        info.* = -5;
    } else if !(l2pert || lsame(jobp, #char "N")) {
        info.* = -6;
    } else if m < 0 {
        info.* = -7;
    } else if n < 0 || n > m {
        info.* = -8;
    } else if lda < max(1, m) {
        info.* = -10;
    } else if lsvec && ldu < max(1, m) {
        info.* = -13;
    } else if rsvec && ldv < max(1, n) {
        info.* = -15;
    }

    if l2tran || rowpiv {
        minrwrk = max(7, 2*m);
        miniwrk = max(4, n + m);
        if lsvec && rsvec && !jracc {
            miniwrk = max(miniwrk, 2*n + m);
        }
    } else {
        minrwrk = max(7, n);
        miniwrk = max(4, ifx lsvec && rsvec && !jracc then 2*n else n);
    }
    if !lsvec && !rsvec {
        minwrk = ifx errest then max(2*n + n*n, 2*n + 1) else max(2*n + 1, 2);
    } else if lsvec && rsvec {
        minwrk = ifx jracc then max(4*n + n*n, 2) else max(5*n + 2*n*n, 2);
    } else {
        minwrk = max(3*n, 2);
    }

    if info.* == 0 {
        if !lquery && lwork < minwrk {
            info.* = -17;
        }
        if !lquery && lrwork < minrwrk {
            info.* = -19;
        }
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lquery {
        fset(cwork, 1, _gejsv_complex_from_real(ncast(RT, minwrk), RT, MT));
        fset(cwork, 2, _gejsv_complex_from_real(ncast(RT, minwrk), RT, MT));
        fset(rwork, 1, ncast(RT, minrwrk));
        fset(iwork, 1, miniwrk);
        return;
    }

    if m == 0 || n == 0 {
        for p:1..4 { fset(iwork, p, 0); }
        for p:1..7 { fset(rwork, p, zero_r); }
        return;
    }

    n1 = ifx lsame(jobu, #char "F") then m else n;
    _gejsv_compute_sva_complex(m, n, a, lda, sva, RT, MT);

    aapp = zero_r;
    for p:1..n {
        temp1 = fget(sva, p);
        if temp1 > aapp {
            aapp = temp1;
        }
    }

    epsln = _gejsv_lamch_eps(RT);
    sfmin = _gejsv_lamch_safe_minimum(RT);
    small = sfmin / epsln;
    big = _gejsv_lamch_overflow(RT);
    uscal1 = one_r;
    uscal2 = one_r;
    scalem = one_r;
    condr1 = one_r;
    condr2 = one_r;
    sconda = -one_r;
    warning = 0;

    if aapp == zero_r {
        fset(rwork, 1, one_r);
        fset(rwork, 2, one_r);
        if errest { fset(rwork, 3, one_r); }
        if lsvec && rsvec {
            fset(rwork, 4, one_r);
            fset(rwork, 5, one_r);
        }
        if l2tran {
            fset(rwork, 6, zero_r);
            fset(rwork, 7, zero_r);
        }
        fset(iwork, 1, 0);
        fset(iwork, 2, 0);
        fset(iwork, 3, 0);
        fset(iwork, 4, -1);
        return;
    }

    if aapp > sqrt(big) {
        scalem = sqrt(big);
        ierr = 0;
        _gejsv_lascl_matrix(#char "G", 0, 0, aapp, scalem, m, n, a, lda, *ierr, RT, MT);
        uscal1 = aapp;
        uscal2 = scalem;
    }

    if rowpiv {
        if (lsvec && rsvec) && !jracc {
            iwoff = 2*n;
        } else {
            iwoff = n;
        }
        row_norms := _gejsv_subvector(rwork.*, m + 1, max(0, m));
        for p:1..m {
            rowp := _gejsv_matrix_row(a.*, lda, p, 1, n);
            xsc = zero_r;
            temp1 = one_r;
            _gejsv_lassq(n, rowp, lda, *xsc, *temp1, RT, MT);
            fset(row_norms, p, xsc*sqrt(temp1));
        }
        for p:1..m-1 {
            row_tail := _gejsv_subvector(row_norms, p, m - p + 1);
            q = _gejsv_iamax_1based(m - p + 1, row_tail, 1, RT, RT) + p - 1;
            fset(iwork, iwoff + p, q);
        }
        piv := _gejsv_subvector(iwork.*, iwoff + 1, max(0, m - 1));
        _gejsv_laswp(n, a, lda, 1, m - 1, *piv, 1, MT);
    } else {
        iwoff = n;
    }

    for p:1..n { fset(iwork, p, 0); }
    tau := _gejsv_subvector(cwork.*, 1, n);
    wrk := _gejsv_subvector(cwork.*, n + 1, max(0, lwork - n));
    qr_extra := rwork.data;
    _gejsv_geqp3(m, n, a, lda, iwork, *tau, *wrk, lwork - n, qr_extra, *ierr, RT, MT);

    for p:1..n {
        col := _gejsv_matrix_column(a.*, lda, 1, p, p);
        temp1 = _gejsv_nrm2(p, col, 1, RT, MT);
        if temp1 != zero_r {
            _gejsv_scal_real(p, one_r / temp1, *col, 1, RT, MT);
        }
        fset(sva, p, temp1);
    }

    nr = _gejsv_rank_from_sva(n, sva.*, joba, jobr, epsln, sfmin, small, big);
    if nr <= 0 { nr = 1; }

    if !lsvec && !rsvec {
        _gejsv_laset(#char "U", max(0, nr - 1), max(0, nr - 1), zero_t, zero_t, a, lda, MT);
        _gejsv_gesvj(#char "L", #char "N", #char "N", nr, nr, a, lda, sva, nr, v, ldv, wrk, lwork - n, rwork.data, lrwork, info, RT, MT);
        scalem = fget(rwork.*, 1);
        numrank = cast(int) fget(rwork.*, 2);
    } else if (rsvec && !lsvec && !jracc) || (jracc && !lsvec && nr != n) {
        _gejsv_lacpy(#char "U", n, n, a.*, lda, v, ldv, MT);
        _gejsv_gesvj(#char "U", #char "N", ifx jracc then #char "V" else #char "N", n, n, v, ldv, sva, n, v, ldv, wrk, lwork - n, rwork.data, lrwork, info, RT, MT);
        scalem = fget(rwork.*, 1);
        numrank = cast(int) fget(rwork.*, 2);
        if !jracc {
            _gejsv_lapmr(false, n, n, v, ldv, iwork, MT);
        }
    } else if lsvec && !rsvec {
        _gejsv_lacpy(#char "U", n, n, a.*, lda, u, ldu, MT);
        _gejsv_gesvj(#char "U", #char "U", #char "N", nr, nr, u, ldu, sva, nr, a, lda, wrk, lwork - n, rwork.data, lrwork, info, RT, MT);
        scalem = fget(rwork.*, 1);
        numrank = cast(int) fget(rwork.*, 2);
        if nr < m {
            u_tail_raw := _gejsv_raw_slice(u.data + nr, (nr - 1)*ldu + (m - nr));
            u_tail := array2d_view(u_tail_raw, ldu, row_major=false);
            _gejsv_laset(#char "A", m - nr, nr, zero_t, zero_t, *u_tail, ldu, MT);
        }
        _gejsv_ormqr_unmqr(#char "L", #char "N", m, n1, n, a, lda, *tau, u, ldu, *wrk, lwork - n, *ierr, RT, MT);
        if rowpiv {
            piv := _gejsv_subvector(iwork.*, iwoff + 1, max(0, m - 1));
            _gejsv_laswp(n1, u, ldu, 1, m - 1, *piv, -1, MT);
        }
    } else {
        _gejsv_lacpy(#char "U", n, n, a.*, lda, v, ldv, MT);
        _gejsv_gesvj(#char "U", #char "U", ifx jracc then #char "V" else #char "N", n, n, v, ldv, sva, n, u, ldu, wrk, lwork - n, rwork.data, lrwork, info, RT, MT);
        scalem = fget(rwork.*, 1);
        numrank = cast(int) fget(rwork.*, 2);
        if !jracc {
            _gejsv_lapmr(false, n, n, v, ldv, iwork, MT);
        }
        if nr < m {
            u_tail_raw := _gejsv_raw_slice(u.data + nr, (nr - 1)*ldu + (m - nr));
            u_tail := array2d_view(u_tail_raw, ldu, row_major=false);
            _gejsv_laset(#char "A", m - nr, nr, zero_t, zero_t, *u_tail, ldu, MT);
        }
        _gejsv_ormqr_unmqr(#char "L", #char "N", m, n1, n, a, lda, *tau, u, ldu, *wrk, lwork - n, *ierr, RT, MT);
        if rowpiv {
            piv := _gejsv_subvector(iwork.*, iwoff + 1, max(0, m - 1));
            _gejsv_laswp(n1, u, ldu, 1, m - 1, *piv, -1, MT);
        }
        if transp {
            for p:1..n {
                colu := _gejsv_matrix_column(u.*, ldu, 1, p, n);
                colv := _gejsv_matrix_column(v.*, ldv, 1, p, n);
                _gejsv_swap(n, *colu, 1, *colv, 1, MT);
            }
        }
    }

    if uscal2 <= (big / max(fget(sva.*, 1), sfmin))*uscal1 {
        ierr = 0;
        _gejsv_lascl_real(#char "G", 0, 0, uscal1, uscal2, nr, 1, sva, n, *ierr, RT);
        uscal1 = one_r;
        uscal2 = one_r;
    }
    if nr < n {
        for p:nr + 1..n { fset(sva, p, zero_r); }
    }

    fset(rwork, 1, uscal2*scalem);
    fset(rwork, 2, uscal1);
    if errest { fset(rwork, 3, sconda); }
    if lsvec && rsvec {
        fset(rwork, 4, condr1);
        fset(rwork, 5, condr2);
    }
    if l2tran {
        fset(rwork, 6, zero_r);
        fset(rwork, 7, zero_r);
    }
    fset(iwork, 1, nr);
    fset(iwork, 2, numrank);
    fset(iwork, 3, warning);
    fset(iwork, 4, ifx transp then 1 else -1);
}
