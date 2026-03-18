// ? LAPACK *gelss()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgelss/cgelss/dgelss/zgelss
// *GELSS computes the minimum-norm solution to a least squares problem
// using the singular value decomposition.

sgelss :: #bake_arguments _gelss_real(RT=float32);
dgelss :: #bake_arguments _gelss_real(RT=float64);
cgelss :: #bake_arguments _gelss_complex(RT=float32, MT=C32);
zgelss :: #bake_arguments _gelss_complex(RT=float64, MT=C64);

#scope_file
_gelss_real :: (
    $$m: int,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    s: *$S/ArrayType1D,
    rcond: $RT,
    rank: *int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           B.data_type == RT &&
           S.data_type == RT &&
           WORK.data_type == RT;
} {
    routine_name :: #run if RT == float32 { "SGELSS " } else { "DGELSS " };
    query_value  :: #run if RT == float32 { sroundup_lwork(ncast(RT, 1)) } else { droundup_lwork(ncast(RT, 1)) };

    // First-pass staging note:
    // SGELSS/DGELSS follow the three LAPACK paths:
    // 1. M >= N, optionally QR first when M >= MNTHR
    // 2a. N >> M with enough workspace, factor through L in WORK(IL)
    // 2. Remaining underdetermined path through bidiagonalization of A
    //
    // The final literal port should keep:
    // - query calls for *GEQRF/*ORMQR/*GEBRD/*ORMBR/*ORGBR/*GELQF/*ORMLQ
    // - scaling via *LANGE/*LASCL/*LASET
    // - SBDSQR/DBDSQR solve
    // - SRSCL/DRSCL thresholding and the GEMM/GEMV right-singular-vector multiply
    //
    // Missing helpers are assumed to exist by project convention.
    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if lda < max(1, m) {
        info.* = -5;
    } else if ldb < max(1, max(m, n)) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lwork == -1 {
        fset(work, 1, query_value);
        return;
    }
    if m == 0 || n == 0 {
        rank.* = 0;
        return;
    }
}

_gelss_complex :: (
    $$m: int,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    s: *$S/ArrayType1D,
    rcond: $RT,
    rank: *int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           S.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CGELSS " } else { "ZGELSS " };

    // First-pass staging note:
    // CGELSS/ZGELSS mirror the real path structure but use:
    // - CLANGE/ZLANGE, CLASCL/ZLASCL, CLASET/ZLASET
    // - CUNMQR/ZUNMQR, CUNMBR/ZUNMBR, CUNGBR/ZUNGBR, CUNMLQ/ZUNMLQ
    // - CBDSQR/ZBDSQR and CSRSCL/ZDRSCL
    // - conjugate GEMM/GEMV for the final right-singular-vector multiply
    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if lda < max(1, m) {
        info.* = -5;
    } else if ldb < max(1, max(m, n)) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lwork == -1 {
        fset(work, 1, cpx(ncast(RT, 1), zero(RT)));
        return;
    }
    if m == 0 || n == 0 {
        rank.* = 0;
        return;
    }
}
