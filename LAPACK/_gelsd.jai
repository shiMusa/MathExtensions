// ? LAPACK *gelsd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgelsd/cgelsd/dgelsd/zgelsd
// *GELSD computes the minimum-norm solution to a least squares problem
// using the SVD and divide-and-conquer on the bidiagonal problem.

sgelsd :: #bake_arguments _gelsd_real(RT=float32);
dgelsd :: #bake_arguments _gelsd_real(RT=float64);
cgelsd :: #bake_arguments _gelsd_complex(RT=float32, MT=C32);
zgelsd :: #bake_arguments _gelsd_complex(RT=float64, MT=C64);

#scope_file
_gelsd_real :: (
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
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           B.data_type == RT &&
           S.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SGELSD " } else { "DGELSD " };
    query_value  :: #run if RT == float32 { sroundup_lwork(ncast(RT, 1)) } else { droundup_lwork(ncast(RT, 1)) };

    // First-pass staging note:
    // SGELSD/DGELSD keep the same three geometric paths as *GELSS, but
    // replace the bidiagonal QR iteration with *LALSD and add LIWORK.
    // The final literal port should preserve:
    // - SMLSIZ/MNTHR/NLVL workspace logic
    // - WLALSD-driven workspace thresholds
    // - *GEBRD + *ORMBR + *LALSD + *ORMBR backtransform
    // - the efficient Path 2a scratch-L factorization in WORK(IL)
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
        fset(iwork, 1, 1);
        return;
    }
    if m == 0 || n == 0 {
        rank.* = 0;
        return;
    }
}

_gelsd_complex :: (
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
    iwork: *$IWORK/ArrayType1D,
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
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CGELSD " } else { "ZGELSD " };

    // First-pass staging note:
    // CGELSD/ZGELSD follow the same path split as the real driver, but
    // use CLALSD/ZLALSD and separate LRWORK/LIWORK formulas.
    // The final literal port should preserve the CGEBRD/ZGEBRD +
    // CUNMBR/ZUNMBR transforms and the explicit RWORK(1)/IWORK(1) query results.
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
        fset(rwork, 1, one(RT));
        fset(iwork, 1, 1);
        return;
    }
    if m == 0 || n == 0 {
        rank.* = 0;
        return;
    }
}
