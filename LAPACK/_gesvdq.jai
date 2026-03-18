// ? LAPACK *gesvdq()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesvdq/cgesvdq/dgesvdq/zgesvdq
// *GESVDQ computes the singular value decomposition using the
// preconditioned QR method.

sgesvdq :: #bake_arguments _gesvdq_real(RT=float32);
dgesvdq :: #bake_arguments _gesvdq_real(RT=float64);
cgesvdq :: #bake_arguments _gesvdq_complex(RT=float32, MT=C32);
zgesvdq :: #bake_arguments _gesvdq_complex(RT=float64, MT=C64);

#scope_file
_gesvdq_real :: (
    $$joba: u8,
    $$jobp: u8,
    $$jobr: u8,
    $$jobu: u8,
    $$jobv: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    s: *$S/ArrayType1D,
    u: *$U/ArrayType2D,
    $$ldu: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    numrank: *int,
    iwork: *$IWORK/ArrayType1D,
    $$liwork: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    $$lrwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           S.data_type == RT &&
           U.data_type == RT &&
           V.data_type == RT &&
           IWORK.data_type == int &&
           WORK.data_type == RT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if RT == float32 { "SGESVDQ" } else { "DGESVDQ" };
    lquery := liwork == -1 || lwork == -1 || lrwork == -1;

    info.* = 0;
    if m < 0 {
        info.* = -6;
    } else if n < 0 || n > m {
        info.* = -7;
    } else if lda < max(1, m) {
        info.* = -9;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if lquery {
        if liwork != 0 { fset(iwork, 1, 1); }
        fset(work, 1, one(RT));
        if lrwork != 0 { fset(rwork, 1, one(RT)); }
        return;
    }

    // First-stage staging note:
    // The full rank-revealing QR / SVD driver body is still pending.
    numrank.* = 0;
}

_gesvdq_complex :: (
    $$joba: u8,
    $$jobp: u8,
    $$jobr: u8,
    $$jobu: u8,
    $$jobv: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    s: *$S/ArrayType1D,
    u: *$U/ArrayType2D,
    $$ldu: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    numrank: *int,
    iwork: *$IWORK/ArrayType1D,
    $$liwork: int,
    cwork: *$CWORK/ArrayType1D,
    $$lcwork: int,
    rwork: *$RWORK/ArrayType1D,
    $$lrwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           S.data_type == RT &&
           U.data_type == MT &&
           V.data_type == MT &&
           IWORK.data_type == int &&
           CWORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CGESVDQ" } else { "ZGESVDQ" };
    lquery := liwork == -1 || lcwork == -1 || lrwork == -1;

    info.* = 0;
    if m < 0 {
        info.* = -6;
    } else if n < 0 || n > m {
        info.* = -7;
    } else if lda < max(1, m) {
        info.* = -9;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if lquery {
        if liwork != 0 { fset(iwork, 1, 1); }
        fset(cwork, 1, cpx(one(RT), zero(RT)));
        if lrwork != 0 { fset(rwork, 1, one(RT)); }
        return;
    }

    // First-stage staging note:
    // The complex *GESVDQ execution body is still pending.
    numrank.* = 0;
}
