// ? LAPACK *gesvj()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesvj/cgesvj/dgesvj/zgesvj
// *GESVJ computes the singular value decomposition of a real or
// complex matrix using the Jacobi SVD algorithm.

sgesvj :: #bake_arguments _gesvj_real(RT=float32);
dgesvj :: #bake_arguments _gesvj_real(RT=float64);
cgesvj :: #bake_arguments _gesvj_complex(RT=float32, MT=C32);
zgesvj :: #bake_arguments _gesvj_complex(RT=float64, MT=C64);

#scope_file
_gesvj_real :: (
    $$joba: u8,
    $$jobu: u8,
    $$jobv: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    sva: *$SVA/ArrayType1D,
    $$mv: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           SVA.data_type == RT &&
           V.data_type == RT &&
           WORK.data_type == RT;
} {
    routine_name :: #run if RT == float32 { "SGESVJ" } else { "DGESVJ" };
    lquery := lwork == -1;

    info.* = 0;
    if m < 0 {
        info.* = -4;
    } else if n < 0 || n > m {
        info.* = -5;
    } else if lda < max(1, m) {
        info.* = -7;
    } else if ldv < 1 {
        info.* = -11;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if lquery {
        fset(work, 1, one(RT));
        return;
    }

    // First-stage staging note:
    // The Jacobi sweep body is still pending.
}

_gesvj_complex :: (
    $$joba: u8,
    $$jobu: u8,
    $$jobv: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    sva: *$SVA/ArrayType1D,
    $$mv: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    cwork: *$CWORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    $$lrwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           SVA.data_type == RT &&
           V.data_type == MT &&
           CWORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CGESVJ" } else { "ZGESVJ" };
    lquery := lwork == -1 || lrwork == -1;

    info.* = 0;
    if m < 0 {
        info.* = -4;
    } else if n < 0 || n > m {
        info.* = -5;
    } else if lda < max(1, m) {
        info.* = -7;
    } else if ldv < 1 {
        info.* = -11;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if lquery {
        fset(cwork, 1, cpx(one(RT), zero(RT)));
        if lrwork != 0 { fset(rwork, 1, one(RT)); }
        return;
    }

    // First-stage staging note:
    // The complex Jacobi sweep body is still pending.
}
