// ? LAPACK *gesvd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesvd/cgesvd/dgesvd/zgesvd
// *GESVD computes the singular value decomposition for GE matrices.

sgesvd :: #bake_arguments _gesvd_real(RT=float32);
dgesvd :: #bake_arguments _gesvd_real(RT=float64);
cgesvd :: #bake_arguments _gesvd_complex(RT=float32, MT=C32);
zgesvd :: #bake_arguments _gesvd_complex(RT=float64, MT=C64);

#scope_file
_gesvd_real :: (
    $$jobu: u8,
    $$jobvt: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    s: *$S/ArrayType1D,
    u: *$U/ArrayType2D,
    $$ldu: int,
    vt: *$VT/ArrayType2D,
    $$ldvt: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           S.data_type == RT &&
           U.data_type == RT &&
           VT.data_type == RT &&
           WORK.data_type == RT;
} {
    routine_name :: #run if RT == float32 { "SGESVD" } else { "DGESVD" };
    lquery := lwork == -1;

    info.* = 0;
    if m < 0 {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if lda < max(1, m) {
        info.* = -6;
    } else if ldu < 1 {
        info.* = -9;
    } else if ldvt < 1 {
        info.* = -11;
    } else if lwork < 1 && !lquery {
        info.* = -13;
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
    // This file establishes the family entry points and interface shape
    // for *GESVD. The large path-dependent QR/LQ bidiagonal SVD driver
    // body is still pending.
    if m == 0 || n == 0 {
        return;
    }
}

_gesvd_complex :: (
    $$jobu: u8,
    $$jobvt: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    s: *$S/ArrayType1D,
    u: *$U/ArrayType2D,
    $$ldu: int,
    vt: *$VT/ArrayType2D,
    $$ldvt: int,
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
           S.data_type == RT &&
           U.data_type == MT &&
           VT.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CGESVD" } else { "ZGESVD" };
    lquery := lwork == -1;

    info.* = 0;
    if m < 0 {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if lda < max(1, m) {
        info.* = -6;
    } else if ldu < 1 {
        info.* = -9;
    } else if ldvt < 1 {
        info.* = -11;
    } else if lwork < 1 && !lquery {
        info.* = -13;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if lquery {
        fset(work, 1, cpx(one(RT), zero(RT)));
        return;
    }

    // First-stage staging note:
    // The complex *GESVD execution body is still pending.
    if m == 0 || n == 0 {
        return;
    }
}
