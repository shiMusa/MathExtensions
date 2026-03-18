// ? LAPACK *gesvdx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesvdx/cgesvdx/dgesvdx/zgesvdx
// *GESVDX computes selected singular values and, optionally, vectors.

sgesvdx :: #bake_arguments _gesvdx_real(RT=float32);
dgesvdx :: #bake_arguments _gesvdx_real(RT=float64);
cgesvdx :: #bake_arguments _gesvdx_complex(RT=float32, MT=C32);
zgesvdx :: #bake_arguments _gesvdx_complex(RT=float64, MT=C64);

#scope_file
_gesvdx_real :: (
    $$jobu: u8,
    $$jobvt: u8,
    $$range: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$vl: $RT,
    $$vu: $RT,
    $$il: int,
    $$iu: int,
    ns: *int,
    s: *$S/ArrayType1D,
    u: *$U/ArrayType2D,
    $$ldu: int,
    vt: *$VT/ArrayType2D,
    $$ldvt: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           S.data_type == RT &&
           U.data_type == RT &&
           VT.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SGESVDX" } else { "DGESVDX" };
    lquery := lwork == -1;

    info.* = 0;
    if m < 0 {
        info.* = -4;
    } else if n < 0 {
        info.* = -5;
    } else if lda < max(1, m) {
        info.* = -7;
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
    // The full bidiagonal reduction + *BDSVDX path is still pending.
    ns.* = 0;
}

_gesvdx_complex :: (
    $$jobu: u8,
    $$jobvt: u8,
    $$range: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$vl: $RT,
    $$vu: $RT,
    $$il: int,
    $$iu: int,
    ns: *int,
    s: *$S/ArrayType1D,
    u: *$U/ArrayType2D,
    $$ldu: int,
    vt: *$VT/ArrayType2D,
    $$ldvt: int,
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
           S.data_type == RT &&
           U.data_type == MT &&
           VT.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CGESVDX" } else { "ZGESVDX" };
    lquery := lwork == -1;

    info.* = 0;
    if m < 0 {
        info.* = -4;
    } else if n < 0 {
        info.* = -5;
    } else if lda < max(1, m) {
        info.* = -7;
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
    // The complex *GESVDX execution body is still pending.
    ns.* = 0;
}
