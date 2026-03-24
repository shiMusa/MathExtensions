// ? LAPACK *laed0()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaed0/claed0/dlaed0/zlaed0
// *LAED0 computes all eigenvalues and corresponding eigenvectors of an
// unreduced symmetric tridiagonal matrix using the divide and conquer
// method.

slaed0 :: #bake_arguments _laed0_real(RT=float32);
dlaed0 :: #bake_arguments _laed0_real(RT=float64);
claed0 :: #bake_arguments _laed0_complex(RT=float32, MT=C32);
zlaed0 :: #bake_arguments _laed0_complex(RT=float64, MT=C64);

#scope_file
_laed0_real :: (
    $$icompq: int,
    $$qsiz: int,
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    qstore: *$QSTORE/ArrayType2D,
    $$ldqs: int,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           E.data_type == RT &&
           Q.data_type == RT &&
           QSTORE.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SLAED0" } else { "DLAED0" };

    info.* = 0;
    if icompq < 0 || icompq > 2 {
        info.* = -1;
    } else if icompq == 1 && qsiz < max(0, n) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ldq < max(1, n) {
        info.* = -7;
    } else if ldqs < max(1, n) {
        info.* = -9;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    // Staged port note:
    // The divide-and-conquer body from *laed0.f is still pending.
    // This file lands the family interface and source-matching
    // validation while the dependent *laed1/*laed7 layer is filled in.
    _ = d; _ = e; _ = q; _ = qstore; _ = work; _ = iwork;
}

_laed0_complex :: (
    $$qsiz: int,
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    qstore: *$QSTORE/ArrayType2D,
    $$ldqs: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT &&
           E.data_type == RT &&
           Q.data_type == MT &&
           QSTORE.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CLAED0" } else { "ZLAED0" };

    info.* = 0;
    if qsiz < max(0, n) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if ldq < max(1, n) {
        info.* = -6;
    } else if ldqs < max(1, n) {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    // Staged port note:
    // The complex divide-and-conquer body is still pending.
    _ = d; _ = e; _ = q; _ = qstore; _ = rwork; _ = iwork;
}
