// ? LAPACK *laqz3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqz3/dlaqz3/claqz3/zlaqz3
// *LAQZ3 performs aggressive early deflation / small-bulge QZ work.

slaqz3 :: #bake_arguments _laqz3_real(RT=float32);
dlaqz3 :: #bake_arguments _laqz3_real(RT=float64);
claqz3 :: #bake_arguments _laqz3_complex(RT=float32, MT=C32);
zlaqz3 :: #bake_arguments _laqz3_complex(RT=float64, MT=C64);

#scope_file
_laqz3_set_work1_real :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && WORK.data_type == RT;
} { fset(work, 1, x); }

_laqz3_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} { fset(work, 1, cpx(x, zero(RT))); }

_laqz3_real :: (
    $$ilschur: bool, $$ilq: bool, $$ilz: bool, $$n: int, $$ilo: int, $$ihi: int, $$nw: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    ns: *int, nd: *int, alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    qc: *$QC/ArrayType2D, $$ldqc: int, zc: *$ZC/ArrayType2D, $$ldzc: int,
    work: *$WORK/ArrayType1D, $$lwork: int, $$rec: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT && Q.data_type == RT && Z.data_type == RT &&
           ALPHAR.data_type == RT && ALPHAI.data_type == RT && BETA.data_type == RT &&
           QC.data_type == RT && ZC.data_type == RT && WORK.data_type == RT;
} {
    ns.* = 0;
    nd.* = 0;
    info.* = 0;
    _laqz3_set_work1_real(work, ncast(RT, max(1, nw)), RT);
    if lwork == -1 || n <= 0 then return;
    _ = ilschur, ilq, ilz, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, alphar, alphai, beta, qc, ldqc, zc, ldzc, rec;
    // Staged port note:
    // The real aggressive-early-deflation body from *laqz3.f is still pending.
}

_laqz3_complex :: (
    $$ilschur: bool, $$ilq: bool, $$ilz: bool, $$n: int, $$ilo: int, $$ihi: int, $$nshifts: int, $$nblock_desired: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    qc: *$QC/ArrayType2D, $$ldqc: int, zc: *$ZC/ArrayType2D, $$ldzc: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           ALPHA.data_type == MT && BETA.data_type == MT &&
           A.data_type == MT && B.data_type == MT && Q.data_type == MT && Z.data_type == MT &&
           QC.data_type == MT && ZC.data_type == MT && WORK.data_type == MT;
} {
    info.* = 0;
    _laqz3_set_work1_complex(work, ncast(RT, max(1, n)), RT, MT);
    if lwork == -1 || n <= 0 then return;
    _ = ilschur, ilq, ilz, ilo, ihi, nshifts, nblock_desired, alpha, beta, a, lda, b, ldb, q, ldq, z, ldz, qc, ldqc, zc, ldzc;
    // Staged port note:
    // The complex blocked QZ sweep / AED body from *laqz3.f is still pending.
}
