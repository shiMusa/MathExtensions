// ? LAPACK *laqz2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqz2/dlaqz2/claqz2/zlaqz2
// *LAQZ2 provides the AED / bulge-chasing helper used by QZ.

slaqz2 :: #bake_arguments _laqz2_real(RT=float32);
dlaqz2 :: #bake_arguments _laqz2_real(RT=float64);
claqz2 :: #bake_arguments _laqz2_complex(RT=float32, MT=C32);
zlaqz2 :: #bake_arguments _laqz2_complex(RT=float64, MT=C64);

#scope_file
_laqz2_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} { fset(work, 1, cpx(x, zero(RT))); }

_laqz2_real :: (
    $$ilq: bool, $$ilz: bool, $$k: int, $$istartm: int, $$istopm: int, $$ihi: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    $$nq: int, $$qstart: int, q: *$Q/ArrayType2D, $$ldq: int,
    $$nz: int, $$zstart: int, z: *$Z/ArrayType2D, $$ldz: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT && Q.data_type == RT && Z.data_type == RT;
} {
    _ = ilq, ilz, k, istartm, istopm, ihi, a, lda, b, ldb, nq, qstart, q, ldq, nz, zstart, z, ldz;
    // Staged port note:
    // The real bulge-chasing helper body from *laqz2.f is still pending.
}

_laqz2_complex :: (
    $$ilschur: bool, $$ilq: bool, $$ilz: bool, $$n: int, $$ilo: int, $$ihi: int, $$nw: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    ns: *int, nd: *int, alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    qc: *$QC/ArrayType2D, $$ldqc: int, zc: *$ZC/ArrayType2D, $$ldzc: int,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$rec: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && B.data_type == MT && Q.data_type == MT && Z.data_type == MT &&
           ALPHA.data_type == MT && BETA.data_type == MT && QC.data_type == MT && ZC.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    ns.* = 0;
    nd.* = 0;
    info.* = 0;
    _laqz2_set_work1_complex(work, ncast(RT, max(1, nw)), RT, MT);
    if lwork == -1 || n <= 0 then return;
    _ = ilschur, ilq, ilz, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, alpha, beta, qc, ldqc, zc, ldzc, rwork, rec;
    // Staged port note:
    // The complex AED helper body from *laqz2.f is still pending.
}
