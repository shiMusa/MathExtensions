// ? LAPACK *laqz0()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqz0/dlaqz0/claqz0/zlaqz0
// *LAQZ0 is the multishift QZ driver.

slaqz0 :: #bake_arguments _laqz0_real(RT=float32);
dlaqz0 :: #bake_arguments _laqz0_real(RT=float64);
claqz0 :: #bake_arguments _laqz0_complex(RT=float32, MT=C32);
zlaqz0 :: #bake_arguments _laqz0_complex(RT=float64, MT=C64);

#scope_file
_laqz0_set_work1_real :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && WORK.data_type == RT;
} { fset(work, 1, x); }

_laqz0_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} { fset(work, 1, cpx(x, zero(RT))); }

_laqz0_real :: (
    $$wants: u8, $$wantq: u8, $$wantz: u8, $$n: int, $$ilo: int, $$ihi: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, $$rec: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT &&
           ALPHAR.data_type == RT && ALPHAI.data_type == RT && BETA.data_type == RT &&
           Q.data_type == RT && Z.data_type == RT && WORK.data_type == RT;
} {
    info.* = 0;
    _laqz0_set_work1_real(work, ncast(RT, max(1, n)), RT);
    if lwork == -1 || n <= 0 then return;
    _ = wants, wantq, wantz, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, rec;
    // Staged port note:
    // The multishift QZ driver body from *laqz0.f is still pending.
}

_laqz0_complex :: (
    $$wants: u8, $$wantq: u8, $$wantz: u8, $$n: int, $$ilo: int, $$ihi: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$rec: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && B.data_type == MT && ALPHA.data_type == MT && BETA.data_type == MT &&
           Q.data_type == MT && Z.data_type == MT && WORK.data_type == MT && RWORK.data_type == RT;
} {
    info.* = 0;
    _laqz0_set_work1_complex(work, ncast(RT, max(1, n)), RT, MT);
    if lwork == -1 || n <= 0 then return;
    _ = wants, wantq, wantz, ilo, ihi, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, rwork, rec;
    // Staged port note:
    // The multishift QZ driver body from *laqz0.f is still pending.
}
