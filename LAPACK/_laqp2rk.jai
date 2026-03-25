// ? LAPACK *laqp2rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqp2rk/dlaqp2rk/claqp2rk/zlaqp2rk
// *LAQP2RK computes a truncated QR factorization with column pivoting
// of a matrix block using Level 2 BLAS.

slaqp2rk :: #bake_arguments _laqp2rk_real(RT=float32);
dlaqp2rk :: #bake_arguments _laqp2rk_real(RT=float64);
claqp2rk :: #bake_arguments _laqp2rk_complex(RT=float32, MT=C32);
zlaqp2rk :: #bake_arguments _laqp2rk_complex(RT=float64, MT=C64);

#scope_file
_laqp2rk_real :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$kmax_: int, $$abstol: $RT, $$reltol: RT,
    $$kp1: int, $$maxc2nrm: RT, a: *$A/ArrayType2D, $$lda: int, k: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, work: *$WORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPIV.data_type == int &&
           TAU.data_type == RT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == RT;
} {
    zero_r :: #run zero(RT);

    info.* = 0;
    k.* = 0;
    maxc2nrmk.* = maxc2nrm;
    relmaxc2nrmk.* = if maxc2nrm == zero_r then zero_r else one(RT);

    // Staged port note:
    // The full rank-truncated Level-2 pivoted QR body from *laqp2rk.f
    // is still pending.
    _ = m; _ = n; _ = nrhs; _ = ioffset; _ = kmax_; _ = abstol; _ = reltol; _ = kp1;
    _ = a; _ = lda; _ = jpiv; _ = tau; _ = vn1; _ = vn2; _ = work;
}

_laqp2rk_complex :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$kmax_: int, $$abstol: $RT, $$reltol: RT,
    $$kp1: int, $$maxc2nrm: RT, a: *$A/ArrayType2D, $$lda: int, k: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           JPIV.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == MT;
} {
    zero_r :: #run zero(RT);

    info.* = 0;
    k.* = 0;
    maxc2nrmk.* = maxc2nrm;
    relmaxc2nrmk.* = if maxc2nrm == zero_r then zero_r else one(RT);

    // Staged port note:
    // The complex rank-truncated Level-2 pivoted QR body is still
    // pending.
    _ = m; _ = n; _ = nrhs; _ = ioffset; _ = kmax_; _ = abstol; _ = reltol; _ = kp1;
    _ = a; _ = lda; _ = jpiv; _ = tau; _ = vn1; _ = vn2; _ = work;
}
