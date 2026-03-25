// ? LAPACK *laqp3rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqp3rk/dlaqp3rk/claqp3rk/zlaqp3rk
// *LAQP3RK computes a truncated QR factorization with column pivoting
// of a matrix block using Level 3 BLAS.

slaqp3rk :: #bake_arguments _laqp3rk_real(RT=float32);
dlaqp3rk :: #bake_arguments _laqp3rk_real(RT=float64);
claqp3rk :: #bake_arguments _laqp3rk_complex(RT=float32, MT=C32);
zlaqp3rk :: #bake_arguments _laqp3rk_complex(RT=float64, MT=C64);

#scope_file
_laqp3rk_real :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$nb: int, $$abstol: $RT, $$reltol: RT,
    $$kp1: int, $$maxc2nrm: RT, a: *$A/ArrayType2D, $$lda: int, done: *bool, kb: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, auxv: *$AUXV/ArrayType1D,
    f: *$F/ArrayType2D, $$ldf: int, iwork: *$IWORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPIV.data_type == int &&
           TAU.data_type == RT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           AUXV.data_type == RT &&
           F.data_type == RT &&
           IWORK.data_type == int;
} {
    zero_r :: #run zero(RT);

    info.* = 0;
    done.* = false;
    kb.* = 0;
    maxc2nrmk.* = maxc2nrm;
    relmaxc2nrmk.* = if maxc2nrm == zero_r then zero_r else one(RT);

    // Staged port note:
    // The full block rank-truncated QR kernel from *laqp3rk.f is still
    // pending.
    _ = m; _ = n; _ = nrhs; _ = ioffset; _ = nb; _ = abstol; _ = reltol; _ = kp1;
    _ = a; _ = lda; _ = jpiv; _ = tau; _ = vn1; _ = vn2; _ = auxv; _ = f; _ = ldf; _ = iwork;
}

_laqp3rk_complex :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$nb: int, $$abstol: $RT, $$reltol: RT,
    $$kp1: int, $$maxc2nrm: RT, a: *$A/ArrayType2D, $$lda: int, done: *bool, kb: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, auxv: *$AUXV/ArrayType1D,
    f: *$F/ArrayType2D, $$ldf: int, iwork: *$IWORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           JPIV.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           AUXV.data_type == MT &&
           F.data_type == MT &&
           IWORK.data_type == int;
} {
    zero_r :: #run zero(RT);

    info.* = 0;
    done.* = false;
    kb.* = 0;
    maxc2nrmk.* = maxc2nrm;
    relmaxc2nrmk.* = if maxc2nrm == zero_r then zero_r else one(RT);

    // Staged port note:
    // The full complex block rank-truncated QR kernel is still pending.
    _ = m; _ = n; _ = nrhs; _ = ioffset; _ = nb; _ = abstol; _ = reltol; _ = kp1;
    _ = a; _ = lda; _ = jpiv; _ = tau; _ = vn1; _ = vn2; _ = auxv; _ = f; _ = ldf; _ = iwork;
}
