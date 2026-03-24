// ? LAPACK *laed8()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaed8/claed8/dlaed8/zlaed8
// *LAED8 merges eigenvalues and deflates the secular equation.

slaed8 :: #bake_arguments _laed8_real(RT=float32);
dlaed8 :: #bake_arguments _laed8_real(RT=float64);
claed8 :: #bake_arguments _laed8_complex(RT=float32, MT=C32);
zlaed8 :: #bake_arguments _laed8_complex(RT=float64, MT=C64);

#scope_file
_laed8_real :: (
    $$icompq: int,
    k: *int,
    $$n: int,
    $$qsiz: int,
    d: *$D/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    indxq: *$INDXQ/ArrayType1D,
    rho: *RT,
    $$cutpnt: int,
    z: *$Z/ArrayType1D,
    dlambda: *$DLAMBDA/ArrayType1D,
    q2: *$Q2/ArrayType2D,
    $$ldq2: int,
    w: *$W/ArrayType1D,
    perm: *$PERM/ArrayType1D,
    givptr: *int,
    givcol: *$GIVCOL/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    indxp: *$INDXP/ArrayType1D,
    indx: *$INDX/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           Q.data_type == RT &&
           INDXQ.data_type == int &&
           Z.data_type == RT &&
           DLAMBDA.data_type == RT &&
           Q2.data_type == RT &&
           W.data_type == RT &&
           PERM.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           INDXP.data_type == int &&
           INDX.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SLAED8" } else { "DLAED8" };

    info.* = 0;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if n < 0 {
        info.* = -3;
    } else if icompq == 1 && qsiz < n {
        info.* = -4;
    } else if ldq < max(1, n) {
        info.* = -7;
    } else if cutpnt < min(1, n) || cutpnt > n {
        info.* = -10;
    } else if ldq2 < max(1, n) {
        info.* = -14;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        k.* = 0;
        return;
    }

    // Staged port note:
    // The deflation/merge body depends on *lamrg and the surrounding
    // divide-and-conquer update layer.
    _ = d; _ = q; _ = indxq; _ = rho; _ = z; _ = dlambda; _ = q2; _ = w;
    _ = perm; _ = givptr; _ = givcol; _ = givnum; _ = indxp; _ = indx;
}

_laed8_complex :: (
    k: *int,
    $$n: int,
    $$qsiz: int,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    d: *$D/ArrayType1D,
    rho: *RT,
    $$cutpnt: int,
    z: *$Z/ArrayType1D,
    dlambda: *$DLAMBDA/ArrayType1D,
    q2: *$Q2/ArrayType2D,
    $$ldq2: int,
    w: *$W/ArrayType1D,
    indxp: *$INDXP/ArrayType1D,
    indx: *$INDX/ArrayType1D,
    indxq: *$INDXQ/ArrayType1D,
    perm: *$PERM/ArrayType1D,
    givptr: *int,
    givcol: *$GIVCOL/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           Q.data_type == MT &&
           D.data_type == RT &&
           Z.data_type == RT &&
           DLAMBDA.data_type == RT &&
           Q2.data_type == MT &&
           W.data_type == RT &&
           INDXP.data_type == int &&
           INDX.data_type == int &&
           INDXQ.data_type == int &&
           PERM.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CLAED8" } else { "ZLAED8" };

    info.* = 0;
    if n < 0 {
        info.* = -2;
    } else if qsiz < n {
        info.* = -3;
    } else if ldq < max(1, n) {
        info.* = -5;
    } else if cutpnt < min(1, n) || cutpnt > n {
        info.* = -8;
    } else if ldq2 < max(1, n) {
        info.* = -12;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        k.* = 0;
        return;
    }

    // Staged port note:
    // The complex deflation/merge body depends on *lamrg and the
    // surrounding divide-and-conquer update layer.
    _ = q; _ = d; _ = rho; _ = z; _ = dlambda; _ = q2; _ = w;
    _ = indxp; _ = indx; _ = indxq; _ = perm; _ = givptr; _ = givcol; _ = givnum;
}
