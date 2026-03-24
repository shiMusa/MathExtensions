// ? LAPACK *laed7()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaed7/claed7/dlaed7/zlaed7
// *LAED7 computes the updated eigensystem of a diagonal matrix after
// modification by a rank-one symmetric matrix.

slaed7 :: #bake_arguments _laed7_real(RT=float32);
dlaed7 :: #bake_arguments _laed7_real(RT=float64);
claed7 :: #bake_arguments _laed7_complex(RT=float32, MT=C32);
zlaed7 :: #bake_arguments _laed7_complex(RT=float64, MT=C64);

#scope_file
_laed7_real :: (
    $$icompq: int,
    $$n: int,
    $$qsiz: int,
    $$tlvls: int,
    $$curlvl: int,
    $$curpbm: int,
    d: *$D/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    indxq: *$INDXQ/ArrayType1D,
    $$rho: $RT,
    $$cutpnt: int,
    qstore: *$QSTORE/ArrayType1D,
    qptr: *$QPTR/ArrayType1D,
    prmptr: *$PRMPTR/ArrayType1D,
    perm: *$PERM/ArrayType1D,
    givptr: *$GIVPTR/ArrayType1D,
    givcol: *$GIVCOL/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           Q.data_type == RT &&
           INDXQ.data_type == int &&
           QSTORE.data_type == RT &&
           QPTR.data_type == int &&
           PRMPTR.data_type == int &&
           PERM.data_type == int &&
           GIVPTR.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SLAED7" } else { "DLAED7" };

    info.* = 0;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if icompq == 1 && qsiz < n {
        info.* = -3;
    } else if ldq < max(1, n) {
        info.* = -9;
    } else if min(1, n) > cutpnt || n < cutpnt {
        info.* = -12;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    // Staged port note:
    // The merge/update body depends on *laed8, *laed9, and *laeda.
    _ = tlvls; _ = curlvl; _ = curpbm; _ = d; _ = q; _ = indxq; _ = rho;
    _ = qstore; _ = qptr; _ = prmptr; _ = perm; _ = givptr; _ = givcol; _ = givnum; _ = work; _ = iwork;
}

_laed7_complex :: (
    $$n: int,
    $$cutpnt: int,
    $$qsiz: int,
    $$tlvls: int,
    $$curlvl: int,
    $$curpbm: int,
    d: *$D/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    $$rho: $RT,
    indxq: *$INDXQ/ArrayType1D,
    qstore: *$QSTORE/ArrayType1D,
    qptr: *$QPTR/ArrayType1D,
    prmptr: *$PRMPTR/ArrayType1D,
    perm: *$PERM/ArrayType1D,
    givptr: *$GIVPTR/ArrayType1D,
    givcol: *$GIVCOL/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT &&
           Q.data_type == MT &&
           INDXQ.data_type == int &&
           QSTORE.data_type == RT &&
           QPTR.data_type == int &&
           PRMPTR.data_type == int &&
           PERM.data_type == int &&
           GIVPTR.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CLAED7" } else { "ZLAED7" };

    info.* = 0;
    if n < 0 {
        info.* = -1;
    } else if min(1, n) > cutpnt || n < cutpnt {
        info.* = -2;
    } else if qsiz < n {
        info.* = -3;
    } else if ldq < max(1, n) {
        info.* = -9;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    // Staged port note:
    // The complex merge/update body depends on *laed8, *laed9, and *laeda.
    _ = tlvls; _ = curlvl; _ = curpbm; _ = d; _ = q; _ = rho; _ = indxq;
    _ = qstore; _ = qptr; _ = prmptr; _ = perm; _ = givptr; _ = givcol; _ = givnum; _ = work; _ = rwork; _ = iwork;
}
