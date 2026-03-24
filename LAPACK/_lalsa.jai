// ? LAPACK *lalsa()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slalsa/dlalsa/clalsa/zlalsa
// *LALSA computes/apply the compact divide-and-conquer SVD factors.

slalsa :: #bake_arguments _lalsa_real(RT=float32);
dlalsa :: #bake_arguments _lalsa_real(RT=float64);
clalsa :: #bake_arguments _lalsa_complex(RT=float32, MT=C32);
zlalsa :: #bake_arguments _lalsa_complex(RT=float64, MT=C64);

#scope_file
_lalsa_real :: (
    $$icompq: int,
    $$smlsiz: int,
    $$n: int,
    $$nrhs: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    bx: *$BX/ArrayType2D,
    $$ldbx: int,
    u: *$U/ArrayType2D,
    $$ldu: int,
    vt: *$VT/ArrayType2D,
    k: *$K/ArrayType1D,
    difl: *$DIFL/ArrayType2D,
    difr: *$DIFR/ArrayType2D,
    z: *$Z/ArrayType2D,
    poles: *$POLES/ArrayType2D,
    givptr: *$GIVPTR/ArrayType1D,
    givcol: *$GIVCOL/ArrayType2D,
    $$ldgcol: int,
    perm: *$PERM/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    c: *$C/ArrayType1D,
    s: *$S/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           B.data_type == RT &&
           BX.data_type == RT &&
           U.data_type == RT &&
           VT.data_type == RT &&
           K.data_type == int &&
           DIFL.data_type == RT &&
           DIFR.data_type == RT &&
           Z.data_type == RT &&
           POLES.data_type == RT &&
           GIVPTR.data_type == int &&
           GIVCOL.data_type == int &&
           PERM.data_type == int &&
           GIVNUM.data_type == RT &&
           C.data_type == RT &&
           S.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SLALSA" } else { "DLALSA" };

    info.* = 0;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if smlsiz < 3 {
        info.* = -2;
    } else if n < smlsiz {
        info.* = -3;
    } else if nrhs < 1 {
        info.* = -4;
    } else if ldb < n {
        info.* = -6;
    } else if ldbx < n {
        info.* = -8;
    } else if ldu < n {
        info.* = -10;
    } else if ldgcol < n {
        info.* = -19;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    // Staged port note:
    // The tree construction and recursive/application body from *lalsa.f
    // is still pending.
    _ = b; _ = bx; _ = u; _ = vt; _ = k; _ = difl; _ = difr; _ = z;
    _ = poles; _ = givptr; _ = givcol; _ = perm; _ = givnum; _ = c; _ = s;
    _ = work; _ = iwork; _ = nrhs;
}

_lalsa_complex :: (
    $$icompq: int,
    $$smlsiz: int,
    $$n: int,
    $$nrhs: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    bx: *$BX/ArrayType2D,
    $$ldbx: int,
    u: *$U/ArrayType2D,
    $$ldu: int,
    vt: *$VT/ArrayType2D,
    k: *$K/ArrayType1D,
    difl: *$DIFL/ArrayType2D,
    difr: *$DIFR/ArrayType2D,
    z: *$Z/ArrayType2D,
    poles: *$POLES/ArrayType2D,
    givptr: *$GIVPTR/ArrayType1D,
    givcol: *$GIVCOL/ArrayType2D,
    $$ldgcol: int,
    perm: *$PERM/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    c: *$C/ArrayType1D,
    s: *$S/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           B.data_type == MT &&
           BX.data_type == MT &&
           U.data_type == RT &&
           VT.data_type == RT &&
           K.data_type == int &&
           DIFL.data_type == RT &&
           DIFR.data_type == RT &&
           Z.data_type == RT &&
           POLES.data_type == RT &&
           GIVPTR.data_type == int &&
           GIVCOL.data_type == int &&
           PERM.data_type == int &&
           GIVNUM.data_type == RT &&
           C.data_type == RT &&
           S.data_type == RT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CLALSA" } else { "ZLALSA" };

    info.* = 0;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if smlsiz < 3 {
        info.* = -2;
    } else if n < smlsiz {
        info.* = -3;
    } else if nrhs < 1 {
        info.* = -4;
    } else if ldb < n {
        info.* = -6;
    } else if ldbx < n {
        info.* = -8;
    } else if ldu < n {
        info.* = -10;
    } else if ldgcol < n {
        info.* = -19;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    // Staged port note:
    // The complex tree construction and application body is still pending.
    _ = b; _ = bx; _ = u; _ = vt; _ = k; _ = difl; _ = difr; _ = z;
    _ = poles; _ = givptr; _ = givcol; _ = perm; _ = givnum; _ = c; _ = s;
    _ = rwork; _ = iwork; _ = nrhs;
}
