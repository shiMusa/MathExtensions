// ? LAPACK *lals0()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slals0/dlals0/clals0/zlals0
// *LALS0 applies back multiplying factors in the divide-and-conquer SVD
// least-squares path.

slals0 :: #bake_arguments _lals0_real(RT=float32);
dlals0 :: #bake_arguments _lals0_real(RT=float64);
clals0 :: #bake_arguments _lals0_complex(RT=float32, MT=C32);
zlals0 :: #bake_arguments _lals0_complex(RT=float64, MT=C64);

#scope_file
_lals0_real :: (
    $$icompq: int,
    $$nl: int,
    $$nr: int,
    $$sqre: int,
    $$nrhs: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    bx: *$BX/ArrayType2D,
    $$ldbx: int,
    perm: *$PERM/ArrayType1D,
    $$givptr: int,
    givcol: *$GIVCOL/ArrayType2D,
    $$ldgcol: int,
    givnum: *$GIVNUM/ArrayType2D,
    $$ldgnum: int,
    poles: *$POLES/ArrayType2D,
    difl: *$DIFL/ArrayType1D,
    difr: *$DIFR/ArrayType2D,
    z: *$Z/ArrayType1D,
    $$k: int,
    $$c: $RT,
    $$s: RT,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           B.data_type == RT &&
           BX.data_type == RT &&
           PERM.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           POLES.data_type == RT &&
           DIFL.data_type == RT &&
           DIFR.data_type == RT &&
           Z.data_type == RT &&
           WORK.data_type == RT;
} {
    routine_name :: #run if RT == float32 { "SLALS0" } else { "DLALS0" };
    n : int = nl + nr + 1;

    info.* = 0;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if nl < 1 {
        info.* = -2;
    } else if nr < 1 {
        info.* = -3;
    } else if sqre < 0 || sqre > 1 {
        info.* = -4;
    } else if nrhs < 1 {
        info.* = -5;
    } else if ldb < n {
        info.* = -7;
    } else if ldbx < n {
        info.* = -9;
    } else if givptr < 0 {
        info.* = -11;
    } else if ldgcol < n {
        info.* = -13;
    } else if ldgnum < n {
        info.* = -15;
    } else if k < 1 {
        info.* = -20;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    // Staged port note:
    // The divide-and-conquer back-transform body from *lals0.f
    // is still pending.
    _ = b; _ = bx; _ = perm; _ = givcol; _ = givnum; _ = poles;
    _ = difl; _ = difr; _ = z; _ = c; _ = s; _ = work; _ = sqre; _ = nrhs;
}

_lals0_complex :: (
    $$icompq: int,
    $$nl: int,
    $$nr: int,
    $$sqre: int,
    $$nrhs: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    bx: *$BX/ArrayType2D,
    $$ldbx: int,
    perm: *$PERM/ArrayType1D,
    $$givptr: int,
    givcol: *$GIVCOL/ArrayType2D,
    $$ldgcol: int,
    givnum: *$GIVNUM/ArrayType2D,
    $$ldgnum: int,
    poles: *$POLES/ArrayType2D,
    difl: *$DIFL/ArrayType1D,
    difr: *$DIFR/ArrayType2D,
    z: *$Z/ArrayType1D,
    $$k: int,
    $$c: $RT,
    $$s: RT,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           B.data_type == MT &&
           BX.data_type == MT &&
           PERM.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           POLES.data_type == RT &&
           DIFL.data_type == RT &&
           DIFR.data_type == RT &&
           Z.data_type == RT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CLALS0" } else { "ZLALS0" };
    n : int = nl + nr + 1;

    info.* = 0;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if nl < 1 {
        info.* = -2;
    } else if nr < 1 {
        info.* = -3;
    } else if sqre < 0 || sqre > 1 {
        info.* = -4;
    } else if nrhs < 1 {
        info.* = -5;
    } else if ldb < n {
        info.* = -7;
    } else if ldbx < n {
        info.* = -9;
    } else if givptr < 0 {
        info.* = -11;
    } else if ldgcol < n {
        info.* = -13;
    } else if ldgnum < n {
        info.* = -15;
    } else if k < 1 {
        info.* = -20;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    // Staged port note:
    // The complex divide-and-conquer back-transform body is still pending.
    _ = b; _ = bx; _ = perm; _ = givcol; _ = givnum; _ = poles;
    _ = difl; _ = difr; _ = z; _ = c; _ = s; _ = rwork; _ = sqre; _ = nrhs;
}
