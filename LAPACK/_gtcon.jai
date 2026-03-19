// ? LAPACK *gtcon()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgtcon/cgtcon/dgtcon/zgtcon
// *GTCON estimates the reciprocal condition number of a tridiagonal
// matrix from its LU factorization computed by *GTTRF.

sgtcon :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    anorm: float32,
    rcond: *float32,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gtcon_real(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info, RT=float32);
}

dgtcon :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    anorm: float64,
    rcond: *float64,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gtcon_real(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info, RT=float64);
}

cgtcon :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    anorm: float32,
    rcond: *float32,
    work: *$WORK/ArrayType1D,
    info: *int
) {
    _gtcon_complex(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info, RT=float32, MT=C32);
}

zgtcon :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    anorm: float64,
    rcond: *float64,
    work: *$WORK/ArrayType1D,
    info: *int
) {
    _gtcon_complex(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info, RT=float64, MT=C64);
}

#scope_file
_gtcon_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gtcon_lacn2_real :: (
    $$n: int,
    v: *$V/ArrayType1D,
    x: *$X/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    est: *$RT,
    kase: *int,
    isave: *$ISAVE/ArrayType1D,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           V.data_type == RT &&
           X.data_type == RT &&
           IWORK.data_type == int &&
           ISAVE.data_type == int;
} {
    #if RT == float32 {
        slacn2(n, v, x, iwork, est, kase, isave);
    } else {
        dlacn2(n, v, x, iwork, est, kase, isave);
    }
}

_gtcon_lacn2_complex :: (
    $$n: int,
    v: *$V/ArrayType1D,
    x: *$X/ArrayType1D,
    est: *$RT,
    kase: *int,
    isave: *$ISAVE/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           V.data_type == MT &&
           X.data_type == MT &&
           ISAVE.data_type == int;
} {
    #if MT == C32 {
        clacn2(n, v, x, est, kase, isave);
    } else {
        zlacn2(n, v, x, est, kase, isave);
    }
}

_gtcon_gttrs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    #if MT == float32 {
        sgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    } else if MT == float64 {
        dgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    } else if MT == C32 {
        cgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    } else {
        zgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    }
}

_gtcon_real :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    anorm: $RT,
    rcond: *RT,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           DL.data_type == RT &&
           D.data_type == RT &&
           DU.data_type == RT &&
           DU2.data_type == RT &&
           IPIV.data_type == int &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    routine_name :: #run if RT == float32 { "SGTCON" } else { "DGTCON" };

    onenrm : bool = ---;
    ainvnm : RT = ---;
    i, kase, kase1 : int = ---;
    isave : [3] int;

    info.* = 0;
    onenrm = norm == #char "1" || lsame(norm, #char "O");
    if !onenrm && !lsame(norm, #char "I") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if anorm < zero_r {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    rcond.* = zero_r;
    if n == 0 {
        rcond.* = one_r;
        return;
    } else if anorm == zero_r {
        return;
    }

    for i:1..n {
        if fget(d, i) == zero_r {
            return;
        }
    }

    ainvnm = zero_r;
    if onenrm {
        kase1 = 1;
    } else {
        kase1 = 2;
    }
    kase = 0;

    x_view := _gtcon_subvector(work.*, 1, n);
    v_view := _gtcon_subvector(work.*, n + 1, n);
    work_mat := array2d_view(x_view, max(1, n), row_major=false);

    while true {
        _gtcon_lacn2_real(n, *v_view, *x_view, iwork, *ainvnm, *kase, *isave, RT);
        if kase == 0 break;

        if kase == kase1 {
            _gtcon_gttrs(#char "N", n, 1, dl, d, du, du2, ipiv, *work_mat, n, info, RT, RT);
        } else {
            _gtcon_gttrs(#char "T", n, 1, dl, d, du, du2, ipiv, *work_mat, n, info, RT, RT);
        }
    }

    if ainvnm != zero_r {
        rcond.* = (one_r / ainvnm) / anorm;
    }
}

_gtcon_complex :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    anorm: $RT,
    rcond: *RT,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           WORK.data_type == MT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    zero_t :: #run cpx(zero(RT), zero(RT));
    routine_name :: #run if MT == C32 { "CGTCON" } else { "ZGTCON" };

    onenrm : bool = ---;
    ainvnm : RT = ---;
    i, kase, kase1 : int = ---;
    isave : [3] int;

    info.* = 0;
    onenrm = norm == #char "1" || lsame(norm, #char "O");
    if !onenrm && !lsame(norm, #char "I") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if anorm < zero_r {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    rcond.* = zero_r;
    if n == 0 {
        rcond.* = one_r;
        return;
    } else if anorm == zero_r {
        return;
    }

    for i:1..n {
        if fget(d, i) == zero_t {
            return;
        }
    }

    ainvnm = zero_r;
    if onenrm {
        kase1 = 1;
    } else {
        kase1 = 2;
    }
    kase = 0;

    x_view := _gtcon_subvector(work.*, 1, n);
    v_view := _gtcon_subvector(work.*, n + 1, n);
    work_mat := array2d_view(x_view, max(1, n), row_major=false);

    while true {
        _gtcon_lacn2_complex(n, *v_view, *x_view, *ainvnm, *kase, *isave, RT, MT);
        if kase == 0 break;

        if kase == kase1 {
            _gtcon_gttrs(#char "N", n, 1, dl, d, du, du2, ipiv, *work_mat, n, info, RT, MT);
        } else {
            _gtcon_gttrs(#char "C", n, 1, dl, d, du, du2, ipiv, *work_mat, n, info, RT, MT);
        }
    }

    if ainvnm != zero_r {
        rcond.* = (one_r / ainvnm) / anorm;
    }
}
