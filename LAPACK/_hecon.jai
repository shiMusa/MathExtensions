// ? LAPACK *hecon()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment checon/zhecon
// *HECON estimates the reciprocal condition number of a Hermitian
// matrix using the factorization computed by *HETRF.

checon :: (
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: float32,
    rcond: *float32,
    work: *$WORK/ArrayType1D,
    info: *int
) {
    _hecon(uplo, n, a, lda, ipiv, anorm, rcond, work, info, RT=float32, MT=C32);
}

zhecon :: (
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: float64,
    rcond: *float64,
    work: *$WORK/ArrayType1D,
    info: *int
) {
    _hecon(uplo, n, a, lda, ipiv, anorm, rcond, work, info, RT=float64, MT=C64);
}

#scope_file
_hecon_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hecon_lacn2 :: (
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

_hecon_hetrs :: (
    $$uplo: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    #if MT == C32 {
        chetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info);
    } else {
        zhetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info);
    }
}

_hecon :: (
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
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
           A.data_type == MT &&
           IPIV.data_type == int &&
           WORK.data_type == MT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == C32 { "CHECON" } else { "ZHECON" };

    upper : bool = ---;
    ainvnm : RT = ---;
    i, kase : int = ---;
    isave : [3] int;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if anorm < zero_r {
        info.* = -6;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    rcond.* = zero_r;
    if n == 0 {
        rcond.* = one_r;
        return;
    } else if anorm <= zero_r {
        return;
    }

    if upper {
        for <i:n..1 {
            if fget(ipiv, i) > 0 && fget(a, i, i) == zero_t {
                return;
            }
        }
    } else {
        for i:1..n {
            if fget(ipiv, i) > 0 && fget(a, i, i) == zero_t {
                return;
            }
        }
    }

    ainvnm = zero_r;
    kase = 0;
    x_view := _hecon_subvector(work.*, 1, n);
    v_view := _hecon_subvector(work.*, n + 1, n);
    work_mat := array2d_view(x_view, max(1, n), row_major=false);

    while true {
        _hecon_lacn2(n, *v_view, *x_view, *ainvnm, *kase, *isave, RT, MT);
        if kase == 0 then break;
        _hecon_hetrs(uplo, n, 1, a, lda, ipiv, *work_mat, n, info, RT, MT);
    }

    if ainvnm != zero_r {
        rcond.* = (one_r / ainvnm) / anorm;
    }
}
