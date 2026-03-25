// ? LAPACK *larfy()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarfy/dlarfy/clarfy/zlarfy
// *LARFY applies an elementary reflector to a symmetric or Hermitian matrix.

slarfy :: #bake_arguments _larfy(RT=float32, MT=float32);
dlarfy :: #bake_arguments _larfy(RT=float64, MT=float64);
clarfy :: #bake_arguments _larfy(RT=float32, MT=C32);
zlarfy :: #bake_arguments _larfy(RT=float64, MT=C64);

#scope_file
_larfy_symv :: (
    $$uplo: u8, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: $X/ArrayType1D, $$incx: int, beta: MT, y: *$Y/ArrayType1D, $$incy: int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == float64 {
        dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == C32 {
        chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_larfy_axpy :: (
    $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        saxpy(n, alpha, x, incx, y, incy);
    } else if MT == float64 {
        daxpy(n, alpha, x, incx, y, incy);
    } else if MT == C32 {
        caxpy(n, alpha, x, incx, y, incy);
    } else {
        zaxpy(n, alpha, x, incx, y, incy);
    }
}

_larfy_syr2_update :: (
    $$uplo: u8, $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int,
    y: $Y/ArrayType1D, $$incy: int, a: *$A/ArrayType2D, $$lda: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT && Y.data_type == MT && A.data_type == MT;
} {
    #if MT == float32 {
        ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda);
    } else if MT == float64 {
        dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda);
    } else if MT == C32 {
        cher2(uplo, n, alpha, x, incx, y, incy, a, lda);
    } else {
        zher2(uplo, n, alpha, x, incx, y, incy, a, lda);
    }
}

_larfy_dot :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && Y.data_type == MT;
} -> MT {
    #if MT == float32 {
        return sdot(n, x, incx, y, incy);
    } else if MT == float64 {
        return ddot(n, x, incx, y, incy);
    } else if MT == C32 {
        return cdotc(n, x, incx, y, incy);
    } else {
        return zdotc(n, x, incx, y, incy);
    }
}

_larfy :: (
    $$uplo: u8, $$n: int, v: $V/ArrayType1D, $$incv: int, tau: $MT, c: *$C/ArrayType2D, $$ldc: int,
    work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    zero_r :: #run zero(RT);
    one_t  :: #run one(MT);
    zero_t :: #run zero(MT);
    half_t : MT = #if MT == float32 || MT == float64 { ncast(MT, 0.5) } else { cpx(ncast(RT, 0.5), zero_r) };

    alpha : MT = ---;

    if tau == zero_t then return;

    _larfy_symv(uplo, n, one_t, c.*, ldc, v, incv, zero_t, work, 1, RT, MT);

    alpha = -half_t*tau*_larfy_dot(n, work.*, 1, v, incv, MT);
    _larfy_axpy(n, alpha, v, incv, work, 1, MT);

    _larfy_syr2_update(uplo, n, -tau, v, incv, work.*, 1, c, ldc, RT, MT);
}
