// ? LAPACK *hetd2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetd2/zhetd2
// *HETD2 reduces a Hermitian matrix to real symmetric tridiagonal form
// using an unblocked algorithm.

chetd2 :: #bake_arguments _hetd2(RT=float32, MT=C32);
zhetd2 :: #bake_arguments _hetd2(RT=float64, MT=C64);

#scope_file
_hetd2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetd2_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _hetd2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetd2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetd2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetd2_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetd2_raw_slice(x.data + (offset - 1), count);
}

_hetd2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetd2_larfg :: (
    $$n: int,
    alpha: *$MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    tau: *$MT,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        clarfg(n, alpha, x, incx, tau);
    } else {
        zlarfg(n, alpha, x, incx, tau);
    }
}

_hetd2_hemv :: (
    $$uplo: u8,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    $$incx: int,
    beta: $MT,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_hetd2_dotc :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) -> MT #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        return cdotc(n, x, incx, y, incy);
    } else {
        return zdotc(n, x, incx, y, incy);
    }
}

_hetd2_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        caxpy(n, alpha, x, incx, y, incy);
    } else {
        zaxpy(n, alpha, x, incx, y, incy);
    }
}

_hetd2_her2 :: (
    $$uplo: u8,
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT &&
           A.data_type == MT;
} {
    #if MT == C32 {
        cher2(uplo, n, alpha, x, incx, y, incy, a, lda);
    } else {
        zher2(uplo, n, alpha, x, incx, y, incy, a, lda);
    }
}

_hetd2 :: (
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           TAU.data_type == MT;
} {
    zero_t :: #run zero(MT);
    one_r  :: #run one(RT);
    half_r :: #run ncast(RT) 0.5;
    one_t  :: #run _hetd2_complex_from_real(one_r, RT, MT);
    half_t :: #run _hetd2_complex_from_real(half_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHETD2" } else { "ZHETD2" };

    upper : bool = ---;
    i : int = ---;
    alpha, taui : MT = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n <= 0 then return;

    if upper {
        fset(a, n, n, _hetd2_complex_from_real(real(fget(a.*, n, n)), RT, MT));
        for <i:1..n - 1 {
            alpha = fget(a.*, i, i + 1);
            col_tail := _hetd2_matrix_column(a.*, lda, 1, i + 1, i - 1);
            _hetd2_larfg(i, *alpha, *col_tail, 1, *taui, MT);
            fset(e, i, real(alpha));

            if taui != zero_t {
                fset(a, i, i + 1, one_t);

                a11_raw := _hetd2_matrix_storage(a.*, lda, 1, 1, i, i);
                a11 := array2d_view(a11_raw, lda, row_major=false);
                v := _hetd2_matrix_column(a.*, lda, 1, i + 1, i);
                tau_head := _hetd2_vector_slice(tau, 1, i);
                _hetd2_hemv(uplo, i, taui, a11, lda, v, 1, zero_t, *tau_head, 1, MT);
                alpha = -half_t*taui*_hetd2_dotc(i, tau_head, 1, v, 1, MT);
                _hetd2_axpy(i, alpha, v, 1, *tau_head, 1, MT);
                _hetd2_her2(uplo, i, -one_t, v, 1, tau_head, 1, a11, lda, MT);
            } else {
                fset(a, i, i, _hetd2_complex_from_real(real(fget(a.*, i, i)), RT, MT));
            }

            fset(a, i, i + 1, _hetd2_complex_from_real(fget(e.*, i), RT, MT));
            fset(d, i + 1, real(fget(a.*, i + 1, i + 1)));
            fset(tau, i, taui);
        }
        fset(d, 1, real(fget(a.*, 1, 1)));
    } else {
        fset(a, 1, 1, _hetd2_complex_from_real(real(fget(a.*, 1, 1)), RT, MT));
        for i:1..n - 1 {
            alpha = fget(a.*, i + 1, i);
            col_tail := _hetd2_matrix_column(a.*, lda, min(i + 2, n), i, n - i - 1);
            _hetd2_larfg(n - i, *alpha, *col_tail, 1, *taui, MT);
            fset(e, i, real(alpha));

            if taui != zero_t {
                fset(a, i + 1, i, one_t);

                a22_raw := _hetd2_matrix_storage(a.*, lda, i + 1, i + 1, n - i, n - i);
                a22 := array2d_view(a22_raw, lda, row_major=false);
                v := _hetd2_matrix_column(a.*, lda, i + 1, i, n - i);
                tau_tail := _hetd2_vector_slice(tau, i, n - i);
                _hetd2_hemv(uplo, n - i, taui, a22, lda, v, 1, zero_t, *tau_tail, 1, MT);
                alpha = -half_t*taui*_hetd2_dotc(n - i, tau_tail, 1, v, 1, MT);
                _hetd2_axpy(n - i, alpha, v, 1, *tau_tail, 1, MT);
                _hetd2_her2(uplo, n - i, -one_t, v, 1, tau_tail, 1, a22, lda, MT);
            } else {
                fset(a, i + 1, i + 1, _hetd2_complex_from_real(real(fget(a.*, i + 1, i + 1)), RT, MT));
            }

            fset(a, i + 1, i, _hetd2_complex_from_real(fget(e.*, i), RT, MT));
            fset(d, i, real(fget(a.*, i, i)));
            fset(tau, i, taui);
        }
        fset(d, n, real(fget(a.*, n, n)));
    }
}
