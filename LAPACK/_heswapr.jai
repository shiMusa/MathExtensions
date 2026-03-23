// ? LAPACK *heswapr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheswapr/zheswapr
// *HESWAPR applies an elementary permutation on the rows and columns of
// a Hermitian matrix.

cheswapr :: #bake_arguments _heswapr(RT=float32, MT=C32);
zheswapr :: #bake_arguments _heswapr(RT=float64, MT=C64);

#scope_file
_heswapr_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_heswapr_vector_view :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _heswapr_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_heswapr_swap :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_heswapr :: (
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$i1: int,
    $$i2: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT;
} {
    upper : bool = ---;
    i : int = ---;
    tmp : MT = ---;

    upper = lsame(uplo, #char "U");
    if upper {
        v1 := _heswapr_vector_view(a.*, lda, 1, i1, i1 - 1);
        v2 := _heswapr_vector_view(a.*, lda, 1, i2, i1 - 1);
        _heswapr_swap(i1 - 1, *v1, 1, *v2, 1, MT);

        tmp = fget(a.*, i1, i1);
        fset(a, i1, i1, fget(a.*, i2, i2));
        fset(a, i2, i2, tmp);

        for i:1..i2 - i1 - 1 {
            tmp = fget(a.*, i1, i1 + i);
            fset(a, i1, i1 + i, conjugate(fget(a.*, i1 + i, i2)));
            fset(a, i1 + i, i2, conjugate(tmp));
        }

        fset(a, i1, i2, conjugate(fget(a.*, i1, i2)));

        for i:i2 + 1..n {
            tmp = fget(a.*, i1, i);
            fset(a, i1, i, fget(a.*, i2, i));
            fset(a, i2, i, tmp);
        }
    } else {
        v1 := _heswapr_vector_view(a.*, lda, i1, 1, 1 + (i1 - 2)*lda);
        v2 := _heswapr_vector_view(a.*, lda, i2, 1, 1 + (i1 - 2)*lda);
        _heswapr_swap(i1 - 1, *v1, lda, *v2, lda, MT);

        tmp = fget(a.*, i1, i1);
        fset(a, i1, i1, fget(a.*, i2, i2));
        fset(a, i2, i2, tmp);

        for i:1..i2 - i1 - 1 {
            tmp = fget(a.*, i1 + i, i1);
            fset(a, i1 + i, i1, conjugate(fget(a.*, i2, i1 + i)));
            fset(a, i2, i1 + i, conjugate(tmp));
        }

        fset(a, i2, i1, conjugate(fget(a.*, i2, i1)));

        for i:i2 + 1..n {
            tmp = fget(a.*, i, i1);
            fset(a, i, i1, fget(a.*, i, i2));
            fset(a, i, i2, tmp);
        }
    }
}
