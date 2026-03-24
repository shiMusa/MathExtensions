// ? LAPACK *lantb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slantb/clantb/dlantb/zlantb
// *LANTB returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the element of largest absolute value of a
// triangular band matrix.

slantb :: #bake_arguments _lantb(RT=float32, MT=float32);
dlantb :: #bake_arguments _lantb(RT=float64, MT=float64);
clantb :: #bake_arguments _lantb(RT=float32, MT=C32);
zlantb :: #bake_arguments _lantb(RT=float64, MT=C64);

#scope_file
_lantb_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lantb_matrix_column :: inline (
    ab: $AB/ArrayType2D,
    $$ldab: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] AB.data_type {
    return _lantb_raw_slice(ab.data + (col - 1)*ldab + (row - 1), count);
}

_lantb_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lantb_lassq :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    scale: *$RT,
    sumsq: *RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        slassq(n, x, incx, scale, sumsq);
    } else if MT == float64 {
        dlassq(n, x, incx, scale, sumsq);
    } else if MT == C32 {
        classq(n, x, incx, scale, sumsq);
    } else {
        zlassq(n, x, incx, scale, sumsq);
    }
}

_lantb :: (
    $$norm: u8,
    $$uplo: u8,
    $$diag: u8,
    $$n: int,
    $$k: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           WORK.data_type == RT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    i, j, l : int = ---;
    scale, sum, value : RT = ---;
    udiag : bool = ---;

    if n == 0 {
        value = zero_r;
    } else if lsame(norm, #char "M") {
        if lsame(diag, #char "U") {
            value = one_r;
            if lsame(uplo, #char "U") {
                for j:1..n {
                    for i:max(k + 2 - j, 1)..k {
                        sum = abs(fget(ab, i, j));
                        if value < sum || _lantb_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                }
            } else {
                for j:1..n {
                    for i:2..min(n + 1 - j, k + 1) {
                        sum = abs(fget(ab, i, j));
                        if value < sum || _lantb_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                }
            }
        } else {
            value = zero_r;
            if lsame(uplo, #char "U") {
                for j:1..n {
                    for i:max(k + 2 - j, 1)..k + 1 {
                        sum = abs(fget(ab, i, j));
                        if value < sum || _lantb_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                }
            } else {
                for j:1..n {
                    for i:1..min(n + 1 - j, k + 1) {
                        sum = abs(fget(ab, i, j));
                        if value < sum || _lantb_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                }
            }
        }
    } else if lsame(norm, #char "O") || norm == #char "1" {
        value = zero_r;
        udiag = lsame(diag, #char "U");
        if lsame(uplo, #char "U") {
            for j:1..n {
                if udiag {
                    sum = one_r;
                    for i:max(k + 2 - j, 1)..k {
                        sum = sum + abs(fget(ab, i, j));
                    }
                } else {
                    sum = zero_r;
                    for i:max(k + 2 - j, 1)..k + 1 {
                        sum = sum + abs(fget(ab, i, j));
                    }
                }
                if value < sum || _lantb_isnan(sum, RT) {
                    value = sum;
                }
            }
        } else {
            for j:1..n {
                if udiag {
                    sum = one_r;
                    for i:2..min(n + 1 - j, k + 1) {
                        sum = sum + abs(fget(ab, i, j));
                    }
                } else {
                    sum = zero_r;
                    for i:1..min(n + 1 - j, k + 1) {
                        sum = sum + abs(fget(ab, i, j));
                    }
                }
                if value < sum || _lantb_isnan(sum, RT) {
                    value = sum;
                }
            }
        }
    } else if lsame(norm, #char "I") {
        value = zero_r;
        if lsame(uplo, #char "U") {
            if lsame(diag, #char "U") {
                for i:1..n {
                    fset(work, i, one_r);
                }
                for j:1..n {
                    l = k + 1 - j;
                    for i:max(1, j-k)..j - 1 {
                        fset(work, i, fget(work.*, i) + abs(fget(ab, l + i, j)));
                    }
                }
            } else {
                for i:1..n {
                    fset(work, i, zero_r);
                }
                for j:1..n {
                    l = k + 1 - j;
                    for i:max(1, j-k)..j {
                        fset(work, i, fget(work.*, i) + abs(fget(ab, l + i, j)));
                    }
                }
            }
        } else {
            if lsame(diag, #char "U") {
                for i:1..n {
                    fset(work, i, one_r);
                }
                for j:1..n {
                    l = 1 - j;
                    for i:j + 1..min(n, j+k) {
                        fset(work, i, fget(work.*, i) + abs(fget(ab, l + i, j)));
                    }
                }
            } else {
                for i:1..n {
                    fset(work, i, zero_r);
                }
                for j:1..n {
                    l = 1 - j;
                    for i:j..min(n, j+k) {
                        fset(work, i, fget(work.*, i) + abs(fget(ab, l + i, j)));
                    }
                }
            }
        }
        for i:1..n {
            sum = fget(work.*, i);
            if value < sum || _lantb_isnan(sum, RT) {
                value = sum;
            }
        }
    } else {
        if lsame(uplo, #char "U") {
            if lsame(diag, #char "U") {
                scale = one_r;
                sum = ncast(RT, n);
                if k > 0 {
                    for j:2..n {
                        _lantb_lassq(min(j - 1, k), _lantb_matrix_column(ab, ldab, max(k + 2 - j, 1), j, min(j - 1, k)), 1, scale, sum, RT, MT);
                    }
                }
            } else {
                scale = zero_r;
                sum = one_r;
                for j:1..n {
                    _lantb_lassq(min(j, k + 1), _lantb_matrix_column(ab, ldab, max(k + 2 - j, 1), j, min(j, k + 1)), 1, scale, sum, RT, MT);
                }
            }
        } else {
            if lsame(diag, #char "U") {
                scale = one_r;
                sum = ncast(RT, n);
                if k > 0 {
                    for j:1..n - 1 {
                        _lantb_lassq(min(n - j, k), _lantb_matrix_column(ab, ldab, 2, j, min(n - j, k)), 1, scale, sum, RT, MT);
                    }
                }
            } else {
                scale = zero_r;
                sum = one_r;
                for j:1..n {
                    _lantb_lassq(min(n - j + 1, k + 1), _lantb_matrix_column(ab, ldab, 1, j, min(n - j + 1, k + 1)), 1, scale, sum, RT, MT);
                }
            }
        }
        value = scale*sqrt(sum);
    }

    return value;
}
