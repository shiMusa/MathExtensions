// ? LAPACK *lansy()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slansy/clansy/dlansy/zlansy
// *LANSY returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// symmetric matrix.

slansy :: #bake_arguments _lansy(RT=float32, MT=float32);
dlansy :: #bake_arguments _lansy(RT=float64, MT=float64);
clansy :: #bake_arguments _lansy(RT=float32, MT=C32);
zlansy :: #bake_arguments _lansy(RT=float64, MT=C64);

#scope_file
_lansy_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lansy_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _lansy_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lansy_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int,
    $$inc: int
) -> [] A.data_type {
    storage_count := ifx count <= 0 then 0 else 1 + (count - 1)*inc;
    return _lansy_raw_slice(a.data + (col - 1)*lda + (row - 1), storage_count);
}

_lansy_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lansy_lassq :: (
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

_lansy :: (
    $$norm: u8,
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           WORK.data_type == RT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    i, j : int = ---;
    absa, scale, sum, value : RT = ---;

    if n == 0 {
        value = zero_r;
    } else if lsame(norm, #char "M") {
        value = zero_r;
        if lsame(uplo, #char "U") {
            for j:1..n {
                for i:1..j {
                    sum = abs(fget(a, i, j));
                    if value < sum || _lansy_isnan(sum, RT) {
                        value = sum;
                    }
                }
            }
        } else {
            for j:1..n {
                for i:j..n {
                    sum = abs(fget(a, i, j));
                    if value < sum || _lansy_isnan(sum, RT) {
                        value = sum;
                    }
                }
            }
        }
    } else if lsame(norm, #char "I") || lsame(norm, #char "O") || norm == #char "1" {
        value = zero_r;
        if lsame(uplo, #char "U") {
            for j:1..n {
                sum = zero_r;
                for i:1..j - 1 {
                    absa = abs(fget(a, i, j));
                    sum = sum + absa;
                    fset(work, i, fget(work.*, i) + absa);
                }
                fset(work, j, sum + abs(fget(a, j, j)));
            }
            for i:1..n {
                sum = fget(work.*, i);
                if value < sum || _lansy_isnan(sum, RT) {
                    value = sum;
                }
            }
        } else {
            for i:1..n {
                fset(work, i, zero_r);
            }
            for j:1..n {
                sum = fget(work.*, j) + abs(fget(a, j, j));
                for i:j + 1..n {
                    absa = abs(fget(a, i, j));
                    sum = sum + absa;
                    fset(work, i, fget(work.*, i) + absa);
                }
                if value < sum || _lansy_isnan(sum, RT) {
                    value = sum;
                }
            }
        }
    } else {
        scale = zero_r;
        sum = one_r;
        if lsame(uplo, #char "U") {
            for j:2..n {
                _lansy_lassq(j - 1, _lansy_matrix_column(a, lda, 1, j, j - 1), 1, scale, sum, RT, MT);
            }
        } else {
            for j:1..n - 1 {
                _lansy_lassq(n - j, _lansy_matrix_column(a, lda, j + 1, j, n - j), 1, scale, sum, RT, MT);
            }
        }
        sum = ncast(RT, 2)*sum;
        _lansy_lassq(n, _lansy_matrix_storage(a, lda, 1, 1, n, lda + 1), lda + 1, scale, sum, RT, MT);
        value = scale*sqrt(sum);
    }

    return value;
}
