// ? LAPACK *lantr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slantr/clantr/dlantr/zlantr
// *LANTR returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the element of largest absolute value of a
// trapezoidal or triangular matrix.

slantr :: #bake_arguments _lantr(RT=float32, MT=float32);
dlantr :: #bake_arguments _lantr(RT=float64, MT=float64);
clantr :: #bake_arguments _lantr(RT=float32, MT=C32);
zlantr :: #bake_arguments _lantr(RT=float64, MT=C64);

#scope_file
_lantr_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lantr_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _lantr_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lantr_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lantr_lassq :: (
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

_lantr :: (
    $$norm: u8,
    $$uplo: u8,
    $$diag: u8,
    $$m: int,
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
    scale, sum, value : RT = ---;
    udiag : bool = ---;

    if min(m, n) == 0 {
        value = zero_r;
    } else if lsame(norm, #char "M") {
        if lsame(diag, #char "U") {
            value = one_r;
            if lsame(uplo, #char "U") {
                for j:1..n {
                    for i:1..min(m, j - 1) {
                        sum = abs(fget(a, i, j));
                        if value < sum || _lantr_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                }
            } else {
                for j:1..n {
                    for i:j + 1..m {
                        sum = abs(fget(a, i, j));
                        if value < sum || _lantr_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                }
            }
        } else {
            value = zero_r;
            if lsame(uplo, #char "U") {
                for j:1..n {
                    for i:1..min(m, j) {
                        sum = abs(fget(a, i, j));
                        if value < sum || _lantr_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                }
            } else {
                for j:1..n {
                    for i:j..m {
                        sum = abs(fget(a, i, j));
                        if value < sum || _lantr_isnan(sum, RT) {
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
                if udiag && j <= m {
                    sum = one_r;
                    for i:1..j - 1 {
                        sum = sum + abs(fget(a, i, j));
                    }
                } else {
                    sum = zero_r;
                    for i:1..min(m, j) {
                        sum = sum + abs(fget(a, i, j));
                    }
                }
                if value < sum || _lantr_isnan(sum, RT) {
                    value = sum;
                }
            }
        } else {
            for j:1..n {
                if udiag {
                    sum = one_r;
                    for i:j + 1..m {
                        sum = sum + abs(fget(a, i, j));
                    }
                } else {
                    sum = zero_r;
                    for i:j..m {
                        sum = sum + abs(fget(a, i, j));
                    }
                }
                if value < sum || _lantr_isnan(sum, RT) {
                    value = sum;
                }
            }
        }
    } else if lsame(norm, #char "I") {
        if lsame(uplo, #char "U") {
            if lsame(diag, #char "U") {
                for i:1..m {
                    fset(work, i, one_r);
                }
                for j:1..n {
                    for i:1..min(m, j - 1) {
                        fset(work, i, fget(work.*, i) + abs(fget(a, i, j)));
                    }
                }
            } else {
                for i:1..m {
                    fset(work, i, zero_r);
                }
                for j:1..n {
                    for i:1..min(m, j) {
                        fset(work, i, fget(work.*, i) + abs(fget(a, i, j)));
                    }
                }
            }
        } else {
            if lsame(diag, #char "U") {
                for i:1..min(m, n) {
                    fset(work, i, one_r);
                }
                for i:n + 1..m {
                    fset(work, i, zero_r);
                }
                for j:1..n {
                    for i:j + 1..m {
                        fset(work, i, fget(work.*, i) + abs(fget(a, i, j)));
                    }
                }
            } else {
                for i:1..m {
                    fset(work, i, zero_r);
                }
                for j:1..n {
                    for i:j..m {
                        fset(work, i, fget(work.*, i) + abs(fget(a, i, j)));
                    }
                }
            }
        }
        value = zero_r;
        for i:1..m {
            sum = fget(work.*, i);
            if value < sum || _lantr_isnan(sum, RT) {
                value = sum;
            }
        }
    } else {
        if lsame(uplo, #char "U") {
            if lsame(diag, #char "U") {
                scale = one_r;
                sum = ncast(RT, min(m, n));
                for j:2..n {
                    _lantr_lassq(min(m, j - 1), _lantr_matrix_column(a, lda, 1, j, min(m, j - 1)), 1, scale, sum, RT, MT);
                }
            } else {
                scale = zero_r;
                sum = one_r;
                for j:1..n {
                    _lantr_lassq(min(m, j), _lantr_matrix_column(a, lda, 1, j, min(m, j)), 1, scale, sum, RT, MT);
                }
            }
        } else {
            if lsame(diag, #char "U") {
                scale = one_r;
                sum = ncast(RT, min(m, n));
                for j:1..n {
                    _lantr_lassq(m - j, _lantr_matrix_column(a, lda, min(m, j + 1), j, m - j), 1, scale, sum, RT, MT);
                }
            } else {
                scale = zero_r;
                sum = one_r;
                for j:1..n {
                    _lantr_lassq(m - j + 1, _lantr_matrix_column(a, lda, j, j, m - j + 1), 1, scale, sum, RT, MT);
                }
            }
        }
        value = scale*sqrt(sum);
    }

    return value;
}
