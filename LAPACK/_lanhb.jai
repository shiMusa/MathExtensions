// ? LAPACK *lanhb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clanhb/zlanhb
// *LANHB returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// Hermitian band matrix.

clanhb :: #bake_arguments _lanhb(RT=float32, MT=C32);
zlanhb :: #bake_arguments _lanhb(RT=float64, MT=C64);

#scope_file
_lanhb_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lanhb_matrix_column :: inline (
    ab: $AB/ArrayType2D,
    $$ldab: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] AB.data_type {
    return _lanhb_raw_slice(ab.data + (col - 1)*ldab + (row - 1), count);
}

_lanhb_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lanhb_lassq :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    scale: *$RT,
    sumsq: *RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        classq(n, x, incx, scale, sumsq);
    } else {
        zlassq(n, x, incx, scale, sumsq);
    }
}

_lanhb :: (
    $$norm: u8,
    $$uplo: u8,
    $$n: int,
    $$k: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           WORK.data_type == RT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    i, j, l : int = ---;
    absa, scale, sum, value : RT = ---;

    if n == 0 {
        value = zero_r;
    } else if lsame(norm, #char "M") {
        value = zero_r;
        if lsame(uplo, #char "U") {
            for j:1..n {
                for i:max(k + 2 - j, 1)..k {
                    sum = abs(fget(ab, i, j));
                    if value < sum || _lanhb_isnan(sum, RT) {
                        value = sum;
                    }
                }
                sum = abs(real(fget(ab, k + 1, j)));
                if value < sum || _lanhb_isnan(sum, RT) {
                    value = sum;
                }
            }
        } else {
            for j:1..n {
                sum = abs(real(fget(ab, 1, j)));
                if value < sum || _lanhb_isnan(sum, RT) {
                    value = sum;
                }
                for i:2..min(n + 1 - j, k + 1) {
                    sum = abs(fget(ab, i, j));
                    if value < sum || _lanhb_isnan(sum, RT) {
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
                l = k + 1 - j;
                for i:max(1, j-k)..j - 1 {
                    absa = abs(fget(ab, l + i, j));
                    sum = sum + absa;
                    fset(work, i, fget(work.*, i) + absa);
                }
                fset(work, j, sum + abs(real(fget(ab, k + 1, j))));
            }
            for i:1..n {
                sum = fget(work.*, i);
                if value < sum || _lanhb_isnan(sum, RT) {
                    value = sum;
                }
            }
        } else {
            for i:1..n {
                fset(work, i, zero_r);
            }
            for j:1..n {
                sum = fget(work.*, j) + abs(real(fget(ab, 1, j)));
                l = 1 - j;
                for i:j + 1..min(n, j + k) {
                    absa = abs(fget(ab, l + i, j));
                    sum = sum + absa;
                    fset(work, i, fget(work.*, i) + absa);
                }
                if value < sum || _lanhb_isnan(sum, RT) {
                    value = sum;
                }
            }
        }
    } else {
        scale = zero_r;
        sum = one_r;
        if k > 0 {
            if lsame(uplo, #char "U") {
                for j:2..n {
                    _lanhb_lassq(min(j - 1, k), _lanhb_matrix_column(ab, ldab, max(k + 2 - j, 1), j, min(j - 1, k)), 1, scale, sum, RT, MT);
                }
                l = k + 1;
            } else {
                for j:1..n - 1 {
                    _lanhb_lassq(min(n - j, k), _lanhb_matrix_column(ab, ldab, 2, j, min(n - j, k)), 1, scale, sum, RT, MT);
                }
                l = 1;
            }
            sum = ncast(RT, 2)*sum;
        } else {
            l = 1;
        }
        for j:1..n {
            if real(fget(ab, l, j)) != zero_r {
                absa = abs(real(fget(ab, l, j)));
                if scale < absa {
                    sum = one_r + sum*(scale / absa)*(scale / absa);
                    scale = absa;
                } else {
                    sum = sum + (absa / scale)*(absa / scale);
                }
            }
        }
        value = scale*sqrt(sum);
    }

    return value;
}
