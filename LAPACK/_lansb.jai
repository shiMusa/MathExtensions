// ? LAPACK *lansb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slansb/dlansb
// *LANSB returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// symmetric band matrix.

slansb :: #bake_arguments _lansb(RT=float32, MT=float32);
dlansb :: #bake_arguments _lansb(RT=float64, MT=float64);
clansb :: #bake_arguments _lansb(RT=float32, MT=C32);
zlansb :: #bake_arguments _lansb(RT=float64, MT=C64);

#scope_file
_lansb_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lansb_matrix_column :: inline (
    ab: $AB/ArrayType2D,
    $$ldab: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] AB.data_type {
    return _lansb_raw_slice(ab.data + (col - 1)*ldab + (row - 1), count);
}

_lansb_matrix_storage :: inline (
    ab: $AB/ArrayType2D,
    $$ldab: int,
    $$row: int,
    $$col: int,
    $$count: int,
    $$inc: int
) -> [] AB.data_type {
    storage_count := ifx count <= 0 then 0 else 1 + (count - 1)*inc;
    return _lansb_raw_slice(ab.data + (col - 1)*ldab + (row - 1), storage_count);
}

_lansb_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lansb_lassq :: (
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

_lansb :: (
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
           (MT == RT || MT == complex_with_inner_type(RT)) &&
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
                for i:max(k + 2 - j, 1)..k + 1 {
                    sum = abs(fget(ab, i, j));
                    if value < sum || _lansb_isnan(sum, RT) {
                        value = sum;
                    }
                }
            }
        } else {
            for j:1..n {
                for i:1..min(n + 1 - j, k + 1) {
                    sum = abs(fget(ab, i, j));
                    if value < sum || _lansb_isnan(sum, RT) {
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
                fset(work, j, sum + abs(fget(ab, k + 1, j)));
            }
            for i:1..n {
                sum = fget(work.*, i);
                if value < sum || _lansb_isnan(sum, RT) {
                    value = sum;
                }
            }
        } else {
            for i:1..n {
                fset(work, i, zero_r);
            }
            for j:1..n {
                sum = fget(work.*, j) + abs(fget(ab, 1, j));
                l = 1 - j;
                for i:j + 1..min(n, j + k) {
                    absa = abs(fget(ab, l + i, j));
                    sum = sum + absa;
                    fset(work, i, fget(work.*, i) + absa);
                }
                if value < sum || _lansb_isnan(sum, RT) {
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
                    _lansb_lassq(min(j - 1, k), _lansb_matrix_column(ab, ldab, max(k + 2 - j, 1), j, min(j - 1, k)), 1, scale, sum, RT, MT);
                }
                l = k + 1;
            } else {
                for j:1..n - 1 {
                    _lansb_lassq(min(n - j, k), _lansb_matrix_column(ab, ldab, 2, j, min(n - j, k)), 1, scale, sum, RT, MT);
                }
                l = 1;
            }
            sum = ncast(RT, 2)*sum;
        } else {
            l = 1;
        }
        _lansb_lassq(n, _lansb_matrix_storage(ab, ldab, l, 1, n, ldab), ldab, scale, sum, RT, MT);
        value = scale*sqrt(sum);
    }

    return value;
}
