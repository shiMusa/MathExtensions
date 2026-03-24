// ? LAPACK *lansp()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slansp/clansp/dlansp/zlansp
// *LANSP returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// symmetric matrix supplied in packed form.

slansp :: #bake_arguments _lansp(RT=float32, MT=float32);
dlansp :: #bake_arguments _lansp(RT=float64, MT=float64);
clansp :: #bake_arguments _lansp(RT=float32, MT=C32);
zlansp :: #bake_arguments _lansp(RT=float64, MT=C64);

#scope_file
_lansp_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lansp_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lansp_lassq :: (
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

_lansp :: (
    $$norm: u8,
    $$uplo: u8,
    $$n: int,
    ap: $AP/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AP.data_type == MT &&
           WORK.data_type == RT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    i, j, k : int = ---;
    absa, scale, sum, value : RT = ---;

    if n == 0 {
        value = zero_r;
    } else if lsame(norm, #char "M") {
        value = zero_r;
        k = 1;
        if lsame(uplo, #char "U") {
            for j:1..n {
                for i:k..k + j - 1 {
                    sum = abs(fget(ap, i));
                    if value < sum || _lansp_isnan(sum, RT) {
                        value = sum;
                    }
                }
                k = k + j;
            }
        } else {
            for j:1..n {
                for i:k..k + n - j {
                    sum = abs(fget(ap, i));
                    if value < sum || _lansp_isnan(sum, RT) {
                        value = sum;
                    }
                }
                k = k + n - j + 1;
            }
        }
    } else if lsame(norm, #char "I") || lsame(norm, #char "O") || norm == #char "1" {
        value = zero_r;
        k = 1;
        if lsame(uplo, #char "U") {
            for j:1..n {
                sum = zero_r;
                for i:1..j - 1 {
                    absa = abs(fget(ap, k));
                    sum = sum + absa;
                    fset(work, i, fget(work.*, i) + absa);
                    k = k + 1;
                }
                fset(work, j, sum + abs(fget(ap, k)));
                k = k + 1;
            }
            for i:1..n {
                sum = fget(work.*, i);
                if value < sum || _lansp_isnan(sum, RT) {
                    value = sum;
                }
            }
        } else {
            for i:1..n {
                fset(work, i, zero_r);
            }
            for j:1..n {
                sum = fget(work.*, j) + abs(fget(ap, k));
                k = k + 1;
                for i:j + 1..n {
                    absa = abs(fget(ap, k));
                    sum = sum + absa;
                    fset(work, i, fget(work.*, i) + absa);
                    k = k + 1;
                }
                if value < sum || _lansp_isnan(sum, RT) {
                    value = sum;
                }
            }
        }
    } else {
        scale = zero_r;
        sum = one_r;
        k = 2;
        if lsame(uplo, #char "U") {
            for j:2..n {
                _lansp_lassq(j - 1, _lansp_raw_slice(ap.data + (k - 1), j - 1), 1, scale, sum, RT, MT);
                k = k + j;
            }
        } else {
            for j:1..n - 1 {
                _lansp_lassq(n - j, _lansp_raw_slice(ap.data + (k - 1), n - j), 1, scale, sum, RT, MT);
                k = k + n - j + 1;
            }
        }
        sum = ncast(RT, 2)*sum;
        k = 1;
        for i:1..n {
            _lansp_lassq(1, _lansp_raw_slice(ap.data + (k - 1), 1), 1, scale, sum, RT, MT);
            if lsame(uplo, #char "U") {
                k = k + i + 1;
            } else {
                k = k + n - i + 1;
            }
        }
        value = scale*sqrt(sum);
    }

    return value;
}
