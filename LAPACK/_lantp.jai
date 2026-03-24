// ? LAPACK *lantp()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slantp/clantp/dlantp/zlantp
// *LANTP returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the element of largest absolute value of a
// triangular matrix supplied in packed form.

slantp :: #bake_arguments _lantp(RT=float32, MT=float32);
dlantp :: #bake_arguments _lantp(RT=float64, MT=float64);
clantp :: #bake_arguments _lantp(RT=float32, MT=C32);
zlantp :: #bake_arguments _lantp(RT=float64, MT=C64);

#scope_file
_lantp_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lantp_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lantp_lassq :: (
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

_lantp :: (
    $$norm: u8,
    $$uplo: u8,
    $$diag: u8,
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
    scale, sum, value : RT = ---;
    udiag : bool = ---;

    if n == 0 {
        value = zero_r;
    } else if lsame(norm, #char "M") {
        k = 1;
        if lsame(diag, #char "U") {
            value = one_r;
            if lsame(uplo, #char "U") {
                for j:1..n {
                    for i:k..k + j - 2 {
                        sum = abs(fget(ap, i));
                        if value < sum || _lantp_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                    k = k + j;
                }
            } else {
                for j:1..n {
                    for i:k + 1..k + n - j {
                        sum = abs(fget(ap, i));
                        if value < sum || _lantp_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                    k = k + n - j + 1;
                }
            }
        } else {
            value = zero_r;
            if lsame(uplo, #char "U") {
                for j:1..n {
                    for i:k..k + j - 1 {
                        sum = abs(fget(ap, i));
                        if value < sum || _lantp_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                    k = k + j;
                }
            } else {
                for j:1..n {
                    for i:k..k + n - j {
                        sum = abs(fget(ap, i));
                        if value < sum || _lantp_isnan(sum, RT) {
                            value = sum;
                        }
                    }
                    k = k + n - j + 1;
                }
            }
        }
    } else if lsame(norm, #char "O") || norm == #char "1" {
        value = zero_r;
        k = 1;
        udiag = lsame(diag, #char "U");
        if lsame(uplo, #char "U") {
            for j:1..n {
                if udiag {
                    sum = one_r;
                    for i:k..k + j - 2 {
                        sum = sum + abs(fget(ap, i));
                    }
                } else {
                    sum = zero_r;
                    for i:k..k + j - 1 {
                        sum = sum + abs(fget(ap, i));
                    }
                }
                k = k + j;
                if value < sum || _lantp_isnan(sum, RT) {
                    value = sum;
                }
            }
        } else {
            for j:1..n {
                if udiag {
                    sum = one_r;
                    for i:k + 1..k + n - j {
                        sum = sum + abs(fget(ap, i));
                    }
                } else {
                    sum = zero_r;
                    for i:k..k + n - j {
                        sum = sum + abs(fget(ap, i));
                    }
                }
                k = k + n - j + 1;
                if value < sum || _lantp_isnan(sum, RT) {
                    value = sum;
                }
            }
        }
    } else if lsame(norm, #char "I") {
        k = 1;
        if lsame(uplo, #char "U") {
            if lsame(diag, #char "U") {
                for i:1..n {
                    fset(work, i, one_r);
                }
                for j:1..n {
                    for i:1..j - 1 {
                        fset(work, i, fget(work.*, i) + abs(fget(ap, k)));
                        k = k + 1;
                    }
                    k = k + 1;
                }
            } else {
                for i:1..n {
                    fset(work, i, zero_r);
                }
                for j:1..n {
                    for i:1..j {
                        fset(work, i, fget(work.*, i) + abs(fget(ap, k)));
                        k = k + 1;
                    }
                }
            }
        } else {
            if lsame(diag, #char "U") {
                for i:1..n {
                    fset(work, i, one_r);
                }
                for j:1..n {
                    k = k + 1;
                    for i:j + 1..n {
                        fset(work, i, fget(work.*, i) + abs(fget(ap, k)));
                        k = k + 1;
                    }
                }
            } else {
                for i:1..n {
                    fset(work, i, zero_r);
                }
                for j:1..n {
                    for i:j..n {
                        fset(work, i, fget(work.*, i) + abs(fget(ap, k)));
                        k = k + 1;
                    }
                }
            }
        }
        value = zero_r;
        for i:1..n {
            sum = fget(work.*, i);
            if value < sum || _lantp_isnan(sum, RT) {
                value = sum;
            }
        }
    } else {
        if lsame(uplo, #char "U") {
            if lsame(diag, #char "U") {
                scale = one_r;
                sum = ncast(RT, n);
                k = 2;
                for j:2..n {
                    _lantp_lassq(j - 1, _lantp_raw_slice(ap.data + (k - 1), j - 1), 1, scale, sum, RT, MT);
                    k = k + j;
                }
            } else {
                scale = zero_r;
                sum = one_r;
                k = 1;
                for j:1..n {
                    _lantp_lassq(j, _lantp_raw_slice(ap.data + (k - 1), j), 1, scale, sum, RT, MT);
                    k = k + j;
                }
            }
        } else {
            if lsame(diag, #char "U") {
                scale = one_r;
                sum = ncast(RT, n);
                k = 2;
                for j:1..n - 1 {
                    _lantp_lassq(n - j, _lantp_raw_slice(ap.data + (k - 1), n - j), 1, scale, sum, RT, MT);
                    k = k + n - j + 1;
                }
            } else {
                scale = zero_r;
                sum = one_r;
                k = 1;
                for j:1..n {
                    _lantp_lassq(n - j + 1, _lantp_raw_slice(ap.data + (k - 1), n - j + 1), 1, scale, sum, RT, MT);
                    k = k + n - j + 1;
                }
            }
        }
        value = scale*sqrt(sum);
    }

    return value;
}
