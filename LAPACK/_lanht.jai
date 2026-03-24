// ? LAPACK *lanht()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clanht/zlanht
// *LANHT returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// Hermitian tridiagonal matrix.

clanht :: #bake_arguments _lanht(RT=float32, MT=C32);
zlanht :: #bake_arguments _lanht(RT=float64, MT=C64);

#scope_file
_lanht_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lanht_complex_lassq :: (
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

_lanht_real_lassq :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    scale: *$RT,
    sumsq: *RT,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        slassq(n, x, incx, scale, sumsq);
    } else {
        dlassq(n, x, incx, scale, sumsq);
    }
}

_lanht :: (
    $$norm: u8,
    $$n: int,
    d: $D/ArrayType1D,
    e: $E/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT &&
           E.data_type == MT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    i : int = ---;
    anorm, scale, sum : RT = ---;

    if n <= 0 {
        anorm = zero_r;
    } else if lsame(norm, #char "M") {
        anorm = abs(fget(d, n));
        for i:1..n - 1 {
            sum = abs(fget(d, i));
            if anorm < sum || _lanht_isnan(sum, RT) {
                anorm = sum;
            }
            sum = abs(fget(e, i));
            if anorm < sum || _lanht_isnan(sum, RT) {
                anorm = sum;
            }
        }
    } else if lsame(norm, #char "O") || norm == #char "1" || lsame(norm, #char "I") {
        if n == 1 {
            anorm = abs(fget(d, 1));
        } else {
            anorm = abs(fget(d, 1)) + abs(fget(e, 1));
            sum = abs(fget(e, n - 1)) + abs(fget(d, n));
            if anorm < sum || _lanht_isnan(sum, RT) {
                anorm = sum;
            }
            for i:2..n - 1 {
                sum = abs(fget(d, i)) + abs(fget(e, i)) + abs(fget(e, i - 1));
                if anorm < sum || _lanht_isnan(sum, RT) {
                    anorm = sum;
                }
            }
        }
    } else {
        scale = zero_r;
        sum = one_r;
        if n > 1 {
            _lanht_complex_lassq(n - 1, e, 1, scale, sum, RT, MT);
            sum = ncast(RT, 2)*sum;
        }
        _lanht_real_lassq(n, d, 1, scale, sum, RT);
        anorm = scale*sqrt(sum);
    }

    return anorm;
}
