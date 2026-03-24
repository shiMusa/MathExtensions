// ? LAPACK *langt()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slangt/clangt/dlangt/zlangt
// *LANGT returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// general tridiagonal matrix.

slangt :: #bake_arguments _langt(RT=float32, MT=float32);
dlangt :: #bake_arguments _langt(RT=float64, MT=float64);
clangt :: #bake_arguments _langt(RT=float32, MT=C32);
zlangt :: #bake_arguments _langt(RT=float64, MT=C64);

#scope_file
_langt_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_langt_lassq :: (
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

_langt :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    i : int = ---;
    anorm, scale, sum, temp : RT = ---;

    if n <= 0 {
        anorm = zero_r;
    } else if lsame(norm, #char "M") {
        anorm = abs(fget(d, n));
        for i:1..n - 1 {
            temp = abs(fget(dl, i));
            if anorm < temp || _langt_isnan(temp, RT) {
                anorm = temp;
            }
            temp = abs(fget(d, i));
            if anorm < temp || _langt_isnan(temp, RT) {
                anorm = temp;
            }
            temp = abs(fget(du, i));
            if anorm < temp || _langt_isnan(temp, RT) {
                anorm = temp;
            }
        }
    } else if lsame(norm, #char "O") || norm == #char "1" {
        if n == 1 {
            anorm = abs(fget(d, 1));
        } else {
            anorm = abs(fget(d, 1)) + abs(fget(dl, 1));
            temp = abs(fget(d, n)) + abs(fget(du, n - 1));
            if anorm < temp || _langt_isnan(temp, RT) {
                anorm = temp;
            }
            for i:2..n - 1 {
                temp = abs(fget(d, i)) + abs(fget(dl, i)) + abs(fget(du, i - 1));
                if anorm < temp || _langt_isnan(temp, RT) {
                    anorm = temp;
                }
            }
        }
    } else if lsame(norm, #char "I") {
        if n == 1 {
            anorm = abs(fget(d, 1));
        } else {
            anorm = abs(fget(d, 1)) + abs(fget(du, 1));
            temp = abs(fget(d, n)) + abs(fget(dl, n - 1));
            if anorm < temp || _langt_isnan(temp, RT) {
                anorm = temp;
            }
            for i:2..n - 1 {
                temp = abs(fget(d, i)) + abs(fget(du, i)) + abs(fget(dl, i - 1));
                if anorm < temp || _langt_isnan(temp, RT) {
                    anorm = temp;
                }
            }
        }
    } else {
        scale = zero_r;
        sum = one_r;
        _langt_lassq(n, d, 1, scale, sum, RT, MT);
        if n > 1 {
            _langt_lassq(n - 1, dl, 1, scale, sum, RT, MT);
            _langt_lassq(n - 1, du, 1, scale, sum, RT, MT);
        }
        anorm = scale*sqrt(sum);
    }

    return anorm;
}
