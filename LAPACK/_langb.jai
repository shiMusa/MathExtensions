// ? LAPACK *langb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slangb/clangb/dlangb/zlangb
// *LANGB returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// general band matrix.

slangb :: #bake_arguments _langb(RT=float32, MT=float32);
dlangb :: #bake_arguments _langb(RT=float64, MT=float64);
clangb :: #bake_arguments _langb(RT=float32, MT=C32);
zlangb :: #bake_arguments _langb(RT=float64, MT=C64);

#scope_file
_langb_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_langb_matrix_column :: inline (
    ab: $AB/ArrayType2D,
    $$ldab: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] AB.data_type {
    return _langb_raw_slice(ab.data + (col - 1)*ldab + (row - 1), count);
}

_langb_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_langb_lassq :: (
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

_langb :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
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

    i, j, k, l : int = ---;
    scale, sum, temp, value : RT = ---;

    if n == 0 {
        value = zero_r;
    } else if lsame(norm, #char "M") {
        value = zero_r;
        for j:1..n {
            for i:max(ku + 2 - j, 1)..min(n + ku + 1 - j, kl + ku + 1) {
                temp = abs(fget(ab, i, j));
                if value < temp || _langb_isnan(temp, RT) {
                    value = temp;
                }
            }
        }
    } else if lsame(norm, #char "O") || norm == #char "1" {
        value = zero_r;
        for j:1..n {
            sum = zero_r;
            for i:max(ku + 2 - j, 1)..min(n + ku + 1 - j, kl + ku + 1) {
                sum = sum + abs(fget(ab, i, j));
            }
            if value < sum || _langb_isnan(sum, RT) {
                value = sum;
            }
        }
    } else if lsame(norm, #char "I") {
        for i:1..n {
            fset(work, i, zero_r);
        }
        for j:1..n {
            k = ku + 1 - j;
            for i:max(1, j-ku)..min(n, j+kl) {
                fset(work, i, fget(work.*, i) + abs(fget(ab, k + i, j)));
            }
        }
        value = zero_r;
        for i:1..n {
            temp = fget(work.*, i);
            if value < temp || _langb_isnan(temp, RT) {
                value = temp;
            }
        }
    } else {
        scale = zero_r;
        sum = one_r;
        for j:1..n {
            l = max(1, j-ku);
            k = ku + 1 - j + l;
            _langb_lassq(min(n, j+kl) - l + 1, _langb_matrix_column(ab, ldab, k, j, min(n, j+kl) - l + 1), 1, scale, sum, RT, MT);
        }
        value = scale*sqrt(sum);
    }

    return value;
}
