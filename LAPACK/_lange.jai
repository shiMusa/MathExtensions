// ? LAPACK *lange()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slange/clange/dlange/zlange
// *LANGE returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// general rectangular matrix.

slange :: #bake_arguments _lange(RT=float32, MT=float32);
dlange :: #bake_arguments _lange(RT=float64, MT=float64);
clange :: #bake_arguments _lange(RT=float32, MT=C32);
zlange :: #bake_arguments _lange(RT=float64, MT=C64);

#scope_file
_lange_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lange_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _lange_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lange_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lange_lassq :: (
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

_lange :: (
    $$norm: u8,
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
    scale, sum, temp, value : RT = ---;

    if min(m, n) == 0 {
        value = zero_r;
    } else if lsame(norm, #char "M") {
        value = zero_r;
        for j:1..n {
            for i:1..m {
                temp = abs(fget(a, i, j));
                if value < temp || _lange_isnan(temp, RT) {
                    value = temp;
                }
            }
        }
    } else if lsame(norm, #char "O") || norm == #char "1" {
        value = zero_r;
        for j:1..n {
            sum = zero_r;
            for i:1..m {
                sum = sum + abs(fget(a, i, j));
            }
            if value < sum || _lange_isnan(sum, RT) {
                value = sum;
            }
        }
    } else if lsame(norm, #char "I") {
        for i:1..m {
            fset(work, i, zero_r);
        }
        for j:1..n {
            for i:1..m {
                fset(work, i, fget(work.*, i) + abs(fget(a, i, j)));
            }
        }
        value = zero_r;
        for i:1..m {
            temp = fget(work.*, i);
            if value < temp || _lange_isnan(temp, RT) {
                value = temp;
            }
        }
    } else {
        scale = zero_r;
        sum = one_r;
        for j:1..n {
            _lange_lassq(m, _lange_matrix_column(a, lda, 1, j, m), 1, scale, sum, RT, MT);
        }
        value = scale*sqrt(sum);
    }

    return value;
}
