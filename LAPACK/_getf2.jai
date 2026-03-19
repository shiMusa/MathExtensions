// ? LAPACK *getf2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgetf2/cgetf2/dgetf2/zgetf2
// *GETF2 computes the LU factorization of a general m-by-n matrix
// using partial pivoting with row interchanges (unblocked algorithm).

sgetf2 :: #bake_arguments _getf2(RT=float32, MT=float32);
dgetf2 :: #bake_arguments _getf2(RT=float64, MT=float64);
cgetf2 :: #bake_arguments _getf2(RT=float32, MT=C32);
zgetf2 :: #bake_arguments _getf2(RT=float64, MT=C64);

#scope_file
_getf2_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_getf2_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $MT: Type
) -> int #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return isamax(n, x, 1) + 1;
    } else if MT == float64 {
        return idamax(n, x, 1) + 1;
    } else if MT == C32 {
        return icamax(n, x, 1) + 1;
    } else {
        return izamax(n, x, 1) + 1;
    }
}

_getf2_scal_real :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        sscal(n, alpha, x, 1);
    } else {
        dscal(n, alpha, x, 1);
    }
}

_getf2_rscl_complex :: (
    $$n: int,
    alpha: $MT,
    x: *$X/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        crscl(n, alpha, x, 1);
    } else {
        zrscl(n, alpha, x, 1);
    }
}

_getf2_swap_rows :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$row1: int,
    $$row2: int
) {
    temp : A.data_type = ---;
    for j:1..n {
        temp = fget(a.*, row1, j);
        fset(a, row1, j, fget(a.*, row2, j));
        fset(a, row2, j, temp);
    }
}

_getf2_matrix_segment :: inline (a: $A/ArrayType2D, $$ld: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    #if CHECKS {
        assert(ld >= 1);
        assert(row >= 1);
        assert(col >= 1);
        assert(count >= 0);
    }

    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld + (row - 1);
    view.count = count;
    return view;
}

_getf2 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SGETF2"
    } else if MT == float64 {
        "DGETF2"
    } else if MT == C32 {
        "CGETF2"
    } else {
        "ZGETF2"
    };

    i, j, jp : int = ---;
    minmn := min(m, n);
    sfmin : RT = ---;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m == 0 || n == 0 then return;

    #if MT == RT {
        sfmin = _getf2_lamch_safe_minimum(RT);
    }

    for j:1..minmn {
        col_view := _getf2_matrix_segment(a.*, lda, j, j, m - j + 1);
        jp = j - 1 + _getf2_iamax_1based(m - j + 1, col_view, MT);
        fset(ipiv, j, jp);

        if fget(a.*, jp, j) != zero_t {
            if jp != j {
                _getf2_swap_rows(n, a, j, jp);
            }

            if j < m {
                multipliers_view := _getf2_matrix_segment(a.*, lda, j + 1, j, m - j);
                #if MT == RT {
                    if abs(fget(a.*, j, j)) >= sfmin {
                        _getf2_scal_real(m - j, one(RT) / fget(a.*, j, j), *multipliers_view, RT);
                    } else {
                        for i:1..m - j {
                            fset(a, j + i, j, fget(a.*, j + i, j) / fget(a.*, j, j));
                        }
                    }
                } else {
                    _getf2_rscl_complex(m - j, fget(a.*, j, j), *multipliers_view, RT, MT);
                }
            }
        } else if info.* == 0 {
            info.* = j;
        }

        if j < minmn {
            for jp:j + 1..n {
                for i:j + 1..m {
                    fset(a, i, jp, fget(a.*, i, jp) - fget(a.*, i, j)*fget(a.*, j, jp));
                }
            }
        }
    }
}
