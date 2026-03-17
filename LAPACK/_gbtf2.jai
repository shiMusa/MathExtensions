// ? LAPACK *gbtf2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbtf2/cgbtf2/dgbtf2/zgbtf2
// *GBTF2 computes the LU factorization of a general band matrix
// using the unblocked version of the algorithm.

sgbtf2 :: #bake_arguments _gbtf2(MT=float32);
dgbtf2 :: #bake_arguments _gbtf2(MT=float64);
cgbtf2 :: #bake_arguments _gbtf2(MT=C32);
zgbtf2 :: #bake_arguments _gbtf2(MT=C64);

#scope_file
_gbtf2_matrix_segment :: inline (a: $A/ArrayType2D, $$ld: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
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

_gbtf2_iamax_1based :: (
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

_gbtf2_scal :: (
    $$n: int,
    alpha: $MT,
    x: *$X/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, 1);
    } else if MT == float64 {
        dscal(n, alpha, x, 1);
    } else if MT == C32 {
        cscal(n, alpha, x, 1);
    } else {
        zscal(n, alpha, x, 1);
    }
}

_gbtf2_swap_band_rows :: (
    ab: *$AB/ArrayType2D,
    $$kv: int,
    $$j: int,
    $$jp: int,
    $$ju: int
) {
    row1 := kv + jp;
    row2 := kv + 1;

    for col:j..ju {
        temp := fget(ab.*, row1, col);
        fset(ab, row1, col, fget(ab.*, row2, col));
        fset(ab, row2, col, temp);
        row1 -= 1;
        row2 -= 1;
    }
}

_gbtf2_rank1_update :: (
    $$km: int,
    $$j: int,
    $$ju: int,
    ab: *$AB/ArrayType2D,
    $$kv: int,
    alpha: $MT
) #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    zero_t :: #run zero(MT);

    if km <= 0 || ju <= j then return;

    temp : MT = ---;
    offset, row, col : int = ---;

    for col:j+1..ju {
        offset = col - (j + 1);
        temp = alpha * fget(ab.*, kv - offset, col);
        if temp != zero_t {
            for i:1..km {
                row = kv + i - offset;
                fset(ab, row, col, fget(ab.*, row, col) + fget(ab.*, kv + 1 + i, j) * temp);
            }
        }
    }
}

_gbtf2 :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT &&
           IPIV.data_type == int;
} {
    one_t  :: #run one(MT);
    zero_t :: #run zero(MT);
    minus_one_t :: #run -one(MT);
    routine_name :: #run if MT == float32 {
        "SGBTF2"
    } else if MT == float64 {
        "DGBTF2"
    } else if MT == C32 {
        "CGBTF2"
    } else {
        "ZGBTF2"
    };

    jp, ju, km, kv : int = ---;
    minmn :: min(m, n);

    kv = ku + kl;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if kl < 0 {
        info.* = -3;
    } else if ku < 0 {
        info.* = -4;
    } else if ldab < kl + kv + 1 {
        info.* = -6;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m == 0 || n == 0 then return;

    for j:ku+2..min(kv, n) {
        for i:kv-j+2..kl {
            fset(ab, i, j, zero_t);
        }
    }

    ju = 1;

    for j:1..minmn {
        if j + kv <= n {
            for i:1..kl {
                fset(ab, i, j + kv, zero_t);
            }
        }

        km = min(kl, m - j);
        column_view := _gbtf2_matrix_segment(ab.*, ldab, kv + 1, j, km + 1);
        jp = _gbtf2_iamax_1based(km + 1, column_view, MT);
        fset(ipiv, j, jp + j - 1);

        if fget(ab.*, kv + jp, j) != zero_t {
            ju = max(ju, min(j + ku + jp - 1, n));

            if jp != 1 {
                _gbtf2_swap_band_rows(ab, kv, j, jp, ju);
            }

            if km > 0 {
                multipliers_view := _gbtf2_matrix_segment(ab.*, ldab, kv + 2, j, km);
                _gbtf2_scal(km, one_t / fget(ab.*, kv + 1, j), *multipliers_view, MT);

                if ju > j {
                    _gbtf2_rank1_update(km, j, ju, ab, kv, minus_one_t);
                }
            }
        } else {
            if info.* == 0 {
                info.* = j;
            }
        }
    }
}
