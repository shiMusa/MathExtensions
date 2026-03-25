// ? LAPACK *laqsy()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqsy/dlaqsy/claqsy/zlaqsy
// *LAQSY scales a symmetric/Hermitian matrix.

slaqsy :: #bake_arguments _laqsy(RT=float32, MT=float32);
dlaqsy :: #bake_arguments _laqsy(RT=float64, MT=float64);
claqsy :: #bake_arguments _laqsy(RT=float32, MT=C32);
zlaqsy :: #bake_arguments _laqsy(RT=float64, MT=C64);

#scope_file
_laqsy_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_laqsy_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_laqsy :: (
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    s: $S/ArrayType1D,
    $$scond: $RT,
    $$amax: RT,
    equed: *u8,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           S.data_type == RT;
} {
    one_r    :: #run one(RT);
    thresh_r :: #run ncast(RT, 0.1);

    i, j : int = ---;
    cj, large, small : RT = ---;

    if n <= 0 {
        equed.* = #char "N";
        return;
    }

    small = _laqsy_lamch_safe_minimum(RT) / _laqsy_lamch_precision(RT);
    large = one_r / small;

    if scond >= thresh_r && amax >= small && amax <= large {
        equed.* = #char "N";
    } else {
        if lsame(uplo, #char "U") {
            for j:1..n {
                cj = fget(s, j);
                for i:1..j {
                    fset(a, i, j, cj*fget(s, i)*fget(a.*, i, j));
                }
            }
        } else {
            for j:1..n {
                cj = fget(s, j);
                for i:j..n {
                    fset(a, i, j, cj*fget(s, i)*fget(a.*, i, j));
                }
            }
        }
        equed.* = #char "Y";
    }
}
