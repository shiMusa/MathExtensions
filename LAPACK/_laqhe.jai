// ? LAPACK *laqhe()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment claqhe/zlaqhe
// *LAQHE scales a Hermitian matrix.

claqhe :: #bake_arguments _laqhe(RT=float32, MT=C32);
zlaqhe :: #bake_arguments _laqhe(RT=float64, MT=C64);

#scope_file
_laqhe_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_laqhe_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_laqhe_real_to_complex :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laqhe :: (
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
           MT == complex_with_inner_type(RT) &&
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

    small = _laqhe_lamch_safe_minimum(RT) / _laqhe_lamch_precision(RT);
    large = one_r / small;

    if scond >= thresh_r && amax >= small && amax <= large {
        equed.* = #char "N";
        return;
    }

    if lsame(uplo, #char "U") {
        for j:1..n {
            cj = fget(s, j);
            for i:1..j - 1 {
                fset(a, i, j, cj*fget(s, i)*fget(a.*, i, j));
            }
            fset(a, j, j, _laqhe_real_to_complex(cj*cj*real(fget(a.*, j, j)), RT, MT));
        }
    } else {
        for j:1..n {
            cj = fget(s, j);
            fset(a, j, j, _laqhe_real_to_complex(cj*cj*real(fget(a.*, j, j)), RT, MT));
            for i:j + 1..n {
                fset(a, i, j, cj*fget(s, i)*fget(a.*, i, j));
            }
        }
    }

    equed.* = #char "Y";
}
