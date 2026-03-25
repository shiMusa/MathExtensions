// ? LAPACK *laqhp()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment claqhp/zlaqhp
// *LAQHP scales a Hermitian matrix stored in packed form.

claqhp :: #bake_arguments _laqhp(RT=float32, MT=C32);
zlaqhp :: #bake_arguments _laqhp(RT=float64, MT=C64);

#scope_file
_laqhp_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_laqhp_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_laqhp_real_to_complex :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laqhp :: (
    $$uplo: u8,
    $$n: int,
    ap: *$AP/ArrayType1D,
    s: $S/ArrayType1D,
    $$scond: $RT,
    $$amax: RT,
    equed: *u8,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT &&
           S.data_type == RT;
} {
    one_r    :: #run one(RT);
    thresh_r :: #run ncast(RT, 0.1);

    i, j, jc : int = ---;
    cj, large, small : RT = ---;

    if n <= 0 {
        equed.* = #char "N";
        return;
    }

    small = _laqhp_lamch_safe_minimum(RT) / _laqhp_lamch_precision(RT);
    large = one_r / small;

    if scond >= thresh_r && amax >= small && amax <= large {
        equed.* = #char "N";
        return;
    }

    if lsame(uplo, #char "U") {
        jc = 1;
        for j:1..n {
            cj = fget(s, j);
            for i:1..j - 1 {
                fset(ap, jc + i - 1, cj*fget(s, i)*fget(ap.*, jc + i - 1));
            }
            fset(ap, jc + j - 1, _laqhp_real_to_complex(cj*cj*real(fget(ap.*, jc + j - 1)), RT, MT));
            jc = jc + j;
        }
    } else {
        jc = 1;
        for j:1..n {
            cj = fget(s, j);
            fset(ap, jc, _laqhp_real_to_complex(cj*cj*real(fget(ap.*, jc)), RT, MT));
            for i:j + 1..n {
                fset(ap, jc + i - j, cj*fget(s, i)*fget(ap.*, jc + i - j));
            }
            jc = jc + n - j + 1;
        }
    }

    equed.* = #char "Y";
}
