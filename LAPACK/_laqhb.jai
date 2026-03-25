// ? LAPACK *laqhb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment claqhb/zlaqhb
// *LAQHB scales a Hermitian band matrix, using scaling factors
// computed by *PBEQU.

claqhb :: #bake_arguments _laqhb(RT=float32, MT=C32);
zlaqhb :: #bake_arguments _laqhb(RT=float64, MT=C64);

#scope_file
_laqhb_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_laqhb_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_laqhb_real_to_complex :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laqhb :: (
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    s: $S/ArrayType1D,
    $$scond: $RT,
    $$amax: RT,
    equed: *u8,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
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

    small = _laqhb_lamch_safe_minimum(RT) / _laqhb_lamch_precision(RT);
    large = one_r / small;

    if scond >= thresh_r && amax >= small && amax <= large {
        equed.* = #char "N";
        return;
    }

    if lsame(uplo, #char "U") {
        for j:1..n {
            cj = fget(s, j);
            for i:max(1, j - kd)..j - 1 {
                fset(ab, kd + 1 + i - j, j, cj*fget(s, i)*fget(ab.*, kd + 1 + i - j, j));
            }
            fset(ab, kd + 1, j, _laqhb_real_to_complex(cj*cj*real(fget(ab.*, kd + 1, j)), RT, MT));
        }
    } else {
        for j:1..n {
            cj = fget(s, j);
            fset(ab, 1, j, _laqhb_real_to_complex(cj*cj*real(fget(ab.*, 1, j)), RT, MT));
            for i:j + 1..min(n, j + kd) {
                fset(ab, 1 + i - j, j, cj*fget(s, i)*fget(ab.*, 1 + i - j, j));
            }
        }
    }

    equed.* = #char "Y";
}
