// ? LAPACK *laqsb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqsb/dlaqsb/claqsb/zlaqsb
// *LAQSB scales a symmetric/Hermitian band matrix.

slaqsb :: #bake_arguments _laqsb(RT=float32, MT=float32);
dlaqsb :: #bake_arguments _laqsb(RT=float64, MT=float64);
claqsb :: #bake_arguments _laqsb(RT=float32, MT=C32);
zlaqsb :: #bake_arguments _laqsb(RT=float64, MT=C64);

#scope_file
_laqsb_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_laqsb_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_laqsb :: (
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
           (MT == RT || MT == complex_with_inner_type(RT)) &&
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

    small = _laqsb_lamch_safe_minimum(RT) / _laqsb_lamch_precision(RT);
    large = one_r / small;

    if scond >= thresh_r && amax >= small && amax <= large {
        equed.* = #char "N";
    } else {
        if lsame(uplo, #char "U") {
            for j:1..n {
                cj = fget(s, j);
                for i:max(1, j - kd)..j {
                    fset(ab, kd + 1 + i - j, j, cj*fget(s, i)*fget(ab.*, kd + 1 + i - j, j));
                }
            }
        } else {
            for j:1..n {
                cj = fget(s, j);
                for i:j..min(n, j + kd) {
                    fset(ab, 1 + i - j, j, cj*fget(s, i)*fget(ab.*, 1 + i - j, j));
                }
            }
        }
        equed.* = #char "Y";
    }
}
