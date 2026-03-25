// ? LAPACK *laqsp()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqsp/dlaqsp/claqsp/zlaqsp
// *LAQSP scales a symmetric/Hermitian matrix in packed storage.

slaqsp :: #bake_arguments _laqsp(RT=float32, MT=float32);
dlaqsp :: #bake_arguments _laqsp(RT=float64, MT=float64);
claqsp :: #bake_arguments _laqsp(RT=float32, MT=C32);
zlaqsp :: #bake_arguments _laqsp(RT=float64, MT=C64);

#scope_file
_laqsp_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_laqsp_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_laqsp :: (
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
           (MT == RT || MT == complex_with_inner_type(RT)) &&
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

    small = _laqsp_lamch_safe_minimum(RT) / _laqsp_lamch_precision(RT);
    large = one_r / small;

    if scond >= thresh_r && amax >= small && amax <= large {
        equed.* = #char "N";
    } else {
        jc = 1;
        if lsame(uplo, #char "U") {
            for j:1..n {
                cj = fget(s, j);
                for i:1..j {
                    fset(ap, jc + i - 1, cj*fget(s, i)*fget(ap.*, jc + i - 1));
                }
                jc = jc + j;
            }
        } else {
            for j:1..n {
                cj = fget(s, j);
                for i:j..n {
                    fset(ap, jc + i - j, cj*fget(s, i)*fget(ap.*, jc + i - j));
                }
                jc = jc + n - j + 1;
            }
        }
        equed.* = #char "Y";
    }
}
