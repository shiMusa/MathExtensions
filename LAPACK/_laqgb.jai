// ? LAPACK *laqgb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqgb/dlaqgb/claqgb/zlaqgb
// *LAQGB scales a general band matrix, using row and column scaling
// factors computed by *GBEQU.

slaqgb :: #bake_arguments _laqgb(RT=float32, MT=float32);
dlaqgb :: #bake_arguments _laqgb(RT=float64, MT=float64);
claqgb :: #bake_arguments _laqgb(RT=float32, MT=C32);
zlaqgb :: #bake_arguments _laqgb(RT=float64, MT=C64);

#scope_file
_laqgb_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_laqgb_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_laqgb :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    r: $R/ArrayType1D,
    c: $C/ArrayType1D,
    $$rowcnd: $RT,
    $$colcnd: RT,
    $$amax: RT,
    equed: *u8,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           R.data_type == RT &&
           C.data_type == RT;
} {
    one_r    :: #run one(RT);
    thresh_r :: #run ncast(RT, 0.1);

    i, j : int = ---;
    cj, large, small : RT = ---;
    row1, row2 : int = ---;

    if m <= 0 || n <= 0 {
        equed.* = #char "N";
        return;
    }

    small = _laqgb_lamch_safe_minimum(RT) / _laqgb_lamch_precision(RT);
    large = one_r / small;

    if rowcnd >= thresh_r && amax >= small && amax <= large {
        if colcnd >= thresh_r {
            equed.* = #char "N";
        } else {
            for j:1..n {
                cj = fget(c, j);
                row1 = max(1, j - ku);
                row2 = min(m, j + kl);
                for i:row1..row2 {
                    fset(ab, ku + 1 + i - j, j, cj*fget(ab.*, ku + 1 + i - j, j));
                }
            }
            equed.* = #char "C";
        }
    } else if colcnd >= thresh_r {
        for j:1..n {
            row1 = max(1, j - ku);
            row2 = min(m, j + kl);
            for i:row1..row2 {
                fset(ab, ku + 1 + i - j, j, fget(r, i)*fget(ab.*, ku + 1 + i - j, j));
            }
        }
        equed.* = #char "R";
    } else {
        for j:1..n {
            cj = fget(c, j);
            row1 = max(1, j - ku);
            row2 = min(m, j + kl);
            for i:row1..row2 {
                fset(ab, ku + 1 + i - j, j, cj*fget(r, i)*fget(ab.*, ku + 1 + i - j, j));
            }
        }
        equed.* = #char "B";
    }
}
