// ? LAPACK *laqge()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqge/dlaqge/claqge/zlaqge
// *LAQGE scales a general rectangular matrix, using row and column
// scaling factors computed by *GEEQU.

slaqge :: #bake_arguments _laqge(RT=float32, MT=float32);
dlaqge :: #bake_arguments _laqge(RT=float64, MT=float64);
claqge :: #bake_arguments _laqge(RT=float32, MT=C32);
zlaqge :: #bake_arguments _laqge(RT=float64, MT=C64);

#scope_file
_laqge_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_laqge_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_laqge :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
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
           A.data_type == MT &&
           R.data_type == RT &&
           C.data_type == RT;
} {
    one_r    :: #run one(RT);
    thresh_r :: #run ncast(RT, 0.1);

    i, j : int = ---;
    cj, large, small : RT = ---;

    if m <= 0 || n <= 0 {
        equed.* = #char "N";
        return;
    }

    small = _laqge_lamch_safe_minimum(RT) / _laqge_lamch_precision(RT);
    large = one_r / small;

    if rowcnd >= thresh_r && amax >= small && amax <= large {
        if colcnd >= thresh_r {
            equed.* = #char "N";
        } else {
            for j:1..n {
                cj = fget(c, j);
                for i:1..m {
                    fset(a, i, j, cj*fget(a.*, i, j));
                }
            }
            equed.* = #char "C";
        }
    } else if colcnd >= thresh_r {
        for j:1..n {
            for i:1..m {
                fset(a, i, j, fget(r, i)*fget(a.*, i, j));
            }
        }
        equed.* = #char "R";
    } else {
        for j:1..n {
            cj = fget(c, j);
            for i:1..m {
                fset(a, i, j, cj*fget(r, i)*fget(a.*, i, j));
            }
        }
        equed.* = #char "B";
    }
}
