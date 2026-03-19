// ? LAPACK *getc2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgetc2/cgetc2/dgetc2/zgetc2
// *GETC2 computes the LU factorization with complete pivoting
// of a general n-by-n matrix.

sgetc2 :: #bake_arguments _getc2(RT=float32, MT=float32);
dgetc2 :: #bake_arguments _getc2(RT=float64, MT=float64);
cgetc2 :: #bake_arguments _getc2(RT=float32, MT=C32);
zgetc2 :: #bake_arguments _getc2(RT=float64, MT=C64);

#scope_file
_getc2_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_getc2_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_getc2_real_to_matrix :: ($x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == RT {
        return x;
    } else {
        return cpx(x, zero(RT));
    }
}

_getc2_swap_rows :: (
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

_getc2_swap_cols :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$col1: int,
    $$col2: int
) {
    temp : A.data_type = ---;
    for i:1..n {
        temp = fget(a.*, i, col1);
        fset(a, i, col1, fget(a.*, i, col2));
        fset(a, i, col2, temp);
    }
}

_getc2 :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    jpiv: *$JPIV/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           IPIV.data_type == int &&
           JPIV.data_type == int;
} {
    zero_r :: #run zero(RT);

    bignum, eps, smin, smlnum, xmax : RT = ---;
    i, ipv, j, jp, jpv : int = ---;

    info.* = 0;
    if n == 0 then return;

    eps = _getc2_lamch_precision(RT);
    smlnum = _getc2_lamch_safe_minimum(RT) / eps;
    bignum = one(RT) / smlnum;
    _ = bignum;

    if n == 1 {
        fset(ipiv, 1, 1);
        fset(jpiv, 1, 1);
        if abs(fget(a.*, 1, 1)) < smlnum {
            info.* = 1;
            fset(a, 1, 1, _getc2_real_to_matrix(smlnum, RT, MT));
        }
        return;
    }

    for i:1..n - 1 {
        xmax = zero_r;
        for jp:i..n {
            for ip:i..n {
                if abs(fget(a.*, ip, jp)) >= xmax {
                    xmax = abs(fget(a.*, ip, jp));
                    ipv = ip;
                    jpv = jp;
                }
            }
        }

        if i == 1 {
            smin = max(eps*xmax, smlnum);
        }

        if ipv != i {
            _getc2_swap_rows(n, a, ipv, i);
        }
        fset(ipiv, i, ipv);

        if jpv != i {
            _getc2_swap_cols(n, a, jpv, i);
        }
        fset(jpiv, i, jpv);

        if abs(fget(a.*, i, i)) < smin {
            info.* = i;
            fset(a, i, i, _getc2_real_to_matrix(smin, RT, MT));
        }

        for j:i + 1..n {
            fset(a, j, i, fget(a.*, j, i) / fget(a.*, i, i));
        }

        if i < n {
            for jp:i + 1..n {
                for j:i + 1..n {
                    fset(a, j, jp, fget(a.*, j, jp) - fget(a.*, j, i)*fget(a.*, i, jp));
                }
            }
        }
    }

    if abs(fget(a.*, n, n)) < smin {
        info.* = n;
        fset(a, n, n, _getc2_real_to_matrix(smin, RT, MT));
    }

    fset(ipiv, n, n);
    fset(jpiv, n, n);
}
