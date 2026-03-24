// ? LAPACK *lagtm()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slagtm/clagtm/dlagtm/zlagtm
// *LAGTM performs a matrix-matrix product with a tridiagonal matrix.

slagtm :: #bake_arguments _lagtm(RT=float32, MT=float32);
dlagtm :: #bake_arguments _lagtm(RT=float64, MT=float64);
clagtm :: #bake_arguments _lagtm(RT=float32, MT=C32);
zlagtm :: #bake_arguments _lagtm(RT=float64, MT=C64);

#scope_file
_lagtm_zero :: inline ($RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == RT {
        return zero(RT);
    } else {
        return cpx(zero(RT), zero(RT));
    }
}

_lagtm_conj_if_needed :: inline ($$trans: u8, x: $MT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == RT {
        _ = RT;
        return x;
    } else {
        if lsame(trans, #char "C") {
            return conjugate(x);
        } else {
            return x;
        }
    }
}

_lagtm :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    $$alpha: $RT,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    x: $X/ArrayType2D,
    $$ldx: int,
    $$beta: RT,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           X.data_type == MT &&
           B.data_type == MT;
} {
    zero_t :: #run _lagtm_zero(RT, MT);
    one_r  :: #run one(RT);

    if n == 0 then return;

    if beta == zero(RT) {
        for j:1..nrhs {
            for i:1..n {
                fset(b, i, j, zero_t);
            }
        }
    } else if beta == -one_r {
        for j:1..nrhs {
            for i:1..n {
                fset(b, i, j, -fget(b.*, i, j));
            }
        }
    }

    if alpha == one_r {
        if lsame(trans, #char "N") {
            for j:1..nrhs {
                if n == 1 {
                    fset(b, 1, j, fget(b.*, 1, j) + fget(d, 1)*fget(x, 1, j));
                } else {
                    fset(b, 1, j, fget(b.*, 1, j) + fget(d, 1)*fget(x, 1, j) + fget(du, 1)*fget(x, 2, j));
                    fset(b, n, j, fget(b.*, n, j) + fget(dl, n - 1)*fget(x, n - 1, j) + fget(d, n)*fget(x, n, j));
                    for i:2..n - 1 {
                        fset(b, i, j, fget(b.*, i, j) + fget(dl, i - 1)*fget(x, i - 1, j) +
                                      fget(d, i)*fget(x, i, j) + fget(du, i)*fget(x, i + 1, j));
                    }
                }
            }
        } else {
            for j:1..nrhs {
                if n == 1 {
                    fset(b, 1, j, fget(b.*, 1, j) + _lagtm_conj_if_needed(trans, fget(d, 1), RT, MT)*fget(x, 1, j));
                } else {
                    fset(b, 1, j, fget(b.*, 1, j) +
                                  _lagtm_conj_if_needed(trans, fget(d, 1), RT, MT)*fget(x, 1, j) +
                                  _lagtm_conj_if_needed(trans, fget(dl, 1), RT, MT)*fget(x, 2, j));
                    fset(b, n, j, fget(b.*, n, j) +
                                  _lagtm_conj_if_needed(trans, fget(du, n - 1), RT, MT)*fget(x, n - 1, j) +
                                  _lagtm_conj_if_needed(trans, fget(d, n), RT, MT)*fget(x, n, j));
                    for i:2..n - 1 {
                        fset(b, i, j, fget(b.*, i, j) +
                                      _lagtm_conj_if_needed(trans, fget(du, i - 1), RT, MT)*fget(x, i - 1, j) +
                                      _lagtm_conj_if_needed(trans, fget(d, i), RT, MT)*fget(x, i, j) +
                                      _lagtm_conj_if_needed(trans, fget(dl, i), RT, MT)*fget(x, i + 1, j));
                    }
                }
            }
        }
    } else if alpha == -one_r {
        if lsame(trans, #char "N") {
            for j:1..nrhs {
                if n == 1 {
                    fset(b, 1, j, fget(b.*, 1, j) - fget(d, 1)*fget(x, 1, j));
                } else {
                    fset(b, 1, j, fget(b.*, 1, j) - fget(d, 1)*fget(x, 1, j) - fget(du, 1)*fget(x, 2, j));
                    fset(b, n, j, fget(b.*, n, j) - fget(dl, n - 1)*fget(x, n - 1, j) - fget(d, n)*fget(x, n, j));
                    for i:2..n - 1 {
                        fset(b, i, j, fget(b.*, i, j) - fget(dl, i - 1)*fget(x, i - 1, j) -
                                      fget(d, i)*fget(x, i, j) - fget(du, i)*fget(x, i + 1, j));
                    }
                }
            }
        } else {
            for j:1..nrhs {
                if n == 1 {
                    fset(b, 1, j, fget(b.*, 1, j) - _lagtm_conj_if_needed(trans, fget(d, 1), RT, MT)*fget(x, 1, j));
                } else {
                    fset(b, 1, j, fget(b.*, 1, j) -
                                  _lagtm_conj_if_needed(trans, fget(d, 1), RT, MT)*fget(x, 1, j) -
                                  _lagtm_conj_if_needed(trans, fget(dl, 1), RT, MT)*fget(x, 2, j));
                    fset(b, n, j, fget(b.*, n, j) -
                                  _lagtm_conj_if_needed(trans, fget(du, n - 1), RT, MT)*fget(x, n - 1, j) -
                                  _lagtm_conj_if_needed(trans, fget(d, n), RT, MT)*fget(x, n, j));
                    for i:2..n - 1 {
                        fset(b, i, j, fget(b.*, i, j) -
                                      _lagtm_conj_if_needed(trans, fget(du, i - 1), RT, MT)*fget(x, i - 1, j) -
                                      _lagtm_conj_if_needed(trans, fget(d, i), RT, MT)*fget(x, i, j) -
                                      _lagtm_conj_if_needed(trans, fget(dl, i), RT, MT)*fget(x, i + 1, j));
                    }
                }
            }
        }
    }
}
