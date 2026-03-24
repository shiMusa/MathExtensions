// ? LAPACK *la_herpvgrw()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cla_herpvgrw/zla_herpvgrw
// *LA_HERPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U)
// for a Hermitian indefinite matrix.

cla_herpvgrw :: #bake_arguments _la_herpvgrw(RT=float32, MT=C32);
zla_herpvgrw :: #bake_arguments _la_herpvgrw(RT=float64, MT=C64);

#scope_file
_la_herpvgrw_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_la_herpvgrw :: (
    $$uplo: u8, $$n: int, $$info_val: int, a: $A/ArrayType2D, $$lda: int,
    af: $AF/ArrayType2D, $$ldaf: int, ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D,
    $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           WORK.data_type == RT;
} {
    i, j, k, kp, ncols : int = ---;
    amax, rpvgrw, tmp, umax : RT = ---;
    upper : bool = ---;

    upper = lsame(uplo, #char "U");
    if info_val == 0 {
        if upper {
            ncols = 1;
        } else {
            ncols = n;
        }
    } else {
        ncols = info_val;
    }

    rpvgrw = 1;
    for i:1..2*n {
        fset(work, i, 0);
    }

    if upper {
        for j:1..n {
            for i:1..j {
                fset(work, n + i, max(_la_herpvgrw_abs1(fget(a, i, j), RT, MT), fget(work.*, n + i)));
                fset(work, n + j, max(_la_herpvgrw_abs1(fget(a, i, j), RT, MT), fget(work.*, n + j)));
            }
        }
    } else {
        for j:1..n {
            for i:j..n {
                fset(work, n + i, max(_la_herpvgrw_abs1(fget(a, i, j), RT, MT), fget(work.*, n + i)));
                fset(work, n + j, max(_la_herpvgrw_abs1(fget(a, i, j), RT, MT), fget(work.*, n + j)));
            }
        }
    }

    if upper {
        k = n;
        // The LAPACK source uses `K .LT. NCOLS` here, which appears inconsistent
        // with the copied pivot walk; use the descending range implied by the algorithm.
        while k >= ncols && k > 0 {
            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    tmp = fget(work.*, n + k);
                    fset(work, n + k, fget(work.*, n + kp));
                    fset(work, n + kp, tmp);
                }
                for i:1..k {
                    fset(work, k, max(_la_herpvgrw_abs1(fget(af, i, k), RT, MT), fget(work.*, k)));
                }
                k -= 1;
            } else {
                kp = -fget(ipiv, k);
                tmp = fget(work.*, n + k - 1);
                fset(work, n + k - 1, fget(work.*, n + kp));
                fset(work, n + kp, tmp);
                for i:1..k - 1 {
                    fset(work, k, max(_la_herpvgrw_abs1(fget(af, i, k), RT, MT), fget(work.*, k)));
                    fset(work, k - 1, max(_la_herpvgrw_abs1(fget(af, i, k - 1), RT, MT), fget(work.*, k - 1)));
                }
                fset(work, k, max(_la_herpvgrw_abs1(fget(af, k, k), RT, MT), fget(work.*, k)));
                k -= 2;
            }
        }

        k = ncols;
        while k <= n {
            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    tmp = fget(work.*, n + k);
                    fset(work, n + k, fget(work.*, n + kp));
                    fset(work, n + kp, tmp);
                }
                k += 1;
            } else {
                kp = -fget(ipiv, k);
                tmp = fget(work.*, n + k);
                fset(work, n + k, fget(work.*, n + kp));
                fset(work, n + kp, tmp);
                k += 2;
            }
        }
    } else {
        k = 1;
        while k <= ncols {
            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    tmp = fget(work.*, n + k);
                    fset(work, n + k, fget(work.*, n + kp));
                    fset(work, n + kp, tmp);
                }
                for i:k..n {
                    fset(work, k, max(_la_herpvgrw_abs1(fget(af, i, k), RT, MT), fget(work.*, k)));
                }
                k += 1;
            } else {
                kp = -fget(ipiv, k);
                tmp = fget(work.*, n + k + 1);
                fset(work, n + k + 1, fget(work.*, n + kp));
                fset(work, n + kp, tmp);
                for i:k + 1..n {
                    fset(work, k, max(_la_herpvgrw_abs1(fget(af, i, k), RT, MT), fget(work.*, k)));
                    fset(work, k + 1, max(_la_herpvgrw_abs1(fget(af, i, k + 1), RT, MT), fget(work.*, k + 1)));
                }
                fset(work, k, max(_la_herpvgrw_abs1(fget(af, k, k), RT, MT), fget(work.*, k)));
                k += 2;
            }
        }

        k = ncols;
        while k >= 1 {
            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    tmp = fget(work.*, n + k);
                    fset(work, n + k, fget(work.*, n + kp));
                    fset(work, n + kp, tmp);
                }
                k -= 1;
            } else {
                kp = -fget(ipiv, k);
                tmp = fget(work.*, n + k);
                fset(work, n + k, fget(work.*, n + kp));
                fset(work, n + kp, tmp);
                k -= 2;
            }
        }
    }

    if upper {
        for i:ncols..n {
            umax = fget(work.*, i);
            amax = fget(work.*, n + i);
            if umax != 0 {
                rpvgrw = min(amax / umax, rpvgrw);
            }
        }
    } else {
        for i:1..ncols {
            umax = fget(work.*, i);
            amax = fget(work.*, n + i);
            if umax != 0 {
                rpvgrw = min(amax / umax, rpvgrw);
            }
        }
    }

    return rpvgrw;
}
