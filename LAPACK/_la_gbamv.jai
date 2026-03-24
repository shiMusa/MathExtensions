// ? LAPACK *la_gbamv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sla_gbamv/cla_gbamv/dla_gbamv/zla_gbamv
// *LA_GBAMV performs a matrix-vector operation to calculate error bounds.

sla_gbamv :: #bake_arguments _la_gbamv(RT=float32, MT=float32);
dla_gbamv :: #bake_arguments _la_gbamv(RT=float64, MT=float64);
cla_gbamv :: #bake_arguments _la_gbamv(RT=float32, MT=C32);
zla_gbamv :: #bake_arguments _la_gbamv(RT=float64, MT=C64);

#scope_file
_la_gbamv_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_la_gbamv_zero_value :: inline ($MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    return zero(MT);
}

_la_gbamv :: (
    $$trans: int, $$m: int, $$n: int, $$kl: int, $$ku: int, alpha: $RT,
    ab: $AB/ArrayType2D, $$ldab: int, x: $X/ArrayType1D, $$incx: int,
    beta: RT, y: *$Y/ArrayType1D, $$incy: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r :: #run one(RT);
    zero_t :: #run _la_gbamv_zero_value(MT);
    routine_name :: #run if MT == float32 {
        "SLA_GBAMV"
    } else if MT == float64 {
        "DLA_GBAMV"
    } else if MT == C32 {
        "CLA_GBAMV"
    } else {
        "ZLA_GBAMV"
    };

    symb_zero : bool = ---;
    temp, safe1 : RT = ---;
    i, info, iy, j, jx, kx, ky, lenx, leny, kd, ke : int = ---;

    info = 0;
    if !((trans == ilatrans(#char "N")) || (trans == ilatrans(#char "T")) || (trans == ilatrans(#char "C"))) {
        info = 1;
    } else if m < 0 {
        info = 2;
    } else if n < 0 {
        info = 3;
    } else if kl < 0 || kl > m - 1 {
        info = 4;
    } else if ku < 0 || ku > n - 1 {
        info = 5;
    } else if ldab < kl + ku + 1 {
        info = 6;
    } else if incx == 0 {
        info = 8;
    } else if incy == 0 {
        info = 11;
    }
    if info != 0 {
        xerbla(routine_name, info);
        return;
    }

    if (m == 0) || (n == 0) || ((alpha == zero_r) && (beta == one_r)) then return;

    if trans == ilatrans(#char "N") {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    if incx > 0 {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1)*incx;
    }
    if incy > 0 {
        ky = 1;
    } else {
        ky = 1 - (leny - 1)*incy;
    }

    safe1 = _la_gbamv_lamch_safe_minimum(RT);
    safe1 = ncast(RT, n + 1) * safe1;

    kd = ku + 1;
    ke = kl + 1;
    iy = ky;

    if incx == 1 {
        if trans == ilatrans(#char "N") {
            for i:1..leny {
                if beta == zero_r {
                    symb_zero = true;
                    fset(y, iy, zero_r);
                } else if fget(y.*, iy) == zero_r {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    fset(y, iy, beta * abs(fget(y.*, iy)));
                }

                if alpha != zero_r {
                    for j:max(i - kl, 1)..min(i + ku, lenx) {
                        temp = _la_gbamv_abs1(fget(ab, kd + i - j, j), RT, MT);
                        symb_zero = symb_zero && (fget(x, j) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_gbamv_abs1(fget(x, j), RT, MT) * temp);
                    }
                }

                if !symb_zero {
                    fset(y, iy, fget(y.*, iy) + sign(safe1, fget(y.*, iy)));
                }
                iy += incy;
            }
        } else {
            for i:1..leny {
                if beta == zero_r {
                    symb_zero = true;
                    fset(y, iy, zero_r);
                } else if fget(y.*, iy) == zero_r {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    fset(y, iy, beta * abs(fget(y.*, iy)));
                }

                if alpha != zero_r {
                    for j:max(i - kl, 1)..min(i + ku, lenx) {
                        temp = _la_gbamv_abs1(fget(ab, ke - i + j, i), RT, MT);
                        symb_zero = symb_zero && (fget(x, j) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_gbamv_abs1(fget(x, j), RT, MT) * temp);
                    }
                }

                if !symb_zero {
                    fset(y, iy, fget(y.*, iy) + sign(safe1, fget(y.*, iy)));
                }
                iy += incy;
            }
        }
    } else {
        if trans == ilatrans(#char "N") {
            for i:1..leny {
                if beta == zero_r {
                    symb_zero = true;
                    fset(y, iy, zero_r);
                } else if fget(y.*, iy) == zero_r {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    fset(y, iy, beta * abs(fget(y.*, iy)));
                }

                if alpha != zero_r {
                    jx = kx;
                    for j:max(i - kl, 1)..min(i + ku, lenx) {
                        temp = _la_gbamv_abs1(fget(ab, kd + i - j, j), RT, MT);
                        symb_zero = symb_zero && (fget(x, jx) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_gbamv_abs1(fget(x, jx), RT, MT) * temp);
                        jx += incx;
                    }
                }

                if !symb_zero {
                    fset(y, iy, fget(y.*, iy) + sign(safe1, fget(y.*, iy)));
                }
                iy += incy;
            }
        } else {
            for i:1..leny {
                if beta == zero_r {
                    symb_zero = true;
                    fset(y, iy, zero_r);
                } else if fget(y.*, iy) == zero_r {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    fset(y, iy, beta * abs(fget(y.*, iy)));
                }

                if alpha != zero_r {
                    jx = kx;
                    for j:max(i - kl, 1)..min(i + ku, lenx) {
                        temp = _la_gbamv_abs1(fget(ab, ke - i + j, i), RT, MT);
                        symb_zero = symb_zero && (fget(x, jx) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_gbamv_abs1(fget(x, jx), RT, MT) * temp);
                        jx += incx;
                    }
                }

                if !symb_zero {
                    fset(y, iy, fget(y.*, iy) + sign(safe1, fget(y.*, iy)));
                }
                iy += incy;
            }
        }
    }
}

_la_gbamv_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "S"); } else { return dlamch(#char "S"); }
}
