// ? LAPACK *la_syamv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sla_syamv/dla_syamv/cla_syamv/zla_syamv
// *LA_SYAMV performs a matrix-vector operation to calculate error bounds.

sla_syamv :: #bake_arguments _la_syamv(RT=float32, MT=float32);
dla_syamv :: #bake_arguments _la_syamv(RT=float64, MT=float64);
cla_syamv :: #bake_arguments _la_syamv(RT=float32, MT=C32);
zla_syamv :: #bake_arguments _la_syamv(RT=float64, MT=C64);

#scope_file
_la_syamv_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_la_syamv_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "S"); } else { return dlamch(#char "S"); }
}

_la_syamv :: (
    $$uplo: int, $$n: int, alpha: $RT,
    a: $A/ArrayType2D, $$lda: int, x: $X/ArrayType1D, $$incx: int,
    beta: RT, y: *$Y/ArrayType1D, $$incy: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SLA_SYAMV"
    } else if MT == float64 {
        "DLA_SYAMV"
    } else if MT == C32 {
        "CLA_SYAMV"
    } else {
        "ZLA_SYAMV"
    };

    symb_zero : bool = ---;
    temp, safe1 : RT = ---;
    i, info, iy, j, jx, kx, ky : int = ---;

    info = 0;
    if uplo != ilauplo(#char "U") && uplo != ilauplo(#char "L") {
        info = 1;
    } else if n < 0 {
        info = 2;
    } else if lda < max(1, n) {
        info = 5;
    } else if incx == 0 {
        info = 7;
    } else if incy == 0 {
        info = 10;
    }
    if info != 0 {
        xerbla(routine_name, info);
        return;
    }

    if (n == 0) || ((alpha == zero_r) && (beta == one_r)) then return;

    if incx > 0 { kx = 1; } else { kx = 1 - (n - 1)*incx; }
    if incy > 0 { ky = 1; } else { ky = 1 - (n - 1)*incy; }

    safe1 = _la_syamv_lamch_safe_minimum(RT);
    safe1 = ncast(RT, n + 1) * safe1;

    iy = ky;
    if incx == 1 {
        if uplo == ilauplo(#char "U") {
            for i:1..n {
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
                    for j:1..i {
                        temp = _la_syamv_abs1(fget(a, j, i), RT, MT);
                        symb_zero = symb_zero && (fget(x, j) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_syamv_abs1(fget(x, j), RT, MT) * temp);
                    }
                    for j:i + 1..n {
                        temp = _la_syamv_abs1(fget(a, i, j), RT, MT);
                        symb_zero = symb_zero && (fget(x, j) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_syamv_abs1(fget(x, j), RT, MT) * temp);
                    }
                }

                if !symb_zero {
                    fset(y, iy, fget(y.*, iy) + sign(safe1, fget(y.*, iy)));
                }
                iy += incy;
            }
        } else {
            for i:1..n {
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
                    for j:1..i {
                        temp = _la_syamv_abs1(fget(a, i, j), RT, MT);
                        symb_zero = symb_zero && (fget(x, j) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_syamv_abs1(fget(x, j), RT, MT) * temp);
                    }
                    for j:i + 1..n {
                        temp = _la_syamv_abs1(fget(a, j, i), RT, MT);
                        symb_zero = symb_zero && (fget(x, j) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_syamv_abs1(fget(x, j), RT, MT) * temp);
                    }
                }

                if !symb_zero {
                    fset(y, iy, fget(y.*, iy) + sign(safe1, fget(y.*, iy)));
                }
                iy += incy;
            }
        }
    } else {
        if uplo == ilauplo(#char "U") {
            for i:1..n {
                if beta == zero_r {
                    symb_zero = true;
                    fset(y, iy, zero_r);
                } else if fget(y.*, iy) == zero_r {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    fset(y, iy, beta * abs(fget(y.*, iy)));
                }

                jx = kx;
                if alpha != zero_r {
                    for j:1..i {
                        temp = _la_syamv_abs1(fget(a, j, i), RT, MT);
                        symb_zero = symb_zero && (fget(x, jx) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_syamv_abs1(fget(x, jx), RT, MT) * temp);
                        jx += incx;
                    }
                    for j:i + 1..n {
                        temp = _la_syamv_abs1(fget(a, i, j), RT, MT);
                        symb_zero = symb_zero && (fget(x, jx) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_syamv_abs1(fget(x, jx), RT, MT) * temp);
                        jx += incx;
                    }
                }

                if !symb_zero {
                    fset(y, iy, fget(y.*, iy) + sign(safe1, fget(y.*, iy)));
                }
                iy += incy;
            }
        } else {
            for i:1..n {
                if beta == zero_r {
                    symb_zero = true;
                    fset(y, iy, zero_r);
                } else if fget(y.*, iy) == zero_r {
                    symb_zero = true;
                } else {
                    symb_zero = false;
                    fset(y, iy, beta * abs(fget(y.*, iy)));
                }

                jx = kx;
                if alpha != zero_r {
                    for j:1..i {
                        temp = _la_syamv_abs1(fget(a, i, j), RT, MT);
                        symb_zero = symb_zero && (fget(x, jx) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_syamv_abs1(fget(x, jx), RT, MT) * temp);
                        jx += incx;
                    }
                    for j:i + 1..n {
                        temp = _la_syamv_abs1(fget(a, j, i), RT, MT);
                        symb_zero = symb_zero && (fget(x, jx) == zero_t || temp == zero_r);
                        fset(y, iy, fget(y.*, iy) + alpha * _la_syamv_abs1(fget(x, jx), RT, MT) * temp);
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
