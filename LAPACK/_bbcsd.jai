// ? LAPACK *bbcsd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sbbcsd/cbbcsd/dbbcsd/zbbcsd
// *BBCSD computes the CS decomposition of a unitary or orthogonal
// matrix in bidiagonal-block form.

sbbcsd :: #bake_arguments _bbcsd(RT=float32, MT=float32);
dbbcsd :: #bake_arguments _bbcsd(RT=float64, MT=float64);
cbbcsd :: #bake_arguments _bbcsd(RT=float32, MT=C32);
zbbcsd :: #bake_arguments _bbcsd(RT=float64, MT=C64);

#scope_file
_bbcsd_clasr_right_forward_variable :: (
    $$m: int,
    $$n: int,
    rwork: $RWORK/ArrayType1D,
    $$c_index: int,
    $$s_index: int,
    a: *$A/ArrayType2D,
    $$row_start: int = 1,
    $$col_start: int = 1
) {
    RT :: #run get_number_inner_type(A.data_type);

    #if CHECKS {
        #run {
            assert(RWORK.data_type == RT);
        }
        assert(m >= 0);
        assert(n >= 0);
    }

    if m <= 0 || n <= 1 then return;

    temp, x, y : A.data_type = ---;
    ctemp, stemp : RT = ---;

    for j:1..n-1 {
        ctemp = fget(rwork, c_index + j - 1);
        stemp = fget(rwork, s_index + j - 1);
        for i:1..m {
            x = fget(a.*, row_start + i - 1, col_start + j - 1);
            y = fget(a.*, row_start + i - 1, col_start + j);
            temp = ctemp*x + stemp*y;
            fset(a, row_start + i - 1, col_start + j, ctemp*y - stemp*x);
            fset(a, row_start + i - 1, col_start + j - 1, temp);
        }
    }
}

_bbcsd_clasr_left_forward_variable :: (
    $$m: int,
    $$n: int,
    rwork: $RWORK/ArrayType1D,
    $$c_index: int,
    $$s_index: int,
    a: *$A/ArrayType2D,
    $$row_start: int = 1,
    $$col_start: int = 1
) {
    RT :: #run get_number_inner_type(A.data_type);

    #if CHECKS {
        #run {
            assert(RWORK.data_type == RT);
        }
        assert(m >= 0);
        assert(n >= 0);
    }

    if m <= 1 || n <= 0 then return;

    temp, x, y : A.data_type = ---;
    ctemp, stemp : RT = ---;

    for i:1..m-1 {
        ctemp = fget(rwork, c_index + i - 1);
        stemp = fget(rwork, s_index + i - 1);
        for j:1..n {
            x = fget(a.*, row_start + i - 1, col_start + j - 1);
            y = fget(a.*, row_start + i,     col_start + j - 1);
            temp = ctemp*x + stemp*y;
            fset(a, row_start + i,     col_start + j - 1, ctemp*y - stemp*x);
            fset(a, row_start + i - 1, col_start + j - 1, temp);
        }
    }
}

_bbcsd_scale_row :: (
    $$n: int,
    alpha: $T,
    a: *$A/ArrayType2D,
    $$row: int,
    $$col_start: int = 1
) {
    #if CHECKS {
        #run {
            assert(T == A.data_type);
        }
        assert(n >= 0);
    }

    if n <= 0 || alpha == one(T) then return;

    for j:1..n {
        fset(a, row, col_start + j - 1, alpha * fget(a.*, row, col_start + j - 1));
    }
}

_bbcsd_scale_col :: (
    $$n: int,
    alpha: $T,
    a: *$A/ArrayType2D,
    $$col: int,
    $$row_start: int = 1
) {
    #if CHECKS {
        #run {
            assert(T == A.data_type);
        }
        assert(n >= 0);
    }

    if n <= 0 || alpha == one(T) then return;

    for i:1..n {
        fset(a, row_start + i - 1, col, alpha * fget(a.*, row_start + i - 1, col));
    }
}

_bbcsd_swap_rows :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$row1: int,
    $$row2: int,
    $$col_start: int = 1
) {
    #if CHECKS {
        assert(n >= 0);
    }

    if n <= 0 || row1 == row2 then return;

    temp : A.data_type = ---;
    for j:1..n {
        temp = fget(a.*, row1, col_start + j - 1);
        fset(a, row1, col_start + j - 1, fget(a.*, row2, col_start + j - 1));
        fset(a, row2, col_start + j - 1, temp);
    }
}

_bbcsd_swap_cols :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$col1: int,
    $$col2: int,
    $$row_start: int = 1
) {
    #if CHECKS {
        assert(n >= 0);
    }

    if n <= 0 || col1 == col2 then return;

    temp : A.data_type = ---;
    for i:1..n {
        temp = fget(a.*, row_start + i - 1, col1);
        fset(a, row_start + i - 1, col1, fget(a.*, row_start + i - 1, col2));
        fset(a, row_start + i - 1, col2, temp);
    }
}

_bbcsd_slas2 :: (f: $RT, g: RT, h: RT, ssmin: *RT, ssmax: *RT) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float64 {
        dlas2(f, g, h, ssmin, ssmax);
    } else {
        slas2(f, g, h, ssmin, ssmax);
    }
}

_bbcsd_slartgs_store :: (
    x: $RT,
    y: RT,
    sigma: RT,
    rwork: *$RWORK/ArrayType1D,
    $$cs_index: int,
    $$sn_index: int
) #modify {
    return RT == float32 || RT == float64;
} {
    #if CHECKS {
        #run assert(RWORK.data_type == RT);
    }

    cs, sn : RT = ---;
    #if RT == float64 {
        dlartgs(x, y, sigma, *cs, *sn);
    } else {
        slartgs(x, y, sigma, *cs, *sn);
    }
    fset(rwork, cs_index, cs);
    fset(rwork, sn_index, sn);
}

_bbcsd_slartgp_store :: (
    f: $RT,
    g: RT,
    rwork: *$RWORK/ArrayType1D,
    $$sn_index: int,
    $$cs_index: int,
    r: *RT
) #modify {
    return RT == float32 || RT == float64;
} {
    #if CHECKS {
        #run assert(RWORK.data_type == RT);
    }

    sn, cs : RT = ---;
    #if RT == float64 {
        dlartgp(f, g, *sn, *cs, r);
    } else {
        slartgp(f, g, *sn, *cs, r);
    }
    fset(rwork, sn_index, sn);
    fset(rwork, cs_index, cs);
}

_bbcsd :: (
    $$jobu1: u8,
    $$jobu2: u8,
    $$jobv1t: u8,
    $$jobv2t: u8,
    $$trans: u8,
    $$m: int,
    $$p: int,
    $$q: int,
    theta: *$Theta/ArrayType1D,
    phi: *$Phi/ArrayType1D,
    u1: *$U1/ArrayType2D,
    $$ldu1: int,
    u2: *$U2/ArrayType2D,
    $$ldu2: int,
    v1t: *$V1T/ArrayType2D,
    $$ldv1t: int,
    v2t: *$V2T/ArrayType2D,
    $$ldv2t: int,
    b11d: *$B11D/ArrayType1D,
    b11e: *$B11E/ArrayType1D,
    b12d: *$B12D/ArrayType1D,
    b12e: *$B12E/ArrayType1D,
    b21d: *$B21D/ArrayType1D,
    b21e: *$B21E/ArrayType1D,
    b22d: *$B22D/ArrayType1D,
    b22e: *$B22E/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    $$lrwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if CHECKS {
        #run {
            assert(Theta.data_type == RT);
            assert(Phi.data_type   == RT);
            assert(U1.data_type    == MT);
            assert(U2.data_type    == MT);
            assert(V1T.data_type   == MT);
            assert(V2T.data_type   == MT);
            assert(B11D.data_type  == RT);
            assert(B11E.data_type  == RT);
            assert(B12D.data_type  == RT);
            assert(B12E.data_type  == RT);
            assert(B21D.data_type  == RT);
            assert(B21E.data_type  == RT);
            assert(B22D.data_type  == RT);
            assert(B22E.data_type  == RT);
            assert(RWORK.data_type == RT);
        }
    }

    maxitr :: 6;
    hundred : RT : 100.0;
    meighth : RT : -0.125;
    one_r   : RT : 1.0;
    ten     : RT : 10.0;
    zero_r  : RT : 0.0;
    negone  : MT : ncast(MT, -1);
    piover2 : RT : 1.57079632679489661923132169163975144210;
    routine_name :: #run if MT == float32 {
        "SBBCSD"
    } else if MT == float64 {
        "DBBCSD"
    } else if MT == C32 {
        "CBBCSD"
    } else {
        "ZBBCSD"
    };

    colmajor  := !char_eq(trans, #char "t");
    lquery    := lrwork == -1;
    wantu1    := char_eq(jobu1,  #char "y");
    wantu2    := char_eq(jobu2,  #char "y");
    wantv1t   := char_eq(jobv1t, #char "y");
    wantv2t   := char_eq(jobv2t, #char "y");

    i, imin, imax, iter, iu1cs, iu1sn, iu2cs, iu2sn, iv1tcs, iv1tsn, iv2tcs, iv2tsn, j : int = ---;
    lrworkmin, lrworkopt, maxit, mini : int = ---;
    restart11, restart12, restart21, restart22 : bool = ---;
    b11bulge, b12bulge, b21bulge, b22bulge, dummy : RT = ---;
    eps, mu, nu, r, sigma11, sigma21 : RT = ---;
    temp, thetamax, thetamin, thresh, tol, tolmul, unfl, x1, x2, y1, y2 : RT = ---;

    info.* = 0;
    if m < 0 {
        info.* = -6;
    } else if p < 0 || p > m {
        info.* = -7;
    } else if q < 0 || q > m {
        info.* = -8;
    } else if q > p || q > m-p || q > m-q {
        info.* = -8;
    } else if wantu1 && ldu1 < p {
        info.* = -12;
    } else if wantu2 && ldu2 < m-p {
        info.* = -14;
    } else if wantv1t && ldv1t < q {
        info.* = -16;
    } else if wantv2t && ldv2t < m-q {
        info.* = -18;
    }

    if info.* == 0 && q == 0 {
        lrworkmin = 1;
        fset(rwork, 1, cast(RT) lrworkmin);
        return;
    }

    if info.* == 0 {
        iu1cs = 1;
        iu1sn = iu1cs + q;
        iu2cs = iu1sn + q;
        iu2sn = iu2cs + q;
        iv1tcs = iu2sn + q;
        iv1tsn = iv1tcs + q;
        iv2tcs = iv1tsn + q;
        iv2tsn = iv2tcs + q;
        lrworkopt = iv2tsn + q - 1;
        lrworkmin = lrworkopt;
        fset(rwork, 1, cast(RT) lrworkopt);
        if lrwork < lrworkmin && !lquery {
            info.* = -28;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    #if RT == float64 {
        eps = EPS_64;
        unfl = SAFMIN_64;
    } else {
        eps = EPS_32;
        unfl = SAFMIN_32;
    }
    tolmul = max(ten, min(hundred, eps**meighth));
    tol = tolmul * eps;
    thresh = max(tol, cast(RT) (maxitr*q*q) * unfl);

    for i:1..q {
        if fget(theta.*, i) < thresh {
            fset(theta, i, zero_r);
        } else if fget(theta.*, i) > piover2 - thresh {
            fset(theta, i, piover2);
        }
    }
    if q > 1 {
        for i:1..q-1 {
            if fget(phi.*, i) < thresh {
                fset(phi, i, zero_r);
            } else if fget(phi.*, i) > piover2 - thresh {
                fset(phi, i, piover2);
            }
        }
    }

    imax = q;
    while imax > 1 {
        if fget(phi.*, imax-1) != zero_r then break;
        imax -= 1;
    }
    imin = imax - 1;
    if imin > 1 {
        while fget(phi.*, imin-1) != zero_r {
            imin -= 1;
            if imin <= 1 then break;
        }
    }

    maxit = maxitr*q*q;
    iter = 0;

    while imax > 1 {
        fset(b11d, imin, cos(fget(theta.*, imin)));
        fset(b21d, imin, -sin(fget(theta.*, imin)));
        for i:imin..imax-1 {
            fset(b11e, i, -sin(fget(theta.*, i)) * sin(fget(phi.*, i)));
            fset(b11d, i+1, cos(fget(theta.*, i+1)) * cos(fget(phi.*, i)));
            fset(b12d, i, sin(fget(theta.*, i)) * cos(fget(phi.*, i)));
            fset(b12e, i, cos(fget(theta.*, i+1)) * sin(fget(phi.*, i)));
            fset(b21e, i, -cos(fget(theta.*, i)) * sin(fget(phi.*, i)));
            fset(b21d, i+1, -sin(fget(theta.*, i+1)) * cos(fget(phi.*, i)));
            fset(b22d, i, cos(fget(theta.*, i)) * cos(fget(phi.*, i)));
            fset(b22e, i, -sin(fget(theta.*, i+1)) * sin(fget(phi.*, i)));
        }
        fset(b12d, imax, sin(fget(theta.*, imax)));
        fset(b22d, imax, cos(fget(theta.*, imax)));

        if iter > maxit {
            info.* = 0;
            if q > 1 {
                for i:1..q-1 {
                    if fget(phi.*, i) != zero_r {
                        info.* += 1;
                    }
                }
            }
            return;
        }

        iter += imax - imin;

        thetamax = fget(theta.*, imin);
        thetamin = fget(theta.*, imin);
        for i:imin+1..imax {
            if fget(theta.*, i) > thetamax then thetamax = fget(theta.*, i);
            if fget(theta.*, i) < thetamin then thetamin = fget(theta.*, i);
        }

        if thetamax > piover2 - thresh {
            mu = zero_r;
            nu = one_r;
        } else if thetamin < thresh {
            mu = one_r;
            nu = zero_r;
        } else {
            _bbcsd_slas2(fget(b11d.*, imax-1), fget(b11e.*, imax-1), fget(b11d.*, imax), *sigma11, *dummy);
            _bbcsd_slas2(fget(b21d.*, imax-1), fget(b21e.*, imax-1), fget(b21d.*, imax), *sigma21, *dummy);

            if sigma11 <= sigma21 {
                mu = sigma11;
                nu = sqrt(one_r - mu**2);
                if mu < thresh {
                    mu = zero_r;
                    nu = one_r;
                }
            } else {
                nu = sigma21;
                mu = sqrt(one_r - nu**2);
                if nu < thresh {
                    mu = one_r;
                    nu = zero_r;
                }
            }
        }

        if mu <= nu {
            _bbcsd_slartgs_store(fget(b11d.*, imin), fget(b11e.*, imin), mu, rwork, iv1tcs+imin-1, iv1tsn+imin-1);
        } else {
            _bbcsd_slartgs_store(fget(b21d.*, imin), fget(b21e.*, imin), nu, rwork, iv1tcs+imin-1, iv1tsn+imin-1);
        }

        temp = fget(rwork.*, iv1tcs+imin-1)*fget(b11d.*, imin) +
               fget(rwork.*, iv1tsn+imin-1)*fget(b11e.*, imin);
        fset(b11e, imin, fget(rwork.*, iv1tcs+imin-1)*fget(b11e.*, imin) -
                         fget(rwork.*, iv1tsn+imin-1)*fget(b11d.*, imin));
        fset(b11d, imin, temp);
        b11bulge = fget(rwork.*, iv1tsn+imin-1)*fget(b11d.*, imin+1);
        fset(b11d, imin+1, fget(rwork.*, iv1tcs+imin-1)*fget(b11d.*, imin+1));
        temp = fget(rwork.*, iv1tcs+imin-1)*fget(b21d.*, imin) +
               fget(rwork.*, iv1tsn+imin-1)*fget(b21e.*, imin);
        fset(b21e, imin, fget(rwork.*, iv1tcs+imin-1)*fget(b21e.*, imin) -
                         fget(rwork.*, iv1tsn+imin-1)*fget(b21d.*, imin));
        fset(b21d, imin, temp);
        b21bulge = fget(rwork.*, iv1tsn+imin-1)*fget(b21d.*, imin+1);
        fset(b21d, imin+1, fget(rwork.*, iv1tcs+imin-1)*fget(b21d.*, imin+1));

        fset(theta, imin, atan2(sqrt(fget(b21d.*, imin)**2 + b21bulge**2),
                                sqrt(fget(b11d.*, imin)**2 + b11bulge**2)));

        if fget(b11d.*, imin)**2 + b11bulge**2 > thresh**2 {
            _bbcsd_slartgp_store(b11bulge, fget(b11d.*, imin), rwork, iu1sn+imin-1, iu1cs+imin-1, *r);
        } else if mu <= nu {
            _bbcsd_slartgs_store(fget(b11e.*, imin), fget(b11d.*, imin+1), mu, rwork, iu1cs+imin-1, iu1sn+imin-1);
        } else {
            _bbcsd_slartgs_store(fget(b12d.*, imin), fget(b12e.*, imin), nu, rwork, iu1cs+imin-1, iu1sn+imin-1);
        }
        if fget(b21d.*, imin)**2 + b21bulge**2 > thresh**2 {
            _bbcsd_slartgp_store(b21bulge, fget(b21d.*, imin), rwork, iu2sn+imin-1, iu2cs+imin-1, *r);
        } else if nu < mu {
            _bbcsd_slartgs_store(fget(b21e.*, imin), fget(b21d.*, imin+1), nu, rwork, iu2cs+imin-1, iu2sn+imin-1);
        } else {
            _bbcsd_slartgs_store(fget(b22d.*, imin), fget(b22e.*, imin), mu, rwork, iu2cs+imin-1, iu2sn+imin-1);
        }
        fset(rwork, iu2cs+imin-1, -fget(rwork.*, iu2cs+imin-1));
        fset(rwork, iu2sn+imin-1, -fget(rwork.*, iu2sn+imin-1));

        temp = fget(rwork.*, iu1cs+imin-1)*fget(b11e.*, imin) +
               fget(rwork.*, iu1sn+imin-1)*fget(b11d.*, imin+1);
        fset(b11d, imin+1, fget(rwork.*, iu1cs+imin-1)*fget(b11d.*, imin+1) -
                           fget(rwork.*, iu1sn+imin-1)*fget(b11e.*, imin));
        fset(b11e, imin, temp);
        if imax > imin+1 {
            b11bulge = fget(rwork.*, iu1sn+imin-1)*fget(b11e.*, imin+1);
            fset(b11e, imin+1, fget(rwork.*, iu1cs+imin-1)*fget(b11e.*, imin+1));
        }
        temp = fget(rwork.*, iu1cs+imin-1)*fget(b12d.*, imin) +
               fget(rwork.*, iu1sn+imin-1)*fget(b12e.*, imin);
        fset(b12e, imin, fget(rwork.*, iu1cs+imin-1)*fget(b12e.*, imin) -
                         fget(rwork.*, iu1sn+imin-1)*fget(b12d.*, imin));
        fset(b12d, imin, temp);
        b12bulge = fget(rwork.*, iu1sn+imin-1)*fget(b12d.*, imin+1);
        fset(b12d, imin+1, fget(rwork.*, iu1cs+imin-1)*fget(b12d.*, imin+1));
        temp = fget(rwork.*, iu2cs+imin-1)*fget(b21e.*, imin) +
               fget(rwork.*, iu2sn+imin-1)*fget(b21d.*, imin+1);
        fset(b21d, imin+1, fget(rwork.*, iu2cs+imin-1)*fget(b21d.*, imin+1) -
                           fget(rwork.*, iu2sn+imin-1)*fget(b21e.*, imin));
        fset(b21e, imin, temp);
        if imax > imin+1 {
            b21bulge = fget(rwork.*, iu2sn+imin-1)*fget(b21e.*, imin+1);
            fset(b21e, imin+1, fget(rwork.*, iu2cs+imin-1)*fget(b21e.*, imin+1));
        }
        temp = fget(rwork.*, iu2cs+imin-1)*fget(b22d.*, imin) +
               fget(rwork.*, iu2sn+imin-1)*fget(b22e.*, imin);
        fset(b22e, imin, fget(rwork.*, iu2cs+imin-1)*fget(b22e.*, imin) -
                         fget(rwork.*, iu2sn+imin-1)*fget(b22d.*, imin));
        fset(b22d, imin, temp);
        b22bulge = fget(rwork.*, iu2sn+imin-1)*fget(b22d.*, imin+1);
        fset(b22d, imin+1, fget(rwork.*, iu2cs+imin-1)*fget(b22d.*, imin+1));

        for i:imin+1..imax-1 {
            x1 = sin(fget(theta.*, i-1))*fget(b11e.*, i-1) + cos(fget(theta.*, i-1))*fget(b21e.*, i-1);
            x2 = sin(fget(theta.*, i-1))*b11bulge + cos(fget(theta.*, i-1))*b21bulge;
            y1 = sin(fget(theta.*, i-1))*fget(b12d.*, i-1) + cos(fget(theta.*, i-1))*fget(b22d.*, i-1);
            y2 = sin(fget(theta.*, i-1))*b12bulge + cos(fget(theta.*, i-1))*b22bulge;

            fset(phi, i-1, atan2(sqrt(x1**2 + x2**2), sqrt(y1**2 + y2**2)));

            restart11 = fget(b11e.*, i-1)**2 + b11bulge**2 <= thresh**2;
            restart21 = fget(b21e.*, i-1)**2 + b21bulge**2 <= thresh**2;
            restart12 = fget(b12d.*, i-1)**2 + b12bulge**2 <= thresh**2;
            restart22 = fget(b22d.*, i-1)**2 + b22bulge**2 <= thresh**2;

            if !restart11 && !restart21 {
                _bbcsd_slartgp_store(x2, x1, rwork, iv1tsn+i-1, iv1tcs+i-1, *r);
            } else if !restart11 && restart21 {
                _bbcsd_slartgp_store(b11bulge, fget(b11e.*, i-1), rwork, iv1tsn+i-1, iv1tcs+i-1, *r);
            } else if restart11 && !restart21 {
                _bbcsd_slartgp_store(b21bulge, fget(b21e.*, i-1), rwork, iv1tsn+i-1, iv1tcs+i-1, *r);
            } else if mu <= nu {
                _bbcsd_slartgs_store(fget(b11d.*, i), fget(b11e.*, i), mu, rwork, iv1tcs+i-1, iv1tsn+i-1);
            } else {
                _bbcsd_slartgs_store(fget(b21d.*, i), fget(b21e.*, i), nu, rwork, iv1tcs+i-1, iv1tsn+i-1);
            }
            fset(rwork, iv1tcs+i-1, -fget(rwork.*, iv1tcs+i-1));
            fset(rwork, iv1tsn+i-1, -fget(rwork.*, iv1tsn+i-1));

            if !restart12 && !restart22 {
                _bbcsd_slartgp_store(y2, y1, rwork, iv2tsn+i-2, iv2tcs+i-2, *r);
            } else if !restart12 && restart22 {
                _bbcsd_slartgp_store(b12bulge, fget(b12d.*, i-1), rwork, iv2tsn+i-2, iv2tcs+i-2, *r);
            } else if restart12 && !restart22 {
                _bbcsd_slartgp_store(b22bulge, fget(b22d.*, i-1), rwork, iv2tsn+i-2, iv2tcs+i-2, *r);
            } else if nu < mu {
                _bbcsd_slartgs_store(fget(b12e.*, i-1), fget(b12d.*, i), nu, rwork, iv2tcs+i-2, iv2tsn+i-2);
            } else {
                _bbcsd_slartgs_store(fget(b22e.*, i-1), fget(b22d.*, i), mu, rwork, iv2tcs+i-2, iv2tsn+i-2);
            }

            temp = fget(rwork.*, iv1tcs+i-1)*fget(b11d.*, i) + fget(rwork.*, iv1tsn+i-1)*fget(b11e.*, i);
            fset(b11e, i, fget(rwork.*, iv1tcs+i-1)*fget(b11e.*, i) - fget(rwork.*, iv1tsn+i-1)*fget(b11d.*, i));
            fset(b11d, i, temp);
            b11bulge = fget(rwork.*, iv1tsn+i-1)*fget(b11d.*, i+1);
            fset(b11d, i+1, fget(rwork.*, iv1tcs+i-1)*fget(b11d.*, i+1));
            temp = fget(rwork.*, iv1tcs+i-1)*fget(b21d.*, i) + fget(rwork.*, iv1tsn+i-1)*fget(b21e.*, i);
            fset(b21e, i, fget(rwork.*, iv1tcs+i-1)*fget(b21e.*, i) - fget(rwork.*, iv1tsn+i-1)*fget(b21d.*, i));
            fset(b21d, i, temp);
            b21bulge = fget(rwork.*, iv1tsn+i-1)*fget(b21d.*, i+1);
            fset(b21d, i+1, fget(rwork.*, iv1tcs+i-1)*fget(b21d.*, i+1));
            temp = fget(rwork.*, iv2tcs+i-2)*fget(b12e.*, i-1) +
                    fget(rwork.*, iv2tsn+i-2)*fget(b12d.*, i);
            fset(b12d, i, fget(rwork.*, iv2tcs+i-2)*fget(b12d.*, i) -
                         fget(rwork.*, iv2tsn+i-2)*fget(b12e.*, i-1));
            fset(b12e, i-1, temp);
            b12bulge = fget(rwork.*, iv2tsn+i-2)*fget(b12e.*, i);
            fset(b12e, i, fget(rwork.*, iv2tcs+i-2)*fget(b12e.*, i));
            temp = fget(rwork.*, iv2tcs+i-2)*fget(b22e.*, i-1) +
                    fget(rwork.*, iv2tsn+i-2)*fget(b22d.*, i);
            fset(b22d, i, fget(rwork.*, iv2tcs+i-2)*fget(b22d.*, i) -
                         fget(rwork.*, iv2tsn+i-2)*fget(b22e.*, i-1));
            fset(b22e, i-1, temp);
            b22bulge = fget(rwork.*, iv2tsn+i-2)*fget(b22e.*, i);
            fset(b22e, i, fget(rwork.*, iv2tcs+i-2)*fget(b22e.*, i));

            x1 = cos(fget(phi.*, i-1))*fget(b11d.*, i) + sin(fget(phi.*, i-1))*fget(b12e.*, i-1);
            x2 = cos(fget(phi.*, i-1))*b11bulge + sin(fget(phi.*, i-1))*b12bulge;
            y1 = cos(fget(phi.*, i-1))*fget(b21d.*, i) + sin(fget(phi.*, i-1))*fget(b22e.*, i-1);
            y2 = cos(fget(phi.*, i-1))*b21bulge + sin(fget(phi.*, i-1))*b22bulge;

            fset(theta, i, atan2(sqrt(y1**2 + y2**2), sqrt(x1**2 + x2**2)));

            restart11 = fget(b11d.*, i)**2 + b11bulge**2 <= thresh**2;
            restart12 = fget(b12e.*, i-1)**2 + b12bulge**2 <= thresh**2;
            restart21 = fget(b21d.*, i)**2 + b21bulge**2 <= thresh**2;
            restart22 = fget(b22e.*, i-1)**2 + b22bulge**2 <= thresh**2;

            if !restart11 && !restart12 {
                _bbcsd_slartgp_store(x2, x1, rwork, iu1sn+i-1, iu1cs+i-1, *r);
            } else if !restart11 && restart12 {
                _bbcsd_slartgp_store(b11bulge, fget(b11d.*, i), rwork, iu1sn+i-1, iu1cs+i-1, *r);
            } else if restart11 && !restart12 {
                _bbcsd_slartgp_store(b12bulge, fget(b12e.*, i-1), rwork, iu1sn+i-1, iu1cs+i-1, *r);
            } else if mu <= nu {
                _bbcsd_slartgs_store(fget(b11e.*, i), fget(b11d.*, i+1), mu, rwork, iu1cs+i-1, iu1sn+i-1);
            } else {
                _bbcsd_slartgs_store(fget(b12d.*, i), fget(b12e.*, i), nu, rwork, iu1cs+i-1, iu1sn+i-1);
            }
            if !restart21 && !restart22 {
                _bbcsd_slartgp_store(y2, y1, rwork, iu2sn+i-1, iu2cs+i-1, *r);
            } else if !restart21 && restart22 {
                _bbcsd_slartgp_store(b21bulge, fget(b21d.*, i), rwork, iu2sn+i-1, iu2cs+i-1, *r);
            } else if restart21 && !restart22 {
                _bbcsd_slartgp_store(b22bulge, fget(b22e.*, i-1), rwork, iu2sn+i-1, iu2cs+i-1, *r);
            } else if nu < mu {
                _bbcsd_slartgs_store(fget(b21e.*, i), fget(b21e.*, i+1), nu, rwork, iu2cs+i-1, iu2sn+i-1);
            } else {
                _bbcsd_slartgs_store(fget(b22d.*, i), fget(b22e.*, i), mu, rwork, iu2cs+i-1, iu2sn+i-1);
            }
            fset(rwork, iu2cs+i-1, -fget(rwork.*, iu2cs+i-1));
            fset(rwork, iu2sn+i-1, -fget(rwork.*, iu2sn+i-1));

            temp = fget(rwork.*, iu1cs+i-1)*fget(b11e.*, i) + fget(rwork.*, iu1sn+i-1)*fget(b11d.*, i+1);
            fset(b11d, i+1, fget(rwork.*, iu1cs+i-1)*fget(b11d.*, i+1) - fget(rwork.*, iu1sn+i-1)*fget(b11e.*, i));
            fset(b11e, i, temp);
            if i < imax - 1 {
                b11bulge = fget(rwork.*, iu1sn+i-1)*fget(b11e.*, i+1);
                fset(b11e, i+1, fget(rwork.*, iu1cs+i-1)*fget(b11e.*, i+1));
            }
            temp = fget(rwork.*, iu2cs+i-1)*fget(b21e.*, i) + fget(rwork.*, iu2sn+i-1)*fget(b21d.*, i+1);
            fset(b21d, i+1, fget(rwork.*, iu2cs+i-1)*fget(b21d.*, i+1) - fget(rwork.*, iu2sn+i-1)*fget(b21e.*, i));
            fset(b21e, i, temp);
            if i < imax - 1 {
                b21bulge = fget(rwork.*, iu2sn+i-1)*fget(b21e.*, i+1);
                fset(b21e, i+1, fget(rwork.*, iu2cs+i-1)*fget(b21e.*, i+1));
            }
            temp = fget(rwork.*, iu1cs+i-1)*fget(b12d.*, i) + fget(rwork.*, iu1sn+i-1)*fget(b12e.*, i);
            fset(b12e, i, fget(rwork.*, iu1cs+i-1)*fget(b12e.*, i) - fget(rwork.*, iu1sn+i-1)*fget(b12d.*, i));
            fset(b12d, i, temp);
            b12bulge = fget(rwork.*, iu1sn+i-1)*fget(b12d.*, i+1);
            fset(b12d, i+1, fget(rwork.*, iu1cs+i-1)*fget(b12d.*, i+1));
            temp = fget(rwork.*, iu2cs+i-1)*fget(b22d.*, i) + fget(rwork.*, iu2sn+i-1)*fget(b22e.*, i);
            fset(b22e, i, fget(rwork.*, iu2cs+i-1)*fget(b22e.*, i) - fget(rwork.*, iu2sn+i-1)*fget(b22d.*, i));
            fset(b22d, i, temp);
            b22bulge = fget(rwork.*, iu2sn+i-1)*fget(b22d.*, i+1);
            fset(b22d, i+1, fget(rwork.*, iu2cs+i-1)*fget(b22d.*, i+1));
        }

        x1 = sin(fget(theta.*, imax-1))*fget(b11e.*, imax-1) +
             cos(fget(theta.*, imax-1))*fget(b21e.*, imax-1);
        y1 = sin(fget(theta.*, imax-1))*fget(b12d.*, imax-1) +
             cos(fget(theta.*, imax-1))*fget(b22d.*, imax-1);
        y2 = sin(fget(theta.*, imax-1))*b12bulge + cos(fget(theta.*, imax-1))*b22bulge;

        fset(phi, imax-1, atan2(abs(x1), sqrt(y1**2 + y2**2)));

        restart12 = fget(b12d.*, imax-1)**2 + b12bulge**2 <= thresh**2;
        restart22 = fget(b22d.*, imax-1)**2 + b22bulge**2 <= thresh**2;

        if !restart12 && !restart22 {
            _bbcsd_slartgp_store(y2, y1, rwork, iv2tsn+imax-2, iv2tcs+imax-2, *r);
        } else if !restart12 && restart22 {
            _bbcsd_slartgp_store(b12bulge, fget(b12d.*, imax-1), rwork, iv2tsn+imax-2, iv2tcs+imax-2, *r);
        } else if restart12 && !restart22 {
            _bbcsd_slartgp_store(b22bulge, fget(b22d.*, imax-1), rwork, iv2tsn+imax-2, iv2tcs+imax-2, *r);
        } else if nu < mu {
            _bbcsd_slartgs_store(fget(b12e.*, imax-1), fget(b12d.*, imax), nu, rwork, iv2tcs+imax-2, iv2tsn+imax-2);
        } else {
            _bbcsd_slartgs_store(fget(b22e.*, imax-1), fget(b22d.*, imax), mu, rwork, iv2tcs+imax-2, iv2tsn+imax-2);
        }

        temp = fget(rwork.*, iv2tcs+imax-2)*fget(b12e.*, imax-1) + fget(rwork.*, iv2tsn+imax-2)*fget(b12d.*, imax);
        fset(b12d, imax, fget(rwork.*, iv2tcs+imax-2)*fget(b12d.*, imax) - fget(rwork.*, iv2tsn+imax-2)*fget(b12e.*, imax-1));
        fset(b12e, imax-1, temp);
        temp = fget(rwork.*, iv2tcs+imax-2)*fget(b22e.*, imax-1) + fget(rwork.*, iv2tsn+imax-2)*fget(b22d.*, imax);
        fset(b22d, imax, fget(rwork.*, iv2tcs+imax-2)*fget(b22d.*, imax) - fget(rwork.*, iv2tsn+imax-2)*fget(b22e.*, imax-1));
        fset(b22e, imax-1, temp);

        if wantu1 {
            if colmajor {
                _bbcsd_clasr_right_forward_variable(p, imax-imin+1, rwork.*, iu1cs+imin-1, iu1sn+imin-1, u1, 1, imin);
            } else {
                _bbcsd_clasr_left_forward_variable(imax-imin+1, p, rwork.*, iu1cs+imin-1, iu1sn+imin-1, u1, imin, 1);
            }
        }
        if wantu2 {
            if colmajor {
                _bbcsd_clasr_right_forward_variable(m-p, imax-imin+1, rwork.*, iu2cs+imin-1, iu2sn+imin-1, u2, 1, imin);
            } else {
                _bbcsd_clasr_left_forward_variable(imax-imin+1, m-p, rwork.*, iu2cs+imin-1, iu2sn+imin-1, u2, imin, 1);
            }
        }
        if wantv1t {
            if colmajor {
                _bbcsd_clasr_left_forward_variable(imax-imin+1, q, rwork.*, iv1tcs+imin-1, iv1tsn+imin-1, v1t, imin, 1);
            } else {
                _bbcsd_clasr_right_forward_variable(q, imax-imin+1, rwork.*, iv1tcs+imin-1, iv1tsn+imin-1, v1t, 1, imin);
            }
        }
        if wantv2t {
            if colmajor {
                _bbcsd_clasr_left_forward_variable(imax-imin+1, m-q, rwork.*, iv2tcs+imin-1, iv2tsn+imin-1, v2t, imin, 1);
            } else {
                _bbcsd_clasr_right_forward_variable(m-q, imax-imin+1, rwork.*, iv2tcs+imin-1, iv2tsn+imin-1, v2t, 1, imin);
            }
        }

        if fget(b11e.*, imax-1) + fget(b21e.*, imax-1) > zero_r {
            fset(b11d, imax, -fget(b11d.*, imax));
            fset(b21d, imax, -fget(b21d.*, imax));
            if wantv1t {
                if colmajor {
                    _bbcsd_scale_row(q, negone, v1t, imax, 1);
                } else {
                    _bbcsd_scale_col(q, negone, v1t, imax, 1);
                }
            }
        }

        x1 = cos(fget(phi.*, imax-1))*fget(b11d.*, imax) +
             sin(fget(phi.*, imax-1))*fget(b12e.*, imax-1);
        y1 = cos(fget(phi.*, imax-1))*fget(b21d.*, imax) +
             sin(fget(phi.*, imax-1))*fget(b22e.*, imax-1);

        fset(theta, imax, atan2(abs(y1), abs(x1)));

        if fget(b11d.*, imax) + fget(b12e.*, imax-1) < zero_r {
            fset(b12d, imax, -fget(b12d.*, imax));
            if wantu1 {
                if colmajor {
                    _bbcsd_scale_col(p, negone, u1, imax, 1);
                } else {
                    _bbcsd_scale_row(p, negone, u1, imax, 1);
                }
            }
        }
        if fget(b21d.*, imax) + fget(b22e.*, imax-1) > zero_r {
            fset(b22d, imax, -fget(b22d.*, imax));
            if wantu2 {
                if colmajor {
                    _bbcsd_scale_col(m-p, negone, u2, imax, 1);
                } else {
                    _bbcsd_scale_row(m-p, negone, u2, imax, 1);
                }
            }
        }

        if fget(b12d.*, imax) + fget(b22d.*, imax) < zero_r {
            if wantv2t {
                if colmajor {
                    _bbcsd_scale_row(m-q, negone, v2t, imax, 1);
                } else {
                    _bbcsd_scale_col(m-q, negone, v2t, imax, 1);
                }
            }
        }

        for i:imin..imax {
            if fget(theta.*, i) < thresh {
                fset(theta, i, zero_r);
            } else if fget(theta.*, i) > piover2 - thresh {
                fset(theta, i, piover2);
            }
        }
        if imax > imin {
            for i:imin..imax-1 {
                if fget(phi.*, i) < thresh {
                    fset(phi, i, zero_r);
                } else if fget(phi.*, i) > piover2 - thresh {
                    fset(phi, i, piover2);
                }
            }
        }

        if imax > 1 {
            while fget(phi.*, imax-1) == zero_r {
                imax -= 1;
                if imax <= 1 then break;
            }
        }
        if imin > imax - 1 then imin = imax - 1;
        if imin > 1 {
            while fget(phi.*, imin-1) != zero_r {
                imin -= 1;
                if imin <= 1 then break;
            }
        }
    }

    for i:1..q {
        mini = i;
        thetamin = fget(theta.*, i);
        for j:i+1..q {
            if fget(theta.*, j) < thetamin {
                mini = j;
                thetamin = fget(theta.*, j);
            }
        }

        if mini != i {
            fset(theta, mini, fget(theta.*, i));
            fset(theta, i, thetamin);
            if colmajor {
                if wantu1 then _bbcsd_swap_cols(p, u1, i, mini, 1);
                if wantu2 then _bbcsd_swap_cols(m-p, u2, i, mini, 1);
                if wantv1t then _bbcsd_swap_rows(q, v1t, i, mini, 1);
                if wantv2t then _bbcsd_swap_rows(m-q, v2t, i, mini, 1);
            } else {
                if wantu1 then _bbcsd_swap_rows(p, u1, i, mini, 1);
                if wantu2 then _bbcsd_swap_rows(m-p, u2, i, mini, 1);
                if wantv1t then _bbcsd_swap_cols(q, v1t, i, mini, 1);
                if wantv2t then _bbcsd_swap_cols(m-q, v2t, i, mini, 1);
            }
        }
    }
}
