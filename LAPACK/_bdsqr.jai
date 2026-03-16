// ? LAPACK *bdsqr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sbdsqr/cbdsqr/dbdsqr/zbdsqr
// *BDSQR computes the singular values and, optionally, the right and/or
// left singular vectors from the singular value decomposition (SVD)
// of a real N-by-N (upper or lower) bidiagonal matrix.

sbdsqr :: #bake_arguments _bdsqr(RT=float32, MT=float32);
dbdsqr :: #bake_arguments _bdsqr(RT=float64, MT=float64);
cbdsqr :: #bake_arguments _bdsqr(RT=float32, MT=C32);
zbdsqr :: #bake_arguments _bdsqr(RT=float64, MT=C64);

#scope_file
_bdsqr_lasr_right_forward_variable :: (
    $$m: int,
    $$n: int,
    work: $WORK/ArrayType1D,
    $$c_index: int,
    $$s_index: int,
    a: *$A/ArrayType2D,
    $$row_start: int = 1,
    $$col_start: int = 1
) {
    RT :: #run get_number_inner_type(A.data_type);

    #if CHECKS {
        #run {
            assert(WORK.data_type == RT);
        }
        assert(m >= 0);
        assert(n >= 0);
    }

    if m <= 0 || n <= 1 then return;

    temp, x, y : A.data_type = ---;
    ctemp, stemp : RT = ---;

    for j:1..n-1 {
        ctemp = fget(work, c_index + j - 1);
        stemp = fget(work, s_index + j - 1);
        for i:1..m {
            x = fget(a.*, row_start + i - 1, col_start + j - 1);
            y = fget(a.*, row_start + i - 1, col_start + j);
            temp = ctemp*x + stemp*y;
            fset(a, row_start + i - 1, col_start + j, ctemp*y - stemp*x);
            fset(a, row_start + i - 1, col_start + j - 1, temp);
        }
    }
}

_bdsqr_lasr_left_forward_variable :: (
    $$m: int,
    $$n: int,
    work: $WORK/ArrayType1D,
    $$c_index: int,
    $$s_index: int,
    a: *$A/ArrayType2D,
    $$row_start: int = 1,
    $$col_start: int = 1
) {
    RT :: #run get_number_inner_type(A.data_type);

    #if CHECKS {
        #run {
            assert(WORK.data_type == RT);
        }
        assert(m >= 0);
        assert(n >= 0);
    }

    if m <= 1 || n <= 0 then return;

    temp, x, y : A.data_type = ---;
    ctemp, stemp : RT = ---;

    for i:1..m-1 {
        ctemp = fget(work, c_index + i - 1);
        stemp = fget(work, s_index + i - 1);
        for j:1..n {
            x = fget(a.*, row_start + i - 1, col_start + j - 1);
            y = fget(a.*, row_start + i,     col_start + j - 1);
            temp = ctemp*x + stemp*y;
            fset(a, row_start + i,     col_start + j - 1, ctemp*y - stemp*x);
            fset(a, row_start + i - 1, col_start + j - 1, temp);
        }
    }
}

_bdsqr_lasr_right_backward_variable :: (
    $$m: int,
    $$n: int,
    work: $WORK/ArrayType1D,
    $$c_index: int,
    $$s_index: int,
    a: *$A/ArrayType2D,
    $$row_start: int = 1,
    $$col_start: int = 1
) {
    RT :: #run get_number_inner_type(A.data_type);

    #if CHECKS {
        #run {
            assert(WORK.data_type == RT);
        }
        assert(m >= 0);
        assert(n >= 0);
    }

    if m <= 0 || n <= 1 then return;

    temp, x, y : A.data_type = ---;
    ctemp, stemp : RT = ---;

    for <j:1..n-1 {
        ctemp = fget(work, c_index + j - 1);
        stemp = fget(work, s_index + j - 1);
        for i:1..m {
            x = fget(a.*, row_start + i - 1, col_start + j - 1);
            y = fget(a.*, row_start + i - 1, col_start + j);
            temp = ctemp*x + stemp*y;
            fset(a, row_start + i - 1, col_start + j, ctemp*y - stemp*x);
            fset(a, row_start + i - 1, col_start + j - 1, temp);
        }
    }
}

_bdsqr_lasr_left_backward_variable :: (
    $$m: int,
    $$n: int,
    work: $WORK/ArrayType1D,
    $$c_index: int,
    $$s_index: int,
    a: *$A/ArrayType2D,
    $$row_start: int = 1,
    $$col_start: int = 1
) {
    RT :: #run get_number_inner_type(A.data_type);

    #if CHECKS {
        #run {
            assert(WORK.data_type == RT);
        }
        assert(m >= 0);
        assert(n >= 0);
    }

    if m <= 1 || n <= 0 then return;

    temp, x, y : A.data_type = ---;
    ctemp, stemp : RT = ---;

    for <i:1..m-1 {
        ctemp = fget(work, c_index + i - 1);
        stemp = fget(work, s_index + i - 1);
        for j:1..n {
            x = fget(a.*, row_start + i - 1, col_start + j - 1);
            y = fget(a.*, row_start + i,     col_start + j - 1);
            temp = ctemp*x + stemp*y;
            fset(a, row_start + i,     col_start + j - 1, ctemp*y - stemp*x);
            fset(a, row_start + i - 1, col_start + j - 1, temp);
        }
    }
}

_bdsqr_rot_rows :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$row1: int,
    $$row2: int,
    c: $RT,
    s: RT,
    $$col_start: int = 1
) {
    #if CHECKS {
        #run assert(RT == get_number_inner_type(A.data_type));
        assert(n >= 0);
    }

    if n <= 0 then return;

    temp, x, y : A.data_type = ---;
    for j:1..n {
        x = fget(a.*, row1, col_start + j - 1);
        y = fget(a.*, row2, col_start + j - 1);
        temp = c*x + s*y;
        fset(a, row2, col_start + j - 1, c*y - s*x);
        fset(a, row1, col_start + j - 1, temp);
    }
}

_bdsqr_rot_cols :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$col1: int,
    $$col2: int,
    c: $RT,
    s: RT,
    $$row_start: int = 1
) {
    #if CHECKS {
        #run assert(RT == get_number_inner_type(A.data_type));
        assert(n >= 0);
    }

    if n <= 0 then return;

    temp, x, y : A.data_type = ---;
    for i:1..n {
        x = fget(a.*, row_start + i - 1, col1);
        y = fget(a.*, row_start + i - 1, col2);
        temp = c*x + s*y;
        fset(a, row_start + i - 1, col2, c*y - s*x);
        fset(a, row_start + i - 1, col1, temp);
    }
}

_bdsqr_scale_row :: (
    $$n: int,
    alpha: $T,
    a: *$A/ArrayType2D,
    $$row: int,
    $$col_start: int = 1
) {
    #if CHECKS {
        #run assert(T == A.data_type);
        assert(n >= 0);
    }

    if n <= 0 || alpha == one(T) then return;

    for j:1..n {
        fset(a, row, col_start + j - 1, alpha * fget(a.*, row, col_start + j - 1));
    }
}

_bdsqr_scale_col :: (
    $$n: int,
    alpha: $T,
    a: *$A/ArrayType2D,
    $$col: int,
    $$row_start: int = 1
) {
    #if CHECKS {
        #run assert(T == A.data_type);
        assert(n >= 0);
    }

    if n <= 0 || alpha == one(T) then return;

    for i:1..n {
        fset(a, row_start + i - 1, col, alpha * fget(a.*, row_start + i - 1, col));
    }
}

_bdsqr_swap_rows :: (
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

_bdsqr_swap_cols :: (
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

_bdsqr_lartg :: (f: $RT, g: RT, cs: *RT, sn: *RT, r: *RT) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float64 {
        dlartg(f, g, cs, sn, r);
    } else {
        slartg(f, g, cs, sn, r);
    }
}

_bdsqr_las2 :: (f: $RT, g: RT, h: RT, ssmin: *RT, ssmax: *RT) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float64 {
        dlas2(f, g, h, ssmin, ssmax);
    } else {
        slas2(f, g, h, ssmin, ssmax);
    }
}

_bdsqr_lasq1 :: (
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return RT == float32 || RT == float64;
} {
    #if CHECKS {
        #run {
            assert(D.data_type    == RT);
            assert(E.data_type    == RT);
            assert(WORK.data_type == RT);
        }
    }

    #if RT == float64 {
        dlasq1(n, d, e, work, info);
    } else {
        slasq1(n, d, e, work, info);
    }
}

_bdsqr_lasv2 :: (
    f: $RT,
    g: RT,
    h: RT,
    ssmin: *RT,
    ssmax: *RT,
    snr: *RT,
    csr: *RT,
    snl: *RT,
    csl: *RT
) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float64 {
        dlasv2(f, g, h, ssmin, ssmax, snr, csr, snl, csl);
    } else {
        slasv2(f, g, h, ssmin, ssmax, snr, csr, snl, csl);
    }
}

_bdsqr :: (
    $$uplo: u8,
    $$n: int,
    $$ncvt: int,
    $$nru: int,
    $$ncc: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    vt: *$VT/ArrayType2D,
    $$ldvt: int,
    u: *$U/ArrayType2D,
    $$ldu: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if CHECKS {
        #run {
            assert(D.data_type    == RT);
            assert(E.data_type    == RT);
            assert(VT.data_type   == MT);
            assert(U.data_type    == MT);
            assert(C.data_type    == MT);
            assert(WORK.data_type == RT);
        }
    }

    maxitr :: 6;
    hndrth : RT : 0.01;
    hundred : RT : 100.0;
    meighth : RT : -0.125;
    one_r : RT : 1.0;
    ten : RT : 10.0;
    zero_r : RT : 0.0;
    negone : MT : ncast(MT, -1);
    routine_name :: #run if MT == float32 {
        "SBDSQR"
    } else if MT == float64 {
        "DBDSQR"
    } else if MT == C32 {
        "CBDSQR"
    } else {
        "ZBDSQR"
    };

    lower := char_eq(uplo, #char "l");
    rotate := (ncvt > 0) || (nru > 0) || (ncc > 0);

    i, idir, isub, iter, iterdivn, j, ll, lll, m, maxitdivn, nm1, nm12, nm13, oldll, oldm : int = ---;
    split : bool = ---;
    abse, abss, cosl, cosr, cs, eps, f, g, h, mu, oldcs, oldsn, r, shift, sigmn, sigmx, sinl, sinr, sll, smax, smin, sminoa, sn, thresh, tol, tolmul, unfl : RT = ---;

    info.* = 0;
    if !char_eq(uplo, #char "u", #char "l") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if ncvt < 0 {
        info.* = -3;
    } else if nru < 0 {
        info.* = -4;
    } else if ncc < 0 {
        info.* = -5;
    } else if (ncvt == 0 && ldvt < 1) || (ncvt > 0 && ldvt < max(1, n)) {
        info.* = -9;
    } else if ldu < max(1, nru) {
        info.* = -11;
    } else if (ncc == 0 && ldc < 1) || (ncc > 0 && ldc < max(1, n)) {
        info.* = -13;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if n == 0 then return;

    if n > 1 {
        if !rotate {
            _bdsqr_lasq1(n, d, e, work, info, RT);
            if info.* != 2 then return;
            info.* = 0;
        }

        nm1 = n - 1;
        nm12 = nm1 + nm1;
        nm13 = nm12 + nm1;
        idir = 0;

        #if RT == float64 {
            eps = EPS_64;
            unfl = SAFMIN_64;
        } else {
            eps = EPS_32;
            unfl = SAFMIN_32;
        }

        if lower {
            for i:1..n-1 {
                _bdsqr_lartg(fget(d.*, i), fget(e.*, i), *cs, *sn, *r);
                fset(d, i, r);
                fset(e, i, sn*fget(d.*, i+1));
                fset(d, i+1, cs*fget(d.*, i+1));
                fset(work, i, cs);
                fset(work, nm1+i, sn);
            }

            if nru > 0 {
                _bdsqr_lasr_right_forward_variable(nru, n, work.*, 1, n, u, 1, 1);
            }
            if ncc > 0 {
                _bdsqr_lasr_left_forward_variable(n, ncc, work.*, 1, n, c, 1, 1);
            }
        }

        tolmul = max(ten, min(hundred, eps**meighth));
        tol = tolmul*eps;

        smax = zero_r;
        for i:1..n {
            smax = max(smax, abs(fget(d.*, i)));
        }
        for i:1..n-1 {
            smax = max(smax, abs(fget(e.*, i)));
        }

        smin = zero_r;
        if tol >= zero_r {
            sminoa = abs(fget(d.*, 1));
            if sminoa != zero_r {
                mu = sminoa;
                for i:2..n {
                    mu = abs(fget(d.*, i))*(mu / (mu + abs(fget(e.*, i-1))));
                    sminoa = min(sminoa, mu);
                    if sminoa == zero_r then break;
                }
            }
            sminoa = sminoa / sqrt(cast(RT) n);
            thresh = max(tol*sminoa, cast(RT) maxitr*(cast(RT) n*(cast(RT) n*unfl)));
        } else {
            thresh = max(abs(tol)*smax, cast(RT) maxitr*(cast(RT) n*(cast(RT) n*unfl)));
        }

        maxitdivn = maxitr*n;
        iterdivn = 0;
        iter = -1;
        oldll = -1;
        oldm = -1;
        m = n;

        while true {
            if m <= 1 then break;

            if iter >= n {
                iter -= n;
                iterdivn += 1;
                if iterdivn >= maxitdivn {
                    info.* = 0;
                    for i:1..n-1 {
                        if fget(e.*, i) != zero_r {
                            info.* += 1;
                        }
                    }
                    return;
                }
            }

            if tol < zero_r && abs(fget(d.*, m)) <= thresh {
                fset(d, m, zero_r);
            }
            smax = abs(fget(d.*, m));

            ll = 0;
            lll = 1;
            while lll <= m - 1 {
                ll = m - lll;
                abss = abs(fget(d.*, ll));
                abse = abs(fget(e.*, ll));
                if tol < zero_r && abss <= thresh {
                    fset(d, ll, zero_r);
                }
                if abse <= thresh {
                    fset(e, ll, zero_r);
                    break;
                }
                smax = max(smax, abss, abse);
                lll += 1;
            }
            if lll > m - 1 {
                ll = 0;
            }

            if ll == m - 1 {
                m -= 1;
                continue;
            }

            ll += 1;

            if ll == m - 1 {
                _bdsqr_lasv2(fget(d.*, m-1), fget(e.*, m-1), fget(d.*, m), *sigmn, *sigmx, *sinr, *cosr, *sinl, *cosl);
                fset(d, m-1, sigmx);
                fset(e, m-1, zero_r);
                fset(d, m, sigmn);

                if ncvt > 0 then _bdsqr_rot_rows(ncvt, vt, m-1, m, cosr, sinr, 1);
                if nru  > 0 then _bdsqr_rot_cols(nru,  u,  m-1, m, cosl, sinl, 1);
                if ncc  > 0 then _bdsqr_rot_rows(ncc, c,  m-1, m, cosl, sinl, 1);

                m -= 2;
                continue;
            }

            if ll > oldm || m < oldll {
                if abs(fget(d.*, ll)) >= abs(fget(d.*, m)) {
                    idir = 1;
                } else {
                    idir = 2;
                }
            }

            if idir == 1 {
                if abs(fget(e.*, m-1)) <= abs(tol)*abs(fget(d.*, m)) ||
                   (tol < zero_r && abs(fget(e.*, m-1)) <= thresh) {
                    fset(e, m-1, zero_r);
                    continue;
                }

                if tol >= zero_r {
                    mu = abs(fget(d.*, ll));
                    smin = mu;
                    split = false;
                    for lll:ll..m-1 {
                        if abs(fget(e.*, lll)) <= tol*mu {
                            fset(e, lll, zero_r);
                            split = true;
                            break;
                        }
                        mu = abs(fget(d.*, lll+1))*(mu / (mu + abs(fget(e.*, lll))));
                        smin = min(smin, mu);
                    }
                    if split then continue;
                }
            } else {
                if abs(fget(e.*, ll)) <= abs(tol)*abs(fget(d.*, ll)) ||
                   (tol < zero_r && abs(fget(e.*, ll)) <= thresh) {
                    fset(e, ll, zero_r);
                    continue;
                }

                if tol >= zero_r {
                    mu = abs(fget(d.*, m));
                    smin = mu;
                    split = false;
                    for <lll:ll..m-1 {
                        if abs(fget(e.*, lll)) <= tol*mu {
                            fset(e, lll, zero_r);
                            split = true;
                            break;
                        }
                        mu = abs(fget(d.*, lll))*(mu / (mu + abs(fget(e.*, lll))));
                        smin = min(smin, mu);
                    }
                    if split then continue;
                }
            }

            oldll = ll;
            oldm = m;

            if tol >= zero_r && cast(RT) n*tol*(smin / smax) <= max(eps, hndrth*tol) {
                shift = zero_r;
            } else {
                if idir == 1 {
                    sll = abs(fget(d.*, ll));
                    _bdsqr_las2(fget(d.*, m-1), fget(e.*, m-1), fget(d.*, m), *shift, *r);
                } else {
                    sll = abs(fget(d.*, m));
                    _bdsqr_las2(fget(d.*, ll), fget(e.*, ll), fget(d.*, ll+1), *shift, *r);
                }

                if sll > zero_r {
                    if (shift / sll)**2 < eps {
                        shift = zero_r;
                    }
                }
            }

            iter += m - ll;

            if shift == zero_r {
                if idir == 1 {
                    cs = one_r;
                    oldcs = one_r;
                    for i:ll..m-1 {
                        _bdsqr_lartg(fget(d.*, i)*cs, fget(e.*, i), *cs, *sn, *r);
                        if i > ll {
                            fset(e, i-1, oldsn*r);
                        }
                        _bdsqr_lartg(oldcs*r, fget(d.*, i+1)*sn, *oldcs, *oldsn, *h);
                        fset(d, i, h);
                        fset(work, i-ll+1, cs);
                        fset(work, i-ll+1+nm1, sn);
                        fset(work, i-ll+1+nm12, oldcs);
                        fset(work, i-ll+1+nm13, oldsn);
                    }
                    h = fget(d.*, m)*cs;
                    fset(d, m, h*oldcs);
                    fset(e, m-1, h*oldsn);

                    if ncvt > 0 {
                        _bdsqr_lasr_left_forward_variable(m-ll+1, ncvt, work.*, 1, n, vt, ll, 1);
                    }
                    if nru > 0 {
                        _bdsqr_lasr_right_forward_variable(nru, m-ll+1, work.*, nm12+1, nm13+1, u, 1, ll);
                    }
                    if ncc > 0 {
                        _bdsqr_lasr_left_forward_variable(m-ll+1, ncc, work.*, nm12+1, nm13+1, c, ll, 1);
                    }

                    if abs(fget(e.*, m-1)) <= thresh {
                        fset(e, m-1, zero_r);
                    }
                } else {
                    cs = one_r;
                    oldcs = one_r;
                    for <i:ll+1..m {
                        _bdsqr_lartg(fget(d.*, i)*cs, fget(e.*, i-1), *cs, *sn, *r);
                        if i < m {
                            fset(e, i, oldsn*r);
                        }
                        _bdsqr_lartg(oldcs*r, fget(d.*, i-1)*sn, *oldcs, *oldsn, *h);
                        fset(d, i, h);
                        fset(work, i-ll, cs);
                        fset(work, i-ll+nm1, -sn);
                        fset(work, i-ll+nm12, oldcs);
                        fset(work, i-ll+nm13, -oldsn);
                    }
                    h = fget(d.*, ll)*cs;
                    fset(d, ll, h*oldcs);
                    fset(e, ll, h*oldsn);

                    if ncvt > 0 {
                        _bdsqr_lasr_left_backward_variable(m-ll+1, ncvt, work.*, nm12+1, nm13+1, vt, ll, 1);
                    }
                    if nru > 0 {
                        _bdsqr_lasr_right_backward_variable(nru, m-ll+1, work.*, 1, n, u, 1, ll);
                    }
                    if ncc > 0 {
                        _bdsqr_lasr_left_backward_variable(m-ll+1, ncc, work.*, 1, n, c, ll, 1);
                    }

                    if abs(fget(e.*, ll)) <= thresh {
                        fset(e, ll, zero_r);
                    }
                }
            } else {
                if idir == 1 {
                    f = (abs(fget(d.*, ll)) - shift) * (sign(one_r, fget(d.*, ll)) + shift / fget(d.*, ll));
                    g = fget(e.*, ll);
                    for i:ll..m-1 {
                        _bdsqr_lartg(f, g, *cosr, *sinr, *r);
                        if i > ll {
                            fset(e, i-1, r);
                        }
                        f = cosr*fget(d.*, i) + sinr*fget(e.*, i);
                        fset(e, i, cosr*fget(e.*, i) - sinr*fget(d.*, i));
                        g = sinr*fget(d.*, i+1);
                        fset(d, i+1, cosr*fget(d.*, i+1));
                        _bdsqr_lartg(f, g, *cosl, *sinl, *r);
                        fset(d, i, r);
                        f = cosl*fget(e.*, i) + sinl*fget(d.*, i+1);
                        fset(d, i+1, cosl*fget(d.*, i+1) - sinl*fget(e.*, i));
                        if i < m - 1 {
                            g = sinl*fget(e.*, i+1);
                            fset(e, i+1, cosl*fget(e.*, i+1));
                        }
                        fset(work, i-ll+1, cosr);
                        fset(work, i-ll+1+nm1, sinr);
                        fset(work, i-ll+1+nm12, cosl);
                        fset(work, i-ll+1+nm13, sinl);
                    }
                    fset(e, m-1, f);

                    if ncvt > 0 {
                        _bdsqr_lasr_left_forward_variable(m-ll+1, ncvt, work.*, 1, n, vt, ll, 1);
                    }
                    if nru > 0 {
                        _bdsqr_lasr_right_forward_variable(nru, m-ll+1, work.*, nm12+1, nm13+1, u, 1, ll);
                    }
                    if ncc > 0 {
                        _bdsqr_lasr_left_forward_variable(m-ll+1, ncc, work.*, nm12+1, nm13+1, c, ll, 1);
                    }

                    if abs(fget(e.*, m-1)) <= thresh {
                        fset(e, m-1, zero_r);
                    }
                } else {
                    f = (abs(fget(d.*, m)) - shift) * (sign(one_r, fget(d.*, m)) + shift / fget(d.*, m));
                    g = fget(e.*, m-1);
                    for <i:ll+1..m {
                        _bdsqr_lartg(f, g, *cosr, *sinr, *r);
                        if i < m {
                            fset(e, i, r);
                        }
                        f = cosr*fget(d.*, i) + sinr*fget(e.*, i-1);
                        fset(e, i-1, cosr*fget(e.*, i-1) - sinr*fget(d.*, i));
                        g = sinr*fget(d.*, i-1);
                        fset(d, i-1, cosr*fget(d.*, i-1));
                        _bdsqr_lartg(f, g, *cosl, *sinl, *r);
                        fset(d, i, r);
                        f = cosl*fget(e.*, i-1) + sinl*fget(d.*, i-1);
                        fset(d, i-1, cosl*fget(d.*, i-1) - sinl*fget(e.*, i-1));
                        if i > ll + 1 {
                            g = sinl*fget(e.*, i-2);
                            fset(e, i-2, cosl*fget(e.*, i-2));
                        }
                        fset(work, i-ll, cosr);
                        fset(work, i-ll+nm1, -sinr);
                        fset(work, i-ll+nm12, cosl);
                        fset(work, i-ll+nm13, -sinl);
                    }
                    fset(e, ll, f);

                    if abs(fget(e.*, ll)) <= thresh {
                        fset(e, ll, zero_r);
                    }

                    if ncvt > 0 {
                        _bdsqr_lasr_left_backward_variable(m-ll+1, ncvt, work.*, nm12+1, nm13+1, vt, ll, 1);
                    }
                    if nru > 0 {
                        _bdsqr_lasr_right_backward_variable(nru, m-ll+1, work.*, 1, n, u, 1, ll);
                    }
                    if ncc > 0 {
                        _bdsqr_lasr_left_backward_variable(m-ll+1, ncc, work.*, 1, n, c, ll, 1);
                    }
                }
            }
        }
    }

    for i:1..n {
        if fget(d.*, i) == zero_r {
            fset(d, i, zero_r);
        }
        if fget(d.*, i) < zero_r {
            fset(d, i, -fget(d.*, i));
            if ncvt > 0 {
                _bdsqr_scale_row(ncvt, negone, vt, i, 1);
            }
        }
    }

    for i:1..n-1 {
        isub = 1;
        smin = fget(d.*, 1);
        for j:2..n+1-i {
            if fget(d.*, j) <= smin {
                isub = j;
                smin = fget(d.*, j);
            }
        }

        if isub != n+1-i {
            fset(d, isub, fget(d.*, n+1-i));
            fset(d, n+1-i, smin);
            if ncvt > 0 then _bdsqr_swap_rows(ncvt, vt, isub, n+1-i, 1);
            if nru  > 0 then _bdsqr_swap_cols(nru,  u,  isub, n+1-i, 1);
            if ncc  > 0 then _bdsqr_swap_rows(ncc, c,  isub, n+1-i, 1);
        }
    }
}
