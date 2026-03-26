// ? LAPACK *hgeqz()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chgeqz/zhgeqz
// *HGEQZ computes the eigenvalues of a complex Hessenberg-triangular pair.

chgeqz :: #bake_arguments _hgeqz(RT=float32, MT=C32);
zhgeqz :: #bake_arguments _hgeqz(RT=float64, MT=C64);

#scope_file
_hgeqz_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hgeqz_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hgeqz_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hgeqz_row_slice :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$n: int
) -> [] A.data_type {
    count := ifx n <= 0 then 0 else (n - 1)*lda + 1;
    return _hgeqz_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hgeqz_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hgeqz_abs1 :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return abs(real(x)) + abs(imag(x));
}

_hgeqz_set_work1 :: inline (work: *$WORK/ArrayType1D, $$x: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} {
    fset(work, 1, _hgeqz_complex_from_real(ncast(RT, x), RT, MT));
}

_hgeqz_laset :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: MT, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_hgeqz_scal :: (
    $$n: int, alpha: $MT, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        cscal(n, alpha, x, incx);
    } else {
        zscal(n, alpha, x, incx);
    }
}

_hgeqz_lartg :: ($$f: $MT, $$g: MT, c: *$RT, s: *MT, r: *MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        clartg(f, g, c, s, r);
    } else {
        zlartg(f, g, c, s, r);
    }
}

_hgeqz_rot :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $$c: $RT, $$s: $MT, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        crot(n, x, incx, y, incy, c, s);
    } else {
        zrot(n, x, incx, y, incy, c, s);
    }
}

_hgeqz_div :: ($$x: $MT, $$y: MT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        return cladiv(x, y);
    } else {
        return zladiv(x, y);
    }
}

_hgeqz_lanhs :: (
    $$norm: u8, $$n: int, a: $A/ArrayType2D, $$lda: int, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && WORK.data_type == RT;
} {
    #if MT == C32 {
        return clanhs(norm, n, a, lda, work);
    } else {
        return zlanhs(norm, n, a, lda, work);
    }
}

_hgeqz_lamch :: ($$cmach: u8, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach);
    } else {
        return dlamch(cmach);
    }
}

_hgeqz :: (
    $$job: u8, $$compq: u8, $$compz: u8, $$n: int, $$ilo: int, $$ihi: int,
    h: *$H/ArrayType2D, $$ldh: int, t: *$T/ArrayType2D, $$ldt: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           H.data_type == MT && T.data_type == MT && ALPHA.data_type == MT && BETA.data_type == MT &&
           Q.data_type == MT && Z.data_type == MT && WORK.data_type == MT && RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    czero  :: #run zero(MT);
    cone   :: #run _hgeqz_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHGEQZ" } else { "ZHGEQZ" };

    do_clear_last, do_qz_step, do_standardize, failed, finished, ilazr2, ilazro, ilq, ilschr, ilz, lquery : bool = ---;
    icompq, icompz, ifirst, ifrstm, iiter, ilast, ilastm, in, ischur, istart, j, jc, jch, jiter, jr, maxit : int = ---;
    absb, anorm, ascale, atol, bnorm, bscale, btol, cval, safmin, temp, temp2, tempr, ulp : RT = ---;
    abi12, abi22, ad11, ad12, ad21, ad22, ctemp, ctemp2, ctemp3, eshift, s, shift, signbc, u12, x, y : MT = ---;

    if lsame(job, #char "E") {
        ilschr = false;
        ischur = 1;
    } else if lsame(job, #char "S") {
        ilschr = true;
        ischur = 2;
    } else {
        ilschr = true;
        ischur = 0;
    }

    if lsame(compq, #char "N") {
        ilq = false;
        icompq = 1;
    } else if lsame(compq, #char "V") {
        ilq = true;
        icompq = 2;
    } else if lsame(compq, #char "I") {
        ilq = true;
        icompq = 3;
    } else {
        ilq = true;
        icompq = 0;
    }

    if lsame(compz, #char "N") {
        ilz = false;
        icompz = 1;
    } else if lsame(compz, #char "V") {
        ilz = true;
        icompz = 2;
    } else if lsame(compz, #char "I") {
        ilz = true;
        icompz = 3;
    } else {
        ilz = true;
        icompz = 0;
    }

    info.* = 0;
    _hgeqz_set_work1(work, max(1, n), RT, MT);
    lquery = lwork == -1;
    if ischur == 0 {
        info.* = -1;
    } else if icompq == 0 {
        info.* = -2;
    } else if icompz == 0 {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if ilo < 1 {
        info.* = -5;
    } else if ihi > n || ihi < ilo - 1 {
        info.* = -6;
    } else if ldh < max(1, n) {
        info.* = -8;
    } else if ldt < max(1, n) {
        info.* = -10;
    } else if ldq < 1 || (ilq && ldq < n) {
        info.* = -14;
    } else if ldz < 1 || (ilz && ldz < n) {
        info.* = -16;
    } else if lwork < max(1, n) && !lquery {
        info.* = -18;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n <= 0 {
        _hgeqz_set_work1(work, 1, RT, MT);
        return;
    }

    if icompq == 3 {
        q_full_raw := _hgeqz_matrix_storage(q.*, ldq, 1, 1, n, n);
        q_full := array2d_view(q_full_raw, ldq, row_major=false);
        _hgeqz_laset(#char "F", n, n, czero, cone, *q_full, ldq, MT);
    }
    if icompz == 3 {
        z_full_raw := _hgeqz_matrix_storage(z.*, ldz, 1, 1, n, n);
        z_full := array2d_view(z_full_raw, ldz, row_major=false);
        _hgeqz_laset(#char "F", n, n, czero, cone, *z_full, ldz, MT);
    }

    in = ihi + 1 - ilo;
    safmin = _hgeqz_lamch(#char "S", RT);
    ulp = _hgeqz_lamch(#char "E", RT) * _hgeqz_lamch(#char "B", RT);
    if in > 0 {
        h_sub_raw := _hgeqz_matrix_storage(h.*, ldh, ilo, ilo, in, in);
        h_sub := array2d_view(h_sub_raw, ldh, row_major=false);
        t_sub_raw := _hgeqz_matrix_storage(t.*, ldt, ilo, ilo, in, in);
        t_sub := array2d_view(t_sub_raw, ldt, row_major=false);
        anorm = _hgeqz_lanhs(#char "F", in, *h_sub, ldh, rwork, RT, MT);
        bnorm = _hgeqz_lanhs(#char "F", in, *t_sub, ldt, rwork, RT, MT);
    } else {
        anorm = zero_r;
        bnorm = zero_r;
    }
    atol = max(safmin, ulp*anorm);
    btol = max(safmin, ulp*bnorm);
    ascale = one_r / max(safmin, anorm);
    bscale = one_r / max(safmin, bnorm);

    for j:ihi + 1..n {
        absb = abs(fget(t.*, j, j));
        if absb > safmin {
            signbc = conjugate(fget(t.*, j, j) / _hgeqz_complex_from_real(absb, RT, MT));
            fset(t, j, j, _hgeqz_complex_from_real(absb, RT, MT));
            if ilschr {
                if j > 1 {
                    t_col := _hgeqz_raw_slice(t.data + (j - 1)*ldt, j - 1);
                    _hgeqz_scal(j - 1, signbc, *t_col, 1, MT);
                }
                h_col := _hgeqz_raw_slice(h.data + (j - 1)*ldh, j);
                _hgeqz_scal(j, signbc, *h_col, 1, MT);
            } else {
                h_jj := _hgeqz_raw_slice(h.data + (j - 1)*ldh + (j - 1), 1);
                _hgeqz_scal(1, signbc, *h_jj, 1, MT);
            }
            if ilz {
                z_col := _hgeqz_raw_slice(z.data + (j - 1)*ldz, n);
                _hgeqz_scal(n, signbc, *z_col, 1, MT);
            }
        } else {
            fset(t, j, j, czero);
        }
        fset(alpha, j, fget(h.*, j, j));
        fset(beta, j, fget(t.*, j, j));
    }

    finished = ihi < ilo;
    failed = false;

    if !finished {
        ilast = ihi;
        if ilschr {
            ifrstm = 1;
            ilastm = n;
        } else {
            ifrstm = ilo;
            ilastm = ihi;
        }
        iiter = 0;
        eshift = czero;
        maxit = 30*(ihi - ilo + 1);

        for jiter:1..maxit {
            do_clear_last = false;
            do_qz_step = false;
            do_standardize = false;

            if ilast == ilo {
                do_standardize = true;
            } else if _hgeqz_abs1(fget(h.*, ilast, ilast - 1), RT, MT) <=
                      max(safmin, ulp*(_hgeqz_abs1(fget(h.*, ilast, ilast), RT, MT) +
                                       _hgeqz_abs1(fget(h.*, ilast - 1, ilast - 1), RT, MT))) {
                fset(h, ilast, ilast - 1, czero);
                do_standardize = true;
            }

            if !do_standardize {
                if abs(fget(t.*, ilast, ilast)) <= btol {
                    fset(t, ilast, ilast, czero);
                    do_clear_last = true;
                } else {
                    failed = true;
                    for <j:ilast - 1..ilo {
                        if j == ilo {
                            ilazro = true;
                        } else if _hgeqz_abs1(fget(h.*, j, j - 1), RT, MT) <=
                                  max(safmin, ulp*(_hgeqz_abs1(fget(h.*, j, j), RT, MT) +
                                                   _hgeqz_abs1(fget(h.*, j - 1, j - 1), RT, MT))) {
                            fset(h, j, j - 1, czero);
                            ilazro = true;
                        } else {
                            ilazro = false;
                        }

                        if abs(fget(t.*, j, j)) < btol {
                            fset(t, j, j, czero);
                            ilazr2 = false;
                            if !ilazro &&
                               _hgeqz_abs1(fget(h.*, j, j - 1), RT, MT)*
                               (ascale*_hgeqz_abs1(fget(h.*, j + 1, j), RT, MT)) <=
                               _hgeqz_abs1(fget(h.*, j, j), RT, MT)*(ascale*atol) {
                                ilazr2 = true;
                            }

                            if ilazro || ilazr2 {
                                failed = false;
                                for jch:j..ilast - 1 {
                                    ctemp = fget(h.*, jch, jch);
                                    _hgeqz_lartg(ctemp, fget(h.*, jch + 1, jch), *cval, *s, *ctemp, RT, MT);
                                    fset(h, jch, jch, ctemp);
                                    fset(h, jch + 1, jch, czero);

                                    h_row_top := _hgeqz_row_slice(h.*, ldh, jch, jch + 1, ilastm - jch);
                                    h_row_bot := _hgeqz_row_slice(h.*, ldh, jch + 1, jch + 1, ilastm - jch);
                                    _hgeqz_rot(ilastm - jch, *h_row_top, ldh, *h_row_bot, ldh, cval, s, RT, MT);

                                    t_row_top := _hgeqz_row_slice(t.*, ldt, jch, jch + 1, ilastm - jch);
                                    t_row_bot := _hgeqz_row_slice(t.*, ldt, jch + 1, jch + 1, ilastm - jch);
                                    _hgeqz_rot(ilastm - jch, *t_row_top, ldt, *t_row_bot, ldt, cval, s, RT, MT);

                                    if ilq {
                                        q_col_jch := _hgeqz_raw_slice(q.data + (jch - 1)*ldq, n);
                                        q_col_jch1 := _hgeqz_raw_slice(q.data + jch*ldq, n);
                                        _hgeqz_rot(n, *q_col_jch, 1, *q_col_jch1, 1, cval, conjugate(s), RT, MT);
                                    }
                                    if ilazr2 {
                                        fset(h, jch, jch - 1, fget(h.*, jch, jch - 1)*_hgeqz_complex_from_real(cval, RT, MT));
                                    }
                                    ilazr2 = false;

                                    if _hgeqz_abs1(fget(t.*, jch + 1, jch + 1), RT, MT) >= btol {
                                        if jch + 1 >= ilast {
                                            do_standardize = true;
                                        } else {
                                            ifirst = jch + 1;
                                            do_qz_step = true;
                                        }
                                        break;
                                    }
                                    fset(t, jch + 1, jch + 1, czero);
                                }

                                if !do_standardize && !do_qz_step {
                                    do_clear_last = true;
                                }
                                break;
                            } else {
                                failed = false;
                                for jch:j..ilast - 1 {
                                    ctemp = fget(t.*, jch, jch + 1);
                                    _hgeqz_lartg(ctemp, fget(t.*, jch + 1, jch + 1), *cval, *s, *ctemp, RT, MT);
                                    fset(t, jch, jch + 1, ctemp);
                                    fset(t, jch + 1, jch + 1, czero);

                                    if jch < ilastm - 1 {
                                        t_row_top := _hgeqz_row_slice(t.*, ldt, jch, jch + 2, ilastm - jch - 1);
                                        t_row_bot := _hgeqz_row_slice(t.*, ldt, jch + 1, jch + 2, ilastm - jch - 1);
                                        _hgeqz_rot(ilastm - jch - 1, *t_row_top, ldt, *t_row_bot, ldt, cval, s, RT, MT);
                                    }

                                    h_row_top := _hgeqz_row_slice(h.*, ldh, jch, jch - 1, ilastm - jch + 2);
                                    h_row_bot := _hgeqz_row_slice(h.*, ldh, jch + 1, jch - 1, ilastm - jch + 2);
                                    _hgeqz_rot(ilastm - jch + 2, *h_row_top, ldh, *h_row_bot, ldh, cval, s, RT, MT);

                                    if ilq {
                                        q_col_jch := _hgeqz_raw_slice(q.data + (jch - 1)*ldq, n);
                                        q_col_jch1 := _hgeqz_raw_slice(q.data + jch*ldq, n);
                                        _hgeqz_rot(n, *q_col_jch, 1, *q_col_jch1, 1, cval, conjugate(s), RT, MT);
                                    }

                                    ctemp = fget(h.*, jch + 1, jch);
                                    _hgeqz_lartg(ctemp, fget(h.*, jch + 1, jch - 1), *cval, *s, *ctemp, RT, MT);
                                    fset(h, jch + 1, jch, ctemp);
                                    fset(h, jch + 1, jch - 1, czero);

                                    h_col_j := _hgeqz_raw_slice(h.data + (jch - 1)*ldh + (ifrstm - 1), jch + 1 - ifrstm);
                                    h_col_jm1 := _hgeqz_raw_slice(h.data + (jch - 2)*ldh + (ifrstm - 1), jch + 1 - ifrstm);
                                    _hgeqz_rot(jch + 1 - ifrstm, *h_col_j, 1, *h_col_jm1, 1, cval, s, RT, MT);

                                    t_col_j := _hgeqz_raw_slice(t.data + (jch - 1)*ldt + (ifrstm - 1), jch - ifrstm);
                                    t_col_jm1 := _hgeqz_raw_slice(t.data + (jch - 2)*ldt + (ifrstm - 1), jch - ifrstm);
                                    _hgeqz_rot(jch - ifrstm, *t_col_j, 1, *t_col_jm1, 1, cval, s, RT, MT);

                                    if ilz {
                                        z_col_j := _hgeqz_raw_slice(z.data + (jch - 1)*ldz, n);
                                        z_col_jm1 := _hgeqz_raw_slice(z.data + (jch - 2)*ldz, n);
                                        _hgeqz_rot(n, *z_col_j, 1, *z_col_jm1, 1, cval, s, RT, MT);
                                    }
                                }
                                do_clear_last = true;
                                break;
                            }
                        } else if ilazro {
                            ifirst = j;
                            do_qz_step = true;
                            failed = false;
                            break;
                        }
                    }

                    if failed {
                        info.* = 2*n + 1;
                        break;
                    }
                }
            }

            if do_clear_last {
                ctemp = fget(h.*, ilast, ilast);
                _hgeqz_lartg(ctemp, fget(h.*, ilast, ilast - 1), *cval, *s, *ctemp, RT, MT);
                fset(h, ilast, ilast, ctemp);
                fset(h, ilast, ilast - 1, czero);

                h_col_last := _hgeqz_raw_slice(h.data + (ilast - 1)*ldh + (ifrstm - 1), ilast - ifrstm);
                h_col_prev := _hgeqz_raw_slice(h.data + (ilast - 2)*ldh + (ifrstm - 1), ilast - ifrstm);
                _hgeqz_rot(ilast - ifrstm, *h_col_last, 1, *h_col_prev, 1, cval, s, RT, MT);

                t_col_last := _hgeqz_raw_slice(t.data + (ilast - 1)*ldt + (ifrstm - 1), ilast - ifrstm);
                t_col_prev := _hgeqz_raw_slice(t.data + (ilast - 2)*ldt + (ifrstm - 1), ilast - ifrstm);
                _hgeqz_rot(ilast - ifrstm, *t_col_last, 1, *t_col_prev, 1, cval, s, RT, MT);

                if ilz {
                    z_col_last := _hgeqz_raw_slice(z.data + (ilast - 1)*ldz, n);
                    z_col_prev := _hgeqz_raw_slice(z.data + (ilast - 2)*ldz, n);
                    _hgeqz_rot(n, *z_col_last, 1, *z_col_prev, 1, cval, s, RT, MT);
                }
                do_standardize = true;
            }

            if do_standardize {
                absb = abs(fget(t.*, ilast, ilast));
                if absb > safmin {
                    signbc = conjugate(fget(t.*, ilast, ilast) / _hgeqz_complex_from_real(absb, RT, MT));
                    fset(t, ilast, ilast, _hgeqz_complex_from_real(absb, RT, MT));
                    if ilschr {
                        t_col := _hgeqz_raw_slice(t.data + (ilast - 1)*ldt + (ifrstm - 1), ilast - ifrstm);
                        _hgeqz_scal(ilast - ifrstm, signbc, *t_col, 1, MT);
                        h_col := _hgeqz_raw_slice(h.data + (ilast - 1)*ldh + (ifrstm - 1), ilast + 1 - ifrstm);
                        _hgeqz_scal(ilast + 1 - ifrstm, signbc, *h_col, 1, MT);
                    } else {
                        h_jj := _hgeqz_raw_slice(h.data + (ilast - 1)*ldh + (ilast - 1), 1);
                        _hgeqz_scal(1, signbc, *h_jj, 1, MT);
                    }
                    if ilz {
                        z_col := _hgeqz_raw_slice(z.data + (ilast - 1)*ldz, n);
                        _hgeqz_scal(n, signbc, *z_col, 1, MT);
                    }
                } else {
                    fset(t, ilast, ilast, czero);
                }
                fset(alpha, ilast, fget(h.*, ilast, ilast));
                fset(beta, ilast, fget(t.*, ilast, ilast));

                ilast -= 1;
                if ilast < ilo {
                    finished = true;
                    break;
                }

                iiter = 0;
                eshift = czero;
                if !ilschr {
                    ilastm = ilast;
                    if ifrstm > ilast {
                        ifrstm = ilo;
                    }
                }
                continue;
            }

            iiter += 1;
            if !ilschr {
                ifrstm = ifirst;
            }

            if ((iiter / 10)*10) != iiter {
                u12 = (bscale*fget(t.*, ilast - 1, ilast)) / (bscale*fget(t.*, ilast, ilast));
                ad11 = (ascale*fget(h.*, ilast - 1, ilast - 1)) / (bscale*fget(t.*, ilast - 1, ilast - 1));
                ad21 = (ascale*fget(h.*, ilast, ilast - 1)) / (bscale*fget(t.*, ilast - 1, ilast - 1));
                ad12 = (ascale*fget(h.*, ilast - 1, ilast)) / (bscale*fget(t.*, ilast, ilast));
                ad22 = (ascale*fget(h.*, ilast, ilast)) / (bscale*fget(t.*, ilast, ilast));
                abi22 = ad22 - u12*ad21;
                abi12 = ad12 - u12*ad11;

                shift = abi22;
                ctemp = sqrt(abi12)*sqrt(ad21);
                temp = _hgeqz_abs1(ctemp, RT, MT);
                if ctemp != czero {
                    x = ncast(RT, 0.5)*(ad11 - shift);
                    temp2 = _hgeqz_abs1(x, RT, MT);
                    temp = max(temp, temp2);
                    ctemp2 = _hgeqz_complex_from_real(temp, RT, MT);
                    y = ctemp2*sqrt((x / ctemp2)**2 + (ctemp / ctemp2)**2);
                    if temp2 > zero_r {
                        ctemp3 = x / _hgeqz_complex_from_real(temp2, RT, MT);
                        if real(ctemp3)*real(y) + imag(ctemp3)*imag(y) < zero_r {
                            y = -y;
                        }
                    }
                    shift -= ctemp*_hgeqz_div(ctemp, x + y, RT, MT);
                }
            } else {
                if ((iiter / 20)*20) == iiter &&
                   bscale*_hgeqz_abs1(fget(t.*, ilast, ilast), RT, MT) > safmin {
                    eshift += (ascale*fget(h.*, ilast, ilast)) / (bscale*fget(t.*, ilast, ilast));
                } else {
                    eshift += (ascale*fget(h.*, ilast, ilast - 1)) / (bscale*fget(t.*, ilast - 1, ilast - 1));
                }
                shift = eshift;
            }

            istart = ifirst;
            for <j:ilast - 1..ifirst + 1 {
                istart = j;
                ctemp = ascale*fget(h.*, j, j) - shift*(bscale*fget(t.*, j, j));
                temp = _hgeqz_abs1(ctemp, RT, MT);
                temp2 = ascale*_hgeqz_abs1(fget(h.*, j + 1, j), RT, MT);
                tempr = max(temp, temp2);
                if tempr < one_r && tempr != zero_r {
                    temp /= tempr;
                    temp2 /= tempr;
                }
                if _hgeqz_abs1(fget(h.*, j, j - 1), RT, MT)*temp2 <= temp*atol {
                    break;
                }
            }

            ctemp = ascale*fget(h.*, istart, istart) - shift*(bscale*fget(t.*, istart, istart));
            ctemp2 = ascale*fget(h.*, istart + 1, istart);
            _hgeqz_lartg(ctemp, ctemp2, *cval, *s, *ctemp3, RT, MT);

            for j:istart..ilast - 1 {
                if j > istart {
                    ctemp = fget(h.*, j, j - 1);
                    _hgeqz_lartg(ctemp, fget(h.*, j + 1, j - 1), *cval, *s, *ctemp, RT, MT);
                    fset(h, j, j - 1, ctemp);
                    fset(h, j + 1, j - 1, czero);
                }

                for jc:j..ilastm {
                    ctemp = _hgeqz_complex_from_real(cval, RT, MT)*fget(h.*, j, jc) + s*fget(h.*, j + 1, jc);
                    fset(h, j + 1, jc, -conjugate(s)*fget(h.*, j, jc) + _hgeqz_complex_from_real(cval, RT, MT)*fget(h.*, j + 1, jc));
                    fset(h, j, jc, ctemp);

                    ctemp2 = _hgeqz_complex_from_real(cval, RT, MT)*fget(t.*, j, jc) + s*fget(t.*, j + 1, jc);
                    fset(t, j + 1, jc, -conjugate(s)*fget(t.*, j, jc) + _hgeqz_complex_from_real(cval, RT, MT)*fget(t.*, j + 1, jc));
                    fset(t, j, jc, ctemp2);
                }

                if ilq {
                    for jr:1..n {
                        ctemp = _hgeqz_complex_from_real(cval, RT, MT)*fget(q.*, jr, j) + conjugate(s)*fget(q.*, jr, j + 1);
                        fset(q, jr, j + 1, -s*fget(q.*, jr, j) + _hgeqz_complex_from_real(cval, RT, MT)*fget(q.*, jr, j + 1));
                        fset(q, jr, j, ctemp);
                    }
                }

                ctemp = fget(t.*, j + 1, j + 1);
                _hgeqz_lartg(ctemp, fget(t.*, j + 1, j), *cval, *s, *ctemp, RT, MT);
                fset(t, j + 1, j + 1, ctemp);
                fset(t, j + 1, j, czero);

                for jr:ifrstm..min(j + 2, ilast) {
                    ctemp = _hgeqz_complex_from_real(cval, RT, MT)*fget(h.*, jr, j + 1) + s*fget(h.*, jr, j);
                    fset(h, jr, j, -conjugate(s)*fget(h.*, jr, j + 1) + _hgeqz_complex_from_real(cval, RT, MT)*fget(h.*, jr, j));
                    fset(h, jr, j + 1, ctemp);
                }
                for jr:ifrstm..j {
                    ctemp = _hgeqz_complex_from_real(cval, RT, MT)*fget(t.*, jr, j + 1) + s*fget(t.*, jr, j);
                    fset(t, jr, j, -conjugate(s)*fget(t.*, jr, j + 1) + _hgeqz_complex_from_real(cval, RT, MT)*fget(t.*, jr, j));
                    fset(t, jr, j + 1, ctemp);
                }
                if ilz {
                    for jr:1..n {
                        ctemp = _hgeqz_complex_from_real(cval, RT, MT)*fget(z.*, jr, j + 1) + s*fget(z.*, jr, j);
                        fset(z, jr, j, -conjugate(s)*fget(z.*, jr, j + 1) + _hgeqz_complex_from_real(cval, RT, MT)*fget(z.*, jr, j));
                        fset(z, jr, j + 1, ctemp);
                    }
                }
            }
        }

        if !finished && !failed {
            info.* = ilast;
            failed = true;
        }
    }

    if finished {
        for j:1..ilo - 1 {
            absb = abs(fget(t.*, j, j));
            if absb > safmin {
                signbc = conjugate(fget(t.*, j, j) / _hgeqz_complex_from_real(absb, RT, MT));
                fset(t, j, j, _hgeqz_complex_from_real(absb, RT, MT));
                if ilschr {
                    if j > 1 {
                        t_col := _hgeqz_raw_slice(t.data + (j - 1)*ldt, j - 1);
                        _hgeqz_scal(j - 1, signbc, *t_col, 1, MT);
                    }
                    h_col := _hgeqz_raw_slice(h.data + (j - 1)*ldh, j);
                    _hgeqz_scal(j, signbc, *h_col, 1, MT);
                } else {
                    h_jj := _hgeqz_raw_slice(h.data + (j - 1)*ldh + (j - 1), 1);
                    _hgeqz_scal(1, signbc, *h_jj, 1, MT);
                }
                if ilz {
                    z_col := _hgeqz_raw_slice(z.data + (j - 1)*ldz, n);
                    _hgeqz_scal(n, signbc, *z_col, 1, MT);
                }
            } else {
                fset(t, j, j, czero);
            }
            fset(alpha, j, fget(h.*, j, j));
            fset(beta, j, fget(t.*, j, j));
        }
        info.* = 0;
    }

    _hgeqz_set_work1(work, n, RT, MT);
}
