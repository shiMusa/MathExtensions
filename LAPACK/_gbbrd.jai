// ? LAPACK *gbbrd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbbrd/cgbbrd/dgbbrd/zgbbrd
// *GBBRD reduces a general band matrix to upper bidiagonal form.

sgbbrd :: (
    $$vect: u8,
    $$m: int,
    $$n: int,
    $$ncc: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    pt: *$PT/ArrayType2D,
    $$ldpt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    info: *int
) {
    _gbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc,
           work, work, 1, max(m, n)+1, info, RT=float32, MT=float32);
}

dgbbrd :: (
    $$vect: u8,
    $$m: int,
    $$n: int,
    $$ncc: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    pt: *$PT/ArrayType2D,
    $$ldpt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    info: *int
) {
    _gbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc,
           work, work, 1, max(m, n)+1, info, RT=float64, MT=float64);
}

cgbbrd :: (
    $$vect: u8,
    $$m: int,
    $$n: int,
    $$ncc: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    pt: *$PT/ArrayType2D,
    $$ldpt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc,
           work, rwork, 1, 1, info, RT=float32, MT=C32);
}

zgbbrd :: (
    $$vect: u8,
    $$m: int,
    $$n: int,
    $$ncc: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    pt: *$PT/ArrayType2D,
    $$ldpt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc,
           work, rwork, 1, 1, info, RT=float64, MT=C64);
}

#scope_file
_gbbrd_set_identity :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$diag: int
) {
    #if CHECKS {
        assert(m >= 0);
        assert(n >= 0);
    }

    zero_t :: #run zero(A.data_type);
    one_t  :: #run one(A.data_type);

    for j:1..n {
        for i:1..m {
            fset(a, i, j, zero_t);
        }
    }
    for i:1..diag {
        fset(a, i, i, one_t);
    }
}

_gbbrd_apply_rotation :: inline (x: $MT, y: MT, c: $RT, s: MT) -> MT, MT {
    return c*x + s*y, c*y - conjugate(s)*x;
}

_gbbrd_rotate_2d_strided :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$x_row: int,
    $$x_col: int,
    $$x_row_inc: int,
    $$x_col_inc: int,
    $$y_row: int,
    $$y_col: int,
    $$y_row_inc: int,
    $$y_col_inc: int,
    c: $RT,
    s: $MT
) {
    #if CHECKS {
        #run {
            assert(A.data_type == MT);
        }
        assert(n >= 0);
    }

    if n <= 0 then return;

    xr := x_row;
    xc := x_col;
    yr := y_row;
    yc := y_col;
    x, y, nx, ny : MT = ---;

    for 1..n {
        x = fget(a.*, xr, xc);
        y = fget(a.*, yr, yc);
        nx, ny = _gbbrd_apply_rotation(x, y, c, s);
        fset(a, xr, xc, nx);
        fset(a, yr, yc, ny);
        xr += x_row_inc;
        xc += x_col_inc;
        yr += y_row_inc;
        yc += y_col_inc;
    }
}

_gbbrd_rotate_vectorized_2d_strided :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$x_row: int,
    $$x_col: int,
    $$x_row_inc: int,
    $$x_col_inc: int,
    $$y_row: int,
    $$y_col: int,
    $$y_row_inc: int,
    $$y_col_inc: int,
    rwork: *$RWORK/ArrayType1D,
    $$c_index: int,
    work: *$WORK/ArrayType1D,
    $$s_index: int,
    $$cs_inc: int,
    $RT: Type,
    $MT: Type
) {
    #if CHECKS {
        #run {
            assert(A.data_type     == MT);
            assert(RWORK.data_type == RT);
            assert(WORK.data_type  == MT);
        }
        assert(n >= 0);
    }

    if n <= 0 then return;

    xr := x_row;
    xc := x_col;
    yr := y_row;
    yc := y_col;
    ci := c_index;
    si := s_index;
    x, y, nx, ny : MT = ---;
    ctemp : RT = ---;
    stemp : MT = ---;

    for 1..n {
        x = fget(a.*, xr, xc);
        y = fget(a.*, yr, yc);
        ctemp = fget(rwork.*, ci);
        stemp = fget(work.*, si);
        nx, ny = _gbbrd_apply_rotation(x, y, ctemp, stemp);
        fset(a, xr, xc, nx);
        fset(a, yr, yc, ny);
        xr += x_row_inc;
        xc += x_col_inc;
        yr += y_row_inc;
        yc += y_col_inc;
        ci += cs_inc;
        si += cs_inc;
    }
}

_gbbrd_rotate_cols :: (
    $$m: int,
    a: *$A/ArrayType2D,
    $$col1: int,
    $$col2: int,
    c: $RT,
    s: $MT
) {
    _gbbrd_rotate_2d_strided(m, a, 1, col1, 1, 0, 1, col2, 1, 0, c, s);
}

_gbbrd_rotate_rows :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$row1: int,
    $$row2: int,
    c: $RT,
    s: $MT
) {
    _gbbrd_rotate_2d_strided(n, a, row1, 1, 0, 1, row2, 1, 0, 1, c, s);
}

_gbbrd_scale_col :: (
    $$m: int,
    alpha: $MT,
    a: *$A/ArrayType2D,
    $$col: int
) {
    #if CHECKS {
        #run {
            assert(A.data_type == MT);
        }
        assert(m >= 0);
    }

    if m <= 0 || alpha == one(MT) then return;

    for i:1..m {
        fset(a, i, col, alpha*fget(a.*, i, col));
    }
}

_gbbrd_scale_row :: (
    $$n: int,
    alpha: $MT,
    a: *$A/ArrayType2D,
    $$row: int
) {
    #if CHECKS {
        #run {
            assert(A.data_type == MT);
        }
        assert(n >= 0);
    }

    if n <= 0 || alpha == one(MT) then return;

    for j:1..n {
        fset(a, row, j, alpha*fget(a.*, row, j));
    }
}

_gbbrd_lartg :: (f: $MT, g: MT, c: *$RT, s: *MT, r: *MT) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 {
        slartg(f, g, c, s, r);
    } else if MT == float64 {
        dlartg(f, g, c, s, r);
    } else if MT == C32 {
        clartg(f, g, c, s, r);
    } else {
        zlartg(f, g, c, s, r);
    }
}

_gbbrd_largv_from_ab_work :: (
    $$n: int,
    ab: *$AB/ArrayType2D,
    $$ab_row: int,
    $$ab_col: int,
    $$ab_col_inc: int,
    work: *$WORK/ArrayType1D,
    $$work_index: int,
    $$work_inc: int,
    rwork: *$RWORK/ArrayType1D,
    $$rwork_index: int,
    $$rwork_inc: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if CHECKS {
        #run {
            assert(AB.data_type    == MT);
            assert(WORK.data_type  == MT);
            assert(RWORK.data_type == RT);
        }
        assert(n >= 0);
    }

    if n <= 0 then return;

    col := ab_col;
    wi := work_index;
    ri := rwork_index;
    ctemp : RT = ---;
    stemp, rtemp : MT = ---;

    for 1..n {
        _gbbrd_lartg(fget(ab.*, ab_row, col), fget(work.*, wi), *ctemp, *stemp, *rtemp);
        fset(ab, ab_row, col, rtemp);
        fset(work, wi, stemp);
        fset(rwork, ri, ctemp);
        col += ab_col_inc;
        wi += work_inc;
        ri += rwork_inc;
    }
}

_gbbrd :: (
    $$vect: u8,
    $$m: int,
    $$n: int,
    $$ncc: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    pt: *$PT/ArrayType2D,
    $$ldpt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    $$work_offset: int,
    $$rwork_offset: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if CHECKS {
        #run {
            assert(AB.data_type    == MT);
            assert(D.data_type     == RT);
            assert(E.data_type     == RT);
            assert(Q.data_type     == MT);
            assert(PT.data_type    == MT);
            assert(C.data_type     == MT);
            assert(WORK.data_type  == MT);
            assert(RWORK.data_type == RT);
        }
    }

    zero_r :: #run zero(RT);
    one_t  :: #run one(MT);

    wantb := char_eq(vect, #char "b");
    wantq := char_eq(vect, #char "q") || wantb;
    wantpt := char_eq(vect, #char "p") || wantb;
    wantc := ncc > 0;
    klu1 := kl + ku + 1;
    routine_name :: #run if MT == float32 {
        "SGBBRD"
    } else if MT == float64 {
        "DGBBRD"
    } else if MT == C32 {
        "CGBBRD"
    } else {
        "ZGBBRD"
    };

    i, j, j1, j2, kb, kb1, kk, klm, kun, l, minmn, ml, ml0, mu, mu0, nr, nrt : int = ---;
    ra, rb, rs, t : MT = ---;
    rc, abst, ctemp : RT = ---;
    stemp : MT = ---;

    info.* = 0;
    if !wantq && !wantpt && !char_eq(vect, #char "n") {
        info.* = -1;
    } else if m < 0 {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ncc < 0 {
        info.* = -4;
    } else if kl < 0 {
        info.* = -5;
    } else if ku < 0 {
        info.* = -6;
    } else if ldab < klu1 {
        info.* = -8;
    } else if ldq < 1 || (wantq && ldq < max(1, m)) {
        info.* = -12;
    } else if ldpt < 1 || (wantpt && ldpt < max(1, n)) {
        info.* = -14;
    } else if ldc < 1 || (wantc && ldc < max(1, m)) {
        info.* = -16;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if wantq then _gbbrd_set_identity(m, m, q, m);
    if wantpt then _gbbrd_set_identity(n, n, pt, n);

    if m == 0 || n == 0 then return;

    minmn = min(m, n);

    if kl + ku > 1 {
        if ku > 0 {
            ml0 = 1;
            mu0 = 2;
        } else {
            ml0 = 2;
            mu0 = 1;
        }

        klm = min(m-1, kl);
        kun = min(n-1, ku);
        kb = klm + kun;
        kb1 = kb + 1;
        nr = 0;
        j1 = klm + 2;
        j2 = 1 - kun;

        for i:1..minmn {
            ml = klm + 1;
            mu = kun + 1;
            for kk:1..kb {
                j1 += kb;
                j2 += kb;

                if nr > 0 {
                    _gbbrd_largv_from_ab_work(
                        nr, ab, klu1, j1-klm-1, kb1,
                        work, work_offset + j1 - 1, kb1,
                        rwork, rwork_offset + j1 - 1, kb1,
                        RT, MT
                    );
                }

                for l:1..kb {
                    if j2-klm+l-1 > n {
                        nrt = nr - 1;
                    } else {
                        nrt = nr;
                    }
                    if nrt > 0 {
                        _gbbrd_rotate_vectorized_2d_strided(
                            nrt,
                            ab,
                            klu1-l,   j1-klm+l-1, 0, kb1,
                            klu1-l+1, j1-klm+l-1, 0, kb1,
                            rwork, rwork_offset + j1 - 1,
                            work,  work_offset + j1 - 1,
                            kb1,
                            RT, MT
                        );
                    }
                }

                if ml > ml0 {
                    if ml <= m-i+1 {
                        _gbbrd_lartg(
                            fget(ab.*, ku+ml-1, i),
                            fget(ab.*, ku+ml, i),
                            *ctemp, *stemp, *ra
                        );
                        fset(rwork, rwork_offset + i + ml - 2, ctemp);
                        fset(work,  work_offset + i + ml - 2, stemp);
                        fset(ab, ku+ml-1, i, ra);
                        if i < n {
                            _gbbrd_rotate_2d_strided(
                                min(ku+ml-2, n-i),
                                ab,
                                ku+ml-2, i+1, -1, 1,
                                ku+ml-1, i+1, -1, 1,
                                fget(rwork.*, rwork_offset + i + ml - 2),
                                fget(work.*,  work_offset + i + ml - 2)
                            );
                        }
                    }
                    nr += 1;
                    j1 -= kb1;
                }

                if wantq {
                    j = j1;
                    while j <= j2 {
                        #if MT == RT {
                            stemp = fget(work.*, work_offset + j - 1);
                        } else {
                            stemp = conjugate(fget(work.*, work_offset + j - 1));
                        }
                        _gbbrd_rotate_cols(
                            m, q, j-1, j,
                            fget(rwork.*, rwork_offset + j - 1),
                            stemp
                        );
                        j += kb1;
                    }
                }

                if wantc {
                    j = j1;
                    while j <= j2 {
                        _gbbrd_rotate_rows(
                            ncc, c, j-1, j,
                            fget(rwork.*, rwork_offset + j - 1),
                            fget(work.*, work_offset + j - 1)
                        );
                        j += kb1;
                    }
                }

                if j2+kun > n {
                    nr -= 1;
                    j2 -= kb1;
                }

                j = j1;
                while j <= j2 {
                    fset(work, work_offset + j + kun - 1,
                         fget(work.*, work_offset + j - 1) * fget(ab.*, 1, j+kun));
                    fset(ab, 1, j+kun,
                         fget(rwork.*, rwork_offset + j - 1) * fget(ab.*, 1, j+kun));
                    j += kb1;
                }

                if nr > 0 {
                    _gbbrd_largv_from_ab_work(
                        nr, ab, 1, j1+kun-1, kb1,
                        work, work_offset + j1 + kun - 1, kb1,
                        rwork, rwork_offset + j1 + kun - 1, kb1,
                        RT, MT
                    );
                }

                for l:1..kb {
                    if j2+l-1 > m {
                        nrt = nr - 1;
                    } else {
                        nrt = nr;
                    }
                    if nrt > 0 {
                        _gbbrd_rotate_vectorized_2d_strided(
                            nrt,
                            ab,
                            l+1, j1+kun-1, 0, kb1,
                            l,   j1+kun,   0, kb1,
                            rwork, rwork_offset + j1 + kun - 1,
                            work,  work_offset + j1 + kun - 1,
                            kb1,
                            RT, MT
                        );
                    }
                }

                if ml == ml0 && mu > mu0 {
                    if mu <= n-i+1 {
                        _gbbrd_lartg(
                            fget(ab.*, ku-mu+3, i+mu-2),
                            fget(ab.*, ku-mu+2, i+mu-1),
                            *ctemp, *stemp, *ra
                        );
                        fset(rwork, rwork_offset + i + mu - 2, ctemp);
                        fset(work,  work_offset + i + mu - 2, stemp);
                        fset(ab, ku-mu+3, i+mu-2, ra);
                        _gbbrd_rotate_2d_strided(
                            min(kl+mu-2, m-i),
                            ab,
                            ku-mu+4, i+mu-2, 1, 0,
                            ku-mu+3, i+mu-1, 1, 0,
                            fget(rwork.*, rwork_offset + i + mu - 2),
                            fget(work.*,  work_offset + i + mu - 2)
                        );
                    }
                    nr += 1;
                    j1 -= kb1;
                }

                if wantpt {
                    j = j1;
                    while j <= j2 {
                        #if MT == RT {
                            stemp = fget(work.*, work_offset + j + kun - 1);
                        } else {
                            stemp = conjugate(fget(work.*, work_offset + j + kun - 1));
                        }
                        _gbbrd_rotate_rows(
                            n, pt, j+kun-1, j+kun,
                            fget(rwork.*, rwork_offset + j + kun - 1),
                            stemp
                        );
                        j += kb1;
                    }
                }

                if j2+kb > m {
                    nr -= 1;
                    j2 -= kb1;
                }

                j = j1;
                while j <= j2 {
                    fset(work, work_offset + j + kb - 1,
                         fget(work.*, work_offset + j + kun - 1) * fget(ab.*, klu1, j+kun));
                    fset(ab, klu1, j+kun,
                         fget(rwork.*, rwork_offset + j + kun - 1) * fget(ab.*, klu1, j+kun));
                    j += kb1;
                }

                if ml > ml0 {
                    ml -= 1;
                } else {
                    mu -= 1;
                }
            }
        }
    }

    #if MT == RT {
        if ku == 0 && kl > 0 {
            for i:1..min(m-1, n) {
                _gbbrd_lartg(fget(ab.*, 1, i), fget(ab.*, 2, i), *rc, *rs, *ra);
                fset(d, i, ra);
                if i < n {
                    fset(e, i, rs * fget(ab.*, 1, i+1));
                    fset(ab, 1, i+1, rc * fget(ab.*, 1, i+1));
                }
                if wantq then _gbbrd_rotate_cols(m, q, i, i+1, rc, rs);
                if wantc then _gbbrd_rotate_rows(ncc, c, i, i+1, rc, rs);
            }
            if m <= n then fset(d, m, fget(ab.*, 1, m));
        } else if ku > 0 {
            if m < n {
                rb = fget(ab.*, ku, m+1);
                for <i:1..m {
                    _gbbrd_lartg(fget(ab.*, ku+1, i), rb, *rc, *rs, *ra);
                    fset(d, i, ra);
                    if i > 1 {
                        rb = -rs * fget(ab.*, ku, i);
                        fset(e, i-1, rc * fget(ab.*, ku, i));
                    }
                    if wantpt then _gbbrd_rotate_rows(n, pt, i, m+1, rc, rs);
                }
            } else {
                for i:1..minmn-1 {
                    fset(e, i, fget(ab.*, ku, i+1));
                }
                for i:1..minmn {
                    fset(d, i, fget(ab.*, ku+1, i));
                }
            }
        } else {
            for i:1..minmn-1 {
                fset(e, i, zero_r);
            }
            for i:1..minmn {
                fset(d, i, fget(ab.*, 1, i));
            }
        }
    } else {
        if ku == 0 && kl > 0 {
            for i:1..min(m-1, n) {
                _gbbrd_lartg(fget(ab.*, 1, i), fget(ab.*, 2, i), *rc, *rs, *ra);
                fset(ab, 1, i, ra);
                if i < n {
                    fset(ab, 2, i, rs * fget(ab.*, 1, i+1));
                    fset(ab, 1, i+1, rc * fget(ab.*, 1, i+1));
                }
                if wantq then _gbbrd_rotate_cols(m, q, i, i+1, rc, conjugate(rs));
                if wantc then _gbbrd_rotate_rows(ncc, c, i, i+1, rc, rs);
            }
        } else {
            if ku > 0 && m < n {
                rb = fget(ab.*, ku, m+1);
                for <i:1..m {
                    _gbbrd_lartg(fget(ab.*, ku+1, i), rb, *rc, *rs, *ra);
                    fset(ab, ku+1, i, ra);
                    if i > 1 {
                        rb = -conjugate(rs) * fget(ab.*, ku, i);
                        fset(ab, ku, i, rc * fget(ab.*, ku, i));
                    }
                    if wantpt then _gbbrd_rotate_rows(n, pt, i, m+1, rc, conjugate(rs));
                }
            }
        }

        t = fget(ab.*, ku+1, 1);
        for i:1..minmn {
            abst = abs(t);
            fset(d, i, abst);
            if abst != zero_r {
                t /= abst;
            } else {
                t = one_t;
            }
            if wantq then _gbbrd_scale_col(m, t, q, i);
            if wantc then _gbbrd_scale_row(ncc, conjugate(t), c, i);

            if i < minmn {
                if ku == 0 && kl == 0 {
                    fset(e, i, zero_r);
                    t = fget(ab.*, 1, i+1);
                } else {
                    if ku == 0 {
                        t = fget(ab.*, 2, i) * conjugate(t);
                    } else {
                        t = fget(ab.*, ku, i+1) * conjugate(t);
                    }
                    abst = abs(t);
                    fset(e, i, abst);
                    if abst != zero_r {
                        t /= abst;
                    } else {
                        t = one_t;
                    }
                    if wantpt then _gbbrd_scale_row(n, t, pt, i+1);
                    t = fget(ab.*, ku+1, i+1) * conjugate(t);
                }
            }
        }
    }
}
