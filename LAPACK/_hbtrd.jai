// ? LAPACK *hbtrd/*sbtrd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment ssbtrd/dsbtrd/chbtrd/zhbtrd
// *HBTRD reduces a symmetric/Hermitian band matrix to tridiagonal form.

ssbtrd :: #bake_arguments _hbtrd(RT=float32, MT=float32);
dsbtrd :: #bake_arguments _hbtrd(RT=float64, MT=float64);
chbtrd :: #bake_arguments _hbtrd(RT=float32, MT=C32);
zhbtrd :: #bake_arguments _hbtrd(RT=float64, MT=C64);

#scope_file
_hbtrd_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_hbtrd_vector_view :: inline (
    x: *$X/ArrayType1D,
    $$index: int,
    $$count: int
) -> [] X.data_type {
    return _hbtrd_raw_slice(x.data + (index - 1), count);
}

_hbtrd_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _hbtrd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hbtrd_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _hbtrd_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_hbtrd_matrix_strided :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int,
    $$inc: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*inc;
    return _hbtrd_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_hbtrd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_hbtrd_real_part :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return x;
    } else {
        return real(x);
    }
}

_hbtrd_conj_if_needed :: inline (x: $MT, $$do_conj: bool, $MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == float32 || MT == float64 {
        return x;
    } else {
        if do_conj then return conjugate(x);
        return x;
    }
}

_hbtrd_make_real_if_needed :: (
    a: *$A/ArrayType2D,
    $$row: int,
    $$col: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    #if MT == C32 || MT == C64 {
        fset(a, row, col, _hbtrd_complex_from_real(real(fget(a.*, row, col)), RT, MT));
    }
}

_hbtrd_laset :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    beta: MT,
    a: *$A/ArrayType2D,
    $$lda: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        slaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == float64 {
        dlaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_hbtrd_lartg :: (
    f: $MT,
    g: MT,
    c: *$RT,
    s: *MT,
    r: *MT,
    $RT: Type,
    $MT: Type
) #modify {
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

_hbtrd_rot :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    c: $RT,
    s: $MT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else if MT == float64 {
        drot(n, x, incx, y, incy, c, s);
    } else if MT == C32 {
        crot(n, x, incx, y, incy, c, s);
    } else {
        zrot(n, x, incx, y, incy, c, s);
    }
}

_hbtrd_largv :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    c: *$C/ArrayType1D,
    $$incc: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT &&
           C.data_type == RT;
} {
    #if MT == float32 {
        slargv(n, x, incx, y, incy, c, incc);
    } else if MT == float64 {
        dlargv(n, x, incx, y, incy, c, incc);
    } else if MT == C32 {
        clargv(n, x, incx, y, incy, c, incc);
    } else {
        zlargv(n, x, incx, y, incy, c, incc);
    }
}

_hbtrd_lartv :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    c: *$C/ArrayType1D,
    $$incc: int,
    s: *$S/ArrayType1D,
    $$incs: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT &&
           C.data_type == RT &&
           S.data_type == MT;
} {
    #if MT == float32 {
        slartv(n, x, incx, y, incy, c, incc, s, incs);
    } else if MT == float64 {
        dlartv(n, x, incx, y, incy, c, incc, s, incs);
    } else if MT == C32 {
        clartv(n, x, incx, y, incy, c, incc, s, incs);
    } else {
        zlartv(n, x, incx, y, incy, c, incc, s, incs);
    }
}

_hbtrd_lar2v :: (
    $$n: int,
    x: *$X/ArrayType1D,
    y: *$Y/ArrayType1D,
    z: *$Z/ArrayType1D,
    $$inc: int,
    c: *$C/ArrayType1D,
    $$incc: int,
    s: *$S/ArrayType1D,
    $$incs: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT &&
           Z.data_type == MT &&
           C.data_type == RT &&
           S.data_type == MT;
} {
    #if MT == float32 {
        slar2v(n, x, y, z, inc, c, incc, s, incs);
    } else if MT == float64 {
        dlar2v(n, x, y, z, inc, c, incc, s, incs);
    } else if MT == C32 {
        clar2v(n, x, y, z, inc, c, incc, s, incs);
    } else {
        zlar2v(n, x, y, z, inc, c, incc, s, incs);
    }
}

_hbtrd_lacgv :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        clacgv(n, x, incx);
    } else if MT == C64 {
        zlacgv(n, x, incx);
    }
}

_hbtrd_scal :: (
    $$n: int,
    alpha: $MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, incx);
    } else if MT == float64 {
        dscal(n, alpha, x, incx);
    } else if MT == C32 {
        cscal(n, alpha, x, incx);
    } else {
        zscal(n, alpha, x, incx);
    }
}

_hbtrd :: (
    $$vect: u8,
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           Q.data_type == MT &&
           WORK.data_type == MT;
} {
    zero_r :: #run zero(RT);
    zero_t :: #run zero(MT);
    one_t  :: #run one(MT);
    routine_name :: #run if MT == float32 {
        "SSBTRD"
    } else if MT == float64 {
        "DSBTRD"
    } else if MT == C32 {
        "CHBTRD"
    } else {
        "ZHBTRD"
    };

    initq, upper, wantq : bool = ---;
    i, i2, ibl, inca, incx, iqaend, iqb, iqend, j, j1, j1end, j1inc, j2, jend,
    jin, jinc, k, kd1, kdm1, kdn, l, last, lend, nq, nr, nrt : int = ---;
    abst, ctemp : RT = ---;
    t, temp, stemp : MT = ---;

    initq = lsame(vect, #char "V");
    wantq = initq || lsame(vect, #char "U");
    upper = lsame(uplo, #char "U");
    kd1 = kd + 1;
    kdm1 = kd - 1;
    incx = ldab - 1;
    iqend = 1;

    info.* = 0;
    if !(wantq || lsame(vect, #char "N")) {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if kd < 0 {
        info.* = -4;
    } else if ldab < kd1 {
        info.* = -6;
    } else if ldq < max(1, n) && wantq {
        info.* = -10;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    if initq {
        _hbtrd_laset(#char "F", n, n, zero_t, one_t, q, ldq, MT);
    }

    inca = kd1*ldab;
    kdn = min(n - 1, kd);

    if upper {
        if kd > 1 {
            nr = 0;
            j1 = kdn + 2;
            j2 = 1;

            _hbtrd_make_real_if_needed(ab, kd1, 1, RT, MT);

            for i:1..n - 2 {
                for <k:kdn + 1..2 {
                    j1 += kdn;
                    j2 += kdn;

                    if nr > 0 {
                        xview := _hbtrd_matrix_column(ab.*, ldab, 1, j1 - 1,
                                                     1 + (nr - 1)*inca);
                        yview := _hbtrd_vector_view(work, j1, 1 + (nr - 1)*kd1);
                        cview := _hbtrd_vector_view(d, j1, 1 + (nr - 1)*kd1);
                        _hbtrd_largv(nr, *xview, inca, *yview, kd1, *cview, kd1, RT, MT);

                        if nr >= 2*kd - 1 {
                            for l:1..kd - 1 {
                                xrot := _hbtrd_matrix_column(ab.*, ldab, l + 1, j1 - 1,
                                                             1 + (nr - 1)*inca);
                                yrot := _hbtrd_matrix_column(ab.*, ldab, l, j1,
                                                             1 + (nr - 1)*inca);
                                cview := _hbtrd_vector_view(d, j1, 1 + (nr - 1)*kd1);
                                sview := _hbtrd_vector_view(work, j1, 1 + (nr - 1)*kd1);
                                _hbtrd_lartv(nr, *xrot, inca, *yrot, inca, *cview, kd1,
                                             *sview, kd1, RT, MT);
                            }
                        } else {
                            jend = j1 + (nr - 1)*kd1;
                            jinc = j1;
                            while jinc <= jend {
                                xrot := _hbtrd_matrix_column(ab.*, ldab, 2, jinc - 1, kdm1);
                                yrot := _hbtrd_matrix_column(ab.*, ldab, 1, jinc, kdm1);
                                _hbtrd_rot(kdm1, *xrot, 1, *yrot, 1,
                                           fget(d.*, jinc), fget(work.*, jinc), RT, MT);
                                jinc += kd1;
                            }
                        }
                    }

                    if k > 2 {
                        if k <= n - i + 1 {
                            _hbtrd_lartg(
                                fget(ab.*, kd - k + 3, i + k - 2),
                                fget(ab.*, kd - k + 2, i + k - 1),
                                *ctemp, *stemp, *temp, RT, MT
                            );
                            fset(d, i + k - 1, ctemp);
                            fset(work, i + k - 1, stemp);
                            fset(ab, kd - k + 3, i + k - 2, temp);

                            xrot := _hbtrd_matrix_column(ab.*, ldab, kd - k + 4, i + k - 2,
                                                         max(0, k - 3));
                            yrot := _hbtrd_matrix_column(ab.*, ldab, kd - k + 3, i + k - 1,
                                                         max(0, k - 3));
                            _hbtrd_rot(k - 3, *xrot, 1, *yrot, 1,
                                       fget(d.*, i + k - 1), fget(work.*, i + k - 1), RT, MT);
                        }
                        nr += 1;
                        j1 -= kdn + 1;
                    }

                    if nr > 0 {
                        xdiag := _hbtrd_matrix_column(ab.*, ldab, kd1, j1 - 1,
                                                      1 + (nr - 1)*inca);
                        ydiag := _hbtrd_matrix_column(ab.*, ldab, kd1, j1,
                                                      1 + (nr - 1)*inca);
                        zdiag := _hbtrd_matrix_column(ab.*, ldab, kd, j1,
                                                      1 + (nr - 1)*inca);
                        cview := _hbtrd_vector_view(d, j1, 1 + (nr - 1)*kd1);
                        sview := _hbtrd_vector_view(work, j1, 1 + (nr - 1)*kd1);
                        _hbtrd_lar2v(nr, *xdiag, *ydiag, *zdiag, inca, *cview, kd1,
                                     *sview, kd1, RT, MT);
                    }

                    if nr > 0 {
                        #if MT == C32 || MT == C64 {
                            sview := _hbtrd_vector_view(work, j1, 1 + (nr - 1)*kd1);
                            _hbtrd_lacgv(nr, *sview, kd1, MT);
                        }

                        if 2*kd - 1 < nr {
                            for l:1..kd - 1 {
                                if j2 + l > n {
                                    nrt = nr - 1;
                                } else {
                                    nrt = nr;
                                }
                                if nrt > 0 {
                                    xrot := _hbtrd_matrix_column(ab.*, ldab, kd - l, j1 + l,
                                                                 1 + (nrt - 1)*inca);
                                    yrot := _hbtrd_matrix_column(ab.*, ldab, kd - l + 1, j1 + l,
                                                                 1 + (nrt - 1)*inca);
                                    cview := _hbtrd_vector_view(d, j1, 1 + (nrt - 1)*kd1);
                                    sview := _hbtrd_vector_view(work, j1, 1 + (nrt - 1)*kd1);
                                    _hbtrd_lartv(nrt, *xrot, inca, *yrot, inca, *cview, kd1,
                                                 *sview, kd1, RT, MT);
                                }
                            }
                        } else {
                            j1end = j1 + kd1*(nr - 2);
                            if j1end >= j1 {
                                jin = j1;
                                while jin <= j1end {
                                    xrot := _hbtrd_matrix_strided(ab.*, ldab, kd - 1, jin + 1, kd - 1,
                                                                  incx);
                                    yrot := _hbtrd_matrix_strided(ab.*, ldab, kd, jin + 1, kd - 1,
                                                                  incx);
                                    _hbtrd_rot(kd - 1, *xrot, incx, *yrot, incx,
                                               fget(d.*, jin), fget(work.*, jin), RT, MT);
                                    jin += kd1;
                                }
                            }
                            lend = min(kdm1, n - j2);
                            last = j1end + kd1;
                            if lend > 0 {
                                xrot := _hbtrd_matrix_strided(ab.*, ldab, kd - 1, last + 1, lend,
                                                              incx);
                                yrot := _hbtrd_matrix_strided(ab.*, ldab, kd, last + 1, lend,
                                                              incx);
                                _hbtrd_rot(lend, *xrot, incx, *yrot, incx,
                                           fget(d.*, last), fget(work.*, last), RT, MT);
                            }
                        }
                    }

                    if wantq {
                        if initq {
                            iqend = max(iqend, j2);
                            i2 = max(0, k - 3);
                            iqaend = 1 + i*kd;
                            if k == 2 then iqaend += kd;
                            iqaend = min(iqaend, iqend);

                            j = j1;
                            while j <= j2 {
                                ibl = i - i2 / kdm1;
                                i2 += 1;
                                iqb = max(1, j - ibl);
                                nq = 1 + iqaend - iqb;
                                iqaend = min(iqaend + kd, iqend);
                                xrot := _hbtrd_matrix_column(q.*, ldq, iqb, j - 1, nq);
                                yrot := _hbtrd_matrix_column(q.*, ldq, iqb, j, nq);
                                _hbtrd_rot(nq, *xrot, 1, *yrot, 1, fget(d.*, j),
                                           _hbtrd_conj_if_needed(fget(work.*, j), upper, MT),
                                           RT, MT);
                                j += kd1;
                            }
                        } else {
                            j = j1;
                            while j <= j2 {
                                xrot := _hbtrd_matrix_column(q.*, ldq, 1, j - 1, n);
                                yrot := _hbtrd_matrix_column(q.*, ldq, 1, j, n);
                                _hbtrd_rot(n, *xrot, 1, *yrot, 1, fget(d.*, j),
                                           _hbtrd_conj_if_needed(fget(work.*, j), upper, MT),
                                           RT, MT);
                                j += kd1;
                            }
                        }
                    }

                    if j2 + kdn > n {
                        nr -= 1;
                        j2 -= kdn + 1;
                    }

                    j = j1;
                    while j <= j2 {
                        fset(work, j + kd, fget(work.*, j) * fget(ab.*, 1, j + kd));
                        fset(ab, 1, j + kd, _hbtrd_complex_from_real(fget(d.*, j), RT, MT) *
                                             fget(ab.*, 1, j + kd));
                        j += kd1;
                    }
                }
            }
        }

        if kd > 0 {
            #if MT == C32 || MT == C64 {
                for i:1..n - 1 {
                    t = fget(ab.*, kd, i + 1);
                    abst = abs(t);
                    fset(ab, kd, i + 1, _hbtrd_complex_from_real(abst, RT, MT));
                    fset(e, i, abst);
                    if abst != zero_r {
                        t /= abst;
                    } else {
                        t = one_t;
                    }
                    if i < n - 1 {
                        fset(ab, kd, i + 2, fget(ab.*, kd, i + 2)*t);
                    }
                    if wantq {
                        qcol := _hbtrd_matrix_column(q.*, ldq, 1, i + 1, n);
                        _hbtrd_scal(n, conjugate(t), *qcol, 1, MT);
                    }
                }
            } else {
                for i:1..n - 1 {
                    fset(e, i, fget(ab.*, kd, i + 1));
                }
            }
        } else {
            for i:1..n - 1 {
                fset(e, i, zero_r);
            }
        }

        for i:1..n {
            fset(d, i, _hbtrd_real_part(fget(ab.*, kd1, i), RT, MT));
        }
    } else {
        if kd > 1 {
            nr = 0;
            j1 = kdn + 2;
            j2 = 1;

            _hbtrd_make_real_if_needed(ab, 1, 1, RT, MT);

            for i:1..n - 2 {
                for <k:kdn + 1..2 {
                    j1 += kdn;
                    j2 += kdn;

                    if nr > 0 {
                        xview := _hbtrd_matrix_column(ab.*, ldab, kd1, j1 - kd1,
                                                     1 + (nr - 1)*inca);
                        yview := _hbtrd_vector_view(work, j1, 1 + (nr - 1)*kd1);
                        cview := _hbtrd_vector_view(d, j1, 1 + (nr - 1)*kd1);
                        _hbtrd_largv(nr, *xview, inca, *yview, kd1, *cview, kd1, RT, MT);

                        if nr > 2*kd - 1 {
                            for l:1..kd - 1 {
                                xrot := _hbtrd_matrix_column(ab.*, ldab, kd1 - l, j1 - kd1 + l,
                                                             1 + (nr - 1)*inca);
                                yrot := _hbtrd_matrix_column(ab.*, ldab, kd1 - l + 1, j1 - kd1 + l,
                                                             1 + (nr - 1)*inca);
                                cview := _hbtrd_vector_view(d, j1, 1 + (nr - 1)*kd1);
                                sview := _hbtrd_vector_view(work, j1, 1 + (nr - 1)*kd1);
                                _hbtrd_lartv(nr, *xrot, inca, *yrot, inca, *cview, kd1,
                                             *sview, kd1, RT, MT);
                            }
                        } else {
                            jend = j1 + kd1*(nr - 1);
                            jinc = j1;
                            while jinc <= jend {
                                xrot := _hbtrd_matrix_strided(ab.*, ldab, kd, jinc - kd, kdm1,
                                                              incx);
                                yrot := _hbtrd_matrix_strided(ab.*, ldab, kd1, jinc - kd, kdm1,
                                                              incx);
                                _hbtrd_rot(kdm1, *xrot, incx, *yrot, incx,
                                           fget(d.*, jinc), fget(work.*, jinc), RT, MT);
                                jinc += kd1;
                            }
                        }
                    }

                    if k > 2 {
                        if k <= n - i + 1 {
                            _hbtrd_lartg(
                                fget(ab.*, k - 1, i),
                                fget(ab.*, k, i),
                                *ctemp, *stemp, *temp, RT, MT
                            );
                            fset(d, i + k - 1, ctemp);
                            fset(work, i + k - 1, stemp);
                            fset(ab, k - 1, i, temp);

                            xrot := _hbtrd_matrix_strided(ab.*, ldab, k - 2, i + 1,
                                                          max(0, k - 3), ldab - 1);
                            yrot := _hbtrd_matrix_strided(ab.*, ldab, k - 1, i + 1,
                                                          max(0, k - 3), ldab - 1);
                            _hbtrd_rot(k - 3, *xrot, ldab - 1, *yrot, ldab - 1,
                                       fget(d.*, i + k - 1), fget(work.*, i + k - 1), RT, MT);
                        }
                        nr += 1;
                        j1 -= kdn + 1;
                    }

                    if nr > 0 {
                        xdiag := _hbtrd_matrix_column(ab.*, ldab, 1, j1 - 1,
                                                      1 + (nr - 1)*inca);
                        ydiag := _hbtrd_matrix_column(ab.*, ldab, 1, j1,
                                                      1 + (nr - 1)*inca);
                        zdiag := _hbtrd_matrix_column(ab.*, ldab, 2, j1 - 1,
                                                      1 + (nr - 1)*inca);
                        cview := _hbtrd_vector_view(d, j1, 1 + (nr - 1)*kd1);
                        sview := _hbtrd_vector_view(work, j1, 1 + (nr - 1)*kd1);
                        _hbtrd_lar2v(nr, *xdiag, *ydiag, *zdiag, inca, *cview, kd1,
                                     *sview, kd1, RT, MT);
                    }

                    if nr > 0 {
                        #if MT == C32 || MT == C64 {
                            sview := _hbtrd_vector_view(work, j1, 1 + (nr - 1)*kd1);
                            _hbtrd_lacgv(nr, *sview, kd1, MT);
                        }

                        if nr > 2*kd - 1 {
                            for l:1..kd - 1 {
                                if j2 + l > n {
                                    nrt = nr - 1;
                                } else {
                                    nrt = nr;
                                }
                                if nrt > 0 {
                                    xrot := _hbtrd_matrix_column(ab.*, ldab, l + 2, j1 - 1,
                                                                 1 + (nrt - 1)*inca);
                                    yrot := _hbtrd_matrix_column(ab.*, ldab, l + 1, j1,
                                                                 1 + (nrt - 1)*inca);
                                    cview := _hbtrd_vector_view(d, j1, 1 + (nrt - 1)*kd1);
                                    sview := _hbtrd_vector_view(work, j1, 1 + (nrt - 1)*kd1);
                                    _hbtrd_lartv(nrt, *xrot, inca, *yrot, inca, *cview, kd1,
                                                 *sview, kd1, RT, MT);
                                }
                            }
                        } else {
                            j1end = j1 + kd1*(nr - 2);
                            if j1end >= j1 {
                                j1inc = j1;
                                while j1inc <= j1end {
                                    xrot := _hbtrd_matrix_column(ab.*, ldab, 3, j1inc - 1, kdm1);
                                    yrot := _hbtrd_matrix_column(ab.*, ldab, 2, j1inc, kdm1);
                                    _hbtrd_rot(kdm1, *xrot, 1, *yrot, 1,
                                               fget(d.*, j1inc), fget(work.*, j1inc), RT, MT);
                                    j1inc += kd1;
                                }
                            }
                            lend = min(kdm1, n - j2);
                            last = j1end + kd1;
                            if lend > 0 {
                                xrot := _hbtrd_matrix_column(ab.*, ldab, 3, last - 1, lend);
                                yrot := _hbtrd_matrix_column(ab.*, ldab, 2, last, lend);
                                _hbtrd_rot(lend, *xrot, 1, *yrot, 1,
                                           fget(d.*, last), fget(work.*, last), RT, MT);
                            }
                        }
                    }

                    if wantq {
                        if initq {
                            iqend = max(iqend, j2);
                            i2 = max(0, k - 3);
                            iqaend = 1 + i*kd;
                            if k == 2 then iqaend += kd;
                            iqaend = min(iqaend, iqend);

                            j = j1;
                            while j <= j2 {
                                ibl = i - i2 / kdm1;
                                i2 += 1;
                                iqb = max(1, j - ibl);
                                nq = 1 + iqaend - iqb;
                                iqaend = min(iqaend + kd, iqend);
                                xrot := _hbtrd_matrix_column(q.*, ldq, iqb, j - 1, nq);
                                yrot := _hbtrd_matrix_column(q.*, ldq, iqb, j, nq);
                                _hbtrd_rot(nq, *xrot, 1, *yrot, 1, fget(d.*, j),
                                           _hbtrd_conj_if_needed(fget(work.*, j), false, MT),
                                           RT, MT);
                                j += kd1;
                            }
                        } else {
                            j = j1;
                            while j <= j2 {
                                xrot := _hbtrd_matrix_column(q.*, ldq, 1, j - 1, n);
                                yrot := _hbtrd_matrix_column(q.*, ldq, 1, j, n);
                                _hbtrd_rot(n, *xrot, 1, *yrot, 1, fget(d.*, j),
                                           _hbtrd_conj_if_needed(fget(work.*, j), false, MT),
                                           RT, MT);
                                j += kd1;
                            }
                        }
                    }

                    if j2 + kdn > n {
                        nr -= 1;
                        j2 -= kdn + 1;
                    }

                    j = j1;
                    while j <= j2 {
                        fset(work, j + kd, fget(work.*, j) * fget(ab.*, kd1, j));
                        fset(ab, kd1, j, _hbtrd_complex_from_real(fget(d.*, j), RT, MT) *
                                           fget(ab.*, kd1, j));
                        j += kd1;
                    }
                }
            }
        }

        if kd > 0 {
            #if MT == C32 || MT == C64 {
                for i:1..n - 1 {
                    t = fget(ab.*, 2, i);
                    abst = abs(t);
                    fset(ab, 2, i, _hbtrd_complex_from_real(abst, RT, MT));
                    fset(e, i, abst);
                    if abst != zero_r {
                        t /= abst;
                    } else {
                        t = one_t;
                    }
                    if i < n - 1 {
                        fset(ab, 2, i + 1, fget(ab.*, 2, i + 1)*t);
                    }
                    if wantq {
                        qcol := _hbtrd_matrix_column(q.*, ldq, 1, i + 1, n);
                        _hbtrd_scal(n, t, *qcol, 1, MT);
                    }
                }
            } else {
                for i:1..n - 1 {
                    fset(e, i, fget(ab.*, 2, i));
                }
            }
        } else {
            for i:1..n - 1 {
                fset(e, i, zero_r);
            }
        }

        for i:1..n {
            fset(d, i, _hbtrd_real_part(fget(ab.*, 1, i), RT, MT));
        }
    }
}
