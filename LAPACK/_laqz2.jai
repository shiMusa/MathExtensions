// ? LAPACK *laqz2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqz2/dlaqz2/claqz2/zlaqz2
// *LAQZ2 provides the AED / bulge-chasing helper used by QZ.

slaqz2 :: #bake_arguments _laqz2_real(RT=float32);
dlaqz2 :: #bake_arguments _laqz2_real(RT=float64);
claqz2 :: #bake_arguments _laqz2_complex(RT=float32, MT=C32);
zlaqz2 :: #bake_arguments _laqz2_complex(RT=float64, MT=C64);

#scope_file
_laqz2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_laqz2_storage :: inline (data: *$T, $$count: int, $$inc: int) -> [] T {
    storage_count := ifx count <= 0 then 0 else 1 + (count - 1)*inc;
    return _laqz2_raw_slice(data, storage_count);
}

_laqz2_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _laqz2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqz2_work_storage :: inline (
    work: $W/ArrayType1D, $$offset: int, $$rows: int, $$cols: int, $$ldw: int
) -> [] W.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else ldw*(cols - 1) + rows;
    return _laqz2_raw_slice(work.data + (offset - 1), count);
}

_laqz2_lartg_real :: (f: $RT, g: RT, c: *RT, s: *RT, r: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        slartg(f, g, c, s, r);
    } else {
        dlartg(f, g, c, s, r);
    }
}

_laqz2_rot_real :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $$c: $RT, $$s: RT, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else {
        drot(n, x, incx, y, incy, c, s);
    }
}

_laqz2_lartg_complex :: (f: $MT, g: MT, c: *$RT, s: *MT, r: *MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        clartg(f, g, c, s, r);
    } else {
        zlartg(f, g, c, s, r);
    }
}

_laqz2_rot_complex :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $$c: $RT, s: $MT, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        crot(n, x, incx, y, incy, c, s);
    } else {
        zrot(n, x, incx, y, incy, c, s);
    }
}

_laqz2_lamch :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch("Safe minimum".data);
    } else {
        return dlamch("Safe minimum".data);
    }
}

_laqz2_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} { fset(work, 1, cpx(x, zero(RT))); }

_laqz2_lacpy_complex :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_laqz2_laset_complex :: (
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

_laqz2_gemm_complex :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int, alpha: $MT,
    a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int, beta: MT,
    c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_laqz2_tgexc_complex :: (
    $$wantq: bool, $$wantz: bool, $$n: int, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int, ifst: *int, ilst: *int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT && Q.data_type == MT && Z.data_type == MT;
} {
    #if MT == C32 {
        ctgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info);
    } else {
        ztgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info);
    }
}

_laqz2_real :: (
    $$ilq: bool, $$ilz: bool, $$k: int, $$istartm: int, $$istopm: int, $$ihi: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    $$nq: int, $$qstart: int, q: *$Q/ArrayType2D, $$ldq: int,
    $$nz: int, $$zstart: int, z: *$Z/ArrayType2D, $$ldz: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT && Q.data_type == RT && Z.data_type == RT;
} {
    zero_r :: #run zero(RT);

    c1, s1, c2, s2, temp : RT = ---;
    h11, h21, h12, h22, h13, h23 : RT = ---;

    if k + 2 == ihi {
        h11 = fget(b.*, ihi - 1, ihi - 2); h21 = fget(b.*, ihi, ihi - 2);
        h12 = fget(b.*, ihi - 1, ihi - 1); h22 = fget(b.*, ihi, ihi - 1);
        h13 = fget(b.*, ihi - 1, ihi);     h23 = fget(b.*, ihi, ihi);

        _laqz2_lartg_real(h11, h21, c1, s1, temp, RT);
        h21 = zero_r;
        h11 = temp;
        temp = c1*h12 + s1*h22; h22 = c1*h22 - s1*h12; h12 = temp;
        temp = c1*h13 + s1*h23; h23 = c1*h23 - s1*h13; h13 = temp;

        _laqz2_lartg_real(h23, h22, c1, s1, temp, RT);
        temp = c1*h13 + s1*h12; h12 = c1*h12 - s1*h13; h13 = temp;
        _laqz2_lartg_real(h12, h11, c2, s2, temp, RT);

        _laqz2_rot_real(ihi - istartm + 1, _laqz2_storage(b.data + (ihi - 1)*ldb + (istartm - 1), ihi - istartm + 1, 1), 1,
                                            _laqz2_storage(b.data + (ihi - 2)*ldb + (istartm - 1), ihi - istartm + 1, 1), 1, c1, s1, RT);
        _laqz2_rot_real(ihi - istartm + 1, _laqz2_storage(b.data + (ihi - 2)*ldb + (istartm - 1), ihi - istartm + 1, 1), 1,
                                            _laqz2_storage(b.data + (ihi - 3)*ldb + (istartm - 1), ihi - istartm + 1, 1), 1, c2, s2, RT);
        fset(b, ihi - 1, ihi - 2, zero_r);
        fset(b, ihi, ihi - 2, zero_r);
        _laqz2_rot_real(ihi - istartm + 1, _laqz2_storage(a.data + (ihi - 1)*lda + (istartm - 1), ihi - istartm + 1, 1), 1,
                                            _laqz2_storage(a.data + (ihi - 2)*lda + (istartm - 1), ihi - istartm + 1, 1), 1, c1, s1, RT);
        _laqz2_rot_real(ihi - istartm + 1, _laqz2_storage(a.data + (ihi - 2)*lda + (istartm - 1), ihi - istartm + 1, 1), 1,
                                            _laqz2_storage(a.data + (ihi - 3)*lda + (istartm - 1), ihi - istartm + 1, 1), 1, c2, s2, RT);
        if ilz {
            _laqz2_rot_real(nz, _laqz2_storage(z.data + (ihi - zstart)*ldz, nz, 1), 1,
                                _laqz2_storage(z.data + (ihi - 1 - zstart)*ldz, nz, 1), 1, c1, s1, RT);
            _laqz2_rot_real(nz, _laqz2_storage(z.data + (ihi - 1 - zstart)*ldz, nz, 1), 1,
                                _laqz2_storage(z.data + (ihi - 2 - zstart)*ldz, nz, 1), 1, c2, s2, RT);
        }

        _laqz2_lartg_real(fget(a.*, ihi - 1, ihi - 2), fget(a.*, ihi, ihi - 2), c1, s1, temp, RT);
        fset(a, ihi - 1, ihi - 2, temp);
        fset(a, ihi, ihi - 2, zero_r);
        _laqz2_rot_real(istopm - ihi + 2, _laqz2_storage(a.data + (ihi - 2)*lda + (ihi - 2), istopm - ihi + 2, lda), lda,
                                            _laqz2_storage(a.data + (ihi - 1)*lda + (ihi - 2), istopm - ihi + 2, lda), lda, c1, s1, RT);
        _laqz2_rot_real(istopm - ihi + 2, _laqz2_storage(b.data + (ihi - 2)*ldb + (ihi - 2), istopm - ihi + 2, ldb), ldb,
                                            _laqz2_storage(b.data + (ihi - 1)*ldb + (ihi - 2), istopm - ihi + 2, ldb), ldb, c1, s1, RT);
        if ilq {
            _laqz2_rot_real(nq, _laqz2_storage(q.data + (ihi - 1 - qstart)*ldq, nq, 1), 1,
                                _laqz2_storage(q.data + (ihi - qstart)*ldq, nq, 1), 1, c1, s1, RT);
        }

        _laqz2_lartg_real(fget(b.*, ihi, ihi), fget(b.*, ihi, ihi - 1), c1, s1, temp, RT);
        fset(b, ihi, ihi, temp);
        fset(b, ihi, ihi - 1, zero_r);
        _laqz2_rot_real(ihi - istartm, _laqz2_storage(b.data + (ihi - 1)*ldb + (istartm - 1), ihi - istartm, 1), 1,
                                        _laqz2_storage(b.data + (ihi - 2)*ldb + (istartm - 1), ihi - istartm, 1), 1, c1, s1, RT);
        _laqz2_rot_real(ihi - istartm + 1, _laqz2_storage(a.data + (ihi - 1)*lda + (istartm - 1), ihi - istartm + 1, 1), 1,
                                            _laqz2_storage(a.data + (ihi - 2)*lda + (istartm - 1), ihi - istartm + 1, 1), 1, c1, s1, RT);
        if ilz {
            _laqz2_rot_real(nz, _laqz2_storage(z.data + (ihi - zstart)*ldz, nz, 1), 1,
                                _laqz2_storage(z.data + (ihi - 1 - zstart)*ldz, nz, 1), 1, c1, s1, RT);
        }
    } else {
        h11 = fget(b.*, k + 1, k); h21 = fget(b.*, k + 2, k);
        h12 = fget(b.*, k + 1, k + 1); h22 = fget(b.*, k + 2, k + 1);
        h13 = fget(b.*, k + 1, k + 2); h23 = fget(b.*, k + 2, k + 2);

        _laqz2_lartg_real(h11, h21, c1, s1, temp, RT);
        h21 = zero_r;
        h11 = temp;
        temp = c1*h12 + s1*h22; h22 = c1*h22 - s1*h12; h12 = temp;
        temp = c1*h13 + s1*h23; h23 = c1*h23 - s1*h13; h13 = temp;

        _laqz2_lartg_real(h23, h22, c1, s1, temp, RT);
        temp = c1*h13 + s1*h12; h12 = c1*h12 - s1*h13; h13 = temp;
        _laqz2_lartg_real(h12, h11, c2, s2, temp, RT);

        _laqz2_rot_real(k + 3 - istartm + 1, _laqz2_storage(a.data + (k + 1)*lda + (istartm - 1), k + 3 - istartm + 1, 1), 1,
                                            _laqz2_storage(a.data + k*lda + (istartm - 1), k + 3 - istartm + 1, 1), 1, c1, s1, RT);
        _laqz2_rot_real(k + 3 - istartm + 1, _laqz2_storage(a.data + k*lda + (istartm - 1), k + 3 - istartm + 1, 1), 1,
                                            _laqz2_storage(a.data + (k - 1)*lda + (istartm - 1), k + 3 - istartm + 1, 1), 1, c2, s2, RT);
        _laqz2_rot_real(k + 2 - istartm + 1, _laqz2_storage(b.data + (k + 1)*ldb + (istartm - 1), k + 2 - istartm + 1, 1), 1,
                                            _laqz2_storage(b.data + k*ldb + (istartm - 1), k + 2 - istartm + 1, 1), 1, c1, s1, RT);
        _laqz2_rot_real(k + 2 - istartm + 1, _laqz2_storage(b.data + k*ldb + (istartm - 1), k + 2 - istartm + 1, 1), 1,
                                            _laqz2_storage(b.data + (k - 1)*ldb + (istartm - 1), k + 2 - istartm + 1, 1), 1, c2, s2, RT);
        if ilz {
            _laqz2_rot_real(nz, _laqz2_storage(z.data + (k + 2 - zstart)*ldz, nz, 1), 1,
                                _laqz2_storage(z.data + (k + 1 - zstart)*ldz, nz, 1), 1, c1, s1, RT);
            _laqz2_rot_real(nz, _laqz2_storage(z.data + (k + 1 - zstart)*ldz, nz, 1), 1,
                                _laqz2_storage(z.data + (k - zstart)*ldz, nz, 1), 1, c2, s2, RT);
        }
        fset(b, k + 1, k, zero_r);
        fset(b, k + 2, k, zero_r);

        _laqz2_lartg_real(fget(a.*, k + 2, k), fget(a.*, k + 3, k), c1, s1, temp, RT);
        fset(a, k + 2, k, temp);
        fset(a, k + 3, k, zero_r);
        _laqz2_lartg_real(fget(a.*, k + 1, k), fget(a.*, k + 2, k), c2, s2, temp, RT);
        fset(a, k + 1, k, temp);
        fset(a, k + 2, k, zero_r);

        _laqz2_rot_real(istopm - k, _laqz2_storage(a.data + k*lda + k, istopm - k, lda), lda,
                                        _laqz2_storage(a.data + (k + 1)*lda + k, istopm - k, lda), lda, c1, s1, RT);
        _laqz2_rot_real(istopm - k, _laqz2_storage(a.data + (k - 1)*lda + k, istopm - k, lda), lda,
                                        _laqz2_storage(a.data + k*lda + k, istopm - k, lda), lda, c2, s2, RT);
        _laqz2_rot_real(istopm - k, _laqz2_storage(b.data + k*ldb + k, istopm - k, ldb), ldb,
                                        _laqz2_storage(b.data + (k + 1)*ldb + k, istopm - k, ldb), ldb, c1, s1, RT);
        _laqz2_rot_real(istopm - k, _laqz2_storage(b.data + (k - 1)*ldb + k, istopm - k, ldb), ldb,
                                        _laqz2_storage(b.data + k*ldb + k, istopm - k, ldb), ldb, c2, s2, RT);
        if ilq {
            _laqz2_rot_real(nq, _laqz2_storage(q.data + (k + 1 - qstart)*ldq, nq, 1), 1,
                                _laqz2_storage(q.data + (k + 2 - qstart)*ldq, nq, 1), 1, c1, s1, RT);
            _laqz2_rot_real(nq, _laqz2_storage(q.data + (k - qstart)*ldq, nq, 1), 1,
                                _laqz2_storage(q.data + (k + 1 - qstart)*ldq, nq, 1), 1, c2, s2, RT);
        }
    }
}

_laqz2_complex :: (
    $$ilschur: bool, $$ilq: bool, $$ilz: bool, $$n: int, $$ilo: int, $$ihi: int, $$nw: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    ns: *int, nd: *int, alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    qc: *$QC/ArrayType2D, $$ldqc: int, zc: *$ZC/ArrayType2D, $$ldzc: int,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$rec: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && B.data_type == MT && Q.data_type == MT && Z.data_type == MT &&
           ALPHA.data_type == MT && BETA.data_type == MT && QC.data_type == MT && ZC.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    zero_t :: #run zero(MT);
    one_t  :: #run cpx(one(RT), zero_r);

    c1 : RT = ---;
    ctgexc_info, ifst, ilst, istartm, istopm, jw, k, k2, kwbot, kwtop, lworkreq, qz_small_info : int = ---;
    s, s1, temp : MT = ---;
    safmax, safmin, smlnum, tempr, ulp : RT = ---;

    ns.* = 0;
    nd.* = 0;
    info.* = 0;
    _laqz2_set_work1_complex(work, ncast(RT, max(1, nw)), RT, MT);
    if n <= 0 {
        return;
    }

    jw = min(nw, ihi - ilo + 1);
    kwtop = ihi - jw + 1;
    if kwtop == ilo {
        s = zero_t;
    } else {
        s = fget(a.*, kwtop, kwtop - 1);
    }

    ifst = 1;
    ilst = jw;
    a_win_query := array2d_view(_laqz2_matrix_storage(a.*, lda, kwtop, kwtop, jw, jw), lda, row_major=false);
    b_win_query := array2d_view(_laqz2_matrix_storage(b.*, ldb, kwtop, kwtop, jw, jw), ldb, row_major=false);
    #if MT == C32 {
        claqz0(#char "S", #char "V", #char "V", jw, 1, jw, *a_win_query, lda, *b_win_query, ldb, alpha, beta, qc, ldqc, zc, ldzc, work, -1, rwork, rec + 1, *qz_small_info);
    } else {
        zlaqz0(#char "S", #char "V", #char "V", jw, 1, jw, *a_win_query, lda, *b_win_query, ldb, alpha, beta, qc, ldqc, zc, ldzc, work, -1, rwork, rec + 1, *qz_small_info);
    }
    lworkreq = ncast(int, real(fget(work.*, 1))) + 2*jw*jw;
    lworkreq = max(lworkreq, n*nw, 2*nw*nw + n);
    if lwork == -1 {
        _laqz2_set_work1_complex(work, ncast(RT, lworkreq), RT, MT);
        return;
    } else if lwork < lworkreq {
        info.* = -26;
    }

    if info.* != 0 {
        xerbla(ifx MT == C32 then "CLAQZ2" else "ZLAQZ2", -info.*);
        return;
    }

    safmin = _laqz2_lamch(RT);
    safmax = one(RT) / safmin;
    ulp = ifx RT == float32 then slamch("Precision".data) else dlamch("Precision".data);
    smlnum = safmin*(ncast(RT, n) / ulp);
    _ = safmax;

    if ihi == kwtop {
        fset(alpha, kwtop, fget(a.*, kwtop, kwtop));
        fset(beta, kwtop, fget(b.*, kwtop, kwtop));
        ns.* = 1;
        nd.* = 0;
        if abs(s) <= max(smlnum, ulp*abs(fget(a.*, kwtop, kwtop))) {
            ns.* = 0;
            nd.* = 1;
            if kwtop > ilo {
                fset(a, kwtop, kwtop - 1, zero_t);
            }
        }
    }

    win_a_copy := array2d_view(_laqz2_work_storage(work, 1, jw, jw, jw), jw, row_major=false);
    _laqz2_lacpy_complex(#char "A", jw, jw, *a_win_query, lda, *win_a_copy, jw, MT);
    win_b_copy := array2d_view(_laqz2_work_storage(work, jw*jw + 1, jw, jw, jw), jw, row_major=false);
    _laqz2_lacpy_complex(#char "A", jw, jw, *b_win_query, ldb, *win_b_copy, jw, MT);

    _laqz2_laset_complex(#char "F", jw, jw, zero_t, one_t, qc, ldqc, MT);
    _laqz2_laset_complex(#char "F", jw, jw, zero_t, one_t, zc, ldzc, MT);
    work_qz := _laqz2_raw_slice(work.data + (2*jw*jw), lwork - 2*jw*jw);
    #if MT == C32 {
        claqz0(#char "S", #char "V", #char "V", jw, 1, jw, *a_win_query, lda, *b_win_query, ldb, alpha, beta, qc, ldqc, zc, ldzc, *work_qz, lwork - 2*jw*jw, rwork, rec + 1, *qz_small_info);
    } else {
        zlaqz0(#char "S", #char "V", #char "V", jw, 1, jw, *a_win_query, lda, *b_win_query, ldb, alpha, beta, qc, ldqc, zc, ldzc, *work_qz, lwork - 2*jw*jw, rwork, rec + 1, *qz_small_info);
    }

    if qz_small_info != 0 {
        nd.* = 0;
        ns.* = jw - qz_small_info;
        _laqz2_lacpy_complex(#char "A", jw, jw, *win_a_copy, jw, *a_win_query, lda, MT);
        _laqz2_lacpy_complex(#char "A", jw, jw, *win_b_copy, jw, *b_win_query, ldb, MT);
        return;
    }

    if kwtop == ilo || s == zero_t {
        kwbot = kwtop - 1;
    } else {
        kwbot = ihi;
        k = 1;
        k2 = 1;
        while k <= jw {
            tempr = abs(fget(a.*, kwbot, kwbot));
            if tempr == zero_r {
                tempr = abs(s);
            }
            if abs(s*fget(qc.*, 1, kwbot - kwtop + 1)) <= max(ulp*tempr, smlnum) {
                kwbot -= 1;
            } else {
                ifst = kwbot - kwtop + 1;
                ilst = k2;
                _laqz2_tgexc_complex(true, true, jw, *a_win_query, lda, *b_win_query, ldb, qc, ldqc, zc, ldzc, *ifst, *ilst, *ctgexc_info, MT);
                k2 += 1;
            }
            k += 1;
        }
    }

    nd.* = ihi - kwbot;
    ns.* = jw - nd.*;
    for kk:kwtop..ihi {
        fset(alpha, kk, fget(a.*, kk, kk));
        fset(beta, kk, fget(b.*, kk, kk));
    }

    if kwtop != ilo && s != zero_t {
        for kk:kwtop..kwbot {
            fset(a, kk, kwtop - 1, fget(a.*, kwtop, kwtop - 1)*conjugate(fget(qc.*, 1, kk - kwtop + 1)));
        }

        for <kk:kwbot - 1..kwtop {
            _laqz2_lartg_complex(fget(a.*, kk, kwtop - 1), fget(a.*, kk + 1, kwtop - 1), c1, s1, temp, RT, MT);
            fset(a, kk, kwtop - 1, temp);
            fset(a, kk + 1, kwtop - 1, zero_t);
            k2 = max(kwtop, kk - 1);
            _laqz2_rot_complex(ihi - k2 + 1, _laqz2_storage(a.data + (kk - 1)*lda + (k2 - 1), ihi - k2 + 1, lda), lda,
                                               _laqz2_storage(a.data + kk*lda + (k2 - 1), ihi - k2 + 1, lda), lda, c1, s1, RT, MT);
            _laqz2_rot_complex(ihi - (kk - 1) + 1, _laqz2_storage(b.data + (kk - 1)*ldb + (kk - 2), ihi - kk + 2, ldb), ldb,
                                                   _laqz2_storage(b.data + kk*ldb + (kk - 2), ihi - kk + 2, ldb), ldb, c1, s1, RT, MT);
            _laqz2_rot_complex(jw, _laqz2_storage(qc.data + (kk - kwtop)*ldqc, jw, 1), 1,
                                   _laqz2_storage(qc.data + (kk + 1 - kwtop)*ldqc, jw, 1), 1, c1, conjugate(s1), RT, MT);
        }

        istartm = kwtop;
        istopm = ihi;
        k = kwbot - 1;
        while k >= kwtop {
            for kk:k..kwbot - 1 {
                #if MT == C32 {
                    claqz1(true, true, kk, kwtop, kwtop + jw - 1, kwbot, a, lda, b, ldb, jw, kwtop, qc, ldqc, jw, kwtop, zc, ldzc);
                } else {
                    zlaqz1(true, true, kk, kwtop, kwtop + jw - 1, kwbot, a, lda, b, ldb, jw, kwtop, qc, ldqc, jw, kwtop, zc, ldzc);
                }
            }
            k -= 1;
        }
    }

    if ilschur {
        istartm = 1;
        istopm = n;
    } else {
        istartm = ilo;
        istopm = ihi;
    }

    if istopm - ihi > 0 {
        a_right := array2d_view(_laqz2_matrix_storage(a.*, lda, kwtop, ihi + 1, jw, istopm - ihi), lda, row_major=false);
        work_right := array2d_view(_laqz2_work_storage(work, 1, jw, istopm - ihi, jw), jw, row_major=false);
        _laqz2_gemm_complex(#char "C", #char "N", jw, istopm - ihi, jw, one_t, qc.*, ldqc, a_right.*, lda, zero_t, *work_right, jw, MT);
        _laqz2_lacpy_complex(#char "A", jw, istopm - ihi, work_right.*, jw, *a_right, lda, MT);
        b_right := array2d_view(_laqz2_matrix_storage(b.*, ldb, kwtop, ihi + 1, jw, istopm - ihi), ldb, row_major=false);
        _laqz2_gemm_complex(#char "C", #char "N", jw, istopm - ihi, jw, one_t, qc.*, ldqc, b_right.*, ldb, zero_t, *work_right, jw, MT);
        _laqz2_lacpy_complex(#char "A", jw, istopm - ihi, work_right.*, jw, *b_right, ldb, MT);
    }
    if ilq {
        q_cols := array2d_view(_laqz2_matrix_storage(q.*, ldq, 1, kwtop, n, jw), ldq, row_major=false);
        work_q := array2d_view(_laqz2_work_storage(work, 1, n, jw, n), n, row_major=false);
        _laqz2_gemm_complex(#char "N", #char "N", n, jw, jw, one_t, q_cols.*, ldq, qc.*, ldqc, zero_t, *work_q, n, MT);
        _laqz2_lacpy_complex(#char "A", n, jw, work_q.*, n, *q_cols, ldq, MT);
    }

    if kwtop - istartm > 0 {
        a_left := array2d_view(_laqz2_matrix_storage(a.*, lda, istartm, kwtop, kwtop - istartm, jw), lda, row_major=false);
        work_left := array2d_view(_laqz2_work_storage(work, 1, kwtop - istartm, jw, kwtop - istartm), kwtop - istartm, row_major=false);
        _laqz2_gemm_complex(#char "N", #char "N", kwtop - istartm, jw, jw, one_t, a_left.*, lda, zc.*, ldzc, zero_t, *work_left, kwtop - istartm, MT);
        _laqz2_lacpy_complex(#char "A", kwtop - istartm, jw, work_left.*, kwtop - istartm, *a_left, lda, MT);
        b_left := array2d_view(_laqz2_matrix_storage(b.*, ldb, istartm, kwtop, kwtop - istartm, jw), ldb, row_major=false);
        _laqz2_gemm_complex(#char "N", #char "N", kwtop - istartm, jw, jw, one_t, b_left.*, ldb, zc.*, ldzc, zero_t, *work_left, kwtop - istartm, MT);
        _laqz2_lacpy_complex(#char "A", kwtop - istartm, jw, work_left.*, kwtop - istartm, *b_left, ldb, MT);
    }
    if ilz {
        z_cols := array2d_view(_laqz2_matrix_storage(z.*, ldz, 1, kwtop, n, jw), ldz, row_major=false);
        work_z := array2d_view(_laqz2_work_storage(work, 1, n, jw, n), n, row_major=false);
        _laqz2_gemm_complex(#char "N", #char "N", n, jw, jw, one_t, z_cols.*, ldz, zc.*, ldzc, zero_t, *work_z, n, MT);
        _laqz2_lacpy_complex(#char "A", n, jw, work_z.*, n, *z_cols, ldz, MT);
    }

    _laqz2_set_work1_complex(work, ncast(RT, lworkreq), RT, MT);
}
