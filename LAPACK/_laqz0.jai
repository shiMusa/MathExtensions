// ? LAPACK *laqz0()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqz0/dlaqz0/claqz0/zlaqz0
// *LAQZ0 is the multishift QZ driver.

slaqz0 :: #bake_arguments _laqz0_real(RT=float32);
dlaqz0 :: #bake_arguments _laqz0_real(RT=float64);
claqz0 :: #bake_arguments _laqz0_complex(RT=float32, MT=C32);
zlaqz0 :: #bake_arguments _laqz0_complex(RT=float64, MT=C64);

#scope_file
_laqz0_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_laqz0_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _laqz0_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqz0_work_storage :: inline (
    work: $W/ArrayType1D, $$offset: int, $$rows: int, $$cols: int, $$ldw: int
) -> [] W.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else ldw*(cols - 1) + rows;
    return _laqz0_raw_slice(work.data + (offset - 1), count);
}

_laqz0_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_laqz0_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_laqz0_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laqz0_set_work1_real :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && WORK.data_type == RT;
} {
    fset(work, 1, x);
}

_laqz0_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} {
    fset(work, 1, cpx(x, zero(RT)));
}

_laqz0_laset_real :: (
    $$uplo: u8, $$m: int, $$n: int, $$alpha: $RT, $$beta: RT, a: *$A/ArrayType2D, $$lda: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT;
} {
    #if RT == float32 {
        slaset(uplo, m, n, alpha, beta, a, lda);
    } else {
        dlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_laqz0_laset_complex :: (
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

_laqz0_lanhs_real :: (
    $$norm: u8, $$n: int, a: $A/ArrayType2D, $$lda: int, work: *$WORK/ArrayType1D, $RT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT && WORK.data_type == RT;
} {
    #if RT == float32 {
        return slanhs(norm, n, a, lda, work);
    } else {
        return dlanhs(norm, n, a, lda, work);
    }
}

_laqz0_lanhs_complex :: (
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

_laqz0_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach.data);
    } else {
        return dlamch(cmach.data);
    }
}

_laqz0_lartg_real :: (f: $RT, g: RT, c: *RT, s: *RT, r: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        slartg(f, g, c, s, r);
    } else {
        dlartg(f, g, c, s, r);
    }
}

_laqz0_lartg_complex :: (f: $MT, g: MT, c: *$RT, s: *MT, r: *MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        clartg(f, g, c, s, r);
    } else {
        zlartg(f, g, c, s, r);
    }
}

_laqz0_rot_real :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $$c: $RT, $$s: RT, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else {
        drot(n, x, incx, y, incy, c, s);
    }
}

_laqz0_rot_complex :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $$c: $RT, $$s: $MT, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        crot(n, x, incx, y, incy, c, s);
    } else {
        zrot(n, x, incx, y, incy, c, s);
    }
}

_laqz0_ilaenv_job :: inline ($$wants: u8, $$wantq: u8, $$wantz: u8) -> string {
    if wants == #char "E" {
        if wantq == #char "N" {
            if wantz == #char "N" return "ENN";
            if wantz == #char "V" return "ENV";
            return "ENI";
        }
        if wantq == #char "V" {
            if wantz == #char "N" return "EVN";
            if wantz == #char "V" return "EVV";
            return "EVI";
        }
        if wantz == #char "N" return "EIN";
        if wantz == #char "V" return "EIV";
        return "EII";
    }

    if wantq == #char "N" {
        if wantz == #char "N" return "SNN";
        if wantz == #char "V" return "SNV";
        return "SNI";
    }
    if wantq == #char "V" {
        if wantz == #char "N" return "SVN";
        if wantz == #char "V" return "SVV";
        return "SVI";
    }
    if wantz == #char "N" return "SIN";
    if wantz == #char "V" return "SIV";
    return "SII";
}

_laqz0_real :: (
    $$wants: u8, $$wantq: u8, $$wantz: u8, $$n: int, $$ilo: int, $$ihi: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, $$rec: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT &&
           ALPHAR.data_type == RT && ALPHAI.data_type == RT && BETA.data_type == RT &&
           Q.data_type == RT && Z.data_type == RT && WORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    routine_name :: #run if RT == float32 { "SLAQZ0" } else { "DLAQZ0" };

    ilschur, ilq, ilz : bool = ---;
    iwants, iwantq, iwantz, istart, istop, iiter, maxit, istart2, k, ld, nshifts,
    nblock, nw, nmin, nibble, n_undeflated, n_deflated, ns, sweep_info, shiftpos, lworkreq,
    k2, istartm, istopm, norm_info, aed_info, nwr, nbr, nsr, itemp1, itemp2, rcost, i : int = ---;
    safmin, safmax, c1, s1, temp, swap, bnorm, btol, eshift, smlnum, ulp : RT = ---;
    jbcmpz := "";

    if lsame(wants, #char "E") {
        ilschur = false;
        iwants = 1;
    } else if lsame(wants, #char "S") {
        ilschur = true;
        iwants = 2;
    } else {
        iwants = 0;
    }

    if lsame(wantq, #char "N") {
        ilq = false;
        iwantq = 1;
    } else if lsame(wantq, #char "V") {
        ilq = true;
        iwantq = 2;
    } else if lsame(wantq, #char "I") {
        ilq = true;
        iwantq = 3;
    } else {
        iwantq = 0;
    }

    if lsame(wantz, #char "N") {
        ilz = false;
        iwantz = 1;
    } else if lsame(wantz, #char "V") {
        ilz = true;
        iwantz = 2;
    } else if lsame(wantz, #char "I") {
        ilz = true;
        iwantz = 3;
    } else {
        iwantz = 0;
    }

    info.* = 0;
    if iwants == 0 {
        info.* = -1;
    } else if iwantq == 0 {
        info.* = -2;
    } else if iwantz == 0 {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if ilo < 1 {
        info.* = -5;
    } else if ihi > n || ihi < ilo - 1 {
        info.* = -6;
    } else if lda < max(1, n) {
        info.* = -8;
    } else if ldb < max(1, n) {
        info.* = -10;
    } else if ldq < 1 || (ilq && ldq < n) {
        info.* = -15;
    } else if ldz < 1 || (ilz && ldz < n) {
        info.* = -17;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n <= 0 {
        _laqz0_set_work1_real(work, one_r, RT);
        return;
    }

    jbcmpz = _laqz0_ilaenv_job(wants, wantq, wantz);
    nmin = ilaenv(12, routine_name, jbcmpz, n, ilo, ihi, lwork);
    nwr = ilaenv(13, routine_name, jbcmpz, n, ilo, ihi, lwork);
    nwr = max(2, nwr);
    nwr = min(ihi - ilo + 1, (n - 1)/3, nwr);
    nibble = ilaenv(14, routine_name, jbcmpz, n, ilo, ihi, lwork);
    nsr = ilaenv(15, routine_name, jbcmpz, n, ilo, ihi, lwork);
    nsr = min(nsr, (n + 6)/9, ihi - ilo);
    nsr = max(2, nsr - mod(nsr, 2));
    rcost = ilaenv(17, routine_name, jbcmpz, n, ilo, ihi, lwork);
    itemp1 = cast(int) (ncast(RT, nsr) / sqrt(one_r + ncast(RT, 2*nsr) /
             ((ncast(RT, rcost)/ncast(RT, 100))*ncast(RT, n))));
    itemp1 = ((itemp1 - 1)/4)*4 + 4;
    nbr = nsr + itemp1;

    if n < nmin || rec >= 2 {
        #if RT == float32 {
            shgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, info);
        } else {
            dhgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, info);
        }
        return;
    }

    nw = max(nwr, nmin);
    qc_query_raw := _laqz0_work_storage(work.*, 1, nw, nw, nw);
    qc_query := array2d_view(qc_query_raw, nw, row_major=false);
    zc_query_raw := _laqz0_work_storage(work.*, 1, nw, nw, nw);
    zc_query := array2d_view(zc_query_raw, nw, row_major=false);
    #if RT == float32 {
        slaqz3(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, *n_undeflated, *n_deflated, alphar, alphai, beta, *qc_query, nw, *zc_query, nw, work, -1, rec, *aed_info);
    } else {
        dlaqz3(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, *n_undeflated, *n_deflated, alphar, alphai, beta, *qc_query, nw, *zc_query, nw, work, -1, rec, *aed_info);
    }
    itemp1 = _laqz0_query_to_int(fget(work.*, 1), RT, RT);

    qc_sweep_query_raw := _laqz0_work_storage(work.*, 1, nbr, nbr, nbr);
    qc_sweep_query := array2d_view(qc_sweep_query_raw, nbr, row_major=false);
    zc_sweep_query_raw := _laqz0_work_storage(work.*, 1 + nbr*nbr, nbr, nbr, nbr);
    zc_sweep_query := array2d_view(zc_sweep_query_raw, nbr, row_major=false);
    #if RT == float32 {
        slaqz4(ilschur, ilq, ilz, n, ilo, ihi, nsr, nbr, alphar, alphai, beta, a, lda, b, ldb, q, ldq, z, ldz, *qc_sweep_query, nbr, *zc_sweep_query, nbr, work, -1, *sweep_info);
    } else {
        dlaqz4(ilschur, ilq, ilz, n, ilo, ihi, nsr, nbr, alphar, alphai, beta, a, lda, b, ldb, q, ldq, z, ldz, *qc_sweep_query, nbr, *zc_sweep_query, nbr, work, -1, *sweep_info);
    }
    itemp2 = _laqz0_query_to_int(fget(work.*, 1), RT, RT);

    lworkreq = max(itemp1 + 2*nw*nw, itemp2 + 2*nbr*nbr);
    if lwork == -1 {
        _laqz0_set_work1_real(work, _laqz0_roundup_lwork(lworkreq, RT), RT);
        return;
    } else if lwork < lworkreq {
        info.* = -19;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if iwantq == 3 {
        q_full_raw := _laqz0_matrix_storage(q.*, ldq, 1, 1, n, n);
        q_full := array2d_view(q_full_raw, ldq, row_major=false);
        _laqz0_laset_real(#char "F", n, n, zero_r, one_r, *q_full, ldq, RT);
    }
    if iwantz == 3 {
        z_full_raw := _laqz0_matrix_storage(z.*, ldz, 1, 1, n, n);
        z_full := array2d_view(z_full_raw, ldz, row_major=false);
        _laqz0_laset_real(#char "F", n, n, zero_r, one_r, *z_full, ldz, RT);
    }

    safmin = _laqz0_lamch("Safe minimum", RT);
    safmax = one_r / safmin;
    ulp = _laqz0_lamch("Precision", RT);
    smlnum = safmin*(ncast(RT, n) / ulp);
    b_sub_raw := _laqz0_matrix_storage(b.*, ldb, ilo, ilo, ihi - ilo + 1, ihi - ilo + 1);
    b_sub := array2d_view(b_sub_raw, ldb, row_major=false);
    bnorm = _laqz0_lanhs_real(#char "F", ihi - ilo + 1, *b_sub, ldb, work, RT);
    btol = max(safmin, ulp*bnorm);
    _ = safmax;

    istart = ilo;
    istop = ihi;
    maxit = 3*(ihi - ilo + 1);
    ld = 0;
    eshift = zero_r;

    for iiter:1..maxit {
        if iiter >= maxit {
            info.* = istop + 1;
            break;
        }
        if istart + 1 >= istop {
            istop = istart;
            break;
        }

        if abs(fget(a.*, istop - 1, istop - 2)) <= max(smlnum, ulp*(abs(fget(a.*, istop - 1, istop - 1)) + abs(fget(a.*, istop - 2, istop - 2)))) {
            fset(a, istop - 1, istop - 2, zero_r);
            istop -= 2;
            ld = 0;
            eshift = zero_r;
        } else if abs(fget(a.*, istop, istop - 1)) <= max(smlnum, ulp*(abs(fget(a.*, istop, istop)) + abs(fget(a.*, istop - 1, istop - 1)))) {
            fset(a, istop, istop - 1, zero_r);
            istop -= 1;
            ld = 0;
            eshift = zero_r;
        }

        if abs(fget(a.*, istart + 2, istart + 1)) <= max(smlnum, ulp*(abs(fget(a.*, istart + 1, istart + 1)) + abs(fget(a.*, istart + 2, istart + 2)))) {
            fset(a, istart + 2, istart + 1, zero_r);
            istart += 2;
            ld = 0;
            eshift = zero_r;
        } else if abs(fget(a.*, istart + 1, istart)) <= max(smlnum, ulp*(abs(fget(a.*, istart, istart)) + abs(fget(a.*, istart + 1, istart + 1)))) {
            fset(a, istart + 1, istart, zero_r);
            istart += 1;
            ld = 0;
            eshift = zero_r;
        }

        if istart + 1 >= istop {
            break;
        }

        istart2 = istart;
        for <k:istop..istart + 1 {
            if abs(fget(a.*, k, k - 1)) <= max(smlnum, ulp*(abs(fget(a.*, k, k)) + abs(fget(a.*, k - 1, k - 1)))) {
                fset(a, k, k - 1, zero_r);
                istart2 = k;
                break;
            }
        }

        if ilschur {
            istartm = 1;
            istopm = n;
        } else {
            istartm = istart2;
            istopm = istop;
        }

        k = istop;
        while k >= istart2 {
            if abs(fget(b.*, k, k)) < btol {
                for <k2:k..istart2 + 1 {
                    _laqz0_lartg_real(fget(b.*, k2 - 1, k2), fget(b.*, k2 - 1, k2 - 1), *c1, *s1, *temp, RT);
                    fset(b, k2 - 1, k2, temp);
                    fset(b, k2 - 1, k2 - 1, zero_r);

                    b_col_k := _laqz0_raw_slice(b.data + (k2 - 1)*ldb + (istartm - 1), k2 - 1 - istartm + 1);
                    b_col_km1 := _laqz0_raw_slice(b.data + (k2 - 2)*ldb + (istartm - 1), k2 - 1 - istartm + 1);
                    _laqz0_rot_real(k2 - istartm, *b_col_k, 1, *b_col_km1, 1, c1, s1, RT);

                    a_col_k := _laqz0_raw_slice(a.data + (k2 - 1)*lda + (istartm - 1), min(k2 + 1, istop) - istartm + 1);
                    a_col_km1 := _laqz0_raw_slice(a.data + (k2 - 2)*lda + (istartm - 1), min(k2 + 1, istop) - istartm + 1);
                    _laqz0_rot_real(min(k2 + 1, istop) - istartm + 1, *a_col_k, 1, *a_col_km1, 1, c1, s1, RT);

                    if ilz {
                        z_col_k := _laqz0_raw_slice(z.data + (k2 - 1)*ldz, n);
                        z_col_km1 := _laqz0_raw_slice(z.data + (k2 - 2)*ldz, n);
                        _laqz0_rot_real(n, *z_col_k, 1, *z_col_km1, 1, c1, s1, RT);
                    }

                    if k2 < istop {
                        _laqz0_lartg_real(fget(a.*, k2, k2 - 1), fget(a.*, k2 + 1, k2 - 1), *c1, *s1, *temp, RT);
                        fset(a, k2, k2 - 1, temp);
                        fset(a, k2 + 1, k2 - 1, zero_r);

                        a_row_k2 := _laqz0_raw_slice(a.data + (k2 - 1)*lda + (k2 - 1), lda*(istopm - k2));
                        a_row_k2p1 := _laqz0_raw_slice(a.data + (k2 - 1)*lda + k2, lda*(istopm - k2));
                        _laqz0_rot_real(istopm - k2 + 1, *a_row_k2, lda, *a_row_k2p1, lda, c1, s1, RT);

                        b_row_k2 := _laqz0_raw_slice(b.data + (k2 - 1)*ldb + (k2 - 1), ldb*(istopm - k2));
                        b_row_k2p1 := _laqz0_raw_slice(b.data + (k2 - 1)*ldb + k2, ldb*(istopm - k2));
                        _laqz0_rot_real(istopm - k2 + 1, *b_row_k2, ldb, *b_row_k2p1, ldb, c1, s1, RT);

                        if ilq {
                            q_col_k2 := _laqz0_raw_slice(q.data + (k2 - 1)*ldq, n);
                            q_col_k2p1 := _laqz0_raw_slice(q.data + k2*ldq, n);
                            _laqz0_rot_real(n, *q_col_k2, 1, *q_col_k2p1, 1, c1, s1, RT);
                        }
                    }
                }

                if istart2 < istop {
                    _laqz0_lartg_real(fget(a.*, istart2, istart2), fget(a.*, istart2 + 1, istart2), *c1, *s1, *temp, RT);
                    fset(a, istart2, istart2, temp);
                    fset(a, istart2 + 1, istart2, zero_r);

                    a_row_top := _laqz0_raw_slice(a.data + istart2*lda + (istart2 - 1), lda*(istopm - (istart2 + 1)));
                    a_row_bot := _laqz0_raw_slice(a.data + istart2*lda + istart2, lda*(istopm - (istart2 + 1)));
                    _laqz0_rot_real(istopm - istart2, *a_row_top, lda, *a_row_bot, lda, c1, s1, RT);

                    b_row_top := _laqz0_raw_slice(b.data + istart2*ldb + (istart2 - 1), ldb*(istopm - (istart2 + 1)));
                    b_row_bot := _laqz0_raw_slice(b.data + istart2*ldb + istart2, ldb*(istopm - (istart2 + 1)));
                    _laqz0_rot_real(istopm - istart2, *b_row_top, ldb, *b_row_bot, ldb, c1, s1, RT);

                    if ilq {
                        q_col_top := _laqz0_raw_slice(q.data + (istart2 - 1)*ldq, n);
                        q_col_bot := _laqz0_raw_slice(q.data + istart2*ldq, n);
                        _laqz0_rot_real(n, *q_col_top, 1, *q_col_bot, 1, c1, s1, RT);
                    }
                }

                istart2 += 1;
            }
            k -= 1;
        }

        if istart2 >= istop {
            istop = istart2 - 1;
            ld = 0;
            eshift = zero_r;
            continue;
        }

        nw = nwr;
        nshifts = nsr;
        nblock = nbr;
        if istop - istart2 + 1 < nmin {
            if istop - istart + 1 < nmin {
                nw = istop - istart + 1;
                istart2 = istart;
            } else {
                nw = istop - istart2 + 1;
            }
        }

        qc_raw := _laqz0_work_storage(work.*, 1, nw, nw, nw);
        qc := array2d_view(qc_raw, nw, row_major=false);
        zc_raw := _laqz0_work_storage(work.*, nw*nw + 1, nw, nw, nw);
        zc := array2d_view(zc_raw, nw, row_major=false);
        work_tail := array1d_view(_laqz0_raw_slice(work.data + 2*nw*nw, lwork - 2*nw*nw));

        #if RT == float32 {
            slaqz3(ilschur, ilq, ilz, n, istart2, istop, nw, a, lda, b, ldb, q, ldq, z, ldz, *n_undeflated, *n_deflated, alphar, alphai, beta, *qc, nw, *zc, nw, *work_tail, lwork - 2*nw*nw, rec, *aed_info);
        } else {
            dlaqz3(ilschur, ilq, ilz, n, istart2, istop, nw, a, lda, b, ldb, q, ldq, z, ldz, *n_undeflated, *n_deflated, alphar, alphai, beta, *qc, nw, *zc, nw, *work_tail, lwork - 2*nw*nw, rec, *aed_info);
        }

        if n_deflated > 0 {
            istop -= n_deflated;
            ld = 0;
            eshift = zero_r;
        }

        if 100*n_deflated > nibble*(n_deflated + n_undeflated) || istop - istart2 + 1 < nmin {
            continue;
        }

        ld += 1;
        ns = min(nshifts, istop - istart2);
        ns = min(ns, n_undeflated);
        shiftpos = istop - n_undeflated + 1;

        i = shiftpos;
        while i <= shiftpos + n_undeflated - 1 {
            if fget(alphai.*, i) != -fget(alphai.*, i + 1) {
                swap = fget(alphar.*, i);
                fset(alphar, i, fget(alphar.*, i + 1));
                fset(alphar, i + 1, fget(alphar.*, i + 2));
                fset(alphar, i + 2, swap);

                swap = fget(alphai.*, i);
                fset(alphai, i, fget(alphai.*, i + 1));
                fset(alphai, i + 1, fget(alphai.*, i + 2));
                fset(alphai, i + 2, swap);

                swap = fget(beta.*, i);
                fset(beta, i, fget(beta.*, i + 1));
                fset(beta, i + 1, fget(beta.*, i + 2));
                fset(beta, i + 2, swap);
            }
            i += 2;
        }

        if mod(ld, 6) == 0 {
            if (ncast(RT, maxit)*safmin)*abs(fget(a.*, istop, istop - 1)) < abs(fget(a.*, istop - 1, istop - 1)) {
                eshift = fget(a.*, istop, istop - 1) / fget(b.*, istop - 1, istop - 1);
            } else {
                eshift += one_r / (safmin*ncast(RT, maxit));
            }
            fset(alphar, shiftpos, one_r);
            fset(alphar, shiftpos + 1, zero_r);
            fset(alphai, shiftpos, zero_r);
            fset(alphai, shiftpos + 1, zero_r);
            fset(beta, shiftpos, eshift);
            fset(beta, shiftpos + 1, eshift);
            ns = 2;
        }

        qc_sweep_raw := _laqz0_work_storage(work.*, 1, nblock, nblock, nblock);
        qc_sweep := array2d_view(qc_sweep_raw, nblock, row_major=false);
        zc_sweep_raw := _laqz0_work_storage(work.*, nblock*nblock + 1, nblock, nblock, nblock);
        zc_sweep := array2d_view(zc_sweep_raw, nblock, row_major=false);
        sweep_work := array1d_view(_laqz0_raw_slice(work.data + 2*nblock*nblock, lwork - 2*nblock*nblock));
        shift_alphar := array1d_view(_laqz0_raw_slice(alphar.data + (shiftpos - 1), n - shiftpos + 1));
        shift_alphai := array1d_view(_laqz0_raw_slice(alphai.data + (shiftpos - 1), n - shiftpos + 1));
        shift_beta := array1d_view(_laqz0_raw_slice(beta.data + (shiftpos - 1), n - shiftpos + 1));
        #if RT == float32 {
            slaqz4(ilschur, ilq, ilz, n, istart2, istop, ns, nblock, *shift_alphar, *shift_alphai, *shift_beta, a, lda, b, ldb, q, ldq, z, ldz, *qc_sweep, nblock, *zc_sweep, nblock, *sweep_work, lwork - 2*nblock*nblock, *sweep_info);
        } else {
            dlaqz4(ilschur, ilq, ilz, n, istart2, istop, ns, nblock, *shift_alphar, *shift_alphai, *shift_beta, a, lda, b, ldb, q, ldq, z, ldz, *qc_sweep, nblock, *zc_sweep, nblock, *sweep_work, lwork - 2*nblock*nblock, *sweep_info);
        }
    }

    #if RT == float32 {
        shgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, *norm_info);
    } else {
        dhgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, *norm_info);
    }
    info.* = norm_info;
}

_laqz0_complex :: (
    $$wants: u8, $$wantq: u8, $$wantz: u8, $$n: int, $$ilo: int, $$ihi: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    q: *$Q/ArrayType2D, $$ldq: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$rec: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && B.data_type == MT && ALPHA.data_type == MT && BETA.data_type == MT &&
           Q.data_type == MT && Z.data_type == MT && WORK.data_type == MT && RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    czero  :: #run zero(MT);
    cone   :: #run _laqz0_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CLAQZ0" } else { "ZLAQZ0" };

    ilschur, ilq, ilz : bool = ---;
    iwants, iwantq, iwantz, istart, istop, iiter, maxit, istart2, k, ld, nshifts,
    nblock, nw, nmin, nibble, n_undeflated, n_deflated, ns, sweep_info, shiftpos, lworkreq,
    k2, istartm, istopm, norm_info, aed_info, nwr, nbr, nsr, itemp1, itemp2, rcost : int = ---;
    safmin, safmax, c1, bnorm, btol, smlnum, ulp : RT = ---;
    eshift, s1, temp : MT = ---;
    jbcmpz := "";

    if lsame(wants, #char "E") {
        ilschur = false;
        iwants = 1;
    } else if lsame(wants, #char "S") {
        ilschur = true;
        iwants = 2;
    } else {
        iwants = 0;
    }

    if lsame(wantq, #char "N") {
        ilq = false;
        iwantq = 1;
    } else if lsame(wantq, #char "V") {
        ilq = true;
        iwantq = 2;
    } else if lsame(wantq, #char "I") {
        ilq = true;
        iwantq = 3;
    } else {
        iwantq = 0;
    }

    if lsame(wantz, #char "N") {
        ilz = false;
        iwantz = 1;
    } else if lsame(wantz, #char "V") {
        ilz = true;
        iwantz = 2;
    } else if lsame(wantz, #char "I") {
        ilz = true;
        iwantz = 3;
    } else {
        iwantz = 0;
    }

    info.* = 0;
    if iwants == 0 {
        info.* = -1;
    } else if iwantq == 0 {
        info.* = -2;
    } else if iwantz == 0 {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if ilo < 1 {
        info.* = -5;
    } else if ihi > n || ihi < ilo - 1 {
        info.* = -6;
    } else if lda < max(1, n) {
        info.* = -8;
    } else if ldb < max(1, n) {
        info.* = -10;
    } else if ldq < 1 || (ilq && ldq < n) {
        info.* = -15;
    } else if ldz < 1 || (ilz && ldz < n) {
        info.* = -17;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n <= 0 {
        _laqz0_set_work1_complex(work, one_r, RT, MT);
        return;
    }

    jbcmpz = _laqz0_ilaenv_job(wants, wantq, wantz);
    nmin = ilaenv(12, routine_name, jbcmpz, n, ilo, ihi, lwork);
    nwr = ilaenv(13, routine_name, jbcmpz, n, ilo, ihi, lwork);
    nwr = max(2, nwr);
    nwr = min(ihi - ilo + 1, (n - 1)/3, nwr);
    nibble = ilaenv(14, routine_name, jbcmpz, n, ilo, ihi, lwork);
    nsr = ilaenv(15, routine_name, jbcmpz, n, ilo, ihi, lwork);
    nsr = min(nsr, (n + 6)/9, ihi - ilo);
    nsr = max(2, nsr - mod(nsr, 2));
    rcost = ilaenv(17, routine_name, jbcmpz, n, ilo, ihi, lwork);
    itemp1 = cast(int) (ncast(RT, nsr) / sqrt(one_r + ncast(RT, 2*nsr) /
             ((ncast(RT, rcost)/ncast(RT, 100))*ncast(RT, n))));
    itemp1 = ((itemp1 - 1)/4)*4 + 4;
    nbr = nsr + itemp1;

    if n < nmin || rec >= 2 {
        #if MT == C32 {
            chgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
        } else {
            zhgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
        }
        return;
    }

    nw = max(nwr, nmin);
    qc_query_raw := _laqz0_work_storage(work.*, 1, nw, nw, nw);
    qc_query := array2d_view(qc_query_raw, nw, row_major=false);
    zc_query_raw := _laqz0_work_storage(work.*, 1, nw, nw, nw);
    zc_query := array2d_view(zc_query_raw, nw, row_major=false);
    #if MT == C32 {
        claqz2(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, *n_undeflated, *n_deflated, alpha, beta, *qc_query, nw, *zc_query, nw, work, -1, rwork, rec, *aed_info);
    } else {
        zlaqz2(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, *n_undeflated, *n_deflated, alpha, beta, *qc_query, nw, *zc_query, nw, work, -1, rwork, rec, *aed_info);
    }
    itemp1 = _laqz0_query_to_int(fget(work.*, 1), RT, MT);

    qc_sweep_query_raw := _laqz0_work_storage(work.*, 1, nbr, nbr, nbr);
    qc_sweep_query := array2d_view(qc_sweep_query_raw, nbr, row_major=false);
    zc_sweep_query_raw := _laqz0_work_storage(work.*, 1 + nbr*nbr, nbr, nbr, nbr);
    zc_sweep_query := array2d_view(zc_sweep_query_raw, nbr, row_major=false);
    #if MT == C32 {
        claqz3(ilschur, ilq, ilz, n, ilo, ihi, nsr, nbr, alpha, beta, a, lda, b, ldb, q, ldq, z, ldz, *qc_sweep_query, nbr, *zc_sweep_query, nbr, work, -1, *sweep_info);
    } else {
        zlaqz3(ilschur, ilq, ilz, n, ilo, ihi, nsr, nbr, alpha, beta, a, lda, b, ldb, q, ldq, z, ldz, *qc_sweep_query, nbr, *zc_sweep_query, nbr, work, -1, *sweep_info);
    }
    itemp2 = _laqz0_query_to_int(fget(work.*, 1), RT, MT);

    lworkreq = max(itemp1 + 2*nw*nw, itemp2 + 2*nbr*nbr);
    if lwork == -1 {
        _laqz0_set_work1_complex(work, ncast(RT, lworkreq), RT, MT);
        return;
    } else if lwork < lworkreq {
        info.* = -19;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if iwantq == 3 {
        q_full_raw := _laqz0_matrix_storage(q.*, ldq, 1, 1, n, n);
        q_full := array2d_view(q_full_raw, ldq, row_major=false);
        _laqz0_laset_complex(#char "F", n, n, czero, cone, *q_full, ldq, MT);
    }
    if iwantz == 3 {
        z_full_raw := _laqz0_matrix_storage(z.*, ldz, 1, 1, n, n);
        z_full := array2d_view(z_full_raw, ldz, row_major=false);
        _laqz0_laset_complex(#char "F", n, n, czero, cone, *z_full, ldz, MT);
    }

    safmin = _laqz0_lamch("Safe minimum", RT);
    safmax = one_r / safmin;
    ulp = _laqz0_lamch("Precision", RT);
    smlnum = safmin*(ncast(RT, n) / ulp);
    b_sub_raw := _laqz0_matrix_storage(b.*, ldb, ilo, ilo, ihi - ilo + 1, ihi - ilo + 1);
    b_sub := array2d_view(b_sub_raw, ldb, row_major=false);
    bnorm = _laqz0_lanhs_complex(#char "F", ihi - ilo + 1, *b_sub, ldb, rwork, RT, MT);
    btol = max(safmin, ulp*bnorm);
    _ = safmax;

    istart = ilo;
    istop = ihi;
    maxit = 30*(ihi - ilo + 1);
    ld = 0;
    eshift = czero;

    for iiter:1..maxit {
        if iiter >= maxit {
            info.* = istop + 1;
            break;
        }
        if istart + 1 >= istop {
            istop = istart;
            break;
        }

        if abs(fget(a.*, istop, istop - 1)) <= max(smlnum, ulp*(abs(fget(a.*, istop, istop)) + abs(fget(a.*, istop - 1, istop - 1)))) {
            fset(a, istop, istop - 1, czero);
            istop -= 1;
            ld = 0;
            eshift = czero;
        }
        if abs(fget(a.*, istart + 1, istart)) <= max(smlnum, ulp*(abs(fget(a.*, istart, istart)) + abs(fget(a.*, istart + 1, istart + 1)))) {
            fset(a, istart + 1, istart, czero);
            istart += 1;
            ld = 0;
            eshift = czero;
        }

        if istart + 1 >= istop {
            break;
        }

        istart2 = istart;
        for <k:istop..istart + 1 {
            if abs(fget(a.*, k, k - 1)) <= max(smlnum, ulp*(abs(fget(a.*, k, k)) + abs(fget(a.*, k - 1, k - 1)))) {
                fset(a, k, k - 1, czero);
                istart2 = k;
                break;
            }
        }

        if ilschur {
            istartm = 1;
            istopm = n;
        } else {
            istartm = istart2;
            istopm = istop;
        }

        k = istop;
        while k >= istart2 {
            if abs(fget(b.*, k, k)) < btol {
                for <k2:k..istart2 + 1 {
                    _laqz0_lartg_complex(fget(b.*, k2 - 1, k2), fget(b.*, k2 - 1, k2 - 1), *c1, *s1, *temp, RT, MT);
                    fset(b, k2 - 1, k2, temp);
                    fset(b, k2 - 1, k2 - 1, czero);

                    b_col_k := _laqz0_raw_slice(b.data + (k2 - 1)*ldb + (istartm - 1), k2 - 1 - istartm + 1);
                    b_col_km1 := _laqz0_raw_slice(b.data + (k2 - 2)*ldb + (istartm - 1), k2 - 1 - istartm + 1);
                    _laqz0_rot_complex(k2 - istartm, *b_col_k, 1, *b_col_km1, 1, c1, s1, RT, MT);

                    a_col_k := _laqz0_raw_slice(a.data + (k2 - 1)*lda + (istartm - 1), min(k2 + 1, istop) - istartm + 1);
                    a_col_km1 := _laqz0_raw_slice(a.data + (k2 - 2)*lda + (istartm - 1), min(k2 + 1, istop) - istartm + 1);
                    _laqz0_rot_complex(min(k2 + 1, istop) - istartm + 1, *a_col_k, 1, *a_col_km1, 1, c1, s1, RT, MT);

                    if ilz {
                        z_col_k := _laqz0_raw_slice(z.data + (k2 - 1)*ldz, n);
                        z_col_km1 := _laqz0_raw_slice(z.data + (k2 - 2)*ldz, n);
                        _laqz0_rot_complex(n, *z_col_k, 1, *z_col_km1, 1, c1, s1, RT, MT);
                    }

                    if k2 < istop {
                        _laqz0_lartg_complex(fget(a.*, k2, k2 - 1), fget(a.*, k2 + 1, k2 - 1), *c1, *s1, *temp, RT, MT);
                        fset(a, k2, k2 - 1, temp);
                        fset(a, k2 + 1, k2 - 1, czero);

                        a_row_k2 := _laqz0_raw_slice(a.data + (k2 - 1)*lda + (k2 - 1), lda*(istopm - k2));
                        a_row_k2p1 := _laqz0_raw_slice(a.data + (k2 - 1)*lda + k2, lda*(istopm - k2));
                        _laqz0_rot_complex(istopm - k2 + 1, *a_row_k2, lda, *a_row_k2p1, lda, c1, s1, RT, MT);

                        b_row_k2 := _laqz0_raw_slice(b.data + (k2 - 1)*ldb + (k2 - 1), ldb*(istopm - k2));
                        b_row_k2p1 := _laqz0_raw_slice(b.data + (k2 - 1)*ldb + k2, ldb*(istopm - k2));
                        _laqz0_rot_complex(istopm - k2 + 1, *b_row_k2, ldb, *b_row_k2p1, ldb, c1, s1, RT, MT);

                        if ilq {
                            q_col_k2 := _laqz0_raw_slice(q.data + (k2 - 1)*ldq, n);
                            q_col_k2p1 := _laqz0_raw_slice(q.data + k2*ldq, n);
                            _laqz0_rot_complex(n, *q_col_k2, 1, *q_col_k2p1, 1, c1, conjugate(s1), RT, MT);
                        }
                    }
                }

                if istart2 < istop {
                    _laqz0_lartg_complex(fget(a.*, istart2, istart2), fget(a.*, istart2 + 1, istart2), *c1, *s1, *temp, RT, MT);
                    fset(a, istart2, istart2, temp);
                    fset(a, istart2 + 1, istart2, czero);

                    a_row_top := _laqz0_raw_slice(a.data + istart2*lda + (istart2 - 1), lda*(istopm - (istart2 + 1)));
                    a_row_bot := _laqz0_raw_slice(a.data + istart2*lda + istart2, lda*(istopm - (istart2 + 1)));
                    _laqz0_rot_complex(istopm - istart2, *a_row_top, lda, *a_row_bot, lda, c1, s1, RT, MT);

                    b_row_top := _laqz0_raw_slice(b.data + istart2*ldb + (istart2 - 1), ldb*(istopm - (istart2 + 1)));
                    b_row_bot := _laqz0_raw_slice(b.data + istart2*ldb + istart2, ldb*(istopm - (istart2 + 1)));
                    _laqz0_rot_complex(istopm - istart2, *b_row_top, ldb, *b_row_bot, ldb, c1, s1, RT, MT);

                    if ilq {
                        q_col_top := _laqz0_raw_slice(q.data + (istart2 - 1)*ldq, n);
                        q_col_bot := _laqz0_raw_slice(q.data + istart2*ldq, n);
                        _laqz0_rot_complex(n, *q_col_top, 1, *q_col_bot, 1, c1, conjugate(s1), RT, MT);
                    }
                }

                istart2 += 1;
            }
            k -= 1;
        }

        if istart2 >= istop {
            istop = istart2 - 1;
            ld = 0;
            eshift = czero;
            continue;
        }

        nw = nwr;
        nshifts = nsr;
        nblock = nbr;
        if istop - istart2 + 1 < nmin {
            if istop - istart + 1 < nmin {
                nw = istop - istart + 1;
                istart2 = istart;
            } else {
                nw = istop - istart2 + 1;
            }
        }

        qc_raw := _laqz0_work_storage(work.*, 1, nw, nw, nw);
        qc := array2d_view(qc_raw, nw, row_major=false);
        zc_raw := _laqz0_work_storage(work.*, nw*nw + 1, nw, nw, nw);
        zc := array2d_view(zc_raw, nw, row_major=false);
        work_tail := array1d_view(_laqz0_raw_slice(work.data + 2*nw*nw, lwork - 2*nw*nw));
        #if MT == C32 {
            claqz2(ilschur, ilq, ilz, n, istart2, istop, nw, a, lda, b, ldb, q, ldq, z, ldz, *n_undeflated, *n_deflated, alpha, beta, *qc, nw, *zc, nw, *work_tail, lwork - 2*nw*nw, rwork, rec, *aed_info);
        } else {
            zlaqz2(ilschur, ilq, ilz, n, istart2, istop, nw, a, lda, b, ldb, q, ldq, z, ldz, *n_undeflated, *n_deflated, alpha, beta, *qc, nw, *zc, nw, *work_tail, lwork - 2*nw*nw, rwork, rec, *aed_info);
        }

        if n_deflated > 0 {
            istop -= n_deflated;
            ld = 0;
            eshift = czero;
        }

        if 100*n_deflated > nibble*(n_deflated + n_undeflated) || istop - istart2 + 1 < nmin {
            continue;
        }

        ld += 1;
        ns = min(nshifts, istop - istart2);
        ns = min(ns, n_undeflated);
        shiftpos = istop - n_undeflated + 1;

        if mod(ld, 6) == 0 {
            if (ncast(RT, maxit)*safmin)*abs(fget(a.*, istop, istop - 1)) < abs(fget(a.*, istop - 1, istop - 1)) {
                eshift = fget(a.*, istop, istop - 1) / fget(b.*, istop - 1, istop - 1);
            } else {
                eshift += cone / _laqz0_complex_from_real(safmin*ncast(RT, maxit), RT, MT);
            }
            fset(alpha, shiftpos, cone);
            fset(beta, shiftpos, eshift);
            ns = 1;
        }

        qc_sweep_raw := _laqz0_work_storage(work.*, 1, nblock, nblock, nblock);
        qc_sweep := array2d_view(qc_sweep_raw, nblock, row_major=false);
        zc_sweep_raw := _laqz0_work_storage(work.*, nblock*nblock + 1, nblock, nblock, nblock);
        zc_sweep := array2d_view(zc_sweep_raw, nblock, row_major=false);
        sweep_work := array1d_view(_laqz0_raw_slice(work.data + 2*nblock*nblock, lwork - 2*nblock*nblock));
        shift_alpha := array1d_view(_laqz0_raw_slice(alpha.data + (shiftpos - 1), n - shiftpos + 1));
        shift_beta := array1d_view(_laqz0_raw_slice(beta.data + (shiftpos - 1), n - shiftpos + 1));
        #if MT == C32 {
            claqz3(ilschur, ilq, ilz, n, istart2, istop, ns, nblock, *shift_alpha, *shift_beta, a, lda, b, ldb, q, ldq, z, ldz, *qc_sweep, nblock, *zc_sweep, nblock, *sweep_work, lwork - 2*nblock*nblock, *sweep_info);
        } else {
            zlaqz3(ilschur, ilq, ilz, n, istart2, istop, ns, nblock, *shift_alpha, *shift_beta, a, lda, b, ldb, q, ldq, z, ldz, *qc_sweep, nblock, *zc_sweep, nblock, *sweep_work, lwork - 2*nblock*nblock, *sweep_info);
        }
    }

    #if MT == C32 {
        chgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, *norm_info);
    } else {
        zhgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, *norm_info);
    }
    info.* = norm_info;
}
