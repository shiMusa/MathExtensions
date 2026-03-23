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

_hgeqz_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
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

    ilazr2, ilazro, ilq, ilschr, ilz, lquery : bool = ---;
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

    if ihi < ilo {
        return;
    }

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

    // Staged port note:
    // The full single-shift QZ iteration from *hgeqz.f is still pending.
    // This file lands the public interface, argument validation,
    // workspace-query behavior, initialization, machine-constant setup,
    // and the trailing-eigenvalue normalization path.
    _ = ilazr2; _ = ilazro; _ = ifirst; _ = istart; _ = jc; _ = jch; _ = jr; _ = jiter;
    _ = cval; _ = temp; _ = temp2; _ = tempr;
    _ = abi12; _ = abi22; _ = ad11; _ = ad12; _ = ad21; _ = ad22;
    _ = ctemp; _ = ctemp2; _ = ctemp3; _ = s; _ = shift; _ = u12; _ = x; _ = y;
    _ = ilast; _ = ifrstm; _ = ilastm; _ = iiter; _ = eshift; _ = maxit;
    _ = atol; _ = btol; _ = ascale; _ = bscale;
}
