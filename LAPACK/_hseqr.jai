// ? LAPACK *hseqr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chseqr/zhseqr
// *HSEQR computes the eigenvalues of a Hessenberg matrix and,
// optionally, the Schur form and Schur vectors.

chseqr :: #bake_arguments _hseqr(RT=float32, MT=C32);
zhseqr :: #bake_arguments _hseqr(RT=float64, MT=C64);

#scope_file
_hseqr_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hseqr_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hseqr_raw_slice(x.data + offset - 1, count);
}

_hseqr_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hseqr_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hseqr_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hseqr_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 { clacpy(uplo, m, n, a, lda, b, ldb); } else { zlacpy(uplo, m, n, a, lda, b, ldb); }
}

_hseqr_laset :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: MT, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 { claset(uplo, m, n, alpha, beta, a, lda); } else { zlaset(uplo, m, n, alpha, beta, a, lda); }
}

_hseqr_lahqr :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int, h: *$H/ArrayType2D, $$ldh: int,
    w: *$W/ArrayType1D, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && H.data_type == MT && W.data_type == MT && Z.data_type == MT;
} {
    #if MT == C32 {
        clahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
    } else {
        zlahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
    }
}

_hseqr_laqr0 :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int, h: *$H/ArrayType2D, $$ldh: int,
    w: *$W/ArrayType1D, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           H.data_type == MT &&
           W.data_type == MT &&
           Z.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        claqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
    } else {
        zlaqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
    }
}

_hseqr :: (
    $$job: u8, $$compz: u8, $$n: int, $$ilo: int, $$ihi: int, h: *$H/ArrayType2D, $$ldh: int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D,
    $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT &&
           W.data_type == MT &&
           Z.data_type == MT &&
           WORK.data_type == MT;
} {
    NTINY :: 15;
    NL :: 49;
    zero_t :: #run zero(MT);
    one_t :: #run _hseqr_complex_from_real(one(RT), RT, MT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == C32 { "CHSEQR" } else { "ZHSEQR" };

    hl_storage: [NL*NL] MT;
    workl_storage: [NL] MT;

    i, kbot, nmin : int = ---;
    opt : RT = ---;
    wantt, initz, lquery, wantz : bool = ---;

    wantt = lsame(job, #char "S");
    initz = lsame(compz, #char "I");
    wantz = initz || lsame(compz, #char "V");
    fset(work, 1, _hseqr_complex_from_real(ncast(RT, max(1, n)), RT, MT));
    lquery = lwork == -1;

    info.* = 0;
    if !(lsame(job, #char "E") || wantt) {
        info.* = -1;
    } else if !(lsame(compz, #char "N") || wantz) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ilo < 1 || ilo > max(1, n) {
        info.* = -4;
    } else if ihi < min(ilo, n) || ihi > n {
        info.* = -5;
    } else if ldh < max(1, n) {
        info.* = -7;
    } else if ldz < 1 || (wantz && ldz < max(1, n)) {
        info.* = -10;
    } else if lwork < max(1, n) && !lquery {
        info.* = -12;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if n == 0 {
        return;
    } else if lquery {
        _hseqr_laqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info, MT);
        opt = max(real(fget(work.*, 1)), ncast(RT, max(1, n)));
        fset(work, 1, _hseqr_complex_from_real(opt, RT, MT));
        return;
    } else {
        for i:1..ilo - 1 {
            fset(w, i, fget(h.*, i, i));
        }
        for i:ihi + 1..n {
            fset(w, i, fget(h.*, i, i));
        }

        if initz {
            _hseqr_laset(#char "A", n, n, zero_t, one_t, z, ldz, MT);
        }

        if ilo == ihi {
            fset(w, ilo, fget(h.*, ilo, ilo));
            return;
        }

        nmin = ilaenv(12, routine_name, " "s, n, ilo, ihi, lwork);
        nmin = max(NTINY, nmin);

        if n > nmin {
            _hseqr_laqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info, MT);
        } else {
            _hseqr_lahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, info, MT);

            if info.* > 0 {
                kbot = info.*;

                if n >= NL {
                    _hseqr_laqr0(wantt, wantz, n, ilo, kbot, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info, MT);
                } else {
                    hl_view := array2d_view(hl_storage[0..], NL, row_major=false);
                    workl_view := workl_storage[0..];

                    _hseqr_lacpy(#char "A", n, n, h.*, ldh, *hl_view, NL, MT);
                    fset(hl_view, n + 1, n, zero_t);
                    if NL - n > 0 {
                        pad_raw := _hseqr_matrix_storage(hl_view, NL, 1, n + 1, NL, NL - n);
                        pad := array2d_view(pad_raw, NL, row_major=false);
                        _hseqr_laset(#char "A", NL, NL - n, zero_t, zero_t, *pad, NL, MT);
                    }
                    _hseqr_laqr0(wantt, wantz, NL, ilo, kbot, *hl_view, NL, w, ilo, ihi, z, ldz, *workl_view, NL, info, MT);
                    if wantt || info.* != 0 {
                        _hseqr_lacpy(#char "A", n, n, hl_view, NL, h, ldh, MT);
                    }
                }
            }
        }

        if (wantt || info.* != 0) && n > 2 {
            lower_raw := _hseqr_matrix_storage(h.*, ldh, 3, 1, n - 2, n - 2);
            lower := array2d_view(lower_raw, ldh, row_major=false);
            _hseqr_laset(#char "L", n - 2, n - 2, zero_t, zero_t, *lower, ldh, MT);
        }

        opt = max(ncast(RT, max(1, n)), real(fget(work.*, 1)));
        fset(work, 1, _hseqr_complex_from_real(opt, RT, MT));
    }
}
