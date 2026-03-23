// ? LAPACK *hpev()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpev/zhpev
// *HPEV computes all eigenvalues and, optionally, eigenvectors of a Hermitian packed matrix.

chpev :: #bake_arguments _hpev(RT=float32, MT=C32);
zhpev :: #bake_arguments _hpev(RT=float64, MT=C64);

#scope_file
_hpev_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hpev_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hpev_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach.data);
    } else {
        return dlamch(cmach.data);
    }
}

_hpev_lanhp :: (
    $$norm: u8, $$uplo: u8, $$n: int, ap: $AP/ArrayType1D, work: *$RWORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 {
        return clanhp(norm, uplo, n, ap, work);
    } else {
        return zlanhp(norm, uplo, n, ap, work);
    }
}

_hpev_ptrd :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, d: *$D/ArrayType1D, e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && D.data_type == RT && E.data_type == RT && TAU.data_type == MT;
} {
    #if MT == C32 {
        chptrd(uplo, n, ap, d, e, tau, info);
    } else {
        zhptrd(uplo, n, ap, d, e, tau, info);
    }
}

_hpev_scal_packed :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hpev_upgtr :: (
    $$uplo: u8, $$n: int, ap: $AP/ArrayType1D, tau: $TAU/ArrayType1D, q: *$Q/ArrayType2D,
    $$ldq: int, work: *$WORK/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           AP.data_type == MT && TAU.data_type == MT && Q.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cupgtr(uplo, n, ap, tau, q, ldq, work, info);
    } else {
        zupgtr(uplo, n, ap, tau, q, ldq, work, info);
    }
}

_hpev_steqr :: (
    $$compz: u8, $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, z: *$Z/ArrayType2D,
    $$ldz: int, work: *$RWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && Z.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 {
        csteqr(compz, n, d, e, z, ldz, work, info);
    } else {
        zsteqr(compz, n, d, e, z, ldz, work, info);
    }
}

_hpev_sterf :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && D.data_type == RT && E.data_type == RT;
} {
    #if RT == float32 { ssterf(n, d, e, info); } else { dsterf(n, d, e, info); }
}

_hpev_scal_real :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 { sscal(n, alpha, x, incx); } else { dscal(n, alpha, x, incx); }
}

_hpev :: (
    $$jobz: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    one_t  :: #run _hpev_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHPEV " } else { "ZHPEV " };

    wantz : bool = ---;
    iinfo, imax, inde, inddrwk, indtau, indwrk, iscale : int = ---;
    anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum : RT = ---;

    wantz = lsame(jobz, #char "V");

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(lsame(uplo, #char "L") || lsame(uplo, #char "U")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ldz < 1 || (wantz && ldz < n) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }
    if n == 1 {
        fset(w, 1, real(fget(ap.*, 1)));
        fset(rwork, 1, one_r);
        if wantz { fset(z, 1, 1, one_t); }
        return;
    }

    safmin = _hpev_lamch("Safe minimum", RT);
    eps = _hpev_lamch("Precision", RT);
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

    anrm = _hpev_lanhp(#char "M", uplo, n, ap.*, rwork, RT, MT);
    iscale = 0;
    if anrm > zero_r && anrm < rmin {
        iscale = 1;
        sigma = rmin / anrm;
    } else if anrm > rmax {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if iscale == 1 {
        _hpev_scal_packed((n*(n + 1)) / 2, sigma, ap, 1, RT, MT);
    }

    inde = 1;
    indtau = 1;
    e_view := _hpev_subvector(rwork.*, inde, n - 1);
    tau_view := _hpev_subvector(work.*, indtau, n - 1);
    _hpev_ptrd(uplo, n, ap, w, *e_view, *tau_view, info, RT, MT);

    if !wantz {
        _hpev_sterf(n, w, *e_view, info, RT);
    } else {
        indwrk = indtau + n;
        work_gen := _hpev_subvector(work.*, indwrk, max(1, n - 1));
        _hpev_upgtr(uplo, n, ap.*, tau_view, z, ldz, *work_gen, *iinfo, MT);
        inddrwk = inde + n;
        rwork_qr := _hpev_subvector(rwork.*, inddrwk, max(1, 2*n - 1));
        _hpev_steqr(jobz, n, w, *e_view, z, ldz, *rwork_qr, info, RT, MT);
    }

    if iscale == 1 {
        if info.* == 0 {
            imax = n;
        } else {
            imax = info.* - 1;
        }
        _hpev_scal_real(imax, one_r / sigma, w, 1, RT);
    }
}
