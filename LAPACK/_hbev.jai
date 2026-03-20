// ? LAPACK *hbev()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chbev/zhbev
// *HBEV computes all eigenvalues and, optionally, eigenvectors of a
// complex Hermitian band matrix.

chbev :: #bake_arguments _hbev(RT=float32, MT=C32);
zhbev :: #bake_arguments _hbev(RT=float64, MT=C64);

#scope_file
_hbev_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hbev_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach.data);
    } else {
        return dlamch(cmach.data);
    }
}

_hbev_lanhb :: (
    $$norm: u8,
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    work: *$RWORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           RWORK.data_type == RT;
} {
    #if MT == C32 {
        return clanhb(norm, uplo, n, kd, ab, ldab, work);
    } else {
        return zlanhb(norm, uplo, n, kd, ab, ldab, work);
    }
}

_hbev_lascl :: (
    $$type: u8,
    $$kl: int,
    $$ku: int,
    cfrom: $RT,
    cto: $RT,
    $$m: int,
    $$n: int,
    a: *$AB/ArrayType2D,
    $$lda: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT;
} {
    #if MT == C32 {
        clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_hbev_hbtrd :: (
    $$vect: u8,
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Z/ArrayType2D,
    $$ldq: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           Z.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        chbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
    } else {
        zhbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
    }
}

_hbev_steqr :: (
    $$compz: u8,
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT &&
           E.data_type == RT &&
           Z.data_type == MT &&
           RWORK.data_type == RT;
} {
    #if MT == C32 {
        csteqr(compz, n, d, e, z, ldz, work, info);
    } else {
        zsteqr(compz, n, d, e, z, ldz, work, info);
    }
}

_hbev_sterf :: (
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           E.data_type == RT;
} {
    #if RT == float32 {
        ssterf(n, d, e, info);
    } else {
        dsterf(n, d, e, info);
    }
}

_hbev_scal :: (
    $$n: int,
    alpha: $RT,
    x: *$W/ArrayType1D,
    $$incx: int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           W.data_type == RT;
} {
    #if RT == float32 {
        sscal(n, alpha, x, incx);
    } else {
        dscal(n, alpha, x, incx);
    }
}

_hbev :: (
    $$jobz: u8,
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           W.data_type == RT &&
           Z.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    one_t  :: #run _hbev_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHBEV " } else { "ZHBEV " };

    lower, wantz : bool = ---;
    iinfo, imax, inde, indrwk, iscale : int = ---;
    anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum : RT = ---;

    wantz = lsame(jobz, #char "V");
    lower = lsame(uplo, #char "L");

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(lower || lsame(uplo, #char "U")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if kd < 0 {
        info.* = -4;
    } else if ldab < kd + 1 {
        info.* = -6;
    } else if ldz < 1 || (wantz && ldz < n) {
        info.* = -9;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    if n == 1 {
        if lower {
            fset(w, 1, real(fget(ab.*, 1, 1)));
        } else {
            fset(w, 1, real(fget(ab.*, kd + 1, 1)));
        }
        if wantz {
            fset(z, 1, 1, one_t);
        }
        return;
    }

    safmin = _hbev_lamch("Safe minimum"s, RT);
    eps    = _hbev_lamch("Precision"s, RT);
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin   = sqrt(smlnum);
    rmax   = sqrt(bignum);

    anrm = _hbev_lanhb(#char "M", uplo, n, kd, ab.*, ldab, rwork, RT, MT);
    iscale = 0;
    if anrm > zero_r && anrm < rmin {
        iscale = 1;
        sigma = rmin / anrm;
    } else if anrm > rmax {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if iscale == 1 {
        if lower {
            _hbev_lascl(#char "B", kd, kd, one_r, sigma, n, n, ab, ldab, info, RT, MT);
        } else {
            _hbev_lascl(#char "Q", kd, kd, one_r, sigma, n, n, ab, ldab, info, RT, MT);
        }
    }

    inde = 1;
    _hbev_hbtrd(jobz, uplo, n, kd, ab, ldab, w, rwork, z, ldz, work, *iinfo, RT, MT);

    if !wantz {
        _hbev_sterf(n, w, rwork, info, RT);
    } else {
        indrwk = inde + n;
        rwork_tail: [] RWORK.data_type = rwork.*;
        rwork_tail.data += indrwk - 1;
        rwork_tail.count -= indrwk - 1;
        _hbev_steqr(jobz, n, w, rwork, z, ldz, *rwork_tail, info, RT, MT);
    }

    if iscale == 1 {
        if info.* == 0 {
            imax = n;
        } else {
            imax = info.* - 1;
        }
        _hbev_scal(imax, one_r / sigma, w, 1, RT);
    }
}
