// ? LAPACK *hbevd_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chbevd_2stage/zhbevd_2stage
// *HBEVD_2STAGE computes all eigenvalues and, optionally, eigenvectors of
// a complex Hermitian band matrix using the 2stage reduction.

chbevd_2stage :: #bake_arguments _hbevd_2stage(RT=float32, MT=C32);
zhbevd_2stage :: #bake_arguments _hbevd_2stage(RT=float64, MT=C64);

#scope_file
_hbevd_2stage_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hbevd_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hbevd_2stage_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach.data);
    } else {
        return dlamch(cmach.data);
    }
}

_hbevd_2stage_lanhb :: (
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

_hbevd_2stage_lascl :: (
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

_hbevd_2stage_scal :: (
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

_hbevd_2stage_sterf :: (
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

_hbevd_2stage_stedc :: (
    $$compz: u8,
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    $$lrwork: int,
    iwork: *$IWORK/ArrayType1D,
    $$liwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT &&
           E.data_type == RT &&
           Z.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    #if MT == C32 {
        cstedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
    } else {
        zstedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
    }
}

_hbevd_2stage_gemm :: (
    $$transa: u8,
    $$transb: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_hbevd_2stage_lacpy :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_hbevd_2stage_ilaenv2stage :: (
    $$ispec: int,
    $$name: string,
    $$opts: u8,
    $$n1: int,
    $$n2: int,
    $$n3: int,
    $$n4: int
) -> int {
    return ilaenv2stage(ispec, name.data, opts, n1, n2, n3, n4);
}

_hbevd_2stage_hetrd_hb2st :: (
    $$stage1: u8,
    $$vect: u8,
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    hous: *$HOUS/ArrayType1D,
    $$lhous: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           HOUS.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        chetrd_hb2st(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
    } else {
        zhetrd_hb2st(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
    }
}

_hbevd_2stage :: (
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
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    $$lrwork: int,
    iwork: *$IWORK/ArrayType1D,
    $$liwork: int,
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
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    czero_t :: #run _hbevd_2stage_complex_from_real(zero_r, RT, MT);
    cone_t  :: #run _hbevd_2stage_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHBEVD_2STAGE" } else { "ZHBEVD_2STAGE" };

    anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum : RT = ---;
    ib, iinfo, imax, inde, indhous, indwk, indwk2, indrwk, iscale, lhtrd, liwmin, llrwk, llwk2, llwork, lrwmin, lwmin, lwtrd : int = ---;
    lower, lquery, wantz : bool = ---;

    wantz = lsame(jobz, #char "V");
    lower = lsame(uplo, #char "L");
    lquery = lwork == -1 || liwork == -1 || lrwork == -1;

    info.* = 0;
    if n <= 1 {
        lwmin = 1;
        lrwmin = 1;
        liwmin = 1;
    } else {
        name := ifx MT == C32 then "CHETRD_HB2ST" else "ZHETRD_HB2ST";
        ib = _hbevd_2stage_ilaenv2stage(2, name, jobz, n, kd, -1, -1);
        lhtrd = _hbevd_2stage_ilaenv2stage(3, name, jobz, n, kd, ib, -1);
        lwtrd = _hbevd_2stage_ilaenv2stage(4, name, jobz, n, kd, ib, -1);
        if wantz {
            lwmin = 2*n*n;
            lrwmin = 1 + 5*n + 2*n*n;
            liwmin = 3 + 5*n;
        } else {
            lwmin = max(n, lhtrd + lwtrd);
            lrwmin = n;
            liwmin = 1;
        }
    }

    if !lsame(jobz, #char "N") {
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

    if info.* == 0 {
        fset(work, 1, _hbevd_2stage_complex_from_real(ncast(RT, lwmin), RT, MT));
        fset(rwork, 1, ncast(RT, lrwmin));
        fset(iwork, 1, liwmin);

        if lwork < lwmin && !lquery {
            info.* = -11;
        } else if lrwork < lrwmin && !lquery {
            info.* = -13;
        } else if liwork < liwmin && !lquery {
            info.* = -15;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        return;
    }

    if n == 1 {
        fset(w, 1, real(fget(ab.*, 1, 1)));
        if wantz {
            fset(z, 1, 1, cone_t);
        }
        return;
    }

    safmin = _hbevd_2stage_lamch("Safe minimum"s, RT);
    eps = _hbevd_2stage_lamch("Precision"s, RT);
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

    anrm = _hbevd_2stage_lanhb(#char "M", uplo, n, kd, ab.*, ldab, rwork, RT, MT);
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
            _hbevd_2stage_lascl(#char "B", kd, kd, one_r, sigma, n, n, ab, ldab, info, RT, MT);
        } else {
            _hbevd_2stage_lascl(#char "Q", kd, kd, one_r, sigma, n, n, ab, ldab, info, RT, MT);
        }
    }

    inde = 1;
    indrwk = inde + n;
    llrwk = lrwork - indrwk + 1;
    indhous = 1;
    indwk = indhous + lhtrd;
    llwork = lwork - indwk + 1;
    indwk2 = indwk + n*n;
    llwk2 = lwork - indwk2 + 1;

    hous_view := _hbevd_2stage_subvector(work.*, indhous, lhtrd);
    work1_view := _hbevd_2stage_subvector(work.*, indwk, llwork);
    _hbevd_2stage_hetrd_hb2st(#char "N", jobz, uplo, n, kd, ab, ldab, w, rwork, *hous_view, lhtrd, *work1_view, llwork, *iinfo, RT, MT);

    if !wantz {
        _hbevd_2stage_sterf(n, w, rwork, info, RT);
    } else {
        work_z := array2d_view(work.*, n, row_major=false);

        work2_view := _hbevd_2stage_subvector(work.*, indwk2, llwk2);
        rwork2_view := _hbevd_2stage_subvector(rwork.*, indrwk, llrwk);
        work2_mat_view := _hbevd_2stage_subvector(work.*, indwk2, n*n);
        work2_mat := array2d_view(work2_mat_view, n, row_major=false);

        _hbevd_2stage_stedc(#char "I", n, w, rwork, *work_z, n, *work2_view, llwk2, *rwork2_view, llrwk, iwork, liwork, info, RT, MT);
        _hbevd_2stage_gemm(#char "N", #char "N", n, n, n, cone_t, z.*, ldz, *work_z, n, czero_t, *work2_mat, n, MT);
        _hbevd_2stage_lacpy(#char "A", n, n, *work2_mat, n, z, ldz, MT);
    }

    if iscale == 1 {
        if info.* == 0 {
            imax = n;
        } else {
            imax = info.* - 1;
        }
        _hbevd_2stage_scal(imax, one_r / sigma, w, 1, RT);
    }

    fset(work, 1, _hbevd_2stage_complex_from_real(ncast(RT, lwmin), RT, MT));
    fset(rwork, 1, ncast(RT, lrwmin));
    fset(iwork, 1, liwmin);
}
