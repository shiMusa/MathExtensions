// ? LAPACK *hbevx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chbevx/zhbevx
// *HBEVX computes selected eigenvalues and, optionally, eigenvectors of
// a complex Hermitian band matrix.

chbevx :: #bake_arguments _hbevx(RT=float32, MT=C32);
zhbevx :: #bake_arguments _hbevx(RT=float64, MT=C64);

#scope_file
_hbevx_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hbevx_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    #if CHECKS {
        assert(ld >= max(1, count));
        assert(col >= 1);
        assert(count >= 0);
    }

    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld;
    view.count = count;
    return view;
}

_hbevx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hbevx_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach.data);
    } else {
        return dlamch(cmach.data);
    }
}

_hbevx_lanhb :: (
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

_hbevx_lascl :: (
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

_hbevx_hbtrd :: (
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
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           Q.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        chbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
    } else {
        zhbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
    }
}

_hbevx_steqr :: (
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

_hbevx_sterf :: (
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

_hbevx_scal :: (
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

_hbevx_lacpy :: (
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

_hbevx_copy_real :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT;
} {
    #if RT == float32 {
        scopy(n, x, incx, y, incy);
    } else {
        dcopy(n, x, incx, y, incy);
    }
}

_hbevx_copy_complex :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_hbevx_gemv :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    beta: MT,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_hbevx_stebz :: (
    $$range: u8,
    $$order: u8,
    $$n: int,
    vl: $RT,
    vu: RT,
    $$il: int,
    $$iu: int,
    abstol: RT,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    m: *int,
    nsplit: *int,
    w: *$W/ArrayType1D,
    iblock: *$IBLOCK/ArrayType1D,
    isplit: *$ISPLIT/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           E.data_type == RT &&
           W.data_type == RT &&
           WORK.data_type == RT &&
           IBLOCK.data_type == int &&
           ISPLIT.data_type == int &&
           IWORK.data_type == int;
} {
    #if RT == float32 {
        sstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info);
    } else {
        dstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info);
    }
}

_hbevx_stein :: (
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    $$m: int,
    w: *$W/ArrayType1D,
    iblock: *$IBLOCK/ArrayType1D,
    isplit: *$ISPLIT/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    ifail: *$IFAIL/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT &&
           E.data_type == RT &&
           W.data_type == RT &&
           WORK.data_type == RT &&
           Z.data_type == MT &&
           IBLOCK.data_type == int &&
           ISPLIT.data_type == int &&
           IWORK.data_type == int &&
           IFAIL.data_type == int;
} {
    #if MT == C32 {
        cstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info);
    } else {
        zstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info);
    }
}

_hbevx_swap :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_hbevx :: (
    $$jobz: u8,
    $$range: u8,
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    vl: $RT,
    vu: RT,
    $$il: int,
    $$iu: int,
    abstol: RT,
    m: *int,
    w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    ifail: *$IFAIL/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           Q.data_type == MT &&
           W.data_type == RT &&
           Z.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int &&
           IFAIL.data_type == int;
} {
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    czero_t :: #run _hbevx_complex_from_real(zero_r, RT, MT);
    cone_t  :: #run _hbevx_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHBEVX" } else { "ZHBEVX" };

    abstll, anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum, tmp1, vll, vuu : RT = ---;
    ctmp1 : MT = ---;
    alleig, indeig, lower, test, use_fallback, valeig, wantz : bool = ---;
    i, iinfo, imax, inde, indee, indd, indibl, indisp, indiwk, indrwk, indwrk, iscale, itmp1, j, jj, nsplit : int = ---;
    order : u8 = ---;

    wantz  = lsame(jobz, #char "V");
    alleig = lsame(range, #char "A");
    valeig = lsame(range, #char "V");
    indeig = lsame(range, #char "I");
    lower  = lsame(uplo, #char "L");

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(alleig || valeig || indeig) {
        info.* = -2;
    } else if !(lower || lsame(uplo, #char "U")) {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if kd < 0 {
        info.* = -5;
    } else if ldab < kd + 1 {
        info.* = -7;
    } else if wantz && ldq < max(1, n) {
        info.* = -9;
    } else {
        if valeig {
            if n > 0 && vu <= vl {
                info.* = -11;
            }
        } else if indeig {
            if il < 1 || il > max(1, n) {
                info.* = -12;
            } else if iu < min(n, il) || iu > n {
                info.* = -13;
            }
        }
    }
    if info.* == 0 {
        if ldz < 1 || (wantz && ldz < n) {
            info.* = -18;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    m.* = 0;
    if n == 0 {
        return;
    }

    if n == 1 {
        m.* = 1;
        if lower {
            ctmp1 = fget(ab.*, 1, 1);
        } else {
            ctmp1 = fget(ab.*, kd + 1, 1);
        }
        tmp1 = real(ctmp1);
        if valeig {
            if !(vl < tmp1 && vu >= tmp1) {
                m.* = 0;
            }
        }
        if m.* == 1 {
            fset(w, 1, tmp1);
            if wantz {
                fset(z, 1, 1, cone_t);
            }
        }
        return;
    }

    safmin = _hbevx_lamch("Safe minimum"s, RT);
    eps    = _hbevx_lamch("Precision"s, RT);
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin   = sqrt(smlnum);
    rmax   = min(sqrt(bignum), one_r / sqrt(sqrt(safmin)));

    iscale = 0;
    abstll = abstol;
    if valeig {
        vll = vl;
        vuu = vu;
    } else {
        vll = zero_r;
        vuu = zero_r;
    }
    anrm = _hbevx_lanhb(#char "M", uplo, n, kd, ab.*, ldab, rwork, RT, MT);
    if anrm > zero_r && anrm < rmin {
        iscale = 1;
        sigma = rmin / anrm;
    } else if anrm > rmax {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if iscale == 1 {
        if lower {
            _hbevx_lascl(#char "B", kd, kd, one_r, sigma, n, n, ab, ldab, info, RT, MT);
        } else {
            _hbevx_lascl(#char "Q", kd, kd, one_r, sigma, n, n, ab, ldab, info, RT, MT);
        }
        if abstol > zero_r {
            abstll = abstol * sigma;
        }
        if valeig {
            vll = vl * sigma;
            vuu = vu * sigma;
        }
    }

    indd   = 1;
    inde   = indd + n;
    indrwk = inde + n;
    indwrk = 1;

    d_view := _hbevx_subvector(rwork.*, indd, n);
    e_view := _hbevx_subvector(rwork.*, inde, n);
    work_view := _hbevx_subvector(work.*, indwrk, n);
    _hbevx_hbtrd(jobz, uplo, n, kd, ab, ldab, *d_view, *e_view, q, ldq, *work_view, *iinfo, RT, MT);

    test = false;
    if indeig {
        if il == 1 && iu == n {
            test = true;
        }
    }

    use_fallback = true;
    if (alleig || test) && (abstol <= zero_r) {
        _hbevx_copy_real(n, *d_view, 1, w, 1, RT);
        indee = indrwk + 2*n;
        ee_view := _hbevx_subvector(rwork.*, indee, max(0, n - 1));
        if !wantz {
            _hbevx_copy_real(n - 1, *e_view, 1, *ee_view, 1, RT);
            _hbevx_sterf(n, w, *ee_view, info, RT);
        } else {
            _hbevx_lacpy(#char "A", n, n, q.*, ldq, z, ldz, MT);
            _hbevx_copy_real(n - 1, *e_view, 1, *ee_view, 1, RT);
            rwork_tail := _hbevx_subvector(rwork.*, indrwk, 7*n - indrwk + 1);
            _hbevx_steqr(jobz, n, w, *ee_view, z, ldz, *rwork_tail, info, RT, MT);
            if info.* == 0 {
                for i: 1..n {
                    fset(ifail, i, 0);
                }
            }
        }
        if info.* == 0 {
            m.* = n;
            use_fallback = false;
        } else {
            info.* = 0;
        }
    }

    if use_fallback {
        if wantz {
            order = #char "B";
        } else {
            order = #char "E";
        }

        indibl = 1;
        indisp = indibl + n;
        indiwk = indisp + n;

        iblock_view := _hbevx_subvector(iwork.*, indibl, n);
        isplit_view := _hbevx_subvector(iwork.*, indisp, n);
        iwork_tail := _hbevx_subvector(iwork.*, indiwk, 5*n - indiwk + 1);
        rwork_tail := _hbevx_subvector(rwork.*, indrwk, 7*n - indrwk + 1);

        _hbevx_stebz(range, order, n, vll, vuu, il, iu, abstll, *d_view, *e_view, m, *nsplit, w, *iblock_view, *isplit_view, *rwork_tail, *iwork_tail, info, RT);

        if wantz {
            _hbevx_stein(n, *d_view, *e_view, m.*, w, *iblock_view, *isplit_view, z, ldz, *rwork_tail, *iwork_tail, ifail, info, RT, MT);

            for j: 1..m.* {
                work_col := _hbevx_subvector(work.*, 1, n);
                z_col := _hbevx_matrix_column(z.*, ldz, j, n);
                _hbevx_copy_complex(n, *z_col, 1, *work_col, 1, MT);
                _hbevx_gemv(#char "N", n, n, cone_t, q.*, ldq, *work_col, 1, czero_t, *z_col, 1, MT);
            }
        }
    }

    if iscale == 1 {
        if info.* == 0 {
            imax = m.*;
        } else {
            imax = info.* - 1;
        }
        _hbevx_scal(imax, one_r / sigma, w, 1, RT);
    }

    if wantz {
        indibl = 1;
        for j: 1..m.* - 1 {
            i = 0;
            tmp1 = fget(w, j);
            for jj: j + 1..m.* {
                if fget(w, jj) < tmp1 {
                    i = jj;
                    tmp1 = fget(w, jj);
                }
            }

            if i != 0 {
                itmp1 = fget(iwork, indibl + i - 1);
                fset(w, i, fget(w, j));
                fset(iwork, indibl + i - 1, fget(iwork, indibl + j - 1));
                fset(w, j, tmp1);
                fset(iwork, indibl + j - 1, itmp1);

                zi_col := _hbevx_matrix_column(z.*, ldz, i, n);
                zj_col := _hbevx_matrix_column(z.*, ldz, j, n);
                _hbevx_swap(n, *zi_col, 1, *zj_col, 1, MT);

                if info.* != 0 {
                    itmp1 = fget(ifail, i);
                    fset(ifail, i, fget(ifail, j));
                    fset(ifail, j, itmp1);
                }
            }
        }
    }
}
