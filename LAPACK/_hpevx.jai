// ? LAPACK *hpevx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpevx/zhpevx
// *HPEVX computes selected eigenvalues and, optionally, eigenvectors of a Hermitian packed matrix.

chpevx :: #bake_arguments _hpevx(RT=float32, MT=C32);
zhpevx :: #bake_arguments _hpevx(RT=float64, MT=C64);

#scope_file
_hpevx_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hpevx_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld;
    view.count = count;
    return view;
}

_hpevx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hpevx_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(cmach.data); } else { return dlamch(cmach.data); }
}

_hpevx_lanhp :: (
    $$norm: u8, $$uplo: u8, $$n: int, ap: $AP/ArrayType1D, work: *$RWORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { return clanhp(norm, uplo, n, ap, work); } else { return zlanhp(norm, uplo, n, ap, work); }
}

_hpevx_ptrd :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, d: *$D/ArrayType1D, e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && D.data_type == RT && E.data_type == RT && TAU.data_type == MT;
} {
    #if MT == C32 { chptrd(uplo, n, ap, d, e, tau, info); } else { zhptrd(uplo, n, ap, d, e, tau, info); }
}

_hpevx_scal_packed :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hpevx_upgtr :: (
    $$uplo: u8, $$n: int, ap: $AP/ArrayType1D, tau: $TAU/ArrayType1D, q: *$Q/ArrayType2D,
    $$ldq: int, work: *$WORK/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           AP.data_type == MT && TAU.data_type == MT && Q.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 { cupgtr(uplo, n, ap, tau, q, ldq, work, info); } else { zupgtr(uplo, n, ap, tau, q, ldq, work, info); }
}

_hpevx_upmtr :: (
    $$side: u8, $$uplo: u8, $$trans: u8, $$m: int, $$n: int, ap: $AP/ArrayType1D, tau: $TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           AP.data_type == MT && TAU.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 { cupmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info); } else { zupmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info); }
}

_hpevx_steqr :: (
    $$compz: u8, $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, z: *$Z/ArrayType2D,
    $$ldz: int, work: *$RWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && Z.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { csteqr(compz, n, d, e, z, ldz, work, info); } else { zsteqr(compz, n, d, e, z, ldz, work, info); }
}

_hpevx_stein :: (
    $$n: int, d: $D/ArrayType1D, e: $E/ArrayType1D, $$m: int, w: $W/ArrayType1D,
    iblock: $IBLOCK/ArrayType1D, isplit: $ISPLIT/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && W.data_type == RT &&
           IBLOCK.data_type == int && ISPLIT.data_type == int && Z.data_type == MT &&
           RWORK.data_type == RT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    #if MT == C32 {
        cstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info);
    } else {
        zstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info);
    }
}

_hpevx_stebz :: (
    $$range: u8, $$order: u8, $$n: int, vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT,
    d: $D/ArrayType1D, e: $E/ArrayType1D, m: *int, nsplit: *int, w: *$W/ArrayType1D,
    iblock: *$IBLOCK/ArrayType1D, isplit: *$ISPLIT/ArrayType1D, work: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT && E.data_type == RT && W.data_type == RT &&
           IBLOCK.data_type == int && ISPLIT.data_type == int && RWORK.data_type == RT && IWORK.data_type == int;
} {
    #if RT == float32 {
        sstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info);
    } else {
        dstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info);
    }
}

_hpevx_sterf :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && D.data_type == RT && E.data_type == RT;
} {
    #if RT == float32 { ssterf(n, d, e, info); } else { dsterf(n, d, e, info); }
}

_hpevx_copy_real :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 { scopy(n, x, incx, y, incy); } else { dcopy(n, x, incx, y, incy); }
}

_hpevx_swap_cols :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hpevx_scal_real :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 { sscal(n, alpha, x, incx); } else { dscal(n, alpha, x, incx); }
}

_hpevx :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, vl: $RT, vu: RT,
    $$il: int, $$iu: int, abstol: RT, m: *int, w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && W.data_type == RT && Z.data_type == MT && WORK.data_type == MT &&
           RWORK.data_type == RT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    cone   :: #run _hpevx_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHPEVX" } else { "ZHPEVX" };

    alleig, indeig, test, valeig, wantz : bool = ---;
    order : u8 = ---;
    i, iinfo, imax, indd, inde, indee, inddisp, indiwrk, inddrwk, indtau, indwrk, iscale, itmp1, j, jj, nsplit : int = ---;
    abstll, anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum, tmp1, vll, vuu : RT = ---;

    wantz = lsame(jobz, #char "V");
    alleig = lsame(range, #char "A");
    valeig = lsame(range, #char "V");
    indeig = lsame(range, #char "I");

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(alleig || valeig || indeig) {
        info.* = -2;
    } else if !(lsame(uplo, #char "L") || lsame(uplo, #char "U")) {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else {
        if valeig {
            if n > 0 && vu <= vl { info.* = -7; }
        } else if indeig {
            if il < 1 || il > max(1, n) {
                info.* = -8;
            } else if iu < min(n, il) || iu > n {
                info.* = -9;
            }
        }
    }
    if info.* == 0 {
        if ldz < 1 || (wantz && ldz < n) { info.* = -14; }
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    m.* = 0;
    if n == 0 { return; }
    if n == 1 {
        if alleig || indeig {
            m.* = 1;
            fset(w, 1, real(fget(ap.*, 1)));
        } else if vl < real(fget(ap.*, 1)) && vu >= real(fget(ap.*, 1)) {
            m.* = 1;
            fset(w, 1, real(fget(ap.*, 1)));
        }
        if wantz { fset(z, 1, 1, cone); }
        return;
    }

    safmin = _hpevx_lamch("Safe minimum", RT);
    eps = _hpevx_lamch("Precision", RT);
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin = sqrt(smlnum);
    rmax = min(sqrt(bignum), one_r / sqrt(sqrt(safmin)));

    iscale = 0;
    abstll = abstol;
    if valeig {
        vll = vl;
        vuu = vu;
    } else {
        vll = zero_r;
        vuu = zero_r;
    }
    anrm = _hpevx_lanhp(#char "M", uplo, n, ap.*, rwork, RT, MT);
    if anrm > zero_r && anrm < rmin {
        iscale = 1;
        sigma = rmin / anrm;
    } else if anrm > rmax {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if iscale == 1 {
        _hpevx_scal_packed((n*(n + 1)) / 2, sigma, ap, 1, RT, MT);
        if abstol > zero_r { abstll = abstol * sigma; }
        if valeig {
            vll = vl * sigma;
            vuu = vu * sigma;
        }
    }

    indd = 1;
    inde = indd + n;
    inddrwk = inde + n;
    indtau = 1;
    indwrk = indtau + n;
    d_view := _hpevx_subvector(rwork.*, indd, n);
    e_view := _hpevx_subvector(rwork.*, inde, n - 1);
    tau_view := _hpevx_subvector(work.*, indtau, n - 1);
    _hpevx_ptrd(uplo, n, ap, *d_view, *e_view, *tau_view, *iinfo, RT, MT);

    test = false;
    if indeig && il == 1 && iu == n {
        test = true;
    }
    if (alleig || test) && (abstol <= zero_r) {
        _hpevx_copy_real(n, d_view, 1, w, 1, RT);
        indee = inddrwk + 2*n;
        e2_view := _hpevx_subvector(rwork.*, indee, n - 1);
        if !wantz {
            _hpevx_copy_real(n - 1, e_view, 1, *e2_view, 1, RT);
            _hpevx_sterf(n, w, *e2_view, info, RT);
        } else {
            gen_view := _hpevx_subvector(work.*, indwrk, max(1, n - 1));
            _hpevx_upgtr(uplo, n, ap.*, tau_view, z, ldz, *gen_view, *iinfo, MT);
            _hpevx_copy_real(n - 1, e_view, 1, *e2_view, 1, RT);
            qr_view := _hpevx_subvector(rwork.*, inddrwk, max(1, 2*n - 1));
            _hpevx_steqr(jobz, n, w, *e2_view, z, ldz, *qr_view, info, RT, MT);
            if info.* == 0 {
                for i:1..n {
                    fset(ifail, i, 0);
                }
            }
        }
        if info.* == 0 {
            m.* = n;
            goto scale_done;
        }
        info.* = 0;
    }

    if wantz {
        order = #char "B";
    } else {
        order = #char "E";
    }
    inddisp = 1 + n;
    indiwrk = inddisp + n;
    iblock := _hpevx_subvector(iwork.*, 1, n);
    isplit := _hpevx_subvector(iwork.*, inddisp, n);
    iwork2 := _hpevx_subvector(iwork.*, indiwrk, 5*n - indiwrk + 1);
    rwork2 := _hpevx_subvector(rwork.*, inddrwk, 7*n - inddrwk + 1);
    _hpevx_stebz(range, order, n, vll, vuu, il, iu, abstll, d_view, e_view, m, *nsplit, w, *iblock, *isplit, *rwork2, *iwork2, info, RT);

    if wantz {
        _hpevx_stein(n, d_view, e_view, m.*, w.*, iblock, isplit, z, ldz, *rwork2, *iwork2, ifail, info, RT, MT);
        gen_view := _hpevx_subvector(work.*, indwrk, max(1, n - 1));
        _hpevx_upmtr(#char "L", uplo, #char "N", n, m.*, ap.*, tau_view, z, ldz, *gen_view, *iinfo, MT);
    }

scale_done:
    if iscale == 1 {
        if info.* == 0 { imax = m.*; } else { imax = info.* - 1; }
        _hpevx_scal_real(imax, one_r / sigma, w, 1, RT);
    }

    if wantz {
        for j:1..m.* - 1 {
            i = 0;
            tmp1 = fget(w.*, j);
            for jj:j + 1..m.* {
                if fget(w.*, jj) < tmp1 {
                    i = jj;
                    tmp1 = fget(w.*, jj);
                }
            }
            if i != 0 {
                itmp1 = fget(iwork.*, i);
                fset(w, i, fget(w.*, j));
                fset(iwork, i, fget(iwork.*, j));
                fset(w, j, tmp1);
                fset(iwork, j, itmp1);
                zi := _hpevx_matrix_column(z.*, ldz, i, n);
                zj := _hpevx_matrix_column(z.*, ldz, j, n);
                _hpevx_swap_cols(n, *zi, 1, *zj, 1, MT);
                if info.* != 0 {
                    itmp1 = fget(ifail.*, i);
                    fset(ifail, i, fget(ifail.*, j));
                    fset(ifail, j, itmp1);
                }
            }
        }
    }
}
