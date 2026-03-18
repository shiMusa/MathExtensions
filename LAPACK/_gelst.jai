// ? LAPACK *gelst()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgelst/cgelst/dgelst/zgelst
// *GELST solves overdetermined or underdetermined systems using blocked QR/LQ.

sgelst :: #bake_arguments _gelst(RT=float32, MT=float32);
dgelst :: #bake_arguments _gelst(RT=float64, MT=float64);
cgelst :: #bake_arguments _gelst(RT=float32, MT=C32);
zgelst :: #bake_arguments _gelst(RT=float64, MT=C64);

#scope_file
_gelst_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gelst_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gelst_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gelst_roundup_lwork_real :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_gelst_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gelst_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_gelst_lange :: (
    $$norm: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           WORK.data_type == RT;
} {
    #if MT == float32 {
        return slange(norm, m, n, a, lda, work);
    } else if MT == float64 {
        return dlange(norm, m, n, a, lda, work);
    } else if MT == C32 {
        return clange(norm, m, n, a, lda, work);
    } else {
        return zlange(norm, m, n, a, lda, work);
    }
}

_gelst_lascl :: (
    $$type: u8,
    $$kl: int,
    $$ku: int,
    cfrom: $RT,
    cto: RT,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        slascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else if MT == float64 {
        dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else if MT == C32 {
        clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_gelst_laset :: (
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

_gelst_geqrt :: (
    $$m: int,
    $$n: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgeqrt(m, n, nb, a, lda, t, ldt, work, info);
    } else if MT == float64 {
        dgeqrt(m, n, nb, a, lda, t, ldt, work, info);
    } else if MT == C32 {
        cgeqrt(m, n, nb, a, lda, t, ldt, work, info);
    } else {
        zgeqrt(m, n, nb, a, lda, t, ldt, work, info);
    }
}

_gelst_gelqt :: (
    $$m: int,
    $$n: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgelqt(m, n, nb, a, lda, t, ldt, work, info);
    } else if MT == float64 {
        dgelqt(m, n, nb, a, lda, t, ldt, work, info);
    } else if MT == C32 {
        cgelqt(m, n, nb, a, lda, t, ldt, work, info);
    } else {
        zgelqt(m, n, nb, a, lda, t, ldt, work, info);
    }
}

_gelst_gemqrt :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    $$nb: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    } else if MT == float64 {
        dgemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    } else if MT == C32 {
        cgemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    } else {
        zgemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    }
}

_gelst_gemlqt :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    $$nb: int,
    v: *$V/ArrayType2D,
    $$ldv: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgemlqt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    } else if MT == float64 {
        dgemlqt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    } else if MT == C32 {
        cgemlqt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    } else {
        zgemlqt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);
    }
}

_gelst_trtrs :: (
    $$uplo: u8,
    $$trans: u8,
    $$diag: u8,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        strtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    } else if MT == float64 {
        dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    } else if MT == C32 {
        ctrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    } else {
        ztrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    }
}

_gelst :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           B.data_type == MT &&
           WORK.data_type == MT;
} {
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SGELST "
    } else if MT == float64 {
        "DGELST "
    } else if MT == C32 {
        "CGELST "
    } else {
        "ZGELST "
    };

    lquery, tpsd : bool = ---;
    brow, i, iascl, ibscl, j, lwopt, mn, mnnrhs, nb, nbmin, scllen : int = ---;
    anrm, bignum, bnrm, smlnum : RT = ---;
    rwork_data: [1] RT;
    t_view : [] MT;
    tail_work : [] MT;
    early_zero : bool = false;

    info.* = 0;
    mn = min(m, n);
    lquery = lwork == -1;
    #if MT == float32 || MT == float64 {
        if !(lsame(trans, #char "N") || lsame(trans, #char "T")) {
            info.* = -1;
        }
    } else {
        if !(lsame(trans, #char "N") || lsame(trans, #char "C")) {
            info.* = -1;
        }
    }
    if info.* == 0 {
        if m < 0 {
            info.* = -2;
        } else if n < 0 {
            info.* = -3;
        } else if nrhs < 0 {
            info.* = -4;
        } else if lda < max(1, m) {
            info.* = -6;
        } else if ldb < max(1, m, n) {
            info.* = -8;
        } else if lwork < max(1, mn + max(mn, nrhs)) && !lquery {
            info.* = -10;
        }
    }

    if info.* == 0 || info.* == -10 {
        tpsd = true;
        if lsame(trans, #char "N") {
            tpsd = false;
        }

        #if MT == float32 {
            nb = ilaenv(1, "SGELST", " "s, m, n, -1, -1);
        } else if MT == float64 {
            nb = ilaenv(1, "DGELST", " "s, m, n, -1, -1);
        } else if MT == C32 {
            nb = ilaenv(1, "CGELST", " "s, m, n, -1, -1);
        } else {
            nb = ilaenv(1, "ZGELST", " "s, m, n, -1, -1);
        }

        mnnrhs = max(mn, nrhs);
        lwopt = max(1, (mn + mnnrhs)*nb);
        #if MT == float32 || MT == float64 {
            fset(work, 1, _gelst_roundup_lwork_real(ncast(RT, lwopt), RT));
        } else {
            fset(work, 1, _gelst_complex_from_real(_gelst_roundup_lwork_real(ncast(RT, lwopt), RT), RT, MT));
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if min(m, n, nrhs) == 0 {
        _gelst_laset(#char "F", max(m, n), nrhs, zero_t, zero_t, b, ldb, MT);
        #if MT == float32 || MT == float64 {
            fset(work, 1, _gelst_roundup_lwork_real(ncast(RT, lwopt), RT));
        } else {
            fset(work, 1, _gelst_complex_from_real(_gelst_roundup_lwork_real(ncast(RT, lwopt), RT), RT, MT));
        }
        return;
    }

    if nb > mn { nb = mn; }
    nb = min(nb, lwork / (mn + mnnrhs));
    #if MT == float32 {
        nbmin = max(2, ilaenv(2, "SGELST", " "s, m, n, -1, -1));
    } else if MT == float64 {
        nbmin = max(2, ilaenv(2, "DGELST", " "s, m, n, -1, -1));
    } else if MT == C32 {
        nbmin = max(2, ilaenv(2, "CGELST", " "s, m, n, -1, -1));
    } else {
        nbmin = max(2, ilaenv(2, "ZGELST", " "s, m, n, -1, -1));
    }
    if nb < nbmin {
        nb = 1;
    }

    smlnum = _gelst_lamch_safe_minimum(RT) / _gelst_lamch_precision(RT);
    bignum = one(RT) / smlnum;

    anrm = _gelst_lange(#char "M", m, n, a.*, lda, *rwork_data, RT, MT);
    iascl = 0;
    if anrm > zero(RT) && anrm < smlnum {
        _gelst_lascl(#char "G", 0, 0, anrm, smlnum, m, n, a, lda, info, RT, MT);
        iascl = 1;
    } else if anrm > bignum {
        _gelst_lascl(#char "G", 0, 0, anrm, bignum, m, n, a, lda, info, RT, MT);
        iascl = 2;
    } else if anrm == zero(RT) {
        _gelst_laset(#char "F", max(m, n), nrhs, zero_t, zero_t, b, ldb, MT);
        early_zero = true;
    }

    if !early_zero {
        brow = m;
        if tpsd {
            brow = n;
        }
        bnrm = _gelst_lange(#char "M", brow, nrhs, b.*, ldb, *rwork_data, RT, MT);
        ibscl = 0;
        if bnrm > zero(RT) && bnrm < smlnum {
            _gelst_lascl(#char "G", 0, 0, bnrm, smlnum, brow, nrhs, b, ldb, info, RT, MT);
            ibscl = 1;
        } else if bnrm > bignum {
            _gelst_lascl(#char "G", 0, 0, bnrm, bignum, brow, nrhs, b, ldb, info, RT, MT);
            ibscl = 2;
        }

        t_view = _gelst_subvector(work.*, 1, mn*nb);
        tail_work = _gelst_subvector(work.*, mn*nb + 1, max(0, lwork - mn*nb));

        if m >= n {
            t_mat := array2d_view(t_view, nb, row_major=false);
            _gelst_geqrt(m, n, nb, a, lda, *t_mat, nb, *tail_work, info, MT);

            if !tpsd {
                #if MT == float32 || MT == float64 {
                    _gelst_gemqrt(#char "L", #char "T", m, nrhs, n, nb, a, lda, *t_mat, nb, b, ldb, *tail_work, info, MT);
                    _gelst_trtrs(#char "U", #char "N", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
                } else {
                    _gelst_gemqrt(#char "L", #char "C", m, nrhs, n, nb, a, lda, *t_mat, nb, b, ldb, *tail_work, info, MT);
                    _gelst_trtrs(#char "U", #char "N", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
                }
                if info.* > 0 { return; }
                scllen = n;
            } else {
                #if MT == float32 || MT == float64 {
                    _gelst_trtrs(#char "U", #char "T", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
                } else {
                    _gelst_trtrs(#char "U", #char "C", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
                }
                if info.* > 0 { return; }
                for j:1..nrhs {
                    for i:n + 1..m {
                        fset(b, i, j, zero_t);
                    }
                }
                _gelst_gemqrt(#char "L", #char "N", m, nrhs, n, nb, a, lda, *t_mat, nb, b, ldb, *tail_work, info, MT);
                scllen = m;
            }
        } else {
            t_mat := array2d_view(t_view, nb, row_major=false);
            _gelst_gelqt(m, n, nb, a, lda, *t_mat, nb, *tail_work, info, MT);

            if !tpsd {
                _gelst_trtrs(#char "L", #char "N", #char "N", m, nrhs, a, lda, b, ldb, info, MT);
                if info.* > 0 { return; }
                for j:1..nrhs {
                    for i:m + 1..n {
                        fset(b, i, j, zero_t);
                    }
                }
                #if MT == float32 || MT == float64 {
                    _gelst_gemlqt(#char "L", #char "T", n, nrhs, m, nb, a, lda, *t_mat, nb, b, ldb, *tail_work, info, MT);
                } else {
                    _gelst_gemlqt(#char "L", #char "C", n, nrhs, m, nb, a, lda, *t_mat, nb, b, ldb, *tail_work, info, MT);
                }
                scllen = n;
            } else {
                _gelst_gemlqt(#char "L", #char "N", n, nrhs, m, nb, a, lda, *t_mat, nb, b, ldb, *tail_work, info, MT);
                #if MT == float32 || MT == float64 {
                    _gelst_trtrs(#char "L", #char "T", #char "N", m, nrhs, a, lda, b, ldb, info, MT);
                } else {
                    _gelst_trtrs(#char "L", #char "C", #char "N", m, nrhs, a, lda, b, ldb, info, MT);
                }
                if info.* > 0 { return; }
                scllen = m;
            }
        }

        if iascl == 1 {
            _gelst_lascl(#char "G", 0, 0, anrm, smlnum, scllen, nrhs, b, ldb, info, RT, MT);
        } else if iascl == 2 {
            _gelst_lascl(#char "G", 0, 0, anrm, bignum, scllen, nrhs, b, ldb, info, RT, MT);
        }
        if ibscl == 1 {
            _gelst_lascl(#char "G", 0, 0, smlnum, bnrm, scllen, nrhs, b, ldb, info, RT, MT);
        } else if ibscl == 2 {
            _gelst_lascl(#char "G", 0, 0, bignum, bnrm, scllen, nrhs, b, ldb, info, RT, MT);
        }
    }

    #if MT == float32 || MT == float64 {
        fset(work, 1, _gelst_roundup_lwork_real(ncast(RT, lwopt), RT));
    } else {
        fset(work, 1, _gelst_complex_from_real(_gelst_roundup_lwork_real(ncast(RT, lwopt), RT), RT, MT));
    }
}
