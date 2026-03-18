// ? LAPACK *gels()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgels/cgels/dgels/zgels
// *GELS solves overdetermined or underdetermined systems using QR or LQ.

sgels :: #bake_arguments _gels(RT=float32, MT=float32);
dgels :: #bake_arguments _gels(RT=float64, MT=float64);
cgels :: #bake_arguments _gels(RT=float32, MT=C32);
zgels :: #bake_arguments _gels(RT=float64, MT=C64);

#scope_file
_gels_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gels_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gels_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gels_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gels_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gels_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gels_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_gels_lange :: (
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

_gels_lascl :: (
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

_gels_laset :: (
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

_gels_geqrf :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == float64 {
        dgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == C32 {
        cgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else {
        zgeqrf(m, n, a, lda, tau, work, lwork, info);
    }
}

_gels_gelqf :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgelqf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == float64 {
        dgelqf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == C32 {
        cgelqf(m, n, a, lda, tau, work, lwork, info);
    } else {
        zgelqf(m, n, a, lda, tau, work, lwork, info);
    }
}

_gels_unmqr :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_gels_unmlq :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sormlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dormlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunmlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_gels_trtrs :: (
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

_gels_roundup_lwork_real :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_gels :: (
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
        "SGELS "
    } else if MT == float64 {
        "DGELS "
    } else if MT == C32 {
        "CGELS "
    } else {
        "ZGELS "
    };

    lquery, tpsd : bool = ---;
    brow, i, iascl, ibscl, j, mn, nb, scllen, wsize : int = ---;
    anrm, bignum, bnrm, smlnum : RT = ---;
    rwork_data: [1] RT;
    tau_view: [] MT;
    tail_work: [] MT;
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

        if m >= n {
            #if MT == float32 {
                nb = ilaenv(1, "SGEQRF", " "s, m, n, -1, -1);
                if tpsd {
                    nb = max(nb, ilaenv(1, "SORMQR", "LN"s, m, nrhs, n, -1));
                } else {
                    nb = max(nb, ilaenv(1, "SORMQR", "LT"s, m, nrhs, n, -1));
                }
            } else if MT == float64 {
                nb = ilaenv(1, "DGEQRF", " "s, m, n, -1, -1);
                if tpsd {
                    nb = max(nb, ilaenv(1, "DORMQR", "LN"s, m, nrhs, n, -1));
                } else {
                    nb = max(nb, ilaenv(1, "DORMQR", "LT"s, m, nrhs, n, -1));
                }
            } else if MT == C32 {
                nb = ilaenv(1, "CGEQRF", " "s, m, n, -1, -1);
                if tpsd {
                    nb = max(nb, ilaenv(1, "CUNMQR", "LN"s, m, nrhs, n, -1));
                } else {
                    nb = max(nb, ilaenv(1, "CUNMQR", "LC"s, m, nrhs, n, -1));
                }
            } else {
                nb = ilaenv(1, "ZGEQRF", " "s, m, n, -1, -1);
                if tpsd {
                    nb = max(nb, ilaenv(1, "ZUNMQR", "LN"s, m, nrhs, n, -1));
                } else {
                    nb = max(nb, ilaenv(1, "ZUNMQR", "LC"s, m, nrhs, n, -1));
                }
            }
        } else {
            #if MT == float32 {
                nb = ilaenv(1, "SGELQF", " "s, m, n, -1, -1);
                if tpsd {
                    nb = max(nb, ilaenv(1, "SORMLQ", "LT"s, n, nrhs, m, -1));
                } else {
                    nb = max(nb, ilaenv(1, "SORMLQ", "LN"s, n, nrhs, m, -1));
                }
            } else if MT == float64 {
                nb = ilaenv(1, "DGELQF", " "s, m, n, -1, -1);
                if tpsd {
                    nb = max(nb, ilaenv(1, "DORMLQ", "LT"s, n, nrhs, m, -1));
                } else {
                    nb = max(nb, ilaenv(1, "DORMLQ", "LN"s, n, nrhs, m, -1));
                }
            } else if MT == C32 {
                nb = ilaenv(1, "CGELQF", " "s, m, n, -1, -1);
                if tpsd {
                    nb = max(nb, ilaenv(1, "CUNMLQ", "LC"s, n, nrhs, m, -1));
                } else {
                    nb = max(nb, ilaenv(1, "CUNMLQ", "LN"s, n, nrhs, m, -1));
                }
            } else {
                nb = ilaenv(1, "ZGELQF", " "s, m, n, -1, -1);
                if tpsd {
                    nb = max(nb, ilaenv(1, "ZUNMLQ", "LC"s, n, nrhs, m, -1));
                } else {
                    nb = max(nb, ilaenv(1, "ZUNMLQ", "LN"s, n, nrhs, m, -1));
                }
            }
        }

        wsize = max(1, mn + max(mn, nrhs)*nb);
        #if MT == float32 || MT == float64 {
            fset(work, 1, _gels_roundup_lwork_real(ncast(RT, wsize), RT));
        } else {
            fset(work, 1, _gels_complex_from_real(ncast(RT, wsize), RT, MT));
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if min(m, n, nrhs) == 0 {
        _gels_laset(#char "F", max(m, n), nrhs, zero_t, zero_t, b, ldb, MT);
        return;
    }

    smlnum = _gels_lamch_safe_minimum(RT) / _gels_lamch_precision(RT);
    bignum = one(RT) / smlnum;

    anrm = _gels_lange(#char "M", m, n, a.*, lda, *rwork_data, RT, MT);
    iascl = 0;
    if anrm > zero(RT) && anrm < smlnum {
        _gels_lascl(#char "G", 0, 0, anrm, smlnum, m, n, a, lda, info, RT, MT);
        iascl = 1;
    } else if anrm > bignum {
        _gels_lascl(#char "G", 0, 0, anrm, bignum, m, n, a, lda, info, RT, MT);
        iascl = 2;
    } else if anrm == zero(RT) {
        _gels_laset(#char "F", max(m, n), nrhs, zero_t, zero_t, b, ldb, MT);
        early_zero = true;
    }

    if !early_zero {
        brow = m;
        if tpsd {
            brow = n;
        }
        bnrm = _gels_lange(#char "M", brow, nrhs, b.*, ldb, *rwork_data, RT, MT);
        ibscl = 0;
        if bnrm > zero(RT) && bnrm < smlnum {
            _gels_lascl(#char "G", 0, 0, bnrm, smlnum, brow, nrhs, b, ldb, info, RT, MT);
            ibscl = 1;
        } else if bnrm > bignum {
            _gels_lascl(#char "G", 0, 0, bnrm, bignum, brow, nrhs, b, ldb, info, RT, MT);
            ibscl = 2;
        }

        tau_view = _gels_subvector(work.*, 1, mn);
        tail_work = _gels_subvector(work.*, mn + 1, max(0, lwork - mn));

        if m >= n {
            _gels_geqrf(m, n, a, lda, *tau_view, *tail_work, lwork - mn, info, MT);

            if !tpsd {
                #if MT == float32 || MT == float64 {
                    _gels_unmqr(#char "L", #char "T", m, nrhs, n, a, lda, *tau_view, b, ldb, *tail_work, lwork - mn, info, RT, MT);
                    _gels_trtrs(#char "U", #char "N", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
                } else {
                    _gels_unmqr(#char "L", #char "C", m, nrhs, n, a, lda, *tau_view, b, ldb, *tail_work, lwork - mn, info, RT, MT);
                    _gels_trtrs(#char "U", #char "N", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
                }
                if info.* > 0 {
                    return;
                }
                scllen = n;
            } else {
                #if MT == float32 || MT == float64 {
                    _gels_trtrs(#char "U", #char "T", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
                } else {
                    _gels_trtrs(#char "U", #char "C", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
                }
                if info.* > 0 {
                    return;
                }

                for j:1..nrhs {
                    for i:n + 1..m {
                        fset(b, i, j, zero_t);
                    }
                }

                _gels_unmqr(#char "L", #char "N", m, nrhs, n, a, lda, *tau_view, b, ldb, *tail_work, lwork - mn, info, RT, MT);
                scllen = m;
            }
        } else {
            _gels_gelqf(m, n, a, lda, *tau_view, *tail_work, lwork - mn, info, MT);

            if !tpsd {
                _gels_trtrs(#char "L", #char "N", #char "N", m, nrhs, a, lda, b, ldb, info, MT);
                if info.* > 0 {
                    return;
                }

                for j:1..nrhs {
                    for i:m + 1..n {
                        fset(b, i, j, zero_t);
                    }
                }

                #if MT == float32 || MT == float64 {
                    _gels_unmlq(#char "L", #char "T", n, nrhs, m, a, lda, *tau_view, b, ldb, *tail_work, lwork - mn, info, RT, MT);
                } else {
                    _gels_unmlq(#char "L", #char "C", n, nrhs, m, a, lda, *tau_view, b, ldb, *tail_work, lwork - mn, info, RT, MT);
                }
                scllen = n;
            } else {
                _gels_unmlq(#char "L", #char "N", n, nrhs, m, a, lda, *tau_view, b, ldb, *tail_work, lwork - mn, info, RT, MT);
                #if MT == float32 || MT == float64 {
                    _gels_trtrs(#char "L", #char "T", #char "N", m, nrhs, a, lda, b, ldb, info, MT);
                } else {
                    _gels_trtrs(#char "L", #char "C", #char "N", m, nrhs, a, lda, b, ldb, info, MT);
                }
                if info.* > 0 {
                    return;
                }
                scllen = m;
            }
        }

        if iascl == 1 {
            _gels_lascl(#char "G", 0, 0, anrm, smlnum, scllen, nrhs, b, ldb, info, RT, MT);
        } else if iascl == 2 {
            _gels_lascl(#char "G", 0, 0, anrm, bignum, scllen, nrhs, b, ldb, info, RT, MT);
        }
        if ibscl == 1 {
            _gels_lascl(#char "G", 0, 0, smlnum, bnrm, scllen, nrhs, b, ldb, info, RT, MT);
        } else if ibscl == 2 {
            _gels_lascl(#char "G", 0, 0, bignum, bnrm, scllen, nrhs, b, ldb, info, RT, MT);
        }
    }

    #if MT == float32 || MT == float64 {
        fset(work, 1, _gels_roundup_lwork_real(ncast(RT, wsize), RT));
    } else {
        fset(work, 1, _gels_complex_from_real(ncast(RT, wsize), RT, MT));
    }
}
