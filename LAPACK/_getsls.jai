// ? LAPACK *getsls()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgetsls/cgetsls/dgetsls/zgetsls
// *GETSLS solves overdetermined or underdetermined systems
// using the modern tall-skinny QR or short-wide LQ drivers.

sgetsls :: #bake_arguments _getsls(RT=float32, MT=float32);
dgetsls :: #bake_arguments _getsls(RT=float64, MT=float64);
cgetsls :: #bake_arguments _getsls(RT=float32, MT=C32);
zgetsls :: #bake_arguments _getsls(RT=float64, MT=C64);

#scope_file
_getsls_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_getsls_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_getsls_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_getsls_query_value_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_getsls_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_getsls_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_getsls_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_getsls_lange :: (
    $$norm: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    work: $WORK/ArrayType1D,
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

_getsls_lascl :: (
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

_getsls_laset :: (
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

_getsls_ge_factor_query :: (
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    $$minimal: bool,
    t_meta: *$TMETA/ArrayType1D,
    tsize_out: *int,
    lwork_out: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           TMETA.data_type == MT;
} {
    workq: [1] MT = ---;
    workq_view := _getsls_raw_slice(workq.data, 1);
    info2: int = 0;

    #if MT == float32 {
        if m >= n {
            sgeqr(m, n, a, lda, t_meta, ifx minimal then -2 else -1, *workq_view, ifx minimal then -2 else -1, *info2);
        } else {
            sgelq(m, n, a, lda, t_meta, ifx minimal then -2 else -1, *workq_view, ifx minimal then -2 else -1, *info2);
        }
    } else if MT == float64 {
        if m >= n {
            dgeqr(m, n, a, lda, t_meta, ifx minimal then -2 else -1, *workq_view, ifx minimal then -2 else -1, *info2);
        } else {
            dgelq(m, n, a, lda, t_meta, ifx minimal then -2 else -1, *workq_view, ifx minimal then -2 else -1, *info2);
        }
    } else if MT == C32 {
        if m >= n {
            cgeqr(m, n, a, lda, t_meta, ifx minimal then -2 else -1, *workq_view, ifx minimal then -2 else -1, *info2);
        } else {
            cgelq(m, n, a, lda, t_meta, ifx minimal then -2 else -1, *workq_view, ifx minimal then -2 else -1, *info2);
        }
    } else {
        if m >= n {
            zgeqr(m, n, a, lda, t_meta, ifx minimal then -2 else -1, *workq_view, ifx minimal then -2 else -1, *info2);
        } else {
            zgelq(m, n, a, lda, t_meta, ifx minimal then -2 else -1, *workq_view, ifx minimal then -2 else -1, *info2);
        }
    }

    tsize_out.* = _getsls_query_value_to_int(fget(t_meta.*, 1), RT, MT);
    lwork_out.* = _getsls_query_value_to_int(workq[0], RT, MT);
}

_getsls_gem_apply_query :: (
    $$qr_path: bool,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    a: $A/ArrayType2D,
    $$lda: int,
    t_meta: $T/ArrayType1D,
    $$tsize: int,
    c: $C/ArrayType2D,
    $$ldc: int,
    lwork_out: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT;
} {
    workq: [1] MT = ---;
    workq_view := _getsls_raw_slice(workq.data, 1);
    info2: int = 0;

    #if MT == float32 {
        if qr_path {
            sgemqr(#char "L", trans, m, n, k, a, lda, t_meta, tsize, c, ldc, *workq_view, -1, *info2);
        } else {
            sgemlq(#char "L", trans, m, n, k, a, lda, t_meta, tsize, c, ldc, *workq_view, -1, *info2);
        }
    } else if MT == float64 {
        if qr_path {
            dgemqr(#char "L", trans, m, n, k, a, lda, t_meta, tsize, c, ldc, *workq_view, -1, *info2);
        } else {
            dgemlq(#char "L", trans, m, n, k, a, lda, t_meta, tsize, c, ldc, *workq_view, -1, *info2);
        }
    } else if MT == C32 {
        if qr_path {
            cgemqr(#char "L", trans, m, n, k, a, lda, t_meta, tsize, c, ldc, *workq_view, -1, *info2);
        } else {
            cgemlq(#char "L", trans, m, n, k, a, lda, t_meta, tsize, c, ldc, *workq_view, -1, *info2);
        }
    } else {
        if qr_path {
            zgemqr(#char "L", trans, m, n, k, a, lda, t_meta, tsize, c, ldc, *workq_view, -1, *info2);
        } else {
            zgemlq(#char "L", trans, m, n, k, a, lda, t_meta, tsize, c, ldc, *workq_view, -1, *info2);
        }
    }

    lwork_out.* = _getsls_query_value_to_int(workq[0], RT, MT);
}

_getsls_ge_factor :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType1D,
    $$tsize: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        if m >= n {
            sgeqr(m, n, a, lda, t, tsize, work, lwork, info);
        } else {
            sgelq(m, n, a, lda, t, tsize, work, lwork, info);
        }
    } else if MT == float64 {
        if m >= n {
            dgeqr(m, n, a, lda, t, tsize, work, lwork, info);
        } else {
            dgelq(m, n, a, lda, t, tsize, work, lwork, info);
        }
    } else if MT == C32 {
        if m >= n {
            cgeqr(m, n, a, lda, t, tsize, work, lwork, info);
        } else {
            cgelq(m, n, a, lda, t, tsize, work, lwork, info);
        }
    } else {
        if m >= n {
            zgeqr(m, n, a, lda, t, tsize, work, lwork, info);
        } else {
            zgelq(m, n, a, lda, t, tsize, work, lwork, info);
        }
    }
}

_getsls_gem_apply :: (
    $$qr_path: bool,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType1D,
    $$tsize: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        if qr_path {
            sgemqr(#char "L", trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
        } else {
            sgemlq(#char "L", trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
        }
    } else if MT == float64 {
        if qr_path {
            dgemqr(#char "L", trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
        } else {
            dgemlq(#char "L", trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
        }
    } else if MT == C32 {
        if qr_path {
            cgemqr(#char "L", trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
        } else {
            cgemlq(#char "L", trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
        }
    } else {
        if qr_path {
            zgemqr(#char "L", trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
        } else {
            zgemlq(#char "L", trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info);
        }
    }
}

_getsls_trtrs :: (
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

_getsls_transposed_char :: ($MT: Type) -> u8 #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == float32 || MT == float64 {
        return #char "T";
    } else {
        return #char "C";
    }
}

_getsls :: (
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
    routine_name :: #run if MT == float32 {
        "SGETSLS"
    } else if MT == float64 {
        "DGETSLS"
    } else if MT == C32 {
        "CGETSLS"
    } else {
        "ZGETSLS"
    };

    zero_t :: #run _getsls_complex_from_real(zero(RT), RT, MT);
    one_r :: #run one(RT);
    transp_char :: #run _getsls_transposed_char(MT);

    lquery, tran : bool = ---;
    anrm, bignum, bnrm, smlnum : RT = ---;
    b_row, iascl, ibscl, lw1, lw2, lwo, lwm, maxmn, scllen, tszo, tszm, wsizeo, wsizem : int = ---;
    norm_work_storage: [1] RT = ---;
    norm_work := _getsls_raw_slice(norm_work_storage.data, 1);
    tq_storage: [5] MT = ---;
    tq_view := _getsls_raw_slice(tq_storage.data, 5);

    info.* = 0;
    maxmn = max(m, n);
    #if MT == float32 || MT == float64 {
        tran = lsame(trans, #char "T");
        if !(lsame(trans, #char "N") || lsame(trans, #char "T")) {
            info.* = -1;
        }
    } else {
        tran = lsame(trans, #char "C");
        if !(lsame(trans, #char "N") || lsame(trans, #char "C")) {
            info.* = -1;
        }
    }
    lquery = (lwork == -1 || lwork == -2);

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
        }
    }
    
    if info.* == 0 {
        if min(m, n, nrhs) == 0 {
            wsizeo = 1;
            wsizem = 1;
            tszo = 1;
            tszm = 1;
            lwo = 1;
            lwm = 1;
        } else if m >= n {
            _getsls_ge_factor_query(m, n, a.*, lda, false, *tq_view, *tszo, *lwo, RT, MT);
            _getsls_gem_apply_query(true, trans, m, nrhs, n, a.*, lda, tq_view, tszo, b.*, ldb, *lwo, RT, MT);
            _getsls_ge_factor_query(m, n, a.*, lda, true, *tq_view, *tszm, *lwm, RT, MT);
            _getsls_gem_apply_query(true, trans, m, nrhs, n, a.*, lda, tq_view, tszm, b.*, ldb, *lwm, RT, MT);
            wsizeo = tszo + lwo;
            wsizem = tszm + lwm;
        } else {
            _getsls_ge_factor_query(m, n, a.*, lda, false, *tq_view, *tszo, *lwo, RT, MT);
            _getsls_gem_apply_query(false, trans, n, nrhs, m, a.*, lda, tq_view, tszo, b.*, ldb, *lwo, RT, MT);
            _getsls_ge_factor_query(m, n, a.*, lda, true, *tq_view, *tszm, *lwm, RT, MT);
            _getsls_gem_apply_query(false, trans, n, nrhs, m, a.*, lda, tq_view, tszm, b.*, ldb, *lwm, RT, MT);
            wsizeo = tszo + lwo;
            wsizem = tszm + lwm;
        }

        if lwork < wsizem && !lquery {
            info.* = -10;
        }

        fset(work, 1, _getsls_complex_from_real(_getsls_roundup_lwork(ncast(RT, wsizeo), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lquery {
        if lwork == -2 {
            fset(work, 1, _getsls_complex_from_real(_getsls_roundup_lwork(ncast(RT, wsizem), RT), RT, MT));
        }
        return;
    }

    if lwork < wsizeo {
        lw1 = tszm;
        lw2 = lwm;
    } else {
        lw1 = tszo;
        lw2 = lwo;
    }

    if min(m, n, nrhs) == 0 {
        _getsls_laset(#char "F", max(m, n), nrhs, zero_t, zero_t, b, ldb, MT);
        return;
    }

    smlnum = _getsls_lamch_safe_minimum(RT) / _getsls_lamch_precision(RT);
    bignum = one_r / smlnum;

    anrm = _getsls_lange(#char "M", m, n, a.*, lda, norm_work, RT, MT);
    iascl = 0;
    if anrm > zero(RT) && anrm < smlnum {
        _getsls_lascl(#char "G", 0, 0, anrm, smlnum, m, n, a, lda, info, RT, MT);
        iascl = 1;
    } else if anrm > bignum {
        _getsls_lascl(#char "G", 0, 0, anrm, bignum, m, n, a, lda, info, RT, MT);
        iascl = 2;
    } else if anrm == zero(RT) {
        _getsls_laset(#char "F", maxmn, nrhs, zero_t, zero_t, b, ldb, MT);
        fset(work, 1, _getsls_complex_from_real(_getsls_roundup_lwork(ncast(RT, tszo + lwo), RT), RT, MT));
        return;
    }

    b_row = m;
    if tran {
        b_row = n;
    }
    bnrm = _getsls_lange(#char "M", b_row, nrhs, b.*, ldb, norm_work, RT, MT);
    ibscl = 0;
    if bnrm > zero(RT) && bnrm < smlnum {
        _getsls_lascl(#char "G", 0, 0, bnrm, smlnum, b_row, nrhs, b, ldb, info, RT, MT);
        ibscl = 1;
    } else if bnrm > bignum {
        _getsls_lascl(#char "G", 0, 0, bnrm, bignum, b_row, nrhs, b, ldb, info, RT, MT);
        ibscl = 2;
    }

    t_work := _getsls_subvector(work.*, lw2 + 1, lw1);
    op_work := _getsls_subvector(work.*, 1, lw2);

    if m >= n {
        _getsls_ge_factor(m, n, a, lda, *t_work, lw1, *op_work, lw2, info, MT);
        if !tran {
            _getsls_gem_apply(true, transp_char, m, nrhs, n, a, lda, *t_work, lw1, b, ldb, *op_work, lw2, info, MT);
            _getsls_trtrs(#char "U", #char "N", #char "N", n, nrhs, a, lda, b, ldb, info, MT);
            if info.* > 0 {
                return;
            }
            scllen = n;
        } else {
            _getsls_trtrs(#char "U", trans, #char "N", n, nrhs, a, lda, b, ldb, info, MT);
            if info.* > 0 {
                return;
            }

            for j:1..nrhs {
                for i:n + 1..m {
                    fset(b, i, j, zero_t);
                }
            }

            _getsls_gem_apply(true, #char "N", m, nrhs, n, a, lda, *t_work, lw1, b, ldb, *op_work, lw2, info, MT);
            scllen = m;
        }
    } else {
        _getsls_ge_factor(m, n, a, lda, *t_work, lw1, *op_work, lw2, info, MT);
        if !tran {
            _getsls_trtrs(#char "L", #char "N", #char "N", m, nrhs, a, lda, b, ldb, info, MT);
            if info.* > 0 {
                return;
            }

            for j:1..nrhs {
                for i:m + 1..n {
                    fset(b, i, j, zero_t);
                }
            }

            _getsls_gem_apply(false, transp_char, n, nrhs, m, a, lda, *t_work, lw1, b, ldb, *op_work, lw2, info, MT);
            scllen = n;
        } else {
            _getsls_gem_apply(false, #char "N", n, nrhs, m, a, lda, *t_work, lw1, b, ldb, *op_work, lw2, info, MT);
            _getsls_trtrs(#char "L", trans, #char "N", m, nrhs, a, lda, b, ldb, info, MT);
            if info.* > 0 {
                return;
            }
            scllen = m;
        }
    }

    if iascl == 1 {
        _getsls_lascl(#char "G", 0, 0, anrm, smlnum, scllen, nrhs, b, ldb, info, RT, MT);
    } else if iascl == 2 {
        _getsls_lascl(#char "G", 0, 0, anrm, bignum, scllen, nrhs, b, ldb, info, RT, MT);
    }
    if ibscl == 1 {
        _getsls_lascl(#char "G", 0, 0, smlnum, bnrm, scllen, nrhs, b, ldb, info, RT, MT);
    } else if ibscl == 2 {
        _getsls_lascl(#char "G", 0, 0, bignum, bnrm, scllen, nrhs, b, ldb, info, RT, MT);
    }
    fset(work, 1, _getsls_complex_from_real(_getsls_roundup_lwork(ncast(RT, tszo + lwo), RT), RT, MT));
}
