// ? LAPACK *gerqf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgerqf/cgerqf/dgerqf/zgerqf
// *GERQF computes an RQ factorization of a general matrix.

sgerqf :: #bake_arguments _gerqf(RT=float32, MT=float32);
dgerqf :: #bake_arguments _gerqf(RT=float64, MT=float64);
cgerqf :: #bake_arguments _gerqf(RT=float32, MT=C32);
zgerqf :: #bake_arguments _gerqf(RT=float64, MT=C64);

#scope_file
_gerqf_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gerqf_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gerqf_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gerqf_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gerqf_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gerqf_set_work1 :: inline (work: *$WORK/ArrayType1D, $$value: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           WORK.data_type == MT;
} {
    #if MT == float32 || MT == float64 {
        fset(work, 1, ncast(MT, value));
    } else {
        fset(work, 1, _gerqf_complex_from_real(ncast(RT, value), RT, MT));
    }
}

_gerqf_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_gerqf_larft :: (
    $$direct: u8,
    $$storev: u8,
    $$n: int,
    $$k: int,
    v: $V/ArrayType2D,
    $$ldv: int,
    tau: $TAU/ArrayType1D,
    t: *$T/ArrayType2D,
    $$ldt: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT &&
           TAU.data_type == MT &&
           T.data_type == MT;
} {
    #if MT == float32 {
        slarft(direct, storev, n, k, v, ldv, tau, t, ldt);
    } else if MT == float64 {
        dlarft(direct, storev, n, k, v, ldv, tau, t, ldt);
    } else if MT == C32 {
        clarft(direct, storev, n, k, v, ldv, tau, t, ldt);
    } else {
        zlarft(direct, storev, n, k, v, ldv, tau, t, ldt);
    }
}

_gerqf_larfb :: (
    $$side: u8,
    $$trans: u8,
    $$direct: u8,
    $$storev: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    v: $V/ArrayType2D,
    $$ldv: int,
    t: $T/ArrayType2D,
    $$ldt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType2D,
    $$ldwork: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else if MT == float64 {
        dlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else if MT == C32 {
        clarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else {
        zlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    }
}

_gerqf_gerq2 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgerq2(m, n, a, lda, tau, work, info);
    } else if MT == float64 {
        dgerq2(m, n, a, lda, tau, work, info);
    } else if MT == C32 {
        cgerq2(m, n, a, lda, tau, work, info);
    } else {
        zgerq2(m, n, a, lda, tau, work, info);
    }
}

_gerqf :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
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
           WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGERQF"
    } else if MT == float64 {
        "DGERQF"
    } else if MT == C32 {
        "CGERQF"
    } else {
        "ZGERQF"
    };

    lquery : bool = ---;
    i, ib, iinfo, iws, k, ki, kk, ldwork, lwkopt, mu, nb, nbmin, nu, nx : int = ---;

    info.* = 0;
    lquery = lwork == -1;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    }

    if info.* == 0 {
        k = min(m, n);
        if k == 0 {
            lwkopt = 1;
        } else {
            nb = ilaenv(1, routine_name, " "s, m, n, -1, -1);
            lwkopt = m*nb;
        }

        #if MT == float32 || MT == float64 {
            fset(work, 1, _gerqf_roundup_lwork(ncast(RT, lwkopt), RT));
        } else {
            fset(work, 1, _gerqf_complex_from_real(_gerqf_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
        }

        if !lquery {
            if lwork <= 0 || (n > 0 && lwork < max(1, m)) {
                info.* = -7;
            }
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if k == 0 {
        _gerqf_set_work1(work, 1, RT, MT);
        return;
    }

    nbmin = 2;
    nx = 1;
    iws = m;
    if nb > 1 && nb < k {
        nx = max(0, ilaenv(3, routine_name, " "s, m, n, -1, -1));
        if nx < k {
            ldwork = m;
            iws = ldwork*nb;
            if lwork < iws {
                nb = lwork / ldwork;
                nbmin = max(2, ilaenv(2, routine_name, " "s, m, n, -1, -1));
            }
        }
    }

    if nb >= nbmin && nb < k && nx < k {
        ki = ((k - nx - 1) / nb) * nb;
        kk = min(k, ki + nb);

        i = k - kk + ki + 1;
        while i >= k - kk + 1 {
            ib = min(k - i + 1, nb);

            a_block_raw := _gerqf_matrix_storage(a.*, lda, m - k + i, 1, ib, n - k + i + ib - 1);
            a_block := array2d_view(a_block_raw, lda, row_major=false);
            tau_block := _gerqf_subvector(tau.*, i, ib);
            _gerqf_gerq2(ib, n - k + i + ib - 1, *a_block, lda, *tau_block, work, *iinfo, RT, MT);

            if m - k + i > 1 {
                t_block_raw := _gerqf_raw_slice(work.data, ldwork*ib);
                t_block := array2d_view(t_block_raw, ldwork, row_major=false);
                _gerqf_larft(#char "B", #char "R", n - k + i + ib - 1, ib, a_block.*, lda, tau_block, *t_block, ldwork, MT);

                c_raw := _gerqf_matrix_storage(a.*, lda, 1, 1, m - k + i - 1, n - k + i + ib - 1);
                c_view := array2d_view(c_raw, lda, row_major=false);
                work2_raw := _gerqf_raw_slice(work.data + ib, ldwork*(m - k + i - 1));
                work2_view := array2d_view(work2_raw, ldwork, row_major=false);
                _gerqf_larfb(#char "R", #char "N", #char "B", #char "R", m - k + i - 1, n - k + i + ib - 1, ib,
                             a_block.*, lda, t_block.*, ldwork, *c_view, lda, *work2_view, ldwork, MT);
            }

            i -= nb;
        }

        mu = m - k + i + nb - 1;
        nu = n - k + i + nb - 1;
    } else {
        mu = m;
        nu = n;
    }

    if mu > 0 && nu > 0 {
        _gerqf_gerq2(mu, nu, a, lda, tau, work, *iinfo, RT, MT);
    }

    #if MT == float32 || MT == float64 {
        fset(work, 1, _gerqf_roundup_lwork(ncast(RT, iws), RT));
    } else {
        fset(work, 1, _gerqf_complex_from_real(_gerqf_roundup_lwork(ncast(RT, iws), RT), RT, MT));
    }
}
