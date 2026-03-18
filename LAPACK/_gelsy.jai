// ? LAPACK *gelsy()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgelsy/cgelsy/dgelsy/zgelsy
// *GELSY computes the minimum-norm solution to a rank-deficient least
// squares problem using complete orthogonal factorization.

sgelsy :: #bake_arguments _gelsy_real(RT=float32);
dgelsy :: #bake_arguments _gelsy_real(RT=float64);
cgelsy :: #bake_arguments _gelsy_complex(RT=float32, MT=C32);
zgelsy :: #bake_arguments _gelsy_complex(RT=float64, MT=C64);

#scope_file
_gelsy_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gelsy_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda;
    view.count = count;
    return view;
}

_gelsy_roundup_lwork_real :: ($RT: Type, x: RT) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_gelsy_work_query_real :: ($RT: Type, x: int) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    return _gelsy_roundup_lwork_real(RT, ncast(RT, x));
}

_gelsy_work_query_complex :: ($RT: Type, $MT: Type, x: int) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(ncast(RT, x), zero(RT));
}

_gelsy_real :: (
    $$m: int,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    jpvt: *$JPVT/ArrayType1D,
    rcond: $RT,
    rank: *int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           B.data_type == RT &&
           WORK.data_type == RT &&
           JPVT.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SGELSY " } else { "DGELSY " };

    IMAX :: 1;
    IMIN :: 2;
    ZERO :: #run zero(RT);
    ONE  :: #run one(RT);

    anrm, bignum, bnrm, c1, c2, s1, s2, smax, smaxpr, smin, sminpr, smlnum, wsize : RT = ---;
    i, ismax, ismin, j, lwkmin, lwkopt, mn, nb, nb1, nb2, nb3, nb4 : int = ---;
    lquery : bool = ---;
    early_zero : bool = false;

    mn = min(m, n);
    ismin = mn + 1;
    ismax = 2*mn + 1;

    info.* = 0;
    lquery = lwork == -1;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if lda < max(1, m) {
        info.* = -5;
    } else if ldb < max(1, m, n) {
        info.* = -7;
    }

    if info.* == 0 {
        if mn == 0 || nrhs == 0 {
            lwkmin = 1;
            lwkopt = 1;
        } else {
            #if RT == float32 {
                nb1 = ilaenv(1, "SGEQRF", " "s, m, n, -1, -1);
                nb2 = ilaenv(1, "SGERQF", " "s, m, n, -1, -1);
                nb3 = ilaenv(1, "SORMQR", " "s, m, n, nrhs, -1);
                nb4 = ilaenv(1, "SORMRQ", " "s, m, n, nrhs, -1);
            } else {
                nb1 = ilaenv(1, "DGEQRF", " "s, m, n, -1, -1);
                nb2 = ilaenv(1, "DGERQF", " "s, m, n, -1, -1);
                nb3 = ilaenv(1, "DORMQR", " "s, m, n, nrhs, -1);
                nb4 = ilaenv(1, "DORMRQ", " "s, m, n, nrhs, -1);
            }
            nb = max(nb1, nb2, nb3, nb4);
            lwkmin = mn + max(2*mn, n + 1, mn + nrhs);
            lwkopt = max(lwkmin, mn + 2*n + nb*(n + 1), 2*mn + nb*nrhs);
        }
        fset(work, 1, _gelsy_work_query_real(RT, lwkopt));

        if lwork < lwkmin && !lquery {
            info.* = -12;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if mn == 0 || nrhs == 0 {
        rank.* = 0;
        return;
    }

    #if RT == float32 {
        smlnum = slamch(#char "s") / slamch(#char "p");
    } else {
        smlnum = dlamch(#char "s") / dlamch(#char "p");
    }
    bignum = ONE / smlnum;

    #if RT == float32 {
        anrm = slange(#char "M", m, n, a.*, lda, work);
    } else {
        anrm = dlange(#char "M", m, n, a.*, lda, work);
    }
    iascl := 0;
    if anrm > ZERO && anrm < smlnum {
        #if RT == float32 {
            slascl(#char "G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        } else {
            dlascl(#char "G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        }
        iascl = 1;
    } else if anrm > bignum {
        #if RT == float32 {
            slascl(#char "G", 0, 0, anrm, bignum, m, n, a, lda, info);
        } else {
            dlascl(#char "G", 0, 0, anrm, bignum, m, n, a, lda, info);
        }
        iascl = 2;
    } else if anrm == ZERO {
        #if RT == float32 {
            slaset(#char "F", max(m, n), nrhs, ZERO, ZERO, b, ldb);
        } else {
            dlaset(#char "F", max(m, n), nrhs, ZERO, ZERO, b, ldb);
        }
        rank.* = 0;
        early_zero = true;
    }

    ibscl := 0;
    if !early_zero {
        #if RT == float32 {
            bnrm = slange(#char "M", m, nrhs, b.*, ldb, work);
        } else {
            bnrm = dlange(#char "M", m, nrhs, b.*, ldb, work);
        }
        if bnrm > ZERO && bnrm < smlnum {
            #if RT == float32 {
                slascl(#char "G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
            } else {
                dlascl(#char "G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
            }
            ibscl = 1;
        } else if bnrm > bignum {
            #if RT == float32 {
                slascl(#char "G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
            } else {
                dlascl(#char "G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
            }
            ibscl = 2;
        }

        tau_q := _gelsy_subvector(work.*, 1, mn);
        work_tail := _gelsy_subvector(work.*, mn + 1, max(0, lwork - mn));

        #if RT == float32 {
            sgeqp3(m, n, a, lda, jpvt, *tau_q, *work_tail, lwork - mn, info);
        } else {
            dgeqp3(m, n, a, lda, jpvt, *tau_q, *work_tail, lwork - mn, info);
        }
        wsize = ncast(RT, mn) + fget(work, mn + 1);

        fset(work, ismin, ONE);
        fset(work, ismax, ONE);
        smax = abs(fget(a.*, 1, 1));
        smin = smax;
        if abs(fget(a.*, 1, 1)) == ZERO {
            rank.* = 0;
            #if RT == float32 {
                slaset(#char "F", max(m, n), nrhs, ZERO, ZERO, b, ldb);
            } else {
                dlaset(#char "F", max(m, n), nrhs, ZERO, ZERO, b, ldb);
            }
            early_zero = true;
        } else {
            rank.* = 1;
        }

        if !early_zero {
            while rank.* < mn {
                i = rank.* + 1;
                #if RT == float32 {
                    slaic1(IMIN, rank.*, _gelsy_subvector(work.*, ismin, rank.*), smin, _gelsy_matrix_column(a.*, lda, i, rank.*), fget(a.*, i, i), *sminpr, *s1, *c1);
                    slaic1(IMAX, rank.*, _gelsy_subvector(work.*, ismax, rank.*), smax, _gelsy_matrix_column(a.*, lda, i, rank.*), fget(a.*, i, i), *smaxpr, *s2, *c2);
                } else {
                    dlaic1(IMIN, rank.*, _gelsy_subvector(work.*, ismin, rank.*), smin, _gelsy_matrix_column(a.*, lda, i, rank.*), fget(a.*, i, i), *sminpr, *s1, *c1);
                    dlaic1(IMAX, rank.*, _gelsy_subvector(work.*, ismax, rank.*), smax, _gelsy_matrix_column(a.*, lda, i, rank.*), fget(a.*, i, i), *smaxpr, *s2, *c2);
                }

                if smaxpr*rcond <= sminpr {
                    for i_:1..rank.* {
                        fset(work, ismin + i_ - 1, s1*fget(work, ismin + i_ - 1));
                        fset(work, ismax + i_ - 1, s2*fget(work, ismax + i_ - 1));
                    }
                    fset(work, ismin + rank.*, c1);
                    fset(work, ismax + rank.*, c2);
                    smin = sminpr;
                    smax = smaxpr;
                    rank.* += 1;
                } else {
                    break;
                }
            }

            tau_rz := _gelsy_subvector(work.*, mn + 1, mn);
            work_apply := _gelsy_subvector(work.*, 2*mn + 1, max(0, lwork - 2*mn));

            if rank.* < n {
                #if RT == float32 {
                    stzrzf(rank.*, n, a, lda, *tau_rz, *work_apply, lwork - 2*mn, info);
                } else {
                    dtzrzf(rank.*, n, a, lda, *tau_rz, *work_apply, lwork - 2*mn, info);
                }
            }

            #if RT == float32 {
                sormqr(#char "L", #char "T", m, nrhs, mn, a, lda, *tau_q, b, ldb, *work_apply, lwork - 2*mn, info);
            } else {
                dormqr(#char "L", #char "T", m, nrhs, mn, a, lda, *tau_q, b, ldb, *work_apply, lwork - 2*mn, info);
            }
            wsize = max(wsize, ncast(RT, 2*mn) + fget(work, 2*mn + 1));

            #if RT == float32 {
                strsm(#char "L", #char "U", #char "N", #char "N", rank.*, nrhs, ONE, a.*, lda, b, ldb);
            } else {
                dtrsm(#char "L", #char "U", #char "N", #char "N", rank.*, nrhs, ONE, a.*, lda, b, ldb);
            }

            for j:1..nrhs {
                for i_:rank.* + 1..n {
                    fset(b, i_, j, ZERO);
                }
            }

            if rank.* < n {
                #if RT == float32 {
                    sormrz(#char "L", #char "T", n, nrhs, rank.*, n - rank.*, a, lda, *tau_rz, b, ldb, *work_apply, lwork - 2*mn, info);
                } else {
                    dormrz(#char "L", #char "T", n, nrhs, rank.*, n - rank.*, a, lda, *tau_rz, b, ldb, *work_apply, lwork - 2*mn, info);
                }
            }

            for j:1..nrhs {
                for i_:1..n {
                    fset(work, fget(jpvt.*, i_), fget(b.*, i_, j));
                }
                for i_:1..n {
                    fset(b, i_, j, fget(work.*, i_));
                }
            }
        }
    }

    if iascl == 1 {
        #if RT == float32 {
            slascl(#char "G", 0, 0, anrm, smlnum, n, nrhs, b, ldb, info);
            slascl(#char "U", 0, 0, smlnum, anrm, rank.*, rank.*, a, lda, info);
        } else {
            dlascl(#char "G", 0, 0, anrm, smlnum, n, nrhs, b, ldb, info);
            dlascl(#char "U", 0, 0, smlnum, anrm, rank.*, rank.*, a, lda, info);
        }
    } else if iascl == 2 {
        #if RT == float32 {
            slascl(#char "G", 0, 0, anrm, bignum, n, nrhs, b, ldb, info);
            slascl(#char "U", 0, 0, bignum, anrm, rank.*, rank.*, a, lda, info);
        } else {
            dlascl(#char "G", 0, 0, anrm, bignum, n, nrhs, b, ldb, info);
            dlascl(#char "U", 0, 0, bignum, anrm, rank.*, rank.*, a, lda, info);
        }
    }
    if ibscl == 1 {
        #if RT == float32 {
            slascl(#char "G", 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info);
        } else {
            dlascl(#char "G", 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info);
        }
    } else if ibscl == 2 {
        #if RT == float32 {
            slascl(#char "G", 0, 0, bignum, bnrm, n, nrhs, b, ldb, info);
        } else {
            dlascl(#char "G", 0, 0, bignum, bnrm, n, nrhs, b, ldb, info);
        }
    }

    fset(work, 1, _gelsy_work_query_real(RT, lwkopt));
}

_gelsy_complex :: (
    $$m: int,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    jpvt: *$JPVT/ArrayType1D,
    rcond: $RT,
    rank: *int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           JPVT.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CGELSY " } else { "ZGELSY " };
    CZERO :: #run cpx(zero(RT), zero(RT));
    CONE  :: #run cpx(one(RT), zero(RT));

    anrm, bignum, bnrm, smax, smaxpr, smin, sminpr, smlnum, wsize : RT = ---;
    c1, c2, s1, s2 : MT = ---;
    i, ismax, ismin, lwkopt, mn, nb, nb1, nb2, nb3, nb4 : int = ---;
    lquery : bool = ---;
    early_zero : bool = false;

    mn = min(m, n);
    ismin = mn + 1;
    ismax = 2*mn + 1;

    info.* = 0;
    #if MT == C32 {
        nb1 = ilaenv(1, "CGEQRF", " "s, m, n, -1, -1);
        nb2 = ilaenv(1, "CGERQF", " "s, m, n, -1, -1);
        nb3 = ilaenv(1, "CUNMQR", " "s, m, n, nrhs, -1);
        nb4 = ilaenv(1, "CUNMRQ", " "s, m, n, nrhs, -1);
    } else {
        nb1 = ilaenv(1, "ZGEQRF", " "s, m, n, -1, -1);
        nb2 = ilaenv(1, "ZGERQF", " "s, m, n, -1, -1);
        nb3 = ilaenv(1, "ZUNMQR", " "s, m, n, nrhs, -1);
        nb4 = ilaenv(1, "ZUNMRQ", " "s, m, n, nrhs, -1);
    }
    nb = max(nb1, nb2, nb3, nb4);
    lwkopt = max(1, mn + 2*n + nb*(n + 1), 2*mn + nb*nrhs);
    fset(work, 1, _gelsy_work_query_complex(RT, MT, lwkopt));
    lquery = lwork == -1;

    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if lda < max(1, m) {
        info.* = -5;
    } else if ldb < max(1, m, n) {
        info.* = -7;
    } else if lwork < mn + max(2*mn, n + 1, mn + nrhs) && !lquery {
        info.* = -12;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if min(m, n, nrhs) == 0 {
        rank.* = 0;
        return;
    }

    #if RT == float32 {
        smlnum = slamch(#char "s") / slamch(#char "p");
    } else {
        smlnum = dlamch(#char "s") / dlamch(#char "p");
    }
    bignum = one(RT) / smlnum;

    #if MT == C32 {
        anrm = clange(#char "M", m, n, a.*, lda, rwork);
    } else {
        anrm = zlange(#char "M", m, n, a.*, lda, rwork);
    }
    iascl := 0;
    if anrm > zero(RT) && anrm < smlnum {
        #if MT == C32 {
            clascl(#char "G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        } else {
            zlascl(#char "G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        }
        iascl = 1;
    } else if anrm > bignum {
        #if MT == C32 {
            clascl(#char "G", 0, 0, anrm, bignum, m, n, a, lda, info);
        } else {
            zlascl(#char "G", 0, 0, anrm, bignum, m, n, a, lda, info);
        }
        iascl = 2;
    } else if anrm == zero(RT) {
        #if MT == C32 {
            claset(#char "F", max(m, n), nrhs, CZERO, CZERO, b, ldb);
        } else {
            zlaset(#char "F", max(m, n), nrhs, CZERO, CZERO, b, ldb);
        }
        rank.* = 0;
        early_zero = true;
    }

    ibscl := 0;
    if !early_zero {
        #if MT == C32 {
            bnrm = clange(#char "M", m, nrhs, b.*, ldb, rwork);
        } else {
            bnrm = zlange(#char "M", m, nrhs, b.*, ldb, rwork);
        }
        if bnrm > zero(RT) && bnrm < smlnum {
            #if MT == C32 {
                clascl(#char "G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
            } else {
                zlascl(#char "G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
            }
            ibscl = 1;
        } else if bnrm > bignum {
            #if MT == C32 {
                clascl(#char "G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
            } else {
                zlascl(#char "G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
            }
            ibscl = 2;
        }

        tau_q := _gelsy_subvector(work.*, 1, mn);
        work_tail := _gelsy_subvector(work.*, mn + 1, max(0, lwork - mn));

        #if MT == C32 {
            cgeqp3(m, n, a, lda, jpvt, *tau_q, *work_tail, lwork - mn, rwork, info);
        } else {
            zgeqp3(m, n, a, lda, jpvt, *tau_q, *work_tail, lwork - mn, rwork, info);
        }
        wsize = ncast(RT, mn) + real(fget(work, mn + 1));

        fset(work, ismin, CONE);
        fset(work, ismax, CONE);
        smax = abs(fget(a.*, 1, 1));
        smin = smax;
        if abs(fget(a.*, 1, 1)) == zero(RT) {
            rank.* = 0;
            #if MT == C32 {
                claset(#char "F", max(m, n), nrhs, CZERO, CZERO, b, ldb);
            } else {
                zlaset(#char "F", max(m, n), nrhs, CZERO, CZERO, b, ldb);
            }
            early_zero = true;
        } else {
            rank.* = 1;
        }

        if !early_zero {
            while rank.* < mn {
                i = rank.* + 1;
                #if MT == C32 {
                    claic1(2, rank.*, _gelsy_subvector(work.*, ismin, rank.*), smin, _gelsy_matrix_column(a.*, lda, i, rank.*), fget(a.*, i, i), *sminpr, *s1, *c1);
                    claic1(1, rank.*, _gelsy_subvector(work.*, ismax, rank.*), smax, _gelsy_matrix_column(a.*, lda, i, rank.*), fget(a.*, i, i), *smaxpr, *s2, *c2);
                } else {
                    zlaic1(2, rank.*, _gelsy_subvector(work.*, ismin, rank.*), smin, _gelsy_matrix_column(a.*, lda, i, rank.*), fget(a.*, i, i), *sminpr, *s1, *c1);
                    zlaic1(1, rank.*, _gelsy_subvector(work.*, ismax, rank.*), smax, _gelsy_matrix_column(a.*, lda, i, rank.*), fget(a.*, i, i), *smaxpr, *s2, *c2);
                }

                if smaxpr*rcond <= sminpr {
                    for i_:1..rank.* {
                        fset(work, ismin + i_ - 1, s1*fget(work, ismin + i_ - 1));
                        fset(work, ismax + i_ - 1, s2*fget(work, ismax + i_ - 1));
                    }
                    fset(work, ismin + rank.*, c1);
                    fset(work, ismax + rank.*, c2);
                    smin = sminpr;
                    smax = smaxpr;
                    rank.* += 1;
                } else {
                    break;
                }
            }

            tau_rz := _gelsy_subvector(work.*, mn + 1, mn);
            work_apply := _gelsy_subvector(work.*, 2*mn + 1, max(0, lwork - 2*mn));

            if rank.* < n {
                #if MT == C32 {
                    ctzrzf(rank.*, n, a, lda, *tau_rz, *work_apply, lwork - 2*mn, info);
                } else {
                    ztzrzf(rank.*, n, a, lda, *tau_rz, *work_apply, lwork - 2*mn, info);
                }
            }

            #if MT == C32 {
                cunmqr(#char "L", #char "C", m, nrhs, mn, a, lda, *tau_q, b, ldb, *work_apply, lwork - 2*mn, info);
            } else {
                zunmqr(#char "L", #char "C", m, nrhs, mn, a, lda, *tau_q, b, ldb, *work_apply, lwork - 2*mn, info);
            }
            wsize = max(wsize, ncast(RT, 2*mn) + real(fget(work, 2*mn + 1)));

            #if MT == C32 {
                ctrsm(#char "L", #char "U", #char "N", #char "N", rank.*, nrhs, CONE, a.*, lda, b, ldb);
            } else {
                ztrsm(#char "L", #char "U", #char "N", #char "N", rank.*, nrhs, CONE, a.*, lda, b, ldb);
            }

            for j:1..nrhs {
                for i_:rank.* + 1..n {
                    fset(b, i_, j, CZERO);
                }
            }

            if rank.* < n {
                #if MT == C32 {
                    cunmrz(#char "L", #char "C", n, nrhs, rank.*, n - rank.*, a, lda, *tau_rz, b, ldb, *work_apply, lwork - 2*mn, info);
                } else {
                    zunmrz(#char "L", #char "C", n, nrhs, rank.*, n - rank.*, a, lda, *tau_rz, b, ldb, *work_apply, lwork - 2*mn, info);
                }
            }

            for j:1..nrhs {
                for i_:1..n {
                    fset(work, fget(jpvt.*, i_), fget(b.*, i_, j));
                }
                for i_:1..n {
                    fset(b, i_, j, fget(work.*, i_));
                }
            }
        }
    }

    if iascl == 1 {
        #if MT == C32 {
            clascl(#char "G", 0, 0, anrm, smlnum, n, nrhs, b, ldb, info);
            clascl(#char "U", 0, 0, smlnum, anrm, rank.*, rank.*, a, lda, info);
        } else {
            zlascl(#char "G", 0, 0, anrm, smlnum, n, nrhs, b, ldb, info);
            zlascl(#char "U", 0, 0, smlnum, anrm, rank.*, rank.*, a, lda, info);
        }
    } else if iascl == 2 {
        #if MT == C32 {
            clascl(#char "G", 0, 0, anrm, bignum, n, nrhs, b, ldb, info);
            clascl(#char "U", 0, 0, bignum, anrm, rank.*, rank.*, a, lda, info);
        } else {
            zlascl(#char "G", 0, 0, anrm, bignum, n, nrhs, b, ldb, info);
            zlascl(#char "U", 0, 0, bignum, anrm, rank.*, rank.*, a, lda, info);
        }
    }
    if ibscl == 1 {
        #if MT == C32 {
            clascl(#char "G", 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info);
        } else {
            zlascl(#char "G", 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info);
        }
    } else if ibscl == 2 {
        #if MT == C32 {
            clascl(#char "G", 0, 0, bignum, bnrm, n, nrhs, b, ldb, info);
        } else {
            zlascl(#char "G", 0, 0, bignum, bnrm, n, nrhs, b, ldb, info);
        }
    }

    fset(work, 1, _gelsy_work_query_complex(RT, MT, lwkopt));
}
