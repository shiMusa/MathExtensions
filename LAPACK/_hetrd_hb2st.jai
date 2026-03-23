// ? LAPACK *hetrd_hb2st()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrd_hb2st/zhetrd_hb2st
// *HETRD_HB2ST reduces a Hermitian band matrix to tridiagonal form.

chetrd_hb2st :: #bake_arguments _hetrd_hb2st(RT=float32, MT=C32);
zhetrd_hb2st :: #bake_arguments _hetrd_hb2st(RT=float64, MT=C64);

#scope_file
_hetrd_hb2st_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrd_hb2st_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetrd_hb2st_raw_slice(x.data + (offset - 1), count);
}

_hetrd_hb2st_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrd_hb2st_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrd_hb2st_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int,
    a: $A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 { clacpy(uplo, m, n, a, lda, b, ldb); } else { zlacpy(uplo, m, n, a, lda, b, ldb); }
}

_hetrd_hb2st_laset :: (
    $$uplo: u8, $$m: int, $$n: int,
    alpha: $MT, beta: $MT,
    a: *$A/ArrayType2D, $$lda: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 { claset(uplo, m, n, alpha, beta, a, lda); } else { zlaset(uplo, m, n, alpha, beta, a, lda); }
}

_hetrd_hb2st_kernel :: (
    $$uplo: u8, $$wantz: bool, $$ttype: int, $$st: int, $$ed: int, $$sweep: int,
    $$n: int, $$nb: int, $$ib: int,
    a: *$A/ArrayType1D, $$lda: int,
    v: *$V/ArrayType1D, tau: *$TAU/ArrayType1D, $$ldvt: int,
    work: *$WORK/ArrayType1D,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           V.data_type == MT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        chb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
    } else {
        zhb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
    }
}

_hetrd_hb2st :: (
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
    rzero :: #run zero(RT);
    zero_t :: #run zero(MT);
    one_t  :: #run _hetrd_hb2st_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHETRD_HB2ST" } else { "ZHETRD_HB2ST" };

    afters1, lquery, upper, wantq : bool = ---;
    i, m, k, ib, myid, sweepid, shift, stt, st, ed, stind, edind, blk_last_ind, colpt, thed, steppercol, grsiz, thgrsiz, thgrnb, thgrid, nbtiles, ttype, tid, nthreads, abdpos, abofdpos, dpos, ofdpos, awpos, inda, indw, apos, sizea, lda, indv, indtau, sizev, sizetau, ldv, lhmin, lwmin : int = ---;
    tmp : MT = ---;
    abstmp : RT = ---;

    info.* = 0;
    afters1 = lsame(stage1, #char "Y");
    wantq = lsame(vect, #char "V");
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1 || lhous == -1;

    ib = ilaenv2stage(2, ifx MT == C32 then "CHETRD_HB2ST" else "ZHETRD_HB2ST", vect, n, kd, -1, -1);
    if n == 0 || kd <= 1 {
        lhmin = 1;
        lwmin = 1;
    } else {
        lhmin = ilaenv2stage(3, ifx MT == C32 then "CHETRD_HB2ST" else "ZHETRD_HB2ST", vect, n, kd, ib, -1);
        lwmin = ilaenv2stage(4, ifx MT == C32 then "CHETRD_HB2ST" else "ZHETRD_HB2ST", vect, n, kd, ib, -1);
    }

    if !afters1 && !lsame(stage1, #char "N") {
        info.* = -1;
    } else if !lsame(vect, #char "N") {
        info.* = -2;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if kd < 0 {
        info.* = -5;
    } else if ldab < kd + 1 {
        info.* = -7;
    } else if lhous < lhmin && !lquery {
        info.* = -11;
    } else if lwork < lwmin && !lquery {
        info.* = -13;
    }

    if info.* == 0 {
        fset(hous, 1, _hetrd_hb2st_complex_from_real(_hetrd_hb2st_roundup_lwork(ncast(RT, lhmin), RT), RT, MT));
        fset(work, 1, _hetrd_hb2st_complex_from_real(_hetrd_hb2st_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        fset(hous, 1, one_t);
        fset(work, 1, one_t);
        return;
    }

    sizetau = 2*n;
    sizev = 2*n;
    indtau = 1;
    indv = indtau + sizetau;
    lda = 2*kd + 1;
    sizea = lda*n;
    inda = 1;
    indw = inda + sizea;
    nthreads = 1;
    tid = 0;
    ldv = kd + ib;

    if upper {
        apos = inda + kd;
        awpos = inda;
        dpos = apos + kd;
        ofdpos = dpos - 1;
        abdpos = kd + 1;
        abofdpos = kd;
    } else {
        apos = inda;
        awpos = inda + kd + 1;
        dpos = apos;
        ofdpos = dpos + 1;
        abdpos = 1;
        abofdpos = 2;
    }

    if kd == 0 {
        for i:1..n {
            fset(d, i, real(fget(ab.*, abdpos, i)));
        }
        for i:1..n - 1 {
            fset(e, i, rzero);
        }
        fset(hous, 1, one_t);
        fset(work, 1, one_t);
        return;
    }

    if kd == 1 {
        for i:1..n {
            fset(d, i, real(fget(ab.*, abdpos, i)));
        }

        if upper {
            for i:1..n - 1 {
                tmp = fget(ab.*, abofdpos, i + 1);
                abstmp = abs(tmp);
                fset(ab, abofdpos, i + 1, _hetrd_hb2st_complex_from_real(abstmp, RT, MT));
                fset(e, i, abstmp);
                if abstmp != rzero {
                    tmp /= abstmp;
                } else {
                    tmp = one_t;
                }
                if i < n - 1 {
                    fset(ab, abofdpos, i + 2, fget(ab.*, abofdpos, i + 2) * tmp);
                }
            }
        } else {
            for i:1..n - 1 {
                tmp = fget(ab.*, abofdpos, i);
                abstmp = abs(tmp);
                fset(ab, abofdpos, i, _hetrd_hb2st_complex_from_real(abstmp, RT, MT));
                fset(e, i, abstmp);
                if abstmp != rzero {
                    tmp /= abstmp;
                } else {
                    tmp = one_t;
                }
                if i < n - 1 {
                    fset(ab, abofdpos, i + 1, fget(ab.*, abofdpos, i + 1) * tmp);
                }
            }
        }

        fset(hous, 1, one_t);
        fset(work, 1, one_t);
        return;
    }

    thgrsiz = n;
    grsiz = 1;
    shift = 3;
    nbtiles = (n + kd - 1) / kd;
    steppercol = (shift + grsiz - 1) / grsiz;
    thgrnb = (n - 1 + thgrsiz - 1) / thgrsiz;
    _ = nbtiles;
    _ = nthreads;
    _ = tid;

    work_a_raw := _hetrd_hb2st_raw_slice(work.data + (apos - 1), sizea - (apos - inda));
    work_a := array2d_view(work_a_raw, lda, row_major=false);
    _hetrd_hb2st_lacpy(#char "A", kd + 1, n, ab.*, ldab, *work_a, lda, MT);

    work_aw_raw := _hetrd_hb2st_raw_slice(work.data + (awpos - 1), max(0, sizea - (awpos - inda)));
    work_aw := array2d_view(work_aw_raw, lda, row_major=false);
    _hetrd_hb2st_laset(#char "A", kd, n, zero_t, zero_t, *work_aw, lda, MT);

    a_view := _hetrd_hb2st_vector_slice(work, inda, sizea);
    v_view := _hetrd_hb2st_vector_slice(hous, indv, sizev);
    tau_view := _hetrd_hb2st_vector_slice(hous, indtau, sizetau);
    scratch_view := _hetrd_hb2st_vector_slice(work, indw, kd);

    for thgrid:1..thgrnb {
        stt = (thgrid - 1)*thgrsiz + 1;
        thed = min(stt + thgrsiz - 1, n - 1);
        for i:stt..n - 1 {
            ed = min(i, thed);
            if stt > ed then break;
            for m:1..steppercol {
                st = stt;
                for sweepid:st..ed {
                    for k:1..grsiz {
                        myid := (i - sweepid)*(steppercol*grsiz) + (m - 1)*grsiz + k;
                        if myid == 1 {
                            ttype = 1;
                        } else {
                            ttype = mod(myid, 2) + 2;
                        }

                        if ttype == 2 {
                            colpt = (myid / 2)*kd + sweepid;
                            stind = colpt - kd + 1;
                            edind = min(colpt, n);
                            blk_last_ind = colpt;
                        } else {
                            colpt = ((myid + 1) / 2)*kd + sweepid;
                            stind = colpt - kd + 1;
                            edind = min(colpt, n);
                            if stind >= edind - 1 && edind == n {
                                blk_last_ind = n;
                            } else {
                                blk_last_ind = 0;
                            }
                        }

                        _hetrd_hb2st_kernel(uplo, wantq, ttype, stind, edind, sweepid, n, kd, ib, *a_view, lda, *v_view, *tau_view, ldv, *scratch_view, RT, MT);
                        if blk_last_ind >= n - 1 {
                            stt += 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    for i:1..n {
        fset(d, i, real(fget(work.*, dpos + (i - 1)*lda)));
    }

    if upper {
        for i:1..n - 1 {
            fset(e, i, real(fget(work.*, ofdpos + i*lda)));
        }
    } else {
        for i:1..n - 1 {
            fset(e, i, real(fget(work.*, ofdpos + (i - 1)*lda)));
        }
    }

    fset(work, 1, _hetrd_hb2st_complex_from_real(_hetrd_hb2st_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
}
