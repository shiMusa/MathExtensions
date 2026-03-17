// ? LAPACK *gedmd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgedmd/cgedmd/dgedmd/zgedmd
// *GEDMD computes the Dynamic Mode Decomposition (DMD) for a pair of
// data snapshot matrices.

sgedmd :: #bake_arguments _gedmd_real(RT=float32);
dgedmd :: #bake_arguments _gedmd_real(RT=float64);
cgedmd :: #bake_arguments _gedmd_complex(RT=float32, MT=C32);
zgedmd :: #bake_arguments _gedmd_complex(RT=float64, MT=C64);

#scope_file
_gedmd_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gedmd_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gedmd_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _gedmd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gedmd_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _gedmd_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_gedmd_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gedmd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gedmd_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_gedmd_is_nonzero :: inline (z: $MT, $RT: Type, $MT: Type) -> bool #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return z != zero(MT);
    } else {
        return real(z) != zero(RT) || imag(z) != zero(RT);
    }
}

_gedmd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gedmd_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_gedmd_lamch_overflow :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "o");
    } else {
        return dlamch(#char "o");
    }
}

_gedmd_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gedmd_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_gedmd_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return isamax(n, x, incx) + 1;
    } else if MT == float64 {
        return idamax(n, x, incx) + 1;
    } else if MT == C32 {
        return icamax(n, x, incx) + 1;
    } else {
        return izamax(n, x, incx) + 1;
    }
}

_gedmd_lassq :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    scale: *$RT,
    sumsq: *RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        slassq(n, x, incx, scale, sumsq);
    } else if MT == float64 {
        dlassq(n, x, incx, scale, sumsq);
    } else if MT == C32 {
        classq(n, x, incx, scale, sumsq);
    } else {
        zlassq(n, x, incx, scale, sumsq);
    }
}

_gedmd_scal :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, incx);
    } else if MT == float64 {
        dscal(n, alpha, x, incx);
    } else if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_gedmd_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        saxpy(n, alpha, x, 1, y, 1);
    } else if MT == float64 {
        daxpy(n, alpha, x, 1, y, 1);
    } else if MT == C32 {
        caxpy(n, alpha, x, 1, y, 1);
    } else {
        zaxpy(n, alpha, x, 1, y, 1);
    }
}

_gedmd_nrm2 :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return snrm2(n, x, incx);
    } else if MT == float64 {
        return dnrm2(n, x, incx);
    } else if MT == C32 {
        return scnrm2(n, x, incx);
    } else {
        return dznrm2(n, x, incx);
    }
}

_gedmd_lange :: (
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

_gedmd_lacpy :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        slacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == float64 {
        dlacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_gedmd_lascl :: (
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

_gedmd_gemm :: (
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
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == float32 {
        sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == float64 {
        dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_gedmd_real :: (
    $$jobs: u8,
    $$jobz: u8,
    $$jobr: u8,
    $$jobf: u8,
    $$whtsvd: int,
    $$m: int,
    $$n: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    y: *$Y/ArrayType2D,
    $$ldy: int,
    $$nrnk: int,
    tol: $RT,
    k: *int,
    reig: *$REIG/ArrayType1D,
    imeig: *$IMEIG/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    res: *$RES/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    w: *$W/ArrayType2D,
    $$ldw: int,
    s: *$S/ArrayType2D,
    $$lds: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    iwork: *$IWORK/ArrayType1D,
    $$liwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT &&
           REIG.data_type == RT &&
           IMEIG.data_type == RT &&
           Z.data_type == RT &&
           RES.data_type == RT &&
           B.data_type == RT &&
           W.data_type == RT &&
           S.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { "SGEDMD" } else { "DGEDMD" };

    ab : [2][2] RT;
    rdummy : [2] RT;
    rdummy2 : [2] RT;

    badxy, lquery, sccolx, sccoly, wntex, wntref, wntres, wntvec : bool = ---;
    jobzl, t_or_n, jsvopt : u8 = ---;
    i, j, iminwr, info1, info2, lwrkev, lwrsdd, lwrsvd, lwrsvq, mlwork, mwrkev, mwrsdd, mwrsvd, mwrsvj, mwrsvq, numrnk, olwork : int = ---;
    ofl, rootsc, scale, small, ssum, xscl1, xscl2 : RT = ---;

    wntres = lsame(jobr, #char "R");
    sccolx = lsame(jobs, #char "S") || lsame(jobs, #char "C");
    sccoly = lsame(jobs, #char "Y");
    wntvec = lsame(jobz, #char "V");
    wntref = lsame(jobf, #char "R");
    wntex  = lsame(jobf, #char "E");
    info.* = 0;
    lquery = (lwork == -1) || (liwork == -1);

    if !(sccolx || sccoly || lsame(jobs, #char "N")) {
        info.* = -1;
    } else if !(wntvec || lsame(jobz, #char "N") || lsame(jobz, #char "F")) {
        info.* = -2;
    } else if !(wntres || lsame(jobr, #char "N")) || (wntres && !wntvec) {
        info.* = -3;
    } else if !(wntref || wntex || lsame(jobf, #char "N")) {
        info.* = -4;
    } else if !((whtsvd == 1) || (whtsvd == 2) || (whtsvd == 3) || (whtsvd == 4)) {
        info.* = -5;
    } else if m < 0 {
        info.* = -6;
    } else if (n < 0) || (n > m) {
        info.* = -7;
    } else if ldx < m {
        info.* = -9;
    } else if ldy < m {
        info.* = -11;
    } else if !((nrnk == -2) || (nrnk == -1) || ((nrnk >= 1) && (nrnk <= n))) {
        info.* = -12;
    } else if (tol < zero_r) || (tol >= one_r) {
        info.* = -13;
    } else if ldz < m {
        info.* = -18;
    } else if (wntref || wntex) && (ldb < m) {
        info.* = -21;
    } else if ldw < n {
        info.* = -23;
    } else if lds < n {
        info.* = -25;
    }

    if info.* == 0 {
        if n == 0 {
            if lquery {
                fset(iwork, 1, 1);
                fset(work, 1, ncast(RT, 2));
                fset(work, 2, ncast(RT, 2));
            } else {
                k.* = 0;
            }
            info.* = 1;
            return;
        }

        mlwork = max(2, n);
        olwork = max(2, n);
        iminwr = 1;

        switch whtsvd {
        case 1:
            mwrsvd = max(1, 3*min(m, n) + max(m, n), 5*min(m, n));
            mlwork = max(mlwork, n + mwrsvd);
            if lquery {
                #if RT == float32 {
                    sgesvd(#char "O", #char "S", m, n, x, ldx, work, b, ldb, w, ldw, rdummy, -1, *info1);
                } else {
                    dgesvd(#char "O", #char "S", m, n, x, ldx, work, b, ldb, w, ldw, rdummy, -1, *info1);
                }
                lwrsvd = max(mwrsvd, cast(int) rdummy[0]);
                olwork = max(olwork, n + lwrsvd);
            }
        case 2:
            mwrsdd = 3*min(m, n)*min(m, n) + max(max(m, n), 5*min(m, n)*min(m, n) + 4*min(m, n));
            mlwork = max(mlwork, n + mwrsdd);
            iminwr = 8*min(m, n);
            if lquery {
                #if RT == float32 {
                    sgesdd(#char "O", m, n, x, ldx, work, b, ldb, w, ldw, rdummy, -1, iwork, *info1);
                } else {
                    dgesdd(#char "O", m, n, x, ldx, work, b, ldb, w, ldw, rdummy, -1, iwork, *info1);
                }
                lwrsdd = max(mwrsdd, cast(int) rdummy[0]);
                olwork = max(olwork, n + lwrsdd);
            }
        case 3:
            numrnk = 0;
            #if RT == float32 {
                sgesvdq(#char "H", #char "P", #char "N", #char "R", #char "R", m, n, x, ldx, work, z, ldz, w, ldw,
                        *numrnk, iwork, -1, rdummy, -1, rdummy2, -1, *info1);
            } else {
                dgesvdq(#char "H", #char "P", #char "N", #char "R", #char "R", m, n, x, ldx, work, z, ldz, w, ldw,
                        *numrnk, iwork, -1, rdummy, -1, rdummy2, -1, *info1);
            }
            iminwr = fget(iwork.*, 1);
            mwrsvq = cast(int) rdummy[1];
            mlwork = max(mlwork, n + mwrsvq + cast(int) rdummy2[0]);
            if lquery {
                lwrsvq = max(mwrsvq, cast(int) rdummy[0]);
                olwork = max(olwork, n + lwrsvq + cast(int) rdummy2[0]);
            }
        case 4:
            jsvopt = #char "J";
            mwrsvj = max(7, 2*m + n, 4*n + n*n, 2*n + n*n + 6);
            mlwork = max(mlwork, n + mwrsvj);
            iminwr = max(3, m + 3*n);
            if lquery {
                olwork = max(olwork, n + mwrsvj);
            }
        }

        if wntvec || wntex || lsame(jobz, #char "F") {
            jobzl = #char "V";
        } else {
            jobzl = #char "N";
        }

        if lsame(jobzl, #char "V") {
            mwrkev = max(1, 4*n);
        } else {
            mwrkev = max(1, 3*n);
        }
        mlwork = max(mlwork, n + mwrkev);
        if lquery {
            #if RT == float32 {
                sgeev(#char "N", jobzl, n, s, lds, reig, imeig, w, ldw, w, ldw, rdummy, -1, *info1);
            } else {
                dgeev(#char "N", jobzl, n, s, lds, reig, imeig, w, ldw, w, ldw, rdummy, -1, *info1);
            }
            lwrkev = max(mwrkev, cast(int) rdummy[0]);
            olwork = max(olwork, n + lwrkev);
        }

        if (liwork < iminwr) && !lquery {
            info.* = -29;
        }
        if (lwork < mlwork) && !lquery {
            info.* = -27;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        fset(iwork, 1, iminwr);
        fset(work, 1, ncast(RT, mlwork));
        fset(work, 2, ncast(RT, olwork));
        return;
    }

    ofl = _gedmd_lamch_overflow(RT);
    small = _gedmd_lamch_safe_minimum(RT);
    badxy = false;

    if sccolx {
        k.* = 0;
        for i:1..n {
            ssum = one_r;
            scale = zero_r;
            x_col := _gedmd_matrix_column(x.*, ldx, 1, i, m);
            _gedmd_lassq(m, x_col, 1, *scale, *ssum, RT, RT);
            if _gedmd_isnan(scale, RT) || _gedmd_isnan(ssum, RT) {
                k.* = 0;
                info.* = -8;
                xerbla(routine_name, -info.*);
            }
            if (scale != zero_r) && (ssum != zero_r) {
                rootsc = sqrt(ssum);
                if scale >= (ofl / rootsc) {
                    x_col_mat := _gedmd_matrix_storage(x.*, ldx, 1, i, m, 1);
                    _gedmd_lascl(#char "G", 0, 0, scale, one_r/rootsc, m, 1, *x_col_mat, m, *info2, RT, RT);
                    fset(work, i, -scale * (rootsc / ncast(RT, m)));
                } else {
                    fset(work, i, scale * rootsc);
                    x_col_mat := _gedmd_matrix_storage(x.*, ldx, 1, i, m, 1);
                    _gedmd_lascl(#char "G", 0, 0, fget(work.*, i), one_r, m, 1, *x_col_mat, m, *info2, RT, RT);
                }
            } else {
                fset(work, i, zero_r);
                k.* += 1;
            }
        }
        if k.* == n {
            k.* = 0;
            info.* = -8;
            xerbla(routine_name, -info.*);
            return;
        }
        for i:1..n {
            yi := _gedmd_matrix_column(y.*, ldy, 1, i, m);
            if fget(work.*, i) > zero_r {
                _gedmd_scal(m, one_r/fget(work.*, i), *yi, 1, RT, RT);
            } else if fget(work.*, i) < zero_r {
                y_col_mat := _gedmd_matrix_storage(y.*, ldy, 1, i, m, 1);
                _gedmd_lascl(#char "G", 0, 0, -fget(work.*, i), one_r/ncast(RT, m), m, 1, *y_col_mat, m, *info2, RT, RT);
            } else if _gedmd_is_nonzero(fget(yi, _gedmd_iamax_1based(m, yi, 1, RT, RT)), RT, RT) {
                badxy = true;
                if lsame(jobs, #char "C") {
                    _gedmd_scal(m, zero_r, *yi, 1, RT, RT);
                }
            }
        }
    }

    if sccoly {
        for i:1..n {
            ssum = one_r;
            scale = zero_r;
            y_col := _gedmd_matrix_column(y.*, ldy, 1, i, m);
            _gedmd_lassq(m, y_col, 1, *scale, *ssum, RT, RT);
            if _gedmd_isnan(scale, RT) || _gedmd_isnan(ssum, RT) {
                k.* = 0;
                info.* = -10;
                xerbla(routine_name, -info.*);
            }
            if (scale != zero_r) && (ssum != zero_r) {
                rootsc = sqrt(ssum);
                if scale >= (ofl / rootsc) {
                    y_col_mat := _gedmd_matrix_storage(y.*, ldy, 1, i, m, 1);
                    _gedmd_lascl(#char "G", 0, 0, scale, one_r/rootsc, m, 1, *y_col_mat, m, *info2, RT, RT);
                    fset(work, i, -scale * (rootsc / ncast(RT, m)));
                } else {
                    fset(work, i, scale * rootsc);
                    y_col_mat := _gedmd_matrix_storage(y.*, ldy, 1, i, m, 1);
                    _gedmd_lascl(#char "G", 0, 0, fget(work.*, i), one_r, m, 1, *y_col_mat, m, *info2, RT, RT);
                }
            } else {
                fset(work, i, zero_r);
            }
        }
        for i:1..n {
            xi := _gedmd_matrix_column(x.*, ldx, 1, i, m);
            if fget(work.*, i) > zero_r {
                _gedmd_scal(m, one_r/fget(work.*, i), *xi, 1, RT, RT);
            } else if fget(work.*, i) < zero_r {
                x_col_mat := _gedmd_matrix_storage(x.*, ldx, 1, i, m, 1);
                _gedmd_lascl(#char "G", 0, 0, -fget(work.*, i), one_r/ncast(RT, m), m, 1, *x_col_mat, m, *info2, RT, RT);
            } else if _gedmd_is_nonzero(fget(xi, _gedmd_iamax_1based(m, xi, 1, RT, RT)), RT, RT) {
                badxy = true;
            }
        }
    }

    numrnk = n;
    switch whtsvd {
    case 1:
        svd_work := _gedmd_subvector(work.*, n + 1, lwork - n);
        #if RT == float32 {
            sgesvd(#char "O", #char "S", m, n, x, ldx, work, b, ldb, w, ldw, *svd_work, lwork - n, *info1);
        } else {
            dgesvd(#char "O", #char "S", m, n, x, ldx, work, b, ldb, w, ldw, *svd_work, lwork - n, *info1);
        }
        t_or_n = #char "T";
    case 2:
        svd_work := _gedmd_subvector(work.*, n + 1, lwork - n);
        #if RT == float32 {
            sgesdd(#char "O", m, n, x, ldx, work, b, ldb, w, ldw, *svd_work, lwork - n, iwork, *info1);
        } else {
            dgesdd(#char "O", m, n, x, ldx, work, b, ldb, w, ldw, *svd_work, lwork - n, iwork, *info1);
        }
        t_or_n = #char "T";
    case 3:
        qwork := _gedmd_subvector(work.*, n + max(2, m) + 1, max(0, lwork - n - max(2, m)));
        qaux  := _gedmd_subvector(work.*, n + 1, max(2, m));
        #if RT == float32 {
            sgesvdq(#char "H", #char "P", #char "N", #char "R", #char "R", m, n, x, ldx, work, z, ldz, w, ldw,
                    *numrnk, iwork, liwork, *qwork, lwork - n - max(2, m), *qaux, max(2, m), *info1);
        } else {
            dgesvdq(#char "H", #char "P", #char "N", #char "R", #char "R", m, n, x, ldx, work, z, ldz, w, ldw,
                    *numrnk, iwork, liwork, *qwork, lwork - n - max(2, m), *qaux, max(2, m), *info1);
        }
        _gedmd_lacpy(#char "A", m, numrnk, z.*, ldz, x, ldx, RT);
        t_or_n = #char "T";
    case 4:
        jsvopt = #char "J";
        svj_work := _gedmd_subvector(work.*, n + 1, lwork - n);
        #if RT == float32 {
            sgejsv(#char "F", #char "U", jsvopt, #char "N", #char "N", #char "P", m, n, x, ldx, work, z, ldz, w, ldw,
                   *svj_work, lwork - n, iwork, *info1);
        } else {
            dgejsv(#char "F", #char "U", jsvopt, #char "N", #char "N", #char "P", m, n, x, ldx, work, z, ldz, w, ldw,
                   *svj_work, lwork - n, iwork, *info1);
        }
        _gedmd_lacpy(#char "A", m, n, z.*, ldz, x, ldx, RT);
        t_or_n = #char "N";
        xscl1 = fget(work.*, n + 1);
        xscl2 = fget(work.*, n + 2);
        if xscl1 != xscl2 {
            _gedmd_lascl(#char "G", 0, 0, xscl1, xscl2, m, n, y, ldy, *info2, RT, RT);
        }
    }

    if info1 > 0 {
        info.* = 2;
        return;
    }
    if fget(work.*, 1) == zero_r {
        k.* = 0;
        info.* = -8;
        xerbla(routine_name, -info.*);
        return;
    }

    switch nrnk {
    case -1:
        k.* = 1;
        for i:2..numrnk {
            if (fget(work.*, i) <= fget(work.*, 1)*tol) || (fget(work.*, i) <= small) break;
            k.* += 1;
        }
    case -2:
        k.* = 1;
        for i:1..numrnk-1 {
            if (fget(work.*, i + 1) <= fget(work.*, i)*tol) || (fget(work.*, i) <= small) break;
            k.* += 1;
        }
    case:
        k.* = 1;
        for i:2..nrnk {
            if fget(work.*, i) <= small break;
            k.* += 1;
        }
    }

    if lsame(t_or_n, #char "N") {
        for i:1..k.* {
            wi := _gedmd_matrix_column(w.*, ldw, 1, i, n);
            _gedmd_scal(n, one_r/fget(work.*, i), *wi, 1, RT, RT);
        }
    } else {
        for i:1..k.* {
            fset(work, n + i, one_r/fget(work.*, i));
        }
        for j:1..n {
            for i:1..k.* {
                fset(w, i, j, fget(work.*, n + i) * fget(w.*, i, j));
            }
        }
    }

    if wntref {
        _gedmd_gemm(#char "N", t_or_n, m, k.*, n, one_r, y.*, ldy, w.*, ldw, zero_r, z, ldz, RT);
        _gedmd_lacpy(#char "A", m, k.*, z.*, ldz, b, ldb, RT);
        _gedmd_gemm(#char "T", #char "N", k.*, k.*, m, one_r, x.*, ldx, z.*, ldz, zero_r, s, lds, RT);
    } else {
        _gedmd_gemm(#char "T", #char "N", k.*, n, m, one_r, x.*, ldx, y.*, ldy, zero_r, z, ldz, RT);
        _gedmd_gemm(#char "N", t_or_n, k.*, k.*, n, one_r, z.*, ldz, w.*, ldw, zero_r, s, lds, RT);
        if wntres || wntex {
            if lsame(t_or_n, #char "N") {
                _gedmd_lacpy(#char "A", n, k.*, w.*, ldw, z, ldz, RT);
            } else {
                _gedmd_lacpy(#char "A", k.*, n, w.*, ldw, z, ldz, RT);
            }
        }
    }

    eig_work := _gedmd_subvector(work.*, n + 1, lwork - n);
    #if RT == float32 {
        sgeev(#char "N", jobzl, k.*, s, lds, reig, imeig, w, ldw, w, ldw, *eig_work, lwork - n, *info1);
    } else {
        dgeev(#char "N", jobzl, k.*, s, lds, reig, imeig, w, ldw, w, ldw, *eig_work, lwork - n, *info1);
    }
    if info1 > 0 {
        info.* = 3;
        return;
    }

    if wntvec || wntex {
        if wntres {
            if wntref {
                _gedmd_gemm(#char "N", #char "N", m, k.*, k.*, one_r, z.*, ldz, w.*, ldw, zero_r, y, ldy, RT);
            } else {
                _gedmd_gemm(t_or_n, #char "N", n, k.*, k.*, one_r, z.*, ldz, w.*, ldw, zero_r, s, lds, RT);
                _gedmd_gemm(#char "N", #char "N", m, k.*, n, one_r, y.*, ldy, s.*, lds, zero_r, z, ldz, RT);
                _gedmd_lacpy(#char "A", m, k.*, z.*, ldz, y, ldy, RT);
                if wntex {
                    _gedmd_lacpy(#char "A", m, k.*, z.*, ldz, b, ldb, RT);
                }
            }
        } else if wntex {
            _gedmd_gemm(t_or_n, #char "N", n, k.*, k.*, one_r, z.*, ldz, w.*, ldw, zero_r, s, lds, RT);
            _gedmd_gemm(#char "N", #char "N", m, k.*, n, one_r, y.*, ldy, s.*, lds, zero_r, b, ldb, RT);
        }

        if wntvec {
            _gedmd_gemm(#char "N", #char "N", m, k.*, k.*, one_r, x.*, ldx, w.*, ldw, zero_r, z, ldz, RT);
        }

        if wntres {
            i = 1;
            while i <= k.* {
                if fget(imeig.*, i) == zero_r {
                    zi := _gedmd_matrix_column(z.*, ldz, 1, i, m);
                    yi := _gedmd_matrix_column(y.*, ldy, 1, i, m);
                    _gedmd_axpy(m, -fget(reig.*, i), zi, *yi, RT);
                    fset(res, i, _gedmd_nrm2(m, yi, 1, RT, RT));
                    i += 1;
                } else {
                    ab[0][0] =  fget(reig.*, i);
                    ab[1][0] = -fget(imeig.*, i);
                    ab[0][1] =  fget(imeig.*, i);
                    ab[1][1] =  fget(reig.*, i);

                    zpair := _gedmd_matrix_storage(z.*, ldz, 1, i, m, 2);
                    ypair := _gedmd_matrix_storage(y.*, ldy, 1, i, m, 2);
                    scratch := _gedmd_subvector(work.*, n + 1, max(1, m));

                    _gedmd_gemm(#char "N", #char "N", m, 2, 2, -one_r, zpair, ldz, ab, 2, one_r, *ypair, ldy, RT);
                    fset(res, i, _gedmd_lange(#char "F", m, 2, ypair, ldy, *scratch, RT, RT));
                    fset(res, i + 1, fget(res.*, i));
                    i += 2;
                }
            }
        }
    }

    if whtsvd == 4 {
        fset(work, n + 1, xscl1);
        fset(work, n + 2, xscl2);
    }

    if !badxy {
        info.* = 0;
    } else {
        info.* = 4;
    }
}

_gedmd_complex :: (
    $$jobs: u8,
    $$jobz: u8,
    $$jobr: u8,
    $$jobf: u8,
    $$whtsvd: int,
    $$m: int,
    $$n: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    y: *$Y/ArrayType2D,
    $$ldy: int,
    $$nrnk: int,
    tol: $RT,
    k: *int,
    eigs: *$EIGS/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    res: *$RES/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    w: *$W/ArrayType2D,
    $$ldw: int,
    s: *$S/ArrayType2D,
    $$lds: int,
    zwork: *$ZWORK/ArrayType1D,
    $$lzwork: int,
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
           X.data_type == MT &&
           Y.data_type == MT &&
           EIGS.data_type == MT &&
           Z.data_type == MT &&
           RES.data_type == RT &&
           B.data_type == MT &&
           W.data_type == MT &&
           S.data_type == MT &&
           ZWORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    one_t :: #run one(MT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == C32 { "CGEDMD" } else { "ZGEDMD" };

    rdummy : [2] RT;

    badxy, lquery, sccolx, sccoly, wntex, wntref, wntres, wntvec : bool = ---;
    jobzl, t_or_n, jsvopt : u8 = ---;
    i, j, iminwr, info1, info2, lwrkev, lwrsdd, lwrsvd, lwrsvj, lwrsvq, mlwork, mlrwrk, mwrkev, mwrsdd, mwrsvd, mwrsvj, mwrsvq, numrnk, olwork : int = ---;
    ofl, rootsc, scale, small, ssum, xscl1, xscl2 : RT = ---;

    wntres = lsame(jobr, #char "R");
    sccolx = lsame(jobs, #char "S") || lsame(jobs, #char "C");
    sccoly = lsame(jobs, #char "Y");
    wntvec = lsame(jobz, #char "V");
    wntref = lsame(jobf, #char "R");
    wntex  = lsame(jobf, #char "E");
    info.* = 0;
    lquery = (lzwork == -1) || (liwork == -1) || (lrwork == -1);

    if !(sccolx || sccoly || lsame(jobs, #char "N")) {
        info.* = -1;
    } else if !(wntvec || lsame(jobz, #char "N") || lsame(jobz, #char "F")) {
        info.* = -2;
    } else if !(wntres || lsame(jobr, #char "N")) || (wntres && !wntvec) {
        info.* = -3;
    } else if !(wntref || wntex || lsame(jobf, #char "N")) {
        info.* = -4;
    } else if !((whtsvd == 1) || (whtsvd == 2) || (whtsvd == 3) || (whtsvd == 4)) {
        info.* = -5;
    } else if m < 0 {
        info.* = -6;
    } else if (n < 0) || (n > m) {
        info.* = -7;
    } else if ldx < m {
        info.* = -9;
    } else if ldy < m {
        info.* = -11;
    } else if !((nrnk == -2) || (nrnk == -1) || ((nrnk >= 1) && (nrnk <= n))) {
        info.* = -12;
    } else if (tol < zero_r) || (tol >= one_r) {
        info.* = -13;
    } else if ldz < m {
        info.* = -17;
    } else if (wntref || wntex) && (ldb < m) {
        info.* = -20;
    } else if ldw < n {
        info.* = -22;
    } else if lds < n {
        info.* = -24;
    }

    if info.* == 0 {
        if n == 0 {
            if lquery {
                fset(iwork, 1, 1);
                fset(rwork, 1, one_r);
                fset(zwork, 1, _gedmd_complex_from_real(ncast(RT, 2), RT, MT));
                fset(zwork, 2, _gedmd_complex_from_real(ncast(RT, 2), RT, MT));
            } else {
                k.* = 0;
            }
            info.* = 1;
            return;
        }

        iminwr = 1;
        mlrwrk = max(1, n);
        mlwork = 2;
        olwork = 2;

        switch whtsvd {
        case 1:
            mwrsvd = max(1, 2*min(m, n) + max(m, n));
            mlwork = max(mlwork, mwrsvd);
            mlrwrk = max(mlrwrk, n + 5*min(m, n));
            if lquery {
                #if MT == C32 {
                    cgesvd(#char "O", #char "S", m, n, x, ldx, rwork, b, ldb, w, ldw, zwork, -1, rdummy, *info1);
                } else {
                    zgesvd(#char "O", #char "S", m, n, x, ldx, rwork, b, ldb, w, ldw, zwork, -1, rdummy, *info1);
                }
                lwrsvd = cast(int) real(fget(zwork.*, 1));
                olwork = max(olwork, lwrsvd);
            }
        case 2:
            mwrsdd = 2*min(m, n)*min(m, n) + 2*min(m, n) + max(m, n);
            mlwork = max(mlwork, mwrsdd);
            iminwr = 8*min(m, n);
            mlrwrk = max(mlrwrk, n + max(5*min(m, n)*min(m, n) + 7*min(m, n),
                                         5*min(m, n)*min(m, n) + 5*min(m, n),
                                         2*max(m, n)*min(m, n) + 2*min(m, n)*min(m, n) + min(m, n)));
            if lquery {
                #if MT == C32 {
                    cgesdd(#char "O", m, n, x, ldx, rwork, b, ldb, w, ldw, zwork, -1, rdummy, iwork, *info1);
                } else {
                    zgesdd(#char "O", m, n, x, ldx, rwork, b, ldb, w, ldw, zwork, -1, rdummy, iwork, *info1);
                }
                lwrsdd = max(mwrsdd, cast(int) real(fget(zwork.*, 1)));
                olwork = max(olwork, lwrsdd);
            }
        case 3:
            numrnk = 0;
            #if MT == C32 {
                cgesvdq(#char "H", #char "P", #char "N", #char "R", #char "R", m, n, x, ldx, rwork, z, ldz, w, ldw,
                        *numrnk, iwork, -1, zwork, -1, rdummy, -1, *info1);
            } else {
                zgesvdq(#char "H", #char "P", #char "N", #char "R", #char "R", m, n, x, ldx, rwork, z, ldz, w, ldw,
                        *numrnk, iwork, -1, zwork, -1, rdummy, -1, *info1);
            }
            iminwr = fget(iwork.*, 1);
            mwrsvq = cast(int) real(fget(zwork.*, 2));
            mlwork = max(mlwork, mwrsvq);
            mlrwrk = max(mlrwrk, n + cast(int) rdummy[0]);
            if lquery {
                lwrsvq = cast(int) real(fget(zwork.*, 1));
                olwork = max(olwork, lwrsvq);
            }
        case 4:
            jsvopt = #char "J";
            #if MT == C32 {
                cgejsv(#char "F", #char "U", jsvopt, #char "N", #char "N", #char "P", m, n, x, ldx, rwork, z, ldz, w, ldw,
                       zwork, -1, rdummy, -1, iwork, *info1);
            } else {
                zgejsv(#char "F", #char "U", jsvopt, #char "R", #char "N", #char "P", m, n, x, ldx, rwork, z, ldz, w, ldw,
                       zwork, -1, rdummy, -1, iwork, *info1);
            }
            iminwr = fget(iwork.*, 1);
            mwrsvj = cast(int) real(fget(zwork.*, 2));
            mlwork = max(mlwork, mwrsvj);
            mlrwrk = max(mlrwrk, n + max(7, cast(int) rdummy[0]));
            if lquery {
                lwrsvj = cast(int) real(fget(zwork.*, 1));
                olwork = max(olwork, lwrsvj);
            }
        }

        if wntvec || wntex || lsame(jobz, #char "F") {
            jobzl = #char "V";
        } else {
            jobzl = #char "N";
        }

        mwrkev = max(1, 2*n);
        mlwork = max(mlwork, mwrkev);
        mlrwrk = max(mlrwrk, n + 2*n);
        if lquery {
            #if MT == C32 {
                cgeev(#char "N", jobzl, n, s, lds, eigs, w, ldw, w, ldw, zwork, -1, rwork, *info1);
            } else {
                zgeev(#char "N", jobzl, n, s, lds, eigs, w, ldw, w, ldw, zwork, -1, rwork, *info1);
            }
            lwrkev = cast(int) real(fget(zwork.*, 1));
            olwork = max(olwork, lwrkev);
            olwork = max(2, olwork);
        }

        if (liwork < iminwr) && !lquery {
            info.* = -30;
        }
        if (lrwork < mlrwrk) && !lquery {
            info.* = -28;
        }
        if (lzwork < mlwork) && !lquery {
            info.* = -26;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        fset(iwork, 1, iminwr);
        fset(rwork, 1, ncast(RT, mlrwrk));
        fset(zwork, 1, _gedmd_complex_from_real(ncast(RT, mlwork), RT, MT));
        fset(zwork, 2, _gedmd_complex_from_real(ncast(RT, olwork), RT, MT));
        return;
    }

    #if MT == C32 {
        ofl = _gedmd_lamch_overflow(RT) * _gedmd_lamch_precision(RT);
    } else {
        ofl = _gedmd_lamch_overflow(RT);
    }
    small = _gedmd_lamch_safe_minimum(RT);
    badxy = false;

    if sccolx {
        k.* = 0;
        for i:1..n {
            ssum = one_r;
            scale = zero_r;
            x_col := _gedmd_matrix_column(x.*, ldx, 1, i, m);
            _gedmd_lassq(m, x_col, 1, *scale, *ssum, RT, MT);
            if _gedmd_isnan(scale, RT) || _gedmd_isnan(ssum, RT) {
                k.* = 0;
                info.* = -8;
                xerbla(routine_name, -info.*);
            }
            if (scale != zero_r) && (ssum != zero_r) {
                rootsc = sqrt(ssum);
                if scale >= (ofl / rootsc) {
                    x_col_mat := _gedmd_matrix_storage(x.*, ldx, 1, i, m, 1);
                    _gedmd_lascl(#char "G", 0, 0, scale, one_r/rootsc, m, 1, *x_col_mat, ldx, *info2, RT, MT);
                    fset(rwork, i, -scale * (rootsc / ncast(RT, m)));
                } else {
                    fset(rwork, i, scale * rootsc);
                    x_col_mat := _gedmd_matrix_storage(x.*, ldx, 1, i, m, 1);
                    _gedmd_lascl(#char "G", 0, 0, fget(rwork.*, i), one_r, m, 1, *x_col_mat, ldx, *info2, RT, MT);
                }
            } else {
                fset(rwork, i, zero_r);
                k.* += 1;
            }
        }
        if k.* == n {
            k.* = 0;
            info.* = -8;
            xerbla(routine_name, -info.*);
            return;
        }
        for i:1..n {
            yi := _gedmd_matrix_column(y.*, ldy, 1, i, m);
            if fget(rwork.*, i) > zero_r {
                _gedmd_scal(m, one_r/fget(rwork.*, i), *yi, 1, RT, MT);
            } else if fget(rwork.*, i) < zero_r {
                y_col_mat := _gedmd_matrix_storage(y.*, ldy, 1, i, m, 1);
                _gedmd_lascl(#char "G", 0, 0, -fget(rwork.*, i), one_r/ncast(RT, m), m, 1, *y_col_mat, ldy, *info2, RT, MT);
            } else if _gedmd_is_nonzero(fget(yi, _gedmd_iamax_1based(m, yi, 1, RT, MT)), RT, MT) {
                badxy = true;
                if lsame(jobs, #char "C") {
                    _gedmd_scal(m, zero_r, *yi, 1, RT, MT);
                }
            }
        }
    }

    if sccoly {
        for i:1..n {
            ssum = one_r;
            scale = zero_r;
            y_col := _gedmd_matrix_column(y.*, ldy, 1, i, m);
            _gedmd_lassq(m, y_col, 1, *scale, *ssum, RT, MT);
            if _gedmd_isnan(scale, RT) || _gedmd_isnan(ssum, RT) {
                k.* = 0;
                info.* = -10;
                xerbla(routine_name, -info.*);
            }
            if (scale != zero_r) && (ssum != zero_r) {
                rootsc = sqrt(ssum);
                if scale >= (ofl / rootsc) {
                    y_col_mat := _gedmd_matrix_storage(y.*, ldy, 1, i, m, 1);
                    _gedmd_lascl(#char "G", 0, 0, scale, one_r/rootsc, m, 1, *y_col_mat, ldy, *info2, RT, MT);
                    fset(rwork, i, -scale * (rootsc / ncast(RT, m)));
                } else {
                    fset(rwork, i, scale * rootsc);
                    y_col_mat := _gedmd_matrix_storage(y.*, ldy, 1, i, m, 1);
                    _gedmd_lascl(#char "G", 0, 0, fget(rwork.*, i), one_r, m, 1, *y_col_mat, ldy, *info2, RT, MT);
                }
            } else {
                fset(rwork, i, zero_r);
            }
        }
        for i:1..n {
            xi := _gedmd_matrix_column(x.*, ldx, 1, i, m);
            if fget(rwork.*, i) > zero_r {
                _gedmd_scal(m, one_r/fget(rwork.*, i), *xi, 1, RT, MT);
            } else if fget(rwork.*, i) < zero_r {
                x_col_mat := _gedmd_matrix_storage(x.*, ldx, 1, i, m, 1);
                _gedmd_lascl(#char "G", 0, 0, -fget(rwork.*, i), one_r/ncast(RT, m), m, 1, *x_col_mat, ldx, *info2, RT, MT);
            } else if _gedmd_is_nonzero(fget(xi, _gedmd_iamax_1based(m, xi, 1, RT, MT)), RT, MT) {
                badxy = true;
            }
        }
    }
