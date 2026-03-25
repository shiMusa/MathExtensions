// ? LAPACK *gghd3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgghd3/cgghd3/dgghd3/zgghd3
// *GGHD3 reduces a matrix pair (A, B) to generalized upper Hessenberg
// form using the blocked interface and workspace query rules.

sgghd3 :: #bake_arguments _gghd3(RT=float32, MT=float32);
dgghd3 :: #bake_arguments _gghd3(RT=float64, MT=float64);
cgghd3 :: #bake_arguments _gghd3(RT=float32, MT=C32);
zgghd3 :: #bake_arguments _gghd3(RT=float64, MT=C64);

#scope_file
_gghd3_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_gghd3_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _gghd3_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gghd3_matrix_col_prefix :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _gghd3_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gghd3_matrix_row_suffix :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _gghd3_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_gghd3_work_storage :: inline (
    work: *$WORK/ArrayType1D, $$offset: int, $$rows: int, $$cols: int, $$ldw: int
) -> [] WORK.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else ldw*(cols - 1) + rows;
    return _gghd3_raw_slice(work.data + (offset - 1), count);
}

_gghd3_mt_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gghd3_conj_if_needed :: inline (x: $MT, $$do_conj: bool, $MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == float32 || MT == float64 {
        return x;
    } else {
        if do_conj then return conj(x);
        return x;
    }
}

_gghd3_trans_h :: inline ($MT: Type) -> u8 {
    #if MT == C32 || MT == C64 {
        return #char "C";
    } else {
        return #char "T";
    }
}

_gghd3_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_gghd3_laset :: (
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

_gghd3_lacpy :: (
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

_gghd3_gemm :: (
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

_gghd3_gemv :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    $$incx: int,
    beta: MT,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == float64 {
        dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_gghd3_trmv :: (
    $$uplo: u8,
    $$trans: u8,
    $$diag: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT;
} {
    #if MT == float32 {
        strmv(uplo, trans, diag, n, a, lda, x, incx);
    } else if MT == float64 {
        dtrmv(uplo, trans, diag, n, a, lda, x, incx);
    } else if MT == C32 {
        ctrmv(uplo, trans, diag, n, a, lda, x, incx);
    } else {
        ztrmv(uplo, trans, diag, n, a, lda, x, incx);
    }
}

_gghd3_rot :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    c: $RT,
    s: $MT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else if MT == float64 {
        drot(n, x, incx, y, incy, c, s);
    } else if MT == C32 {
        crot(n, x, incx, y, incy, c, s);
    } else {
        zrot(n, x, incx, y, incy, c, s);
    }
}

_gghd3_lartg :: (
    f: $MT,
    g: MT,
    c: *$RT,
    s: *MT,
    r: *MT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 {
        slartg(f, g, c, s, r);
    } else if MT == float64 {
        dlartg(f, g, c, s, r);
    } else if MT == C32 {
        clartg(f, g, c, s, r);
    } else {
        zlartg(f, g, c, s, r);
    }
}

_gghd3_unm22 :: (
    $$side: u8,
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$n1: int,
    $$n2: int,
    q: $Q/ArrayType2D,
    $$ldq: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           Q.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sorm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dorm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
    } else {
        zunm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
    }
}

_gghd3 :: (
    $$compq: u8,
    $$compz: u8,
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    z: *$Z/ArrayType2D,
    $$ldz: int,
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
           Q.data_type == MT &&
           Z.data_type == MT &&
           WORK.data_type == MT;
} {
    zero_t :: #run _gghd3_mt_from_real(zero(RT), RT, MT);
    one_t  :: #run _gghd3_mt_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == float32 { "SGGHD3" } else if MT == float64 { "DGGHD3" } else if MT == C32 { "CGGHD3" } else { "ZGGHD3" };
    trans_h :: u8 = _gghd3_trans_h(MT);

    initq, initz, lquery, wantq, wantz : bool = ---;
    blk22 : bool = ---;
    c : RT = ---;
    ctemp, c1, c2, s, s1, s2, temp, temp1, temp2, temp3 : MT = ---;
    cola, i, ierr, j, j0, jcol, jj, jrow, k, kacc22, len, lwkopt, n2nb, nb, nblst, nbmin, nh, nnb, nx, ppw, ppwo, pw, top, topq : int = ---;
    compq2, compz2 : u8 = ---;

    info.* = 0;
    nb = ilaenv(1, routine_name, " "s, n, ilo, ihi, -1);
    nh = ihi - ilo + 1;
    if nh <= 1 {
        lwkopt = 1;
    } else {
        lwkopt = 6*n*nb;
    }
    fset(work, 1, _gghd3_mt_from_real(_gghd3_roundup_lwork(lwkopt, RT), RT, MT));

    initq = lsame(compq, #char "I");
    wantq = initq || lsame(compq, #char "V");
    initz = lsame(compz, #char "I");
    wantz = initz || lsame(compz, #char "V");
    lquery = lwork == -1;

    if !lsame(compq, #char "N") && !wantq {
        info.* = -1;
    } else if !lsame(compz, #char "N") && !wantz {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ilo < 1 {
        info.* = -4;
    } else if ihi > n || ihi < ilo - 1 {
        info.* = -5;
    } else if lda < max(1, n) {
        info.* = -7;
    } else if ldb < max(1, n) {
        info.* = -9;
    } else if (wantq && ldq < n) || ldq < 1 {
        info.* = -11;
    } else if (wantz && ldz < n) || ldz < 1 {
        info.* = -13;
    } else if lwork < 1 && !lquery {
        info.* = -15;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if initq {
        _gghd3_laset(#char "A", n, n, zero_t, one_t, q, ldq, MT);
    }
    if initz {
        _gghd3_laset(#char "A", n, n, zero_t, one_t, z, ldz, MT);
    }

    if n > 1 {
        for j:1..n - 1 {
            for i:j + 1..n {
                fset(b, i, j, zero_t);
            }
        }
    }

    if nh <= 1 {
        fset(work, 1, one_t);
        return;
    }

    nbmin = ilaenv(2, routine_name, " "s, n, ilo, ihi, -1);
    if nb > 1 && nb < nh {
        nx = max(nb, ilaenv(3, routine_name, " "s, n, ilo, ihi, -1));
        if nx < nh {
            if lwork < lwkopt {
                nbmin = max(2, ilaenv(2, routine_name, " "s, n, ilo, ihi, -1));
                if lwork >= 6*n*nbmin {
                    nb = lwork / (6*n);
                } else {
                    nb = 1;
                }
            }
        }
    }

    jcol = ilo;
    if !(nb < nbmin || nb >= nh) {
        kacc22 = ilaenv(16, routine_name, " "s, n, ilo, ihi, -1);
        blk22 = kacc22 == 2;

        while jcol <= ihi - 2 {
            nnb = min(nb, ihi - jcol - 1);
            n2nb = (ihi - jcol - 1) / nnb - 1;
            nblst = ihi - jcol - n2nb*nnb;

            work_nblst_raw := _gghd3_work_storage(work, 1, nblst, nblst, nblst);
            work_nblst := array2d_view(work_nblst_raw, nblst, row_major=false);
            _gghd3_laset(#char "A", nblst, nblst, zero_t, one_t, *work_nblst, nblst, MT);
            pw = nblst*nblst + 1;
            for ii:1..n2nb {
                work_blk_raw := _gghd3_work_storage(work, pw, 2*nnb, 2*nnb, 2*nnb);
                work_blk := array2d_view(work_blk_raw, 2*nnb, row_major=false);
                _gghd3_laset(#char "A", 2*nnb, 2*nnb, zero_t, one_t, *work_blk, 2*nnb, MT);
                pw += 4*nnb*nnb;
            }

            for j:jcol..jcol + nnb - 1 {
                for <i:ihi..j + 2 {
                    temp = fget(a.*, i - 1, j);
                    _gghd3_lartg(temp, fget(a.*, i, j), *c, *s, *temp, RT, MT);
                    fset(a, i - 1, j, temp);
                    fset(a, i, j, _gghd3_mt_from_real(c, RT, MT));
                    fset(b, i, j, s);
                }

                ppw = (nblst + 1)*(nblst - 2) - j + jcol + 1;
                len = 2 + j - jcol;
                jrow = j + n2nb*nnb + 2;
                for <i:ihi..jrow {
                    ctemp = fget(a.*, i, j);
                    s = fget(b.*, i, j);
                    for jj:ppw..ppw + len - 1 {
                        temp = fget(work.*, jj + nblst);
                        fset(work, jj + nblst, ctemp*temp - s*fget(work.*, jj));
                        fset(work, jj, _gghd3_conj_if_needed(s, MT == C32 || MT == C64, MT)*temp + ctemp*fget(work.*, jj));
                    }
                    len += 1;
                    ppw -= nblst + 1;
                }

                ppwo = nblst*nblst + (nnb + j - jcol - 1)*2*nnb + nnb;
                j0 = jrow - nnb;
                for <jrow:j0..j + 2 {
                    ppw = ppwo;
                    len = 2 + j - jcol;
                    for <i:jrow + nnb - 1..jrow {
                        ctemp = fget(a.*, i, j);
                        s = fget(b.*, i, j);
                        for jj:ppw..ppw + len - 1 {
                            temp = fget(work.*, jj + 2*nnb);
                            fset(work, jj + 2*nnb, ctemp*temp - s*fget(work.*, jj));
                            fset(work, jj, _gghd3_conj_if_needed(s, MT == C32 || MT == C64, MT)*temp + ctemp*fget(work.*, jj));
                        }
                        len += 1;
                        ppw -= 2*nnb + 1;
                    }
                    ppwo += 4*nnb*nnb;
                }

                if jcol <= 2 { top = 0; } else { top = jcol; }

                for <jj:n..j + 1 {
                    for <i:min(jj + 1, ihi)..j + 2 {
                        ctemp = fget(a.*, i, j);
                        s = fget(b.*, i, j);
                        temp = fget(b.*, i, jj);
                        fset(b, i, jj, ctemp*temp - _gghd3_conj_if_needed(s, MT == C32 || MT == C64, MT)*fget(b.*, i - 1, jj));
                        fset(b, i - 1, jj, s*temp + ctemp*fget(b.*, i - 1, jj));
                    }

                    if jj < ihi {
                        temp = fget(b.*, jj + 1, jj + 1);
                        _gghd3_lartg(temp, fget(b.*, jj + 1, jj), *c, *s, *temp, RT, MT);
                        fset(b, jj + 1, jj + 1, temp);
                        fset(b, jj + 1, jj, zero_t);
                        brow1 := _gghd3_matrix_col_prefix(b.*, ldb, top + 1, jj + 1, jj - top);
                        brow2 := _gghd3_matrix_col_prefix(b.*, ldb, top + 1, jj, jj - top);
                        _gghd3_rot(jj - top, *brow1, 1, *brow2, 1, c, s, RT, MT);
                        fset(a, jj + 1, j, _gghd3_mt_from_real(c, RT, MT));
                        fset(b, jj + 1, j, -_gghd3_conj_if_needed(s, MT == C32 || MT == C64, MT));
                    }
                }

                jj = (ihi - j - 1) % 3;
                for <i:ihi - j - 3..jj + 1 {
                    ctemp = fget(a.*, j + 1 + i, j);
                    s = -fget(b.*, j + 1 + i, j);
                    c1 = fget(a.*, j + 2 + i, j);
                    s1 = -fget(b.*, j + 2 + i, j);
                    c2 = fget(a.*, j + 3 + i, j);
                    s2 = -fget(b.*, j + 3 + i, j);

                    for k:top + 1..ihi {
                        temp = fget(a.*, k, j + i);
                        temp1 = fget(a.*, k, j + i + 1);
                        temp2 = fget(a.*, k, j + i + 2);
                        temp3 = fget(a.*, k, j + i + 3);
                        fset(a, k, j + i + 3, c2*temp3 + _gghd3_conj_if_needed(s2, MT == C32 || MT == C64, MT)*temp2);
                        temp2 = -s2*temp3 + c2*temp2;
                        fset(a, k, j + i + 2, c1*temp2 + _gghd3_conj_if_needed(s1, MT == C32 || MT == C64, MT)*temp1);
                        temp1 = -s1*temp2 + c1*temp1;
                        fset(a, k, j + i + 1, ctemp*temp1 + _gghd3_conj_if_needed(s, MT == C32 || MT == C64, MT)*temp);
                        fset(a, k, j + i, -s*temp1 + ctemp*temp);
                    }
                }

                if jj > 0 {
                    for <i:jj..1 {
                        c = real(fget(a.*, j + 1 + i, j));
                        arow1 := _gghd3_matrix_row_suffix(a.*, lda, top + 1, j + i + 1, ihi - top);
                        arow2 := _gghd3_matrix_row_suffix(a.*, lda, top + 1, j + i, ihi - top);
                        _gghd3_rot(ihi - top, *arow1, 1, *arow2, 1, c, -_gghd3_conj_if_needed(fget(b.*, j + 1 + i, j), MT == C32 || MT == C64, MT), RT, MT);
                    }
                }

                if j < jcol + nnb - 1 {
                    len = 1 + j - jcol;
                    jrow = ihi - nblst + 1;

                    work_vec := _gghd3_raw_slice(work.data + (pw - 1), nblst);
                    work_trail := _gghd3_raw_slice(work.data + (pw + len - 1), nblst - len);

                    work_u := array2d_view(_gghd3_work_storage(work, 1, nblst, nblst, nblst), nblst, row_major=false);
                    a_tail := _gghd3_matrix_col_prefix(a.*, lda, jrow, j + 1, nblst);
                    _gghd3_gemv(trans_h, nblst, len, one_t, work_u.*, nblst, a_tail, 1, zero_t, *work_vec, 1, MT);

                    ppw = pw + len;
                    for ii:jrow..jrow + nblst - len - 1 {
                        fset(work, ppw, fget(a.*, ii, j + 1));
                        ppw += 1;
                    }

                    work_lower := array2d_view(_gghd3_work_storage(work, len*nblst + 1, nblst - len, nblst - len, nblst), nblst, row_major=false);
                    _gghd3_trmv(#char "L", trans_h, #char "N", nblst - len, work_lower.*, nblst, *work_trail, 1, MT);

                    work_rect := array2d_view(_gghd3_work_storage(work, (len + 1)*nblst - len + 1, len, nblst - len, nblst), nblst, row_major=false);
                    a_tail2 := _gghd3_matrix_col_prefix(a.*, lda, jrow + nblst - len, j + 1, nblst - len);
                    _gghd3_gemv(trans_h, len, nblst - len, one_t, work_rect.*, nblst, a_tail2, 1, one_t, *work_trail, 1, MT);

                    ppw = pw;
                    for ii:jrow..jrow + nblst - 1 {
                        fset(a, ii, j + 1, fget(work.*, ppw));
                        ppw += 1;
                    }

                    ppwo = 1 + nblst*nblst;
                    j0 = jrow - nnb;
                    for <jrow:j0..jcol + 1 {
                        work_tail := _gghd3_raw_slice(work.data + (pw + len - 1), nnb);
                        for ii:jrow..jrow + nnb - 1 {
                            fset(work, pw + len + ii - jrow, fget(a.*, ii, j + 1));
                        }
                        for ii:jrow + nnb..jrow + nnb + len - 1 {
                            fset(work, pw + ii - (jrow + nnb), fget(a.*, ii, j + 1));
                        }

                        upper_u := array2d_view(_gghd3_work_storage(work, ppwo + nnb, len, len, 2*nnb), 2*nnb, row_major=false);
                        work_head := _gghd3_raw_slice(work.data + (pw - 1), len);
                        _gghd3_trmv(#char "U", trans_h, #char "N", len, upper_u.*, 2*nnb, *work_head, 1, MT);

                        lower_u := array2d_view(_gghd3_work_storage(work, ppwo + 2*len*nnb, nnb, nnb, 2*nnb), 2*nnb, row_major=false);
                        _gghd3_trmv(#char "L", trans_h, #char "N", nnb, lower_u.*, 2*nnb, *work_tail, 1, MT);

                        u12 := array2d_view(_gghd3_work_storage(work, ppwo, nnb, len, 2*nnb), 2*nnb, row_major=false);
                        avec1 := _gghd3_matrix_col_prefix(a.*, lda, jrow, j + 1, nnb);
                        _gghd3_gemv(trans_h, nnb, len, one_t, u12.*, 2*nnb, avec1, 1, one_t, *work_head, 1, MT);

                        u21 := array2d_view(_gghd3_work_storage(work, ppwo + 2*len*nnb + nnb, len, nnb, 2*nnb), 2*nnb, row_major=false);
                        avec2 := _gghd3_matrix_col_prefix(a.*, lda, jrow + nnb, j + 1, len);
                        _gghd3_gemv(trans_h, len, nnb, one_t, u21.*, 2*nnb, avec2, 1, one_t, *work_tail, 1, MT);

                        ppw = pw;
                        for ii:jrow..jrow + len + nnb - 1 {
                            fset(a, ii, j + 1, fget(work.*, ppw));
                            ppw += 1;
                        }
                        ppwo += 4*nnb*nnb;
                    }
                }
            }

            cola = n - jcol - nnb + 1;
            j = ihi - nblst + 1;
            u_last := array2d_view(_gghd3_work_storage(work, 1, nblst, nblst, nblst), nblst, row_major=false);
            a_right_last := array2d_view(_gghd3_matrix_storage(a.*, lda, j, jcol + nnb, nblst, cola), lda, row_major=false);
            work_out_last := array2d_view(_gghd3_work_storage(work, pw, nblst, cola, nblst), nblst, row_major=false);
            _gghd3_gemm(trans_h, #char "N", nblst, cola, nblst, one_t, u_last.*, nblst, a_right_last.*, lda, zero_t, *work_out_last, nblst, MT);
            _gghd3_lacpy(#char "A", nblst, cola, work_out_last.*, nblst, *a_right_last, lda, MT);

            ppwo = nblst*nblst + 1;
            j0 = j - nnb;
            for <j:j0..jcol + 1 {
                ublk := array2d_view(_gghd3_work_storage(work, ppwo, 2*nnb, 2*nnb, 2*nnb), 2*nnb, row_major=false);
                a_right := array2d_view(_gghd3_matrix_storage(a.*, lda, j, jcol + nnb, 2*nnb, cola), lda, row_major=false);
                if blk22 {
                    work_unm := _gghd3_raw_slice(work.data + (pw - 1), lwork - pw + 1);
                    _gghd3_unm22(#char "L", trans_h, 2*nnb, cola, nnb, nnb, ublk.*, 2*nnb, *a_right, lda, *work_unm, lwork - pw + 1, *ierr, MT);
                } else {
                    work_out := array2d_view(_gghd3_work_storage(work, pw, 2*nnb, cola, 2*nnb), 2*nnb, row_major=false);
                    _gghd3_gemm(trans_h, #char "N", 2*nnb, cola, 2*nnb, one_t, ublk.*, 2*nnb, a_right.*, lda, zero_t, *work_out, 2*nnb, MT);
                    _gghd3_lacpy(#char "A", 2*nnb, cola, work_out.*, 2*nnb, *a_right, lda, MT);
                }
                ppwo += 4*nnb*nnb;
            }

            if wantq {
                j = ihi - nblst + 1;
                if initq {
                    topq = max(2, j - jcol + 1);
                    nh = ihi - topq + 1;
                } else {
                    topq = 1;
                    nh = n;
                }
                q_last := array2d_view(_gghd3_matrix_storage(q.*, ldq, topq, j, nh, nblst), ldq, row_major=false);
                work_out_q := array2d_view(_gghd3_work_storage(work, pw, nh, nblst, nh), nh, row_major=false);
                _gghd3_gemm(#char "N", #char "N", nh, nblst, nblst, one_t, q_last.*, ldq, u_last.*, nblst, zero_t, *work_out_q, nh, MT);
                _gghd3_lacpy(#char "A", nh, nblst, work_out_q.*, nh, *q_last, ldq, MT);

                ppwo = nblst*nblst + 1;
                j0 = j - nnb;
                for <j:j0..jcol + 1 {
                    if initq {
                        topq = max(2, j - jcol + 1);
                        nh = ihi - topq + 1;
                    }
                    q_blk := array2d_view(_gghd3_matrix_storage(q.*, ldq, topq, j, nh, 2*nnb), ldq, row_major=false);
                    ublk := array2d_view(_gghd3_work_storage(work, ppwo, 2*nnb, 2*nnb, 2*nnb), 2*nnb, row_major=false);
                    if blk22 {
                        work_unm := _gghd3_raw_slice(work.data + (pw - 1), lwork - pw + 1);
                        _gghd3_unm22(#char "R", #char "N", nh, 2*nnb, nnb, nnb, ublk.*, 2*nnb, *q_blk, ldq, *work_unm, lwork - pw + 1, *ierr, MT);
                    } else {
                        work_out := array2d_view(_gghd3_work_storage(work, pw, nh, 2*nnb, nh), nh, row_major=false);
                        _gghd3_gemm(#char "N", #char "N", nh, 2*nnb, 2*nnb, one_t, q_blk.*, ldq, ublk.*, 2*nnb, zero_t, *work_out, nh, MT);
                        _gghd3_lacpy(#char "A", nh, 2*nnb, work_out.*, nh, *q_blk, ldq, MT);
                    }
                    ppwo += 4*nnb*nnb;
                }
            }

            if wantz || top > 0 {
                work_nblst := array2d_view(_gghd3_work_storage(work, 1, nblst, nblst, nblst), nblst, row_major=false);
                _gghd3_laset(#char "A", nblst, nblst, zero_t, one_t, *work_nblst, nblst, MT);
                pw = nblst*nblst + 1;
                for ii:1..n2nb {
                    work_blk_raw := _gghd3_work_storage(work, pw, 2*nnb, 2*nnb, 2*nnb);
                    work_blk := array2d_view(work_blk_raw, 2*nnb, row_major=false);
                    _gghd3_laset(#char "A", 2*nnb, 2*nnb, zero_t, one_t, *work_blk, 2*nnb, MT);
                    pw += 4*nnb*nnb;
                }

                for j:jcol..jcol + nnb - 1 {
                    ppw = (nblst + 1)*(nblst - 2) - j + jcol + 1;
                    len = 2 + j - jcol;
                    jrow = j + n2nb*nnb + 2;
                    for <i:ihi..jrow {
                        ctemp = fget(a.*, i, j);
                        fset(a, i, j, zero_t);
                        s = fget(b.*, i, j);
                        fset(b, i, j, zero_t);
                        for jj:ppw..ppw + len - 1 {
                            temp = fget(work.*, jj + nblst);
                            fset(work, jj + nblst, ctemp*temp - _gghd3_conj_if_needed(s, MT == C32 || MT == C64, MT)*fget(work.*, jj));
                            fset(work, jj, s*temp + ctemp*fget(work.*, jj));
                        }
                        len += 1;
                        ppw -= nblst + 1;
                    }

                    ppwo = nblst*nblst + (nnb + j - jcol - 1)*2*nnb + nnb;
                    j0 = jrow - nnb;
                    for <jrow:j0..j + 2 {
                        ppw = ppwo;
                        len = 2 + j - jcol;
                        for <i:jrow + nnb - 1..jrow {
                            ctemp = fget(a.*, i, j);
                            fset(a, i, j, zero_t);
                            s = fget(b.*, i, j);
                            fset(b, i, j, zero_t);
                            for jj:ppw..ppw + len - 1 {
                                temp = fget(work.*, jj + 2*nnb);
                                fset(work, jj + 2*nnb, ctemp*temp - _gghd3_conj_if_needed(s, MT == C32 || MT == C64, MT)*fget(work.*, jj));
                                fset(work, jj, s*temp + ctemp*fget(work.*, jj));
                            }
                            len += 1;
                            ppw -= 2*nnb + 1;
                        }
                        ppwo += 4*nnb*nnb;
                    }
                }
            } else {
                if ihi - jcol - 1 > 0 {
                    lower_a := array2d_view(_gghd3_matrix_storage(a.*, lda, jcol + 2, jcol, ihi - jcol - 1, nnb), lda, row_major=false);
                    _gghd3_laset(#char "L", ihi - jcol - 1, nnb, zero_t, zero_t, *lower_a, lda, MT);
                    lower_b := array2d_view(_gghd3_matrix_storage(b.*, ldb, jcol + 2, jcol, ihi - jcol - 1, nnb), ldb, row_major=false);
                    _gghd3_laset(#char "L", ihi - jcol - 1, nnb, zero_t, zero_t, *lower_b, ldb, MT);
                }
            }

            if top > 0 {
                j = ihi - nblst + 1;
                a_top_last := array2d_view(_gghd3_matrix_storage(a.*, lda, 1, j, top, nblst), lda, row_major=false);
                work_out_top := array2d_view(_gghd3_work_storage(work, pw, top, nblst, top), top, row_major=false);
                _gghd3_gemm(#char "N", #char "N", top, nblst, nblst, one_t, a_top_last.*, lda, u_last.*, nblst, zero_t, *work_out_top, top, MT);
                _gghd3_lacpy(#char "A", top, nblst, work_out_top.*, top, *a_top_last, lda, MT);
                ppwo = nblst*nblst + 1;
                j0 = j - nnb;
                for <j:j0..jcol + 1 {
                    ublk := array2d_view(_gghd3_work_storage(work, ppwo, 2*nnb, 2*nnb, 2*nnb), 2*nnb, row_major=false);
                    a_top_blk := array2d_view(_gghd3_matrix_storage(a.*, lda, 1, j, top, 2*nnb), lda, row_major=false);
                    if blk22 {
                        work_unm := _gghd3_raw_slice(work.data + (pw - 1), lwork - pw + 1);
                        _gghd3_unm22(#char "R", #char "N", top, 2*nnb, nnb, nnb, ublk.*, 2*nnb, *a_top_blk, lda, *work_unm, lwork - pw + 1, *ierr, MT);
                    } else {
                        work_out := array2d_view(_gghd3_work_storage(work, pw, top, 2*nnb, top), top, row_major=false);
                        _gghd3_gemm(#char "N", #char "N", top, 2*nnb, 2*nnb, one_t, a_top_blk.*, lda, ublk.*, 2*nnb, zero_t, *work_out, top, MT);
                        _gghd3_lacpy(#char "A", top, 2*nnb, work_out.*, top, *a_top_blk, lda, MT);
                    }
                    ppwo += 4*nnb*nnb;
                }

                b_top_last := array2d_view(_gghd3_matrix_storage(b.*, ldb, 1, j, top, nblst), ldb, row_major=false);
                _gghd3_gemm(#char "N", #char "N", top, nblst, nblst, one_t, b_top_last.*, ldb, u_last.*, nblst, zero_t, *work_out_top, top, MT);
                _gghd3_lacpy(#char "A", top, nblst, work_out_top.*, top, *b_top_last, ldb, MT);
                ppwo = nblst*nblst + 1;
                j0 = j - nnb;
                for <j:j0..jcol + 1 {
                    ublk := array2d_view(_gghd3_work_storage(work, ppwo, 2*nnb, 2*nnb, 2*nnb), 2*nnb, row_major=false);
                    b_top_blk := array2d_view(_gghd3_matrix_storage(b.*, ldb, 1, j, top, 2*nnb), ldb, row_major=false);
                    if blk22 {
                        work_unm := _gghd3_raw_slice(work.data + (pw - 1), lwork - pw + 1);
                        _gghd3_unm22(#char "R", #char "N", top, 2*nnb, nnb, nnb, ublk.*, 2*nnb, *b_top_blk, ldb, *work_unm, lwork - pw + 1, *ierr, MT);
                    } else {
                        work_out := array2d_view(_gghd3_work_storage(work, pw, top, 2*nnb, top), top, row_major=false);
                        _gghd3_gemm(#char "N", #char "N", top, 2*nnb, 2*nnb, one_t, b_top_blk.*, ldb, ublk.*, 2*nnb, zero_t, *work_out, top, MT);
                        _gghd3_lacpy(#char "A", top, 2*nnb, work_out.*, top, *b_top_blk, ldb, MT);
                    }
                    ppwo += 4*nnb*nnb;
                }
            }

            if wantz {
                j = ihi - nblst + 1;
                if initq {
                    topq = max(2, j - jcol + 1);
                    nh = ihi - topq + 1;
                } else {
                    topq = 1;
                    nh = n;
                }
                z_last := array2d_view(_gghd3_matrix_storage(z.*, ldz, topq, j, nh, nblst), ldz, row_major=false);
                work_out_z := array2d_view(_gghd3_work_storage(work, pw, nh, nblst, nh), nh, row_major=false);
                _gghd3_gemm(#char "N", #char "N", nh, nblst, nblst, one_t, z_last.*, ldz, u_last.*, nblst, zero_t, *work_out_z, nh, MT);
                _gghd3_lacpy(#char "A", nh, nblst, work_out_z.*, nh, *z_last, ldz, MT);
                ppwo = nblst*nblst + 1;
                j0 = j - nnb;
                for <j:j0..jcol + 1 {
                    if initq {
                        topq = max(2, j - jcol + 1);
                        nh = ihi - topq + 1;
                    }
                    z_blk := array2d_view(_gghd3_matrix_storage(z.*, ldz, topq, j, nh, 2*nnb), ldz, row_major=false);
                    ublk := array2d_view(_gghd3_work_storage(work, ppwo, 2*nnb, 2*nnb, 2*nnb), 2*nnb, row_major=false);
                    if blk22 {
                        work_unm := _gghd3_raw_slice(work.data + (pw - 1), lwork - pw + 1);
                        _gghd3_unm22(#char "R", #char "N", nh, 2*nnb, nnb, nnb, ublk.*, 2*nnb, *z_blk, ldz, *work_unm, lwork - pw + 1, *ierr, MT);
                    } else {
                        work_out := array2d_view(_gghd3_work_storage(work, pw, nh, 2*nnb, nh), nh, row_major=false);
                        _gghd3_gemm(#char "N", #char "N", nh, 2*nnb, 2*nnb, one_t, z_blk.*, ldz, ublk.*, 2*nnb, zero_t, *work_out, nh, MT);
                        _gghd3_lacpy(#char "A", nh, 2*nnb, work_out.*, nh, *z_blk, ldz, MT);
                    }
                    ppwo += 4*nnb*nnb;
                }
            }

            jcol += nnb;
        }
    }

    compq2 = compq;
    compz2 = compz;
    if jcol != ilo {
        if wantq then compq2 = #char "V";
        if wantz then compz2 = #char "V";
    }

    if jcol < ihi {
        #if MT == float32 {
            sgghrd(compq2, compz2, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, *ierr);
        } else if MT == float64 {
            dgghrd(compq2, compz2, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, *ierr);
        } else if MT == C32 {
            cgghrd(compq2, compz2, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, *ierr);
        } else {
            zgghrd(compq2, compz2, n, jcol, ihi, a, lda, b, ldb, q, ldq, z, ldz, *ierr);
        }
    }

    fset(work, 1, _gghd3_mt_from_real(_gghd3_roundup_lwork(lwkopt, RT), RT, MT));
}
