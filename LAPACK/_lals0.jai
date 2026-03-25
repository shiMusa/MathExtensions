// ? LAPACK *lals0()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slals0/dlals0/clals0/zlals0
// *LALS0 applies back multiplying factors in the divide-and-conquer SVD
// least-squares path.

slals0 :: #bake_arguments _lals0_real(RT=float32);
dlals0 :: #bake_arguments _lals0_real(RT=float64);
clals0 :: #bake_arguments _lals0_complex(RT=float32, MT=C32);
zlals0 :: #bake_arguments _lals0_complex(RT=float64, MT=C64);

#scope_file
_lals0_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_lals0_matrix_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _lals0_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lals0_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _lals0_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lals0_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _lals0_raw_slice(x.data + (offset - 1), count);
}

_lals0_real_from_complex :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return real(z);
}

_lals0_imag_from_complex :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return imag(z);
}

_lals0_complex_from_parts :: inline (x: $RT, y: RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, y);
}

_lals0_copy_real :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 {
        scopy(n, x, incx, y, incy);
    } else {
        dcopy(n, x, incx, y, incy);
    }
}

_lals0_copy_complex :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_lals0_lacpy_real :: (
    $$uplo: u8, $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT && B.data_type == RT;
} {
    #if RT == float32 {
        slacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        dlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_lals0_lacpy_complex :: (
    $$uplo: u8, $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_lals0_lascl_real :: (
    $$type: u8, $$kl: int, $$ku: int, $$cfrom: $RT, $$cto: RT,
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT;
} {
    #if RT == float32 {
        slascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_lals0_lascl_complex :: (
    $$type: u8, $$kl: int, $$ku: int, $$cfrom: $RT, $$cto: RT,
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && A.data_type == MT;
} {
    #if MT == C32 {
        clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_lals0_scal_real :: (
    $$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        sscal(n, alpha, x, incx);
    } else {
        dscal(n, alpha, x, incx);
    }
}

_lals0_scal_complex_real :: (
    $$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_lals0_rot_real :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $$c: $RT, $$s: RT, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else {
        drot(n, x, incx, y, incy, c, s);
    }
}

_lals0_rot_complex :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $$c: $RT, $$s: RT, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        csrot(n, x, incx, y, incy, c, s);
    } else {
        zdrot(n, x, incx, y, incy, c, s);
    }
}

_lals0_gemv_real :: (
    $$trans: u8, $$m: int, $$n: int, $$alpha: $RT,
    a: *$A/ArrayType2D, $$lda: int, x: *$X/ArrayType1D, $$incx: int,
    $$beta: RT, y: *$Y/ArrayType1D, $$incy: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 {
        sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_lals0_lamc3 :: ($x: $RT, $y: RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamc3(x, y);
    } else {
        return dlamc3(x, y);
    }
}

_lals0_nrm2 :: ($$n: int, x: *$X/ArrayType1D, $$incx: int, $RT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return snrm2(n, x, incx);
    } else {
        return dnrm2(n, x, incx);
    }
}

_lals0_real :: (
    $$icompq: int,
    $$nl: int,
    $$nr: int,
    $$sqre: int,
    $$nrhs: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    bx: *$BX/ArrayType2D,
    $$ldbx: int,
    perm: *$PERM/ArrayType1D,
    $$givptr: int,
    givcol: *$GIVCOL/ArrayType2D,
    $$ldgcol: int,
    givnum: *$GIVNUM/ArrayType2D,
    $$ldgnum: int,
    poles: *$POLES/ArrayType2D,
    difl: *$DIFL/ArrayType1D,
    difr: *$DIFR/ArrayType2D,
    z: *$Z/ArrayType1D,
    $$k: int,
    $$c: $RT,
    $$s: RT,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           B.data_type == RT &&
           BX.data_type == RT &&
           PERM.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           POLES.data_type == RT &&
           DIFL.data_type == RT &&
           DIFR.data_type == RT &&
           Z.data_type == RT &&
           WORK.data_type == RT;
} {
    routine_name :: #run if RT == float32 { "SLALS0" } else { "DLALS0" };
    negone :: #run -one(RT);
    one    :: #run one(RT);
    zero   :: #run zero(RT);
    diflj, difrj, dj, dsigj, dsigjp, temp : RT = ---;
    i, j, m, n, nlp1 : int = ---;

    info.* = 0;
    n = nl + nr + 1;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if nl < 1 {
        info.* = -2;
    } else if nr < 1 {
        info.* = -3;
    } else if sqre < 0 || sqre > 1 {
        info.* = -4;
    } else if nrhs < 1 {
        info.* = -5;
    } else if ldb < n {
        info.* = -7;
    } else if ldbx < n {
        info.* = -9;
    } else if givptr < 0 {
        info.* = -11;
    } else if ldgcol < n {
        info.* = -13;
    } else if ldgnum < n {
        info.* = -15;
    } else if k < 1 {
        info.* = -20;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    m = n + sqre;
    nlp1 = nl + 1;

    if icompq == 0 {
        for i:1..givptr {
            row2 := _lals0_matrix_column(b.*, ldb, fget(givcol, i, 2), 1, nrhs);
            row1 := _lals0_matrix_column(b.*, ldb, fget(givcol, i, 1), 1, nrhs);
            _lals0_rot_real(nrhs, *row2, ldb, *row1, ldb, fget(givnum, i, 2), fget(givnum, i, 1), RT);
        }

        bx_first := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
        b_mid := _lals0_matrix_column(b.*, ldb, nlp1, 1, nrhs);
        _lals0_copy_real(nrhs, *b_mid, ldb, *bx_first, ldbx, RT);
        for i:2..n {
            b_src := _lals0_matrix_column(b.*, ldb, fget(perm, i), 1, nrhs);
            bx_dst := _lals0_matrix_column(bx.*, ldbx, i, 1, nrhs);
            _lals0_copy_real(nrhs, *b_src, ldb, *bx_dst, ldbx, RT);
        }

        if k == 1 {
            b_top := _lals0_matrix_column(b.*, ldb, 1, 1, nrhs);
            bx_top := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
            _lals0_copy_real(nrhs, *bx_top, ldbx, *b_top, ldb, RT);
            if fget(z, 1) < zero {
                _lals0_scal_real(nrhs, negone, *b_top, ldb, RT);
            }
        } else {
            for j:1..k {
                diflj = fget(difl, j);
                dj = fget(poles, j, 1);
                dsigj = -fget(poles, j, 2);
                if j < k {
                    difrj = -fget(difr, j, 1);
                    dsigjp = -fget(poles, j + 1, 2);
                }
                if fget(z, j) == zero || fget(poles, j, 2) == zero {
                    fset(work, j, zero);
                } else {
                    fset(work, j, -fget(poles, j, 2)*fget(z, j) / diflj / (fget(poles, j, 2) + dj));
                }
                for i:1..j - 1 {
                    if fget(z, i) == zero || fget(poles, i, 2) == zero {
                        fset(work, i, zero);
                    } else {
                        fset(work, i, fget(poles, i, 2)*fget(z, i) /
                            (_lals0_lamc3(fget(poles, i, 2), dsigj, RT) - diflj) /
                            (fget(poles, i, 2) + dj));
                    }
                }
                for i:j + 1..k {
                    if fget(z, i) == zero || fget(poles, i, 2) == zero {
                        fset(work, i, zero);
                    } else {
                        fset(work, i, fget(poles, i, 2)*fget(z, i) /
                            (_lals0_lamc3(fget(poles, i, 2), dsigjp, RT) + difrj) /
                            (fget(poles, i, 2) + dj));
                    }
                }
                fset(work, 1, negone);
                temp = _lals0_nrm2(k, work, 1, RT);
                bx_top_raw := _lals0_matrix_view(bx.*, ldbx, 1, 1, k, nrhs);
                bx_top := array2d_view(bx_top_raw, ldbx, row_major=false);
                b_row := _lals0_matrix_column(b.*, ldb, j, 1, nrhs);
                _lals0_gemv_real(#char "T", k, nrhs, one, *bx_top, ldbx, work, 1, zero, *b_row, ldb, RT);
                b_row_mat_raw := _lals0_matrix_view(b.*, ldb, j, 1, 1, nrhs);
                b_row_mat := array2d_view(b_row_mat_raw, ldb, row_major=false);
                _lals0_lascl_real(#char "G", 0, 0, temp, one, 1, nrhs, *b_row_mat, ldb, info, RT);
            }
        }

        if k < max(m, n) {
            bx_tail_raw := _lals0_matrix_view(bx.*, ldbx, k + 1, 1, n - k, nrhs);
            bx_tail := array2d_view(bx_tail_raw, ldbx, row_major=false);
            b_tail_raw := _lals0_matrix_view(b.*, ldb, k + 1, 1, n - k, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _lals0_lacpy_real(#char "A", n - k, nrhs, *bx_tail, ldbx, *b_tail, ldb, RT);
        }
    } else {
        if k == 1 {
            bx_top := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
            b_top := _lals0_matrix_column(b.*, ldb, 1, 1, nrhs);
            _lals0_copy_real(nrhs, *b_top, ldb, *bx_top, ldbx, RT);
        } else {
            for j:1..k {
                dsigj = fget(poles, j, 2);
                if fget(z, j) == zero {
                    fset(work, j, zero);
                } else {
                    fset(work, j, -fget(z, j) / fget(difl, j) / (dsigj + fget(poles, j, 1)) / fget(difr, j, 2));
                }
                for i:1..j - 1 {
                    if fget(z, j) == zero {
                        fset(work, i, zero);
                    } else {
                        fset(work, i, fget(z, j) /
                            (_lals0_lamc3(dsigj, -fget(poles, i + 1, 2), RT) - fget(difr, i, 1)) /
                            (dsigj + fget(poles, i, 1)) / fget(difr, i, 2));
                    }
                }
                for i:j + 1..k {
                    if fget(z, j) == zero {
                        fset(work, i, zero);
                    } else {
                        fset(work, i, fget(z, j) /
                            (_lals0_lamc3(dsigj, -fget(poles, i, 2), RT) - fget(difl, i)) /
                            (dsigj + fget(poles, i, 1)) / fget(difr, i, 2));
                    }
                }
                b_top_raw := _lals0_matrix_view(b.*, ldb, 1, 1, k, nrhs);
                b_top := array2d_view(b_top_raw, ldb, row_major=false);
                bx_row := _lals0_matrix_column(bx.*, ldbx, j, 1, nrhs);
                _lals0_gemv_real(#char "T", k, nrhs, one, *b_top, ldb, work, 1, zero, *bx_row, ldbx, RT);
            }
        }

        if sqre == 1 {
            bx_m := _lals0_matrix_column(bx.*, ldbx, m, 1, nrhs);
            b_m := _lals0_matrix_column(b.*, ldb, m, 1, nrhs);
            _lals0_copy_real(nrhs, *b_m, ldb, *bx_m, ldbx, RT);
            bx_1 := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
            _lals0_rot_real(nrhs, *bx_1, ldbx, *bx_m, ldbx, c, s, RT);
        }
        if k < max(m, n) {
            b_tail_raw := _lals0_matrix_view(b.*, ldb, k + 1, 1, n - k, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            bx_tail_raw := _lals0_matrix_view(bx.*, ldbx, k + 1, 1, n - k, nrhs);
            bx_tail := array2d_view(bx_tail_raw, ldbx, row_major=false);
            _lals0_lacpy_real(#char "A", n - k, nrhs, *b_tail, ldb, *bx_tail, ldbx, RT);
        }

        bx_1 := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
        b_nlp1 := _lals0_matrix_column(b.*, ldb, nlp1, 1, nrhs);
        _lals0_copy_real(nrhs, *bx_1, ldbx, *b_nlp1, ldb, RT);
        if sqre == 1 {
            bx_m := _lals0_matrix_column(bx.*, ldbx, m, 1, nrhs);
            b_m := _lals0_matrix_column(b.*, ldb, m, 1, nrhs);
            _lals0_copy_real(nrhs, *bx_m, ldbx, *b_m, ldb, RT);
        }
        for i:2..n {
            bx_src := _lals0_matrix_column(bx.*, ldbx, i, 1, nrhs);
            b_dst := _lals0_matrix_column(b.*, ldb, fget(perm, i), 1, nrhs);
            _lals0_copy_real(nrhs, *bx_src, ldbx, *b_dst, ldb, RT);
        }

        for <i:1..givptr {
            row2 := _lals0_matrix_column(b.*, ldb, fget(givcol, i, 2), 1, nrhs);
            row1 := _lals0_matrix_column(b.*, ldb, fget(givcol, i, 1), 1, nrhs);
            _lals0_rot_real(nrhs, *row2, ldb, *row1, ldb, fget(givnum, i, 2), -fget(givnum, i, 1), RT);
        }
    }
}

_lals0_complex :: (
    $$icompq: int,
    $$nl: int,
    $$nr: int,
    $$sqre: int,
    $$nrhs: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    bx: *$BX/ArrayType2D,
    $$ldbx: int,
    perm: *$PERM/ArrayType1D,
    $$givptr: int,
    givcol: *$GIVCOL/ArrayType2D,
    $$ldgcol: int,
    givnum: *$GIVNUM/ArrayType2D,
    $$ldgnum: int,
    poles: *$POLES/ArrayType2D,
    difl: *$DIFL/ArrayType1D,
    difr: *$DIFR/ArrayType2D,
    z: *$Z/ArrayType1D,
    $$k: int,
    $$c: $RT,
    $$s: RT,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           B.data_type == MT &&
           BX.data_type == MT &&
           PERM.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           POLES.data_type == RT &&
           DIFL.data_type == RT &&
           DIFR.data_type == RT &&
           Z.data_type == RT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CLALS0" } else { "ZLALS0" };
    negone :: #run -one(RT);
    one    :: #run one(RT);
    zero   :: #run zero(RT);
    diflj, difrj, dj, dsigj, dsigjp, temp : RT = ---;
    i, j, jcol, jrow, m, n, nlp1 : int = ---;

    info.* = 0;
    n = nl + nr + 1;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if nl < 1 {
        info.* = -2;
    } else if nr < 1 {
        info.* = -3;
    } else if sqre < 0 || sqre > 1 {
        info.* = -4;
    } else if nrhs < 1 {
        info.* = -5;
    } else if ldb < n {
        info.* = -7;
    } else if ldbx < n {
        info.* = -9;
    } else if givptr < 0 {
        info.* = -11;
    } else if ldgcol < n {
        info.* = -13;
    } else if ldgnum < n {
        info.* = -15;
    } else if k < 1 {
        info.* = -20;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    m = n + sqre;
    nlp1 = nl + 1;

    if icompq == 0 {
        for i:1..givptr {
            row2 := _lals0_matrix_column(b.*, ldb, fget(givcol, i, 2), 1, nrhs);
            row1 := _lals0_matrix_column(b.*, ldb, fget(givcol, i, 1), 1, nrhs);
            _lals0_rot_complex(nrhs, *row2, ldb, *row1, ldb, fget(givnum, i, 2), fget(givnum, i, 1), RT, MT);
        }

        bx_first := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
        b_mid := _lals0_matrix_column(b.*, ldb, nlp1, 1, nrhs);
        _lals0_copy_complex(nrhs, *b_mid, ldb, *bx_first, ldbx, MT);
        for i:2..n {
            b_src := _lals0_matrix_column(b.*, ldb, fget(perm, i), 1, nrhs);
            bx_dst := _lals0_matrix_column(bx.*, ldbx, i, 1, nrhs);
            _lals0_copy_complex(nrhs, *b_src, ldb, *bx_dst, ldbx, MT);
        }

        if k == 1 {
            bx_top := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
            b_top := _lals0_matrix_column(b.*, ldb, 1, 1, nrhs);
            _lals0_copy_complex(nrhs, *bx_top, ldbx, *b_top, ldb, MT);
            if fget(z, 1) < zero {
                _lals0_scal_complex_real(nrhs, negone, *b_top, ldb, RT, MT);
            }
        } else {
            for j:1..k {
                diflj = fget(difl, j);
                dj = fget(poles, j, 1);
                dsigj = -fget(poles, j, 2);
                if j < k {
                    difrj = -fget(difr, j, 1);
                    dsigjp = -fget(poles, j + 1, 2);
                }
                if fget(z, j) == zero || fget(poles, j, 2) == zero {
                    fset(rwork, j, zero);
                } else {
                    fset(rwork, j, -fget(poles, j, 2)*fget(z, j) / diflj / (fget(poles, j, 2) + dj));
                }
                for i:1..j - 1 {
                    if fget(z, i) == zero || fget(poles, i, 2) == zero {
                        fset(rwork, i, zero);
                    } else {
                        fset(rwork, i, fget(poles, i, 2)*fget(z, i) /
                            (_lals0_lamc3(fget(poles, i, 2), dsigj, RT) - diflj) /
                            (fget(poles, i, 2) + dj));
                    }
                }
                for i:j + 1..k {
                    if fget(z, i) == zero || fget(poles, i, 2) == zero {
                        fset(rwork, i, zero);
                    } else {
                        fset(rwork, i, fget(poles, i, 2)*fget(z, i) /
                            (_lals0_lamc3(fget(poles, i, 2), dsigjp, RT) + difrj) /
                            (fget(poles, i, 2) + dj));
                    }
                }
                fset(rwork, 1, negone);
                temp = _lals0_nrm2(k, rwork, 1, RT);

                i = k + nrhs*2;
                for jcol:1..nrhs {
                    for jrow:1..k {
                        i += 1;
                        fset(rwork, i, _lals0_real_from_complex(fget(bx, jrow, jcol), RT, MT));
                    }
                }
                real_mat_raw := _lals0_raw_slice(rwork.data + (k + nrhs*2), k*nrhs);
                real_mat := array2d_view(real_mat_raw, k, row_major=false);
                weights := _lals0_vector_slice(rwork, 1, k);
                out_real := _lals0_vector_slice(rwork, 1 + k, nrhs);
                _lals0_gemv_real(#char "T", k, nrhs, one, *real_mat, k, *weights, 1, zero, *out_real, 1, RT);

                i = k + nrhs*2;
                for jcol:1..nrhs {
                    for jrow:1..k {
                        i += 1;
                        fset(rwork, i, _lals0_imag_from_complex(fget(bx, jrow, jcol), RT, MT));
                    }
                }
                imag_mat_raw := _lals0_raw_slice(rwork.data + (k + nrhs*2), k*nrhs);
                imag_mat := array2d_view(imag_mat_raw, k, row_major=false);
                out_imag := _lals0_vector_slice(rwork, 1 + k + nrhs, nrhs);
                _lals0_gemv_real(#char "T", k, nrhs, one, *imag_mat, k, *weights, 1, zero, *out_imag, 1, RT);

                for jcol:1..nrhs {
                    fset(b, j, jcol, _lals0_complex_from_parts(fget(rwork, jcol + k), fget(rwork, jcol + k + nrhs), RT, MT));
                }
                b_row_mat_raw := _lals0_matrix_view(b.*, ldb, j, 1, 1, nrhs);
                b_row_mat := array2d_view(b_row_mat_raw, ldb, row_major=false);
                _lals0_lascl_complex(#char "G", 0, 0, temp, one, 1, nrhs, *b_row_mat, ldb, info, RT, MT);
            }
        }

        if k < max(m, n) {
            bx_tail_raw := _lals0_matrix_view(bx.*, ldbx, k + 1, 1, n - k, nrhs);
            bx_tail := array2d_view(bx_tail_raw, ldbx, row_major=false);
            b_tail_raw := _lals0_matrix_view(b.*, ldb, k + 1, 1, n - k, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _lals0_lacpy_complex(#char "A", n - k, nrhs, *bx_tail, ldbx, *b_tail, ldb, MT);
        }
    } else {
        if k == 1 {
            bx_top := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
            b_top := _lals0_matrix_column(b.*, ldb, 1, 1, nrhs);
            _lals0_copy_complex(nrhs, *b_top, ldb, *bx_top, ldbx, MT);
        } else {
            for j:1..k {
                dsigj = fget(poles, j, 2);
                if fget(z, j) == zero {
                    fset(rwork, j, zero);
                } else {
                    fset(rwork, j, -fget(z, j) / fget(difl, j) / (dsigj + fget(poles, j, 1)) / fget(difr, j, 2));
                }
                for i:1..j - 1 {
                    if fget(z, j) == zero {
                        fset(rwork, i, zero);
                    } else {
                        fset(rwork, i, fget(z, j) /
                            (_lals0_lamc3(dsigj, -fget(poles, i + 1, 2), RT) - fget(difr, i, 1)) /
                            (dsigj + fget(poles, i, 1)) / fget(difr, i, 2));
                    }
                }
                for i:j + 1..k {
                    if fget(z, j) == zero {
                        fset(rwork, i, zero);
                    } else {
                        fset(rwork, i, fget(z, j) /
                            (_lals0_lamc3(dsigj, -fget(poles, i, 2), RT) - fget(difl, i)) /
                            (dsigj + fget(poles, i, 1)) / fget(difr, i, 2));
                    }
                }

                i = k + nrhs*2;
                for jcol:1..nrhs {
                    for jrow:1..k {
                        i += 1;
                        fset(rwork, i, _lals0_real_from_complex(fget(b, jrow, jcol), RT, MT));
                    }
                }
                real_mat_raw := _lals0_raw_slice(rwork.data + (k + nrhs*2), k*nrhs);
                real_mat := array2d_view(real_mat_raw, k, row_major=false);
                weights := _lals0_vector_slice(rwork, 1, k);
                out_real := _lals0_vector_slice(rwork, 1 + k, nrhs);
                _lals0_gemv_real(#char "T", k, nrhs, one, *real_mat, k, *weights, 1, zero, *out_real, 1, RT);

                i = k + nrhs*2;
                for jcol:1..nrhs {
                    for jrow:1..k {
                        i += 1;
                        fset(rwork, i, _lals0_imag_from_complex(fget(b, jrow, jcol), RT, MT));
                    }
                }
                imag_mat_raw := _lals0_raw_slice(rwork.data + (k + nrhs*2), k*nrhs);
                imag_mat := array2d_view(imag_mat_raw, k, row_major=false);
                out_imag := _lals0_vector_slice(rwork, 1 + k + nrhs, nrhs);
                _lals0_gemv_real(#char "T", k, nrhs, one, *imag_mat, k, *weights, 1, zero, *out_imag, 1, RT);

                for jcol:1..nrhs {
                    fset(bx, j, jcol, _lals0_complex_from_parts(fget(rwork, jcol + k), fget(rwork, jcol + k + nrhs), RT, MT));
                }
            }
        }

        if sqre == 1 {
            bx_m := _lals0_matrix_column(bx.*, ldbx, m, 1, nrhs);
            b_m := _lals0_matrix_column(b.*, ldb, m, 1, nrhs);
            _lals0_copy_complex(nrhs, *b_m, ldb, *bx_m, ldbx, MT);
            bx_1 := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
            _lals0_rot_complex(nrhs, *bx_1, ldbx, *bx_m, ldbx, c, s, RT, MT);
        }
        if k < max(m, n) {
            b_tail_raw := _lals0_matrix_view(b.*, ldb, k + 1, 1, n - k, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            bx_tail_raw := _lals0_matrix_view(bx.*, ldbx, k + 1, 1, n - k, nrhs);
            bx_tail := array2d_view(bx_tail_raw, ldbx, row_major=false);
            _lals0_lacpy_complex(#char "A", n - k, nrhs, *b_tail, ldb, *bx_tail, ldbx, MT);
        }

        bx_1 := _lals0_matrix_column(bx.*, ldbx, 1, 1, nrhs);
        b_nlp1 := _lals0_matrix_column(b.*, ldb, nlp1, 1, nrhs);
        _lals0_copy_complex(nrhs, *bx_1, ldbx, *b_nlp1, ldb, MT);
        if sqre == 1 {
            bx_m := _lals0_matrix_column(bx.*, ldbx, m, 1, nrhs);
            b_m := _lals0_matrix_column(b.*, ldb, m, 1, nrhs);
            _lals0_copy_complex(nrhs, *bx_m, ldbx, *b_m, ldb, MT);
        }
        for i:2..n {
            bx_src := _lals0_matrix_column(bx.*, ldbx, i, 1, nrhs);
            b_dst := _lals0_matrix_column(b.*, ldb, fget(perm, i), 1, nrhs);
            _lals0_copy_complex(nrhs, *bx_src, ldbx, *b_dst, ldb, MT);
        }

        for <i:1..givptr {
            row2 := _lals0_matrix_column(b.*, ldb, fget(givcol, i, 2), 1, nrhs);
            row1 := _lals0_matrix_column(b.*, ldb, fget(givcol, i, 1), 1, nrhs);
            _lals0_rot_complex(nrhs, *row2, ldb, *row1, ldb, fget(givnum, i, 2), -fget(givnum, i, 1), RT, MT);
        }
    }
}
