// ? LAPACK *getri()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgetri/cgetri/dgetri/zgetri
// *GETRI computes the inverse of a matrix using the LU factorization
// computed by *GETRF.

sgetri :: #bake_arguments _getri(RT=float32, MT=float32);
dgetri :: #bake_arguments _getri(RT=float64, MT=float64);
cgetri :: #bake_arguments _getri(RT=float32, MT=C32);
zgetri :: #bake_arguments _getri(RT=float64, MT=C64);

#scope_file
_getri_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda + (row - 1);
    view.count = count;
    return view;
}

_getri_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_getri_work_storage :: inline (
    work: $WORK/ArrayType1D, $$ldwork: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] WORK.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*ldwork + m;
    view: [] WORK.data_type;
    view.data = work.data + (col - 1)*ldwork + (row - 1);
    view.count = count;
    return view;
}

_getri_query_value :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == RT {
        #if RT == float32 {
            return sroundup_lwork(x);
        } else {
            return x;
        }
    } else {
        #if RT == float32 {
            return cpx(sroundup_lwork(x), zero(RT));
        } else {
            return cpx(x, zero(RT));
        }
    }
}

_getri_trtri :: (
    $$uplo: u8,
    $$diag: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        strtri(uplo, diag, n, a, lda, info);
    } else if MT == float64 {
        dtrtri(uplo, diag, n, a, lda, info);
    } else if MT == C32 {
        ctrtri(uplo, diag, n, a, lda, info);
    } else {
        ztrtri(uplo, diag, n, a, lda, info);
    }
}

_getri_gemv :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    $$incx: int,
    beta: $MT,
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

_getri_gemm :: (
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

_getri_trsm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
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
        strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == float64 {
        dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_getri_swap :: (
    $$n: int,
    x: *$X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sswap(n, x, 1, y, 1);
    } else if MT == float64 {
        dswap(n, x, 1, y, 1);
    } else if MT == C32 {
        cswap(n, x, 1, y, 1);
    } else {
        zswap(n, x, 1, y, 1);
    }
}

_getri :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: $IPIV/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           IPIV.data_type == int &&
           WORK.data_type == MT;
} {
    one_t  :: #run ifx MT == RT then ncast(MT, one(RT)) else cpx(one(RT), zero(RT));
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SGETRI"
    } else if MT == float64 {
        "DGETRI"
    } else if MT == C32 {
        "CGETRI"
    } else {
        "ZGETRI"
    };

    i, iws, j, jb, jj, jp, ldwork, lwkopt, nb, nbmin, nn : int = ---;
    lquery := lwork == -1;

    info.* = 0;
    nb = ilaenv(1, routine_name, #char " ", n, -1, -1, -1);
    lwkopt = max(1, n*nb);
    fset(work, 1, _getri_query_value(ncast(RT, lwkopt), RT, MT));

    if n < 0 {
        info.* = -1;
    } else if lda < max(1, n) {
        info.* = -3;
    } else if lwork < max(1, n) && !lquery {
        info.* = -6;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;

    _getri_trtri(#char "U", #char "N", n, a, lda, info, MT);
    if info.* > 0 then return;

    nbmin = 2;
    ldwork = n;
    if nb > 1 && nb < n {
        iws = max(ldwork*nb, 1);
        if lwork < iws {
            nb = lwork / ldwork;
            nbmin = max(2, ilaenv(2, routine_name, #char " ", n, -1, -1, -1));
        }
    } else {
        iws = n;
    }

    if nb < nbmin || nb >= n {
        for <j:1..n {
            for i:j + 1..n {
                fset(work, i, fget(a.*, i, j));
                fset(a, i, j, zero_t);
            }

            if j < n {
                a_right_raw := _getri_matrix_storage(a.*, lda, 1, j + 1, n, n - j);
                a_right_view := array2d_view(a_right_raw, lda, row_major=false);
                work_view := _getri_subvector(work.*, j + 1, n - j);
                y_view: [] MT;
                y_view.data = a.data + (j - 1)*lda;
                y_view.count = n;
                _getri_gemv(#char "N", n, n - j, -one_t, *a_right_view, lda, work_view, 1, one_t, *y_view, 1, MT);
            }
        }
    } else {
        nn = ((n - 1) / nb)*nb + 1;
        j = nn;
        while j >= 1 {
            jb = min(nb, n - j + 1);

            for jj:j..j + jb - 1 {
                for i:jj + 1..n {
                    fset(work, i + (jj - j)*ldwork, fget(a.*, i, jj));
                    fset(a, i, jj, zero_t);
                }
            }

            if j + jb <= n {
                a_right_raw := _getri_matrix_storage(a.*, lda, 1, j + jb, n, n - j - jb + 1);
                a_right_view := array2d_view(a_right_raw, lda, row_major=false);
                work_block_raw := _getri_work_storage(work.*, ldwork, j + jb, 1, n - j - jb + 1, jb);
                work_block_view := array2d_view(work_block_raw, ldwork, row_major=false);
                a_block_raw := _getri_matrix_storage(a.*, lda, 1, j, n, jb);
                a_block_view := array2d_view(a_block_raw, lda, row_major=false);
                _getri_gemm(#char "N", #char "N", n, jb, n - j - jb + 1, -one_t, *a_right_view, lda, *work_block_view, ldwork, one_t, *a_block_view, lda, MT);
            }

            work_j_raw := _getri_work_storage(work.*, ldwork, j, 1, n, jb);
            work_j_view := array2d_view(work_j_raw, ldwork, row_major=false);
            a_j_raw := _getri_matrix_storage(a.*, lda, 1, j, n, jb);
            a_j_view := array2d_view(a_j_raw, lda, row_major=false);
            _getri_trsm(#char "R", #char "L", #char "N", #char "U", n, jb, one_t, *work_j_view, ldwork, *a_j_view, lda, MT);
            j -= nb;
        }
    }

    for <j:1..n - 1 {
        jp = fget(ipiv, j);
        if jp != j {
            col_j: [] MT;
            col_j.data = a.data + (j - 1)*lda;
            col_j.count = n;
            col_jp: [] MT;
            col_jp.data = a.data + (jp - 1)*lda;
            col_jp.count = n;
            _getri_swap(n, *col_j, *col_jp, MT);
        }
    }

    fset(work, 1, _getri_query_value(ncast(RT, iws), RT, MT));
}
