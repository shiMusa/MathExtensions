// ? LAPACK *hegs2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chegs2/zhegs2
// *HEGS2 reduces a complex Hermitian-definite generalized eigenproblem
// to standard form using an unblocked algorithm.

chegs2 :: #bake_arguments _hegs2(RT=float32, MT=C32);
zhegs2 :: #bake_arguments _hegs2(RT=float64, MT=C64);

#scope_file
_hegs2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hegs2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hegs2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hegs2_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _hegs2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hegs2_matrix_strided :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int,
    $$inc: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*inc;
    return _hegs2_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_hegs2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hegs2_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        caxpy(n, alpha, x, incx, y, incy);
    } else {
        zaxpy(n, alpha, x, incx, y, incy);
    }
}

_hegs2_her2 :: (
    $$uplo: u8,
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT &&
           A.data_type == MT;
} {
    #if MT == C32 {
        cher2(uplo, n, alpha, x, incx, y, incy, a, lda);
    } else {
        zher2(uplo, n, alpha, x, incx, y, incy, a, lda);
    }
}

_hegs2_lacgv :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        clacgv(n, x, incx);
    } else {
        zlacgv(n, x, incx);
    }
}

_hegs2_scal_real :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_hegs2_trmv :: (
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
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT;
} {
    #if MT == C32 {
        ctrmv(uplo, trans, diag, n, a, lda, x, incx);
    } else {
        ztrmv(uplo, trans, diag, n, a, lda, x, incx);
    }
}

_hegs2_trsv :: (
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
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT;
} {
    #if MT == C32 {
        ctrsv(uplo, trans, diag, n, a, lda, x, incx);
    } else {
        ztrsv(uplo, trans, diag, n, a, lda, x, incx);
    }
}

_hegs2 :: (
    $$itype: int,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    one_r  :: #run one(RT);
    half_r :: #run cast(RT) 0.5;
    cone_t :: #run _hegs2_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHEGS2" } else { "ZHEGS2" };

    upper : bool = ---;
    k : int = ---;
    akk, bkk : RT = ---;
    ct : MT = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldb < max(1, n) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if itype == 1 {
        if upper {
            for k:1..n {
                akk = real(fget(a.*, k, k));
                bkk = real(fget(b.*, k, k));
                akk /= bkk*bkk;
                fset(a, k, k, _hegs2_complex_from_real(akk, RT, MT));
                if k < n {
                    ntrail := n - k;
                    arow := _hegs2_matrix_strided(a.*, lda, k, k + 1, ntrail, lda);
                    brow := _hegs2_matrix_strided(b.*, ldb, k, k + 1, ntrail, ldb);
                    a22_raw := _hegs2_matrix_storage(a.*, lda, k + 1, k + 1, ntrail, ntrail);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    b22_raw := _hegs2_matrix_storage(b.*, ldb, k + 1, k + 1, ntrail, ntrail);
                    b22 := array2d_view(b22_raw, ldb, row_major=false);

                    _hegs2_scal_real(ntrail, one_r / bkk, *arow, lda, RT, MT);
                    ct = _hegs2_complex_from_real(-half_r*akk, RT, MT);
                    _hegs2_lacgv(ntrail, *arow, lda, MT);
                    _hegs2_lacgv(ntrail, *brow, ldb, MT);
                    _hegs2_axpy(ntrail, ct, brow, ldb, *arow, lda, MT);
                    _hegs2_her2(uplo, ntrail, -cone_t, arow, lda, brow, ldb, *a22, lda, MT);
                    _hegs2_axpy(ntrail, ct, brow, ldb, *arow, lda, MT);
                    _hegs2_lacgv(ntrail, *brow, ldb, MT);
                    _hegs2_trsv(uplo, #char "C", #char "N", ntrail, b22, ldb, *arow, lda, MT);
                    _hegs2_lacgv(ntrail, *arow, lda, MT);
                }
            }
        } else {
            for k:1..n {
                akk = real(fget(a.*, k, k));
                bkk = real(fget(b.*, k, k));
                akk /= bkk*bkk;
                fset(a, k, k, _hegs2_complex_from_real(akk, RT, MT));
                if k < n {
                    ntrail := n - k;
                    acol := _hegs2_matrix_column(a.*, lda, k + 1, k, ntrail);
                    bcol := _hegs2_matrix_column(b.*, ldb, k + 1, k, ntrail);
                    a22_raw := _hegs2_matrix_storage(a.*, lda, k + 1, k + 1, ntrail, ntrail);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    b22_raw := _hegs2_matrix_storage(b.*, ldb, k + 1, k + 1, ntrail, ntrail);
                    b22 := array2d_view(b22_raw, ldb, row_major=false);

                    _hegs2_scal_real(ntrail, one_r / bkk, *acol, 1, RT, MT);
                    ct = _hegs2_complex_from_real(-half_r*akk, RT, MT);
                    _hegs2_axpy(ntrail, ct, bcol, 1, *acol, 1, MT);
                    _hegs2_her2(uplo, ntrail, -cone_t, acol, 1, bcol, 1, *a22, lda, MT);
                    _hegs2_axpy(ntrail, ct, bcol, 1, *acol, 1, MT);
                    _hegs2_trsv(uplo, #char "N", #char "N", ntrail, b22, ldb, *acol, 1, MT);
                }
            }
        }
    } else {
        if upper {
            for k:1..n {
                akk = real(fget(a.*, k, k));
                bkk = real(fget(b.*, k, k));
                if k > 1 {
                    kprev := k - 1;
                    acol := _hegs2_matrix_column(a.*, lda, 1, k, kprev);
                    bcol := _hegs2_matrix_column(b.*, ldb, 1, k, kprev);
                    a11_raw := _hegs2_matrix_storage(a.*, lda, 1, 1, kprev, kprev);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    b11_raw := _hegs2_matrix_storage(b.*, ldb, 1, 1, kprev, kprev);
                    b11 := array2d_view(b11_raw, ldb, row_major=false);

                    _hegs2_trmv(uplo, #char "N", #char "N", kprev, b11, ldb, *acol, 1, MT);
                    ct = _hegs2_complex_from_real(half_r*akk, RT, MT);
                    _hegs2_axpy(kprev, ct, bcol, 1, *acol, 1, MT);
                    _hegs2_her2(uplo, kprev, cone_t, acol, 1, bcol, 1, *a11, lda, MT);
                    _hegs2_axpy(kprev, ct, bcol, 1, *acol, 1, MT);
                    _hegs2_scal_real(kprev, bkk, *acol, 1, RT, MT);
                }
                fset(a, k, k, _hegs2_complex_from_real(akk*bkk*bkk, RT, MT));
            }
        } else {
            for k:1..n {
                akk = real(fget(a.*, k, k));
                bkk = real(fget(b.*, k, k));
                if k > 1 {
                    kprev := k - 1;
                    arow := _hegs2_matrix_strided(a.*, lda, k, 1, kprev, lda);
                    brow := _hegs2_matrix_strided(b.*, ldb, k, 1, kprev, ldb);
                    a11_raw := _hegs2_matrix_storage(a.*, lda, 1, 1, kprev, kprev);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    b11_raw := _hegs2_matrix_storage(b.*, ldb, 1, 1, kprev, kprev);
                    b11 := array2d_view(b11_raw, ldb, row_major=false);

                    _hegs2_lacgv(kprev, *arow, lda, MT);
                    _hegs2_trmv(uplo, #char "C", #char "N", kprev, b11, ldb, *arow, lda, MT);
                    ct = _hegs2_complex_from_real(half_r*akk, RT, MT);
                    _hegs2_lacgv(kprev, *brow, ldb, MT);
                    _hegs2_axpy(kprev, ct, brow, ldb, *arow, lda, MT);
                    _hegs2_her2(uplo, kprev, cone_t, arow, lda, brow, ldb, *a11, lda, MT);
                    _hegs2_axpy(kprev, ct, brow, ldb, *arow, lda, MT);
                    _hegs2_lacgv(kprev, *brow, ldb, MT);
                    _hegs2_scal_real(kprev, bkk, *arow, lda, RT, MT);
                    _hegs2_lacgv(kprev, *arow, lda, MT);
                }
                fset(a, k, k, _hegs2_complex_from_real(akk*bkk*bkk, RT, MT));
            }
        }
    }
}
