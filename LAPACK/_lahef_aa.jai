// ? LAPACK *lahef_aa()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clahef_aa/zlahef_aa
// *LAHEF_AA factorizes a panel of a complex Hermitian matrix
// using Aasen's algorithm.

clahef_aa :: #bake_arguments _lahef_aa(MT=C32, RT=float32);
zlahef_aa :: #bake_arguments _lahef_aa(MT=C64, RT=float64);

#scope_file
_lahef_aa_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_lahef_aa_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _lahef_aa_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lahef_aa_vector_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _lahef_aa_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lahef_aa_work_view :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _lahef_aa_raw_slice(x.data + (offset - 1), count);
}

_lahef_aa_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_lahef_aa_copy :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { ccopy(n, x, incx, y, incy); } else { zcopy(n, x, incx, y, incy); }
}

_lahef_aa_axpy :: (
    $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { caxpy(n, alpha, x, incx, y, incy); } else { zaxpy(n, alpha, x, incx, y, incy); }
}

_lahef_aa_scal :: (
    $$n: int, alpha: $MT, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { cscal(n, alpha, x, incx); } else { zscal(n, alpha, x, incx); }
}

_lahef_aa_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_lahef_aa_gemv :: (
    $$trans: u8, $$m: int, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: $X/ArrayType1D, $$incx: int, beta: $MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }
    else { zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }
}

_lahef_aa_lacgv :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { clacgv(n, x, incx); } else { zlacgv(n, x, incx); }
}

_lahef_aa_amax :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, $MT: Type
) -> int #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { return icamax(n, x, incx) + 1; } else { return izamax(n, x, incx) + 1; }
}

_lahef_aa_laset :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: $MT, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 { claset(uplo, m, n, alpha, beta, a, lda); } else { zlaset(uplo, m, n, alpha, beta, a, lda); }
}

_lahef_aa :: (
    $$uplo: u8,
    $$j1: int,
    $$m: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    h: *$H/ArrayType2D,
    $$ldh: int,
    work: *$WORK/ArrayType1D,
    $MT: Type,
    $RT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           RT == inner_type_of_complex(MT) &&
           A.data_type == MT &&
           IPIV.data_type == int &&
           H.data_type == MT &&
           WORK.data_type == MT;
} {
    zero_t :: #run zero(MT);
    one_t  :: #run one(MT);

    alpha, piv : MT = ---;
    i1, i2, j, k, k1, mj, panel_n : int = ---;

    panel_n = min(m, nb);
    k1 = (2 - j1) + 1;
    j = 1;

    if lsame(uplo, #char "U") {
        while j <= panel_n {
            k = j1 + j - 1;
            mj = ifx j == m then 1 else m - j + 1;

            if k > 2 {
                avec := _lahef_aa_vector_view(a.*, lda, 1, j, j - k1);
                _lahef_aa_lacgv(j - k1, *avec, 1, MT);

                hblk_raw := _lahef_aa_matrix_storage(h.*, ldh, j, k1, mj, j - k1);
                hblk := array2d_view(hblk_raw, ldh, row_major=false);
                hcol := _lahef_aa_vector_view(h.*, ldh, j, j, mj);
                _lahef_aa_gemv(#char "N", mj, j - k1, -one_t, hblk.*, ldh, avec, 1, one_t, *hcol, 1, MT);

                _lahef_aa_lacgv(j - k1, *avec, 1, MT);
            }

            hcol_copy := _lahef_aa_vector_view(h.*, ldh, j, j, mj);
            workv := _lahef_aa_work_view(work, 1, mj);
            _lahef_aa_copy(mj, hcol_copy, 1, *workv, 1, MT);

            if j > k1 {
                alpha = -conjugate(fget(a.*, k - 1, j));
                ax := _lahef_aa_vector_view(a.*, lda, k - 2, j, mj);
                _lahef_aa_axpy(mj, alpha, ax, lda, *workv, 1, MT);
            }

            fset(a, k, j, _lahef_aa_complex_from_real(real(fget(work, 1)), RT, MT));

            if j < m {
                if k > 1 {
                    alpha = -fget(a.*, k, j);
                    ax := _lahef_aa_vector_view(a.*, lda, k - 1, j + 1, m - j);
                    work2 := _lahef_aa_work_view(work, 2, m - j);
                    _lahef_aa_axpy(m - j, alpha, ax, lda, *work2, 1, MT);
                }

                work2 := _lahef_aa_work_view(work, 2, m - j);
                i2 = _lahef_aa_amax(m - j, work2, 1, MT) + 1;
                piv = fget(work, i2);

                if i2 != 2 && piv != zero_t {
                    i1 = 2;
                    fset(work, i2, fget(work, i1));
                    fset(work, i1, piv);

                    i1 = i1 + j - 1;
                    i2 = i2 + j - 1;

                    x1 := _lahef_aa_vector_view(a.*, lda, j1 + i1 - 1, i1 + 1, i2 - i1 - 1);
                    x2 := _lahef_aa_vector_view(a.*, lda, j1 + i1, i2, i2 - i1 - 1);
                    _lahef_aa_swap(i2 - i1 - 1, *x1, lda, *x2, 1, MT);

                    x1c := _lahef_aa_vector_view(a.*, lda, j1 + i1 - 1, i1 + 1, i2 - i1);
                    _lahef_aa_lacgv(i2 - i1, *x1c, lda, MT);
                    x2c := _lahef_aa_vector_view(a.*, lda, j1 + i1, i2, i2 - i1 - 1);
                    _lahef_aa_lacgv(i2 - i1 - 1, *x2c, 1, MT);

                    if i2 < m {
                        y1 := _lahef_aa_vector_view(a.*, lda, j1 + i1 - 1, i2 + 1, m - i2);
                        y2 := _lahef_aa_vector_view(a.*, lda, j1 + i2 - 1, i2 + 1, m - i2);
                        _lahef_aa_swap(m - i2, *y1, lda, *y2, lda, MT);
                    }

                    piv = fget(a.*, i1 + j1 - 1, i1);
                    fset(a, i1 + j1 - 1, i1, fget(a.*, i1 + j1 - 1 + (i2 - i1), i2));
                    fset(a, i2 + j1 - 1, i2, piv);

                    hx1 := _lahef_aa_vector_view(h.*, ldh, i1, 1, i1 - 1);
                    hx2 := _lahef_aa_vector_view(h.*, ldh, i2, 1, i1 - 1);
                    _lahef_aa_swap(i1 - 1, *hx1, ldh, *hx2, ldh, MT);
                    fset(ipiv, i1, i2);

                    if i1 > k1 - 1 {
                        lx1 := _lahef_aa_vector_view(a.*, lda, 1, i1, i1 - k1 + 1);
                        lx2 := _lahef_aa_vector_view(a.*, lda, 1, i2, i1 - k1 + 1);
                        _lahef_aa_swap(i1 - k1 + 1, *lx1, 1, *lx2, 1, MT);
                    }
                } else {
                    fset(ipiv, j + 1, j + 1);
                }

                fset(a, k, j + 1, fget(work, 2));

                if j < nb {
                    src := _lahef_aa_vector_view(a.*, lda, k + 1, j + 1, m - j);
                    dst := _lahef_aa_vector_view(h.*, ldh, j + 1, j + 1, m - j);
                    _lahef_aa_copy(m - j, src, lda, *dst, 1, MT);
                }

                if j < m - 1 {
                    if fget(a.*, k, j + 1) != zero_t {
                        alpha = one_t / fget(a.*, k, j + 1);
                        src := _lahef_aa_work_view(work, 3, m - j - 1);
                        dst := _lahef_aa_vector_view(a.*, lda, k, j + 2, m - j - 1);
                        _lahef_aa_copy(m - j - 1, src, 1, *dst, lda, MT);
                        _lahef_aa_scal(m - j - 1, alpha, *dst, lda, MT);
                    } else {
                        block_raw := _lahef_aa_matrix_storage(a.*, lda, k, j + 2, 1, m - j - 1);
                        block := array2d_view(block_raw, lda, row_major=false);
                        _lahef_aa_laset(#char "F", 1, m - j - 1, zero_t, zero_t, *block, lda, MT);
                    }
                }
            }

            j += 1;
        }
    } else {
        while j <= panel_n {
            k = j1 + j - 1;
            mj = ifx j == m then 1 else m - j + 1;

            if k > 2 {
                avec := _lahef_aa_vector_view(a.*, lda, j, 1, j - k1);
                _lahef_aa_lacgv(j - k1, *avec, lda, MT);

                hblk_raw := _lahef_aa_matrix_storage(h.*, ldh, j, k1, mj, j - k1);
                hblk := array2d_view(hblk_raw, ldh, row_major=false);
                hcol := _lahef_aa_vector_view(h.*, ldh, j, j, mj);
                _lahef_aa_gemv(#char "N", mj, j - k1, -one_t, hblk.*, ldh, avec, lda, one_t, *hcol, 1, MT);

                _lahef_aa_lacgv(j - k1, *avec, lda, MT);
            }

            hcol_copy := _lahef_aa_vector_view(h.*, ldh, j, j, mj);
            workv := _lahef_aa_work_view(work, 1, mj);
            _lahef_aa_copy(mj, hcol_copy, 1, *workv, 1, MT);

            if j > k1 {
                alpha = -conjugate(fget(a.*, j, k - 1));
                ax := _lahef_aa_vector_view(a.*, lda, j, k - 2, mj);
                _lahef_aa_axpy(mj, alpha, ax, 1, *workv, 1, MT);
            }

            fset(a, j, k, _lahef_aa_complex_from_real(real(fget(work, 1)), RT, MT));

            if j < m {
                if k > 1 {
                    alpha = -fget(a.*, j, k);
                    ax := _lahef_aa_vector_view(a.*, lda, j + 1, k - 1, m - j);
                    work2 := _lahef_aa_work_view(work, 2, m - j);
                    _lahef_aa_axpy(m - j, alpha, ax, 1, *work2, 1, MT);
                }

                work2 := _lahef_aa_work_view(work, 2, m - j);
                i2 = _lahef_aa_amax(m - j, work2, 1, MT) + 1;
                piv = fget(work, i2);

                if i2 != 2 && piv != zero_t {
                    i1 = 2;
                    fset(work, i2, fget(work, i1));
                    fset(work, i1, piv);

                    i1 = i1 + j - 1;
                    i2 = i2 + j - 1;

                    x1 := _lahef_aa_vector_view(a.*, lda, i1 + 1, j1 + i1 - 1, i2 - i1 - 1);
                    x2 := _lahef_aa_vector_view(a.*, lda, i2, j1 + i1, i2 - i1 - 1);
                    _lahef_aa_swap(i2 - i1 - 1, *x1, 1, *x2, lda, MT);

                    x1c := _lahef_aa_vector_view(a.*, lda, i1 + 1, j1 + i1 - 1, i2 - i1);
                    _lahef_aa_lacgv(i2 - i1, *x1c, 1, MT);
                    x2c := _lahef_aa_vector_view(a.*, lda, i2, j1 + i1, i2 - i1 - 1);
                    _lahef_aa_lacgv(i2 - i1 - 1, *x2c, lda, MT);

                    if i2 < m {
                        y1 := _lahef_aa_vector_view(a.*, lda, i2 + 1, j1 + i1 - 1, m - i2);
                        y2 := _lahef_aa_vector_view(a.*, lda, i2 + 1, j1 + i2 - 1, m - i2);
                        _lahef_aa_swap(m - i2, *y1, 1, *y2, 1, MT);
                    }

                    piv = fget(a.*, i1, j1 + i1 - 1);
                    fset(a, i1, j1 + i1 - 1, fget(a.*, i2, j1 + i2 - 1));
                    fset(a, i2, j1 + i2 - 1, piv);

                    hx1 := _lahef_aa_vector_view(h.*, ldh, i1, 1, i1 - 1);
                    hx2 := _lahef_aa_vector_view(h.*, ldh, i2, 1, i1 - 1);
                    _lahef_aa_swap(i1 - 1, *hx1, ldh, *hx2, ldh, MT);
                    fset(ipiv, i1, i2);

                    if i1 > k1 - 1 {
                        lx1 := _lahef_aa_vector_view(a.*, lda, i1, 1, i1 - k1 + 1);
                        lx2 := _lahef_aa_vector_view(a.*, lda, i2, 1, i1 - k1 + 1);
                        _lahef_aa_swap(i1 - k1 + 1, *lx1, lda, *lx2, lda, MT);
                    }
                } else {
                    fset(ipiv, j + 1, j + 1);
                }

                fset(a, j + 1, k, fget(work, 2));

                if j < nb {
                    src := _lahef_aa_vector_view(a.*, lda, j + 1, k + 1, m - j);
                    dst := _lahef_aa_vector_view(h.*, ldh, j + 1, j + 1, m - j);
                    _lahef_aa_copy(m - j, src, 1, *dst, 1, MT);
                }

                if j < m - 1 {
                    if fget(a.*, j + 1, k) != zero_t {
                        alpha = one_t / fget(a.*, j + 1, k);
                        src := _lahef_aa_work_view(work, 3, m - j - 1);
                        dst := _lahef_aa_vector_view(a.*, lda, j + 2, k, m - j - 1);
                        _lahef_aa_copy(m - j - 1, src, 1, *dst, 1, MT);
                        _lahef_aa_scal(m - j - 1, alpha, *dst, 1, MT);
                    } else {
                        block_raw := _lahef_aa_matrix_storage(a.*, lda, j + 2, k, m - j - 1, 1);
                        block := array2d_view(block_raw, lda, row_major=false);
                        _lahef_aa_laset(#char "F", m - j - 1, 1, zero_t, zero_t, *block, lda, MT);
                    }
                }
            }

            j += 1;
        }
    }
}
