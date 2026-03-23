// ? LAPACK *hetri()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetri/zhetri
// *HETRI computes the inverse of a complex Hermitian indefinite matrix
// using the factorization computed by *HETRF.

chetri :: #bake_arguments _hetri(RT=float32, MT=C32);
zhetri :: #bake_arguments _hetri(RT=float64, MT=C64);

#scope_file
_hetri_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetri_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetri_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetri_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _hetri_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetri_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetri_raw_slice(x.data + (offset - 1), count);
}

_hetri_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetri_copy :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { ccopy(n, x, incx, y, incy); } else { zcopy(n, x, incx, y, incy); }
}

_hetri_hemv :: (
    $$uplo: u8, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: $X/ArrayType1D, $$incx: int, beta: $MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy); } else { zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy); }
}

_hetri_dotc :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $MT: Type
) -> MT #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { return cdotc(n, x, incx, y, incy); } else { return zdotc(n, x, incx, y, incy); }
}

_hetri_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hetri :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    one_r  :: #run one(RT);
    one_t  :: #run _hetri_complex_from_real(one_r, RT, MT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == C32 { "CHETRI" } else { "ZHETRI" };

    upper : bool = ---;
    j, k, kp, kstep : int = ---;
    ak, akp1, d, t : RT = ---;
    akkp1, temp : MT = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    if upper {
        for <info.*:1..n {
            if fget(ipiv, info.*) > 0 && fget(a, info.*, info.*) == zero_t {
                return;
            }
        }
    } else {
        for info.*:1..n {
            if fget(ipiv, info.*) > 0 && fget(a, info.*, info.*) == zero_t {
                return;
            }
        }
    }
    info.* = 0;

    if upper {
        k = 1;
        while true {
            if k > n then break;

            if fget(ipiv, k) > 0 {
                fset(a, k, k, _hetri_complex_from_real(one_r / real(fget(a.*, k, k)), RT, MT));
                if k > 1 {
                    work_head := _hetri_vector_slice(work, 1, k - 1);
                    col_k := _hetri_matrix_column(a.*, lda, 1, k, k - 1);
                    _hetri_copy(k - 1, col_k, 1, *work_head, 1, MT);
                    a11_raw := _hetri_matrix_storage(a.*, lda, 1, 1, k - 1, k - 1);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    _hetri_hemv(uplo, k - 1, -one_t, *a11, lda, work_head, 1, zero_t, *col_k, 1, MT);
                    fset(a, k, k, fget(a.*, k, k) - _hetri_complex_from_real(real(_hetri_dotc(k - 1, work_head, 1, col_k, 1, MT)), RT, MT));
                }
                kstep = 1;
            } else {
                t = abs(fget(a.*, k, k + 1));
                ak = real(fget(a.*, k, k)) / t;
                akp1 = real(fget(a.*, k + 1, k + 1)) / t;
                akkp1 = fget(a.*, k, k + 1) / t;
                d = t*(ak*akp1 - one_r);
                fset(a, k, k, _hetri_complex_from_real(akp1 / d, RT, MT));
                fset(a, k + 1, k + 1, _hetri_complex_from_real(ak / d, RT, MT));
                fset(a, k, k + 1, -akkp1 / d);
                if k > 1 {
                    work_head := _hetri_vector_slice(work, 1, k - 1);
                    col_k := _hetri_matrix_column(a.*, lda, 1, k, k - 1);
                    col_kp1 := _hetri_matrix_column(a.*, lda, 1, k + 1, k - 1);
                    a11_raw := _hetri_matrix_storage(a.*, lda, 1, 1, k - 1, k - 1);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    _hetri_copy(k - 1, col_k, 1, *work_head, 1, MT);
                    _hetri_hemv(uplo, k - 1, -one_t, *a11, lda, work_head, 1, zero_t, *col_k, 1, MT);
                    fset(a, k, k, fget(a.*, k, k) - _hetri_complex_from_real(real(_hetri_dotc(k - 1, work_head, 1, col_k, 1, MT)), RT, MT));
                    fset(a, k, k + 1, fget(a.*, k, k + 1) - _hetri_dotc(k - 1, col_k, 1, col_kp1, 1, MT));
                    _hetri_copy(k - 1, col_kp1, 1, *work_head, 1, MT);
                    _hetri_hemv(uplo, k - 1, -one_t, *a11, lda, work_head, 1, zero_t, *col_kp1, 1, MT);
                    fset(a, k + 1, k + 1, fget(a.*, k + 1, k + 1) - _hetri_complex_from_real(real(_hetri_dotc(k - 1, work_head, 1, col_kp1, 1, MT)), RT, MT));
                }
                kstep = 2;
            }

            kp = abs(fget(ipiv, k));
            if kp != k {
                if kp > 1 {
                    col_k_pre := _hetri_matrix_column(a.*, lda, 1, k, kp - 1);
                    col_kp_pre := _hetri_matrix_column(a.*, lda, 1, kp, kp - 1);
                    _hetri_swap(kp - 1, *col_k_pre, 1, *col_kp_pre, 1, MT);
                }
                for j:kp + 1..k - 1 {
                    temp = conjugate(fget(a.*, j, k));
                    fset(a, j, k, conjugate(fget(a.*, kp, j)));
                    fset(a, kp, j, temp);
                }
                fset(a, kp, k, conjugate(fget(a.*, kp, k)));
                temp = fget(a.*, k, k);
                fset(a, k, k, fget(a.*, kp, kp));
                fset(a, kp, kp, temp);
                if kstep == 2 {
                    temp = fget(a.*, k, k + 1);
                    fset(a, k, k + 1, fget(a.*, kp, k + 1));
                    fset(a, kp, k + 1, temp);
                }
            }

            k += kstep;
        }
    } else {
        k = n;
        while true {
            if k < 1 then break;

            if fget(ipiv, k) > 0 {
                fset(a, k, k, _hetri_complex_from_real(one_r / real(fget(a.*, k, k)), RT, MT));
                if k < n {
                    work_tail := _hetri_vector_slice(work, 1, n - k);
                    col_k_tail := _hetri_matrix_column(a.*, lda, k + 1, k, n - k);
                    _hetri_copy(n - k, col_k_tail, 1, *work_tail, 1, MT);
                    a22_raw := _hetri_matrix_storage(a.*, lda, k + 1, k + 1, n - k, n - k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    _hetri_hemv(uplo, n - k, -one_t, *a22, lda, work_tail, 1, zero_t, *col_k_tail, 1, MT);
                    fset(a, k, k, fget(a.*, k, k) - _hetri_complex_from_real(real(_hetri_dotc(n - k, work_tail, 1, col_k_tail, 1, MT)), RT, MT));
                }
                kstep = 1;
            } else {
                t = abs(fget(a.*, k, k - 1));
                ak = real(fget(a.*, k - 1, k - 1)) / t;
                akp1 = real(fget(a.*, k, k)) / t;
                akkp1 = fget(a.*, k, k - 1) / t;
                d = t*(ak*akp1 - one_r);
                fset(a, k - 1, k - 1, _hetri_complex_from_real(akp1 / d, RT, MT));
                fset(a, k, k, _hetri_complex_from_real(ak / d, RT, MT));
                fset(a, k, k - 1, -akkp1 / d);
                if k < n {
                    work_tail := _hetri_vector_slice(work, 1, n - k);
                    col_k_tail := _hetri_matrix_column(a.*, lda, k + 1, k, n - k);
                    col_km1_tail := _hetri_matrix_column(a.*, lda, k + 1, k - 1, n - k);
                    a22_raw := _hetri_matrix_storage(a.*, lda, k + 1, k + 1, n - k, n - k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    _hetri_copy(n - k, col_k_tail, 1, *work_tail, 1, MT);
                    _hetri_hemv(uplo, n - k, -one_t, *a22, lda, work_tail, 1, zero_t, *col_k_tail, 1, MT);
                    fset(a, k, k, fget(a.*, k, k) - _hetri_complex_from_real(real(_hetri_dotc(n - k, work_tail, 1, col_k_tail, 1, MT)), RT, MT));
                    fset(a, k, k - 1, fget(a.*, k, k - 1) - _hetri_dotc(n - k, col_k_tail, 1, col_km1_tail, 1, MT));
                    _hetri_copy(n - k, col_km1_tail, 1, *work_tail, 1, MT);
                    _hetri_hemv(uplo, n - k, -one_t, *a22, lda, work_tail, 1, zero_t, *col_km1_tail, 1, MT);
                    fset(a, k - 1, k - 1, fget(a.*, k - 1, k - 1) - _hetri_complex_from_real(real(_hetri_dotc(n - k, work_tail, 1, col_km1_tail, 1, MT)), RT, MT));
                }
                kstep = 2;
            }

            kp = abs(fget(ipiv, k));
            if kp != k {
                if kp < n {
                    col_k_tail := _hetri_matrix_column(a.*, lda, kp + 1, k, n - kp);
                    col_kp_tail := _hetri_matrix_column(a.*, lda, kp + 1, kp, n - kp);
                    _hetri_swap(n - kp, *col_k_tail, 1, *col_kp_tail, 1, MT);
                }
                for j:k + 1..kp - 1 {
                    temp = conjugate(fget(a.*, j, k));
                    fset(a, j, k, conjugate(fget(a.*, kp, j)));
                    fset(a, kp, j, temp);
                }
                fset(a, kp, k, conjugate(fget(a.*, kp, k)));
                temp = fget(a.*, k, k);
                fset(a, k, k, fget(a.*, kp, kp));
                fset(a, kp, kp, temp);
                if kstep == 2 {
                    temp = fget(a.*, k, k - 1);
                    fset(a, k, k - 1, fget(a.*, kp, k - 1));
                    fset(a, kp, k - 1, temp);
                }
            }

            k -= kstep;
        }
    }
}
