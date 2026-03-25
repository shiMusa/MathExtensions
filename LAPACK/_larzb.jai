// ? LAPACK *larzb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarzb/dlarzb/clarzb/zlarzb
// *LARZB applies a block reflector to a matrix.

slarzb :: #bake_arguments _larzb(RT=float32, MT=float32);
dlarzb :: #bake_arguments _larzb(RT=float64, MT=float64);
clarzb :: #bake_arguments _larzb(RT=float32, MT=C32);
zlarzb :: #bake_arguments _larzb(RT=float64, MT=C64);

#scope_file
_larzb_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_larzb_matrix_storage :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larzb_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larzb_matrix_storage_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larzb_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larzb_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int, alpha: $MT,
    a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int, beta: $MT,
    c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT && C.data_type == MT;
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

_larzb_trmm :: (
    $$side: u8, $$uplo: u8, $$transa: u8, $$diag: u8, $$m: int, $$n: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT;
} {
    #if MT == float32 {
        strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == float64 {
        dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == C32 {
        ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_larzb :: (
    $$side: u8, $$trans: u8, $$direct: u8, $$storev: u8, $$m: int, $$n: int, $$k: int, $$l: int,
    v: $V/ArrayType2D, $$ldv: int, t: $T/ArrayType2D, $$ldt: int, c: *$C/ArrayType2D, $$ldc: int,
    work: *$WORK/ArrayType2D, $$ldwork: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && T.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    one_t :: #run one(MT);
    neg_one_t : MT = -one_t;
    routine_name :: #run if MT == float32 {
        "SLARZB"
    } else if MT == float64 {
        "DLARZB"
    } else if MT == C32 {
        "CLARZB"
    } else {
        "ZLARZB"
    };
    transt : u8 = if lsame(trans, #char "N") {
        #if MT == C32 || MT == C64 { #char "C" } else { #char "T" }
    } else {
        #char "N"
    };
    info : int = 0;

    if m <= 0 || n <= 0 then return;

    if !lsame(direct, #char "B") {
        info = -3;
    } else if !lsame(storev, #char "R") {
        info = -4;
    }
    if info != 0 {
        xerbla(routine_name, -info);
        return;
    }

    if lsame(side, #char "L") {
        for j:1..k {
            for i:1..n {
                fset(work, i, j, fget(c.*, j, i));
            }
        }

        if l > 0 {
            ctail_raw := _larzb_matrix_storage(c, ldc, m - l + 1, 1, l, n);
            ctail := array2d_view(ctail_raw, ldc, row_major=false);
            vblk_raw := _larzb_matrix_storage_view(v, ldv, 1, 1, k, l);
            vblk := array2d_view(vblk_raw, ldv, row_major=false);
            wblk_raw := _larzb_matrix_storage(work, ldwork, 1, 1, n, k);
            wblk := array2d_view(wblk_raw, ldwork, row_major=false);
            transb :: u8 = if MT == C32 || MT == C64 { #char "C" } else { #char "T" };
            _larzb_gemm(#char "T", transb, n, k, l, one_t, ctail.*, ldc, vblk.*, ldv, one_t, *wblk, ldwork, MT);
        }

        tblk_raw := _larzb_matrix_storage_view(t, ldt, 1, 1, k, k);
        tblk := array2d_view(tblk_raw, ldt, row_major=false);
        wblk_raw := _larzb_matrix_storage(work, ldwork, 1, 1, n, k);
        wblk := array2d_view(wblk_raw, ldwork, row_major=false);
        _larzb_trmm(#char "R", #char "L", transt, #char "N", n, k, one_t, tblk.*, ldt, *wblk, ldwork, MT);

        for j:1..n {
            for i:1..k {
                fset(c, i, j, fget(c.*, i, j) - fget(work.*, j, i));
            }
        }

        if l > 0 {
            ctail_raw := _larzb_matrix_storage(c, ldc, m - l + 1, 1, l, n);
            ctail := array2d_view(ctail_raw, ldc, row_major=false);
            vblk_raw := _larzb_matrix_storage_view(v, ldv, 1, 1, k, l);
            vblk := array2d_view(vblk_raw, ldv, row_major=false);
            wblk2_raw := _larzb_matrix_storage(work, ldwork, 1, 1, n, k);
            wblk2 := array2d_view(wblk2_raw, ldwork, row_major=false);
            _larzb_gemm(#char "T", #char "T", l, n, k, neg_one_t, vblk.*, ldv, wblk2.*, ldwork, one_t, *ctail, ldc, MT);
        }
    } else if lsame(side, #char "R") {
        for j:1..k {
            for i:1..m {
                fset(work, i, j, fget(c.*, i, j));
            }
        }

        if l > 0 {
            ctail_raw := _larzb_matrix_storage(c, ldc, 1, n - l + 1, m, l);
            ctail := array2d_view(ctail_raw, ldc, row_major=false);
            vblk_raw := _larzb_matrix_storage_view(v, ldv, 1, 1, k, l);
            vblk := array2d_view(vblk_raw, ldv, row_major=false);
            wblk_raw := _larzb_matrix_storage(work, ldwork, 1, 1, m, k);
            wblk := array2d_view(wblk_raw, ldwork, row_major=false);
            _larzb_gemm(#char "N", #char "T", m, k, l, one_t, ctail.*, ldc, vblk.*, ldv, one_t, *wblk, ldwork, MT);
        }

        #if MT == float32 || MT == float64 {
            tblk_raw := _larzb_matrix_storage_view(t, ldt, 1, 1, k, k);
            tblk := array2d_view(tblk_raw, ldt, row_major=false);
            wblk_raw := _larzb_matrix_storage(work, ldwork, 1, 1, m, k);
            wblk := array2d_view(wblk_raw, ldwork, row_major=false);
            _larzb_trmm(#char "R", #char "L", trans, #char "N", m, k, one_t, tblk.*, ldt, *wblk, ldwork, MT);
        } else {
            sum : MT = ---;
            if lsame(trans, #char "N") {
                for row:1..m {
                    for j:1..k {
                        sum = zero(MT);
                        for p:j..k {
                            sum += fget(work.*, row, p) * conjugate(fget(t, p, j));
                        }
                        fset(work, row, j, sum);
                    }
                }
            } else {
                for row:1..m {
                    for <j:1..k {
                        sum = zero(MT);
                        for p:1..j {
                            sum += fget(work.*, row, p) * fget(t, j, p);
                        }
                        fset(work, row, j, sum);
                    }
                }
            }
        }

        for j:1..k {
            for i:1..m {
                fset(c, i, j, fget(c.*, i, j) - fget(work.*, i, j));
            }
        }

        if l > 0 {
            ctail_raw := _larzb_matrix_storage(c, ldc, 1, n - l + 1, m, l);
            ctail := array2d_view(ctail_raw, ldc, row_major=false);
            wblk2_raw := _larzb_matrix_storage(work, ldwork, 1, 1, m, k);
            wblk2 := array2d_view(wblk2_raw, ldwork, row_major=false);
            vblk2_raw := _larzb_matrix_storage_view(v, ldv, 1, 1, k, l);
            vblk2 := array2d_view(vblk2_raw, ldv, row_major=false);
            #if MT == float32 || MT == float64 {
                _larzb_gemm(#char "N", #char "N", m, l, k, neg_one_t, wblk2.*, ldwork, vblk2.*, ldv, one_t, *ctail, ldc, MT);
            } else {
                sum2 : MT = ---;
                for j:1..l {
                    for i:1..m {
                        sum2 = zero(MT);
                        for p:1..k {
                            sum2 += fget(work.*, i, p) * conjugate(fget(v, p, j));
                        }
                        fset(c, i, n - l + j, fget(c.*, i, n - l + j) - sum2);
                    }
                }
            }
        }
    }
}
