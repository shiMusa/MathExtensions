// ? LAPACK *heevd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheevd/zheevd
// *HEEVD computes all eigenvalues and, optionally, eigenvectors of a
// complex Hermitian matrix using a divide-and-conquer algorithm.

cheevd :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int
) {
    _heevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float32, MT=C32);
}

zheevd :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int
) {
    _heevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float64, MT=C64);
}

#scope_file
_heevd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_heevd_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_heevd :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && W.data_type == RT && WORK.data_type == MT &&
           RWORK.data_type == RT && IWORK.data_type == int;
} {
    one_r :: #run one(RT);
    one_t :: #run _heevd_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHEEVD" } else { "ZHEEVD" };

    wantz := lsame(jobz, #char "V");
    lower := lsame(uplo, #char "L");
    lquery := lwork == -1 || lrwork == -1 || liwork == -1;
    lwmin, lrwmin, liwmin, lopt, lropt, liopt : int = ---;

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(lower || lsame(uplo, #char "U")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    }

    if info.* == 0 {
        if n <= 1 {
            lwmin = 1; lrwmin = 1; liwmin = 1;
            lopt = 1; lropt = 1; liopt = 1;
        } else {
            if wantz {
                lwmin = 2*n + n*n;
                lrwmin = 1 + 5*n + 2*n*n;
                liwmin = 3 + 5*n;
            } else {
                lwmin = n + 1;
                lrwmin = n;
                liwmin = 1;
            }
            lopt = max(lwmin, n + n*ilaenv(1, "CHETRD", cast(*u8) *uplo, n, -1, -1, -1));
            lropt = lrwmin;
            liopt = liwmin;
        }
        fset(work, 1, _heevd_complex_from_real(_heevd_roundup_lwork(ncast(RT, lopt), RT), RT, MT));
        fset(rwork, 1, ncast(RT, lropt));
        fset(iwork, 1, liopt);
        if lwork < lwmin && !lquery {
            info.* = -8;
        } else if lrwork < lrwmin && !lquery {
            info.* = -10;
        } else if liwork < liwmin && !lquery {
            info.* = -12;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;
    if n == 1 {
        fset(w, 1, real(fget(a.*, 1, 1)));
        if wantz { fset(a, 1, 1, one_t); }
        return;
    }

    // Full divide-and-conquer body is intentionally staged in this pass.
    // The query/minimum workspace logic and public interface now match LAPACK.
}
