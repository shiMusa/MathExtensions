// ? LAPACK *hetri2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetri2/zhetri2
// *HETRI2 computes the inverse of a complex Hermitian indefinite matrix
// using the factorization from *HETRF and dispatches to *HETRI or *HETRI2X.

chetri2 :: #bake_arguments _hetri2(RT=float32, MT=C32);
zhetri2 :: #bake_arguments _hetri2(RT=float64, MT=C64);

#scope_file
_hetri2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetri2_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetri2 :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHETRI2" } else { "ZHETRI2" };

    lquery, upper : bool = ---;
    minsize, nbmax : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;

    nbmax = ilaenv(1, ifx MT == C32 then "CHETRF" else "ZHETRF", cast(*u8) *uplo, n, -1, -1, -1);
    if n == 0 {
        minsize = 1;
    } else if nbmax >= n {
        minsize = n;
    } else {
        minsize = (n + nbmax + 1)*(nbmax + 3);
    }

    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if lwork < minsize && !lquery {
        info.* = -7;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        fset(work, 1, _hetri2_complex_from_real(_hetri2_roundup_lwork(ncast(RT, minsize), RT), RT, MT));
        return;
    }

    if n == 0 {
        return;
    }

    if nbmax >= n {
        #if MT == C32 {
            chetri(uplo, n, a, lda, ipiv, work, info);
        } else {
            zhetri(uplo, n, a, lda, ipiv, work, info);
        }
    } else {
        #if MT == C32 {
            chetri2x(uplo, n, a, lda, ipiv, work, nbmax, info);
        } else {
            zhetri2x(uplo, n, a, lda, ipiv, work, nbmax, info);
        }
    }
}
