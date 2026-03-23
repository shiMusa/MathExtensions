// ? LAPACK *hetrf_aa_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrf_aa_2stage/zhetrf_aa_2stage
// *HETRF_AA_2STAGE computes the two-stage Aasen factorization of a
// complex Hermitian matrix.

chetrf_aa_2stage :: #bake_arguments _hetrf_aa_2stage(RT=float32, MT=C32);
zhetrf_aa_2stage :: #bake_arguments _hetrf_aa_2stage(RT=float64, MT=C64);

#scope_file
_hetrf_aa_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrf_aa_2stage_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrf_aa_2stage :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tb: *$TB/ArrayType1D, $$ltb: int, ipiv: *$IPIV/ArrayType1D,
    ipiv2: *$IPIV2/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && TB.data_type == MT &&
           IPIV.data_type == int && IPIV2.data_type == int &&
           WORK.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHETRF_AA_2STAGE" } else { "ZHETRF_AA_2STAGE" };

    tquery, upper, wquery : bool = ---;
    nb : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    wquery = lwork == -1;
    tquery = ltb == -1;
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if ltb < max(1, 4*n) && !tquery {
        info.* = -6;
    } else if lwork < max(1, n) && !wquery {
        info.* = -10;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    nb = ilaenv(1, ifx MT == C32 then "CHETRF_AA_2STAGE" else "ZHETRF_AA_2STAGE", cast(*u8) *uplo, n, -1, -1, -1);
    if tquery {
        fset(tb, 1, _hetrf_aa_2stage_complex_from_real(_hetrf_aa_2stage_roundup_lwork(ncast(RT, max(1, (3*nb + 1)*n)), RT), RT, MT));
    }
    if wquery {
        fset(work, 1, _hetrf_aa_2stage_complex_from_real(_hetrf_aa_2stage_roundup_lwork(ncast(RT, max(1, n*nb)), RT), RT, MT));
    }
    if tquery || wquery {
        return;
    }

    if n == 0 {
        return;
    }

    // Staged port note:
    // This file establishes the family entry points, argument checking,
    // and workspace-query behavior for *HETRF_AA_2STAGE. The large
    // two-stage Aasen execution body is still pending.
    fset(tb, 1, _hetrf_aa_2stage_complex_from_real(ncast(RT, nb), RT, MT));
    fset(ipiv, 1, 1);
    fset(ipiv2, 1, 1);
}
