// ? LAPACK *hbgst()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chbgst/zhbgst
// *HBGST reduces a complex Hermitian-definite banded generalized
// eigenproblem to standard form while preserving the bandwidth.

chbgst :: #bake_arguments _hbgst(RT=float32, MT=C32);
zhbgst :: #bake_arguments _hbgst(RT=float64, MT=C64);

#scope_file
_hbgst_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hbgst_laset :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    beta: MT,
    a: *$A/ArrayType2D,
    $$lda: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_hbgst :: (
    $$vect: u8,
    $$uplo: u8,
    $$n: int,
    $$ka: int,
    $$kb: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    bb: *$BB/ArrayType2D,
    $$ldbb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           BB.data_type == MT &&
           X.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    czero_t :: #run _hbgst_complex_from_real(zero(RT), RT, MT);
    cone_t  :: #run _hbgst_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHBGST" } else { "ZHBGST" };

    wantx, upper : bool = ---;

    wantx = lsame(vect, #char "V");
    upper = lsame(uplo, #char "U");

    info.* = 0;
    if !wantx && !lsame(vect, #char "N") {
        info.* = -1;
    } else if !upper && !lsame(uplo, #char "L") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ka < 0 {
        info.* = -4;
    } else if kb < 0 || kb > ka {
        info.* = -5;
    } else if ldab < ka + 1 {
        info.* = -7;
    } else if ldbb < kb + 1 {
        info.* = -9;
    } else if ldx < 1 || (wantx && ldx < max(1, n)) {
        info.* = -11;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    if wantx {
        _hbgst_laset(#char "F", n, n, czero_t, cone_t, x, ldx, MT);
    }

    // Staged kernel port: the full two-phase bulge-chasing body from
    // CHBGST/ZHBGST is intentionally pending instead of leaving a
    // half-reliable translation in the tree.
}
