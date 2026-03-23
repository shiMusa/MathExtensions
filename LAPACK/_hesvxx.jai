// ? LAPACK *hesvxx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chesvxx/zhesvxx
// *HESVXX uses the diagonal pivoting factorization to compute the
// solution to a system of linear equations for Hermitian matrices.

chesvxx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    a: *$A/ArrayType2D, $$lda: int, af: *$AF/ArrayType2D, $$ldaf: int,
    ipiv: *$IPIV/ArrayType1D, equed: *u8, s: *$S/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    rcond: *float32, rpvgrw: *float32, berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int, err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D, $$nparams: int,
    params: *$PARAMS/ArrayType1D, work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D, info: *int
) {
    _hesvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float32, MT=C32);
}

zhesvxx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    a: *$A/ArrayType2D, $$lda: int, af: *$AF/ArrayType2D, $$ldaf: int,
    ipiv: *$IPIV/ArrayType1D, equed: *u8, s: *$S/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    rcond: *float64, rpvgrw: *float64, berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int, err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D, $$nparams: int,
    params: *$PARAMS/ArrayType1D, work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D, info: *int
) {
    _hesvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_hesvxx_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "s"); } else { return dlamch(#char "s"); }
}

_hesvxx_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 { clacpy(uplo, m, n, a, lda, b, ldb); } else { zlacpy(uplo, m, n, a, lda, b, ldb); }
}

_hesvxx_heequb :: (
    $$uplo: u8, $$n: int, a: $A/ArrayType2D, $$lda: int, s: *$S/ArrayType1D,
    scond: *$RT, amax: *$RT, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && S.data_type == RT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cheequb(uplo, n, a, lda, s, scond, amax, work, info);
    } else {
        zheequb(uplo, n, a, lda, s, scond, amax, work, info);
    }
}

_hesvxx_laqhe :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, s: $S/ArrayType1D,
    scond: $RT, amax: $RT, equed: *u8, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && S.data_type == RT;
} {
    #if MT == C32 {
        claqhe(uplo, n, a, lda, s, scond, amax, equed);
    } else {
        zlaqhe(uplo, n, a, lda, s, scond, amax, equed);
    }
}

_hesvxx_hetrf :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, ipiv: *$IPIV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 { chetrf(uplo, n, a, lda, ipiv, work, lwork, info); } else { zhetrf(uplo, n, a, lda, ipiv, work, lwork, info); }
}

_hesvxx_hetrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, af: $AF/ArrayType2D, $$ldaf: int,
    ipiv: $IPIV/ArrayType1D, x: *$X/ArrayType2D, $$ldx: int,
    info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           AF.data_type == MT && IPIV.data_type == int && X.data_type == MT;
} {
    #if MT == C32 { chetrs(uplo, n, nrhs, af, ldaf, ipiv, x, ldx, info); } else { zhetrs(uplo, n, nrhs, af, ldaf, ipiv, x, ldx, info); }
}

_hesvxx_lascl2 :: (
    $$m: int, $$n: int, d: $D/ArrayType1D, x: *$X/ArrayType2D, $$ldx: int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT && X.data_type == MT;
} {
    #if MT == C32 { clascl2(m, n, d, x, ldx); } else { zlascl2(m, n, d, x, ldx); }
}

_hesvxx_herpvgrw :: (
    $$uplo: u8, $$n: int, $$info_val: int, a: $A/ArrayType2D, $$lda: int,
    af: $AF/ArrayType2D, $$ldaf: int, ipiv: $IPIV/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && AF.data_type == MT &&
           IPIV.data_type == int && RWORK.data_type == RT;
} {
    #if MT == C32 {
        return cla_herpvgrw(uplo, n, info_val, a, lda, af, ldaf, ipiv, rwork);
    } else {
        return zla_herpvgrw(uplo, n, info_val, a, lda, af, ldaf, ipiv, rwork);
    }
}

_hesvxx_herfsx :: (
    $$uplo: u8, $$equed: u8, $$n: int, $$nrhs: int, a: $A/ArrayType2D, $$lda: int,
    af: $AF/ArrayType2D, $$ldaf: int, ipiv: $IPIV/ArrayType1D, s: $S/ArrayType1D,
    b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    rcond: *$RT, berr: *$BERR/ArrayType1D, $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D, err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int, params: *$PARAMS/ArrayType1D, work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && AF.data_type == MT && IPIV.data_type == int &&
           S.data_type == RT && B.data_type == MT && X.data_type == MT &&
           BERR.data_type == RT && ERR_BNDS_NORM.data_type == RT &&
           ERR_BNDS_COMP.data_type == RT && PARAMS.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 {
        cherfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
    } else {
        zherfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
    }
}

_hesvxx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    a: *$A/ArrayType2D, $$lda: int, af: *$AF/ArrayType2D, $$ldaf: int,
    ipiv: *$IPIV/ArrayType1D, equed: *u8, s: *$S/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    rcond: *$RT, rpvgrw: *$RT, berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int, err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D, $$nparams: int,
    params: *$PARAMS/ArrayType1D, work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && AF.data_type == MT && IPIV.data_type == int &&
           S.data_type == RT && B.data_type == MT && X.data_type == MT &&
           BERR.data_type == RT && ERR_BNDS_NORM.data_type == RT &&
           ERR_BNDS_COMP.data_type == RT && PARAMS.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    routine_name :: #run if MT == C32 { "CHESVXX" } else { "ZHESVXX" };

    amax, bignum, scond, smax, smin, smlnum : RT = ---;
    equil, nofact, rcequ, upper : bool = ---;
    infequ, j : int = ---;

    info.* = 0;
    nofact = lsame(fact, #char "N");
    equil = lsame(fact, #char "E");
    upper = lsame(uplo, #char "U");
    smlnum = _hesvxx_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    if nofact || equil {
        equed.* = #char "N";
        rcequ = false;
    } else {
        rcequ = lsame(equed.*, #char "Y");
    }

    rpvgrw.* = zero_r;

    if !nofact && !equil && !lsame(fact, #char "F") {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if nrhs < 0 {
        info.* = -4;
    } else if lda < max(1, n) {
        info.* = -6;
    } else if ldaf < max(1, n) {
        info.* = -8;
    } else if lsame(fact, #char "F") && !(rcequ || lsame(equed.*, #char "N")) {
        info.* = -9;
    } else {
        if rcequ {
            smin = bignum;
            smax = zero_r;
            for j:1..n {
                smin = min(smin, fget(s.*, j));
                smax = max(smax, fget(s.*, j));
            }
            if smin <= zero_r {
                info.* = -10;
            } else if n > 0 {
                scond = max(smin, smlnum) / min(smax, bignum);
            } else {
                scond = one_r;
            }
        }
        if info.* == 0 {
            if ldb < max(1, n) {
                info.* = -12;
            } else if ldx < max(1, n) {
                info.* = -14;
            }
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if equil {
        _hesvxx_heequb(uplo, n, a.*, lda, s, *scond, *amax, work, *infequ, RT, MT);
        if infequ == 0 {
            _hesvxx_laqhe(uplo, n, a, lda, s.*, scond, amax, equed, RT, MT);
            rcequ = lsame(equed.*, #char "Y");
        }
    }

    if rcequ {
        _hesvxx_lascl2(n, nrhs, s.*, b, ldb, RT, MT);
    }

    if nofact || equil {
        _hesvxx_lacpy(uplo, n, n, a.*, lda, af, ldaf, MT);
        _hesvxx_hetrf(uplo, n, af, ldaf, ipiv, work, 5*max(1, n), info, MT);
        if info.* > 0 {
            if n > 0 {
                rpvgrw.* = _hesvxx_herpvgrw(uplo, n, info.*, a.*, lda, af.*, ldaf, ipiv.*, rwork, RT, MT);
            }
            return;
        }
    }

    if n > 0 {
        rpvgrw.* = _hesvxx_herpvgrw(uplo, n, info.*, a.*, lda, af.*, ldaf, ipiv.*, rwork, RT, MT);
    }

    _hesvxx_lacpy(#char "F", n, nrhs, b.*, ldb, x, ldx, MT);
    _hesvxx_hetrs(uplo, n, nrhs, af.*, ldaf, ipiv.*, x, ldx, info, MT);
    _hesvxx_herfsx(uplo, equed.*, n, nrhs, a.*, lda, af.*, ldaf, ipiv.*, s.*, b.*, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT, MT);

    if rcequ {
        _hesvxx_lascl2(n, nrhs, s.*, x, ldx, RT, MT);
    }
}
