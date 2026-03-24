// ? LAPACK *lalsd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slalsd/dlalsd/clalsd/zlalsd
// *LALSD uses the singular value decomposition of a bidiagonal matrix to
// solve a least-squares problem.

slalsd :: #bake_arguments _lalsd_real(RT=float32);
dlalsd :: #bake_arguments _lalsd_real(RT=float64);
clalsd :: #bake_arguments _lalsd_complex(RT=float32, MT=C32);
zlalsd :: #bake_arguments _lalsd_complex(RT=float64, MT=C64);

#scope_file
_lalsd_lamch :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "E");
    } else {
        return dlamch(#char "E");
    }
}

_lalsd_laset_real :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $RT, beta: RT, a: *$A/ArrayType2D, $$lda: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT;
} {
    #if RT == float32 {
        slaset(uplo, m, n, alpha, beta, a, lda);
    } else {
        dlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_lalsd_laset_complex :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: MT, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_lalsd_lascl_real :: (
    $$type: u8, $$kl: int, $$ku: int, cfrom: $RT, cto: RT,
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT;
} {
    #if RT == float32 {
        slascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_lalsd_lascl_complex :: (
    $$type: u8, $$kl: int, $$ku: int, cfrom: $RT, cto: RT,
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT;
} {
    #if MT == C32 {
        clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_lalsd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_lalsd_real :: (
    $$uplo: u8,
    $$smlsiz: int,
    $$n: int,
    $$nrhs: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $$rcond: $RT,
    rank: *int,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           E.data_type == RT &&
           B.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SLALSD" } else { "DLALSD" };
    zero_r :: #run zero(RT);
    one_r :: #run one(RT);

    eps, rcnd : RT = ---;

    info.* = 0;
    if n < 0 {
        info.* = -3;
    } else if nrhs < 1 {
        info.* = -4;
    } else if ldb < 1 || ldb < n {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    eps = _lalsd_lamch(RT);
    if rcond <= zero_r || rcond >= one_r {
        rcnd = eps;
    } else {
        rcnd = rcond;
    }

    rank.* = 0;
    if n == 0 then return;

    if n == 1 {
        if fget(d.*, 1) == zero_r {
            _lalsd_laset_real(#char "A", 1, nrhs, zero_r, zero_r, b, ldb, RT);
        } else {
            rank.* = 1;
            _lalsd_lascl_real(#char "G", 0, 0, fget(d.*, 1), one_r, 1, nrhs, b, ldb, info, RT);
            fset(d, 1, abs(fget(d.*, 1)));
        }
        return;
    }

    // Staged port note:
    // The full divide-and-conquer solve tree from *lalsd.f is still pending.
    _ = uplo; _ = smlsiz; _ = e; _ = work; _ = iwork; _ = rcnd;
}

_lalsd_complex :: (
    $$uplo: u8,
    $$smlsiz: int,
    $$n: int,
    $$nrhs: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $$rcond: $RT,
    rank: *int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT &&
           E.data_type == RT &&
           B.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CLALSD" } else { "ZLALSD" };
    zero_r :: #run zero(RT);
    one_r :: #run one(RT);
    zero_t :: #run _lalsd_complex_from_real(zero_r, RT, MT);

    eps, rcnd : RT = ---;

    info.* = 0;
    if n < 0 {
        info.* = -3;
    } else if nrhs < 1 {
        info.* = -4;
    } else if ldb < 1 || ldb < n {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    eps = _lalsd_lamch(RT);
    if rcond <= zero_r || rcond >= one_r {
        rcnd = eps;
    } else {
        rcnd = rcond;
    }

    rank.* = 0;
    if n == 0 then return;

    if n == 1 {
        if fget(d.*, 1) == zero_r {
            _lalsd_laset_complex(#char "A", 1, nrhs, zero_t, zero_t, b, ldb, MT);
        } else {
            rank.* = 1;
            _lalsd_lascl_complex(#char "G", 0, 0, fget(d.*, 1), one_r, 1, nrhs, b, ldb, info, RT, MT);
            fset(d, 1, abs(fget(d.*, 1)));
        }
        return;
    }

    // Staged port note:
    // The complex divide-and-conquer solve tree is still pending.
    _ = uplo; _ = smlsiz; _ = e; _ = work; _ = rwork; _ = iwork; _ = rcnd;
}
