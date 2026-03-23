// ? LAPACK *hbgv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chbgv/zhbgv
// *HBGV computes all eigenvalues and, optionally, eigenvectors of a
// complex generalized Hermitian-definite banded eigenproblem.

chbgv :: #bake_arguments _hbgv(RT=float32, MT=C32);
zhbgv :: #bake_arguments _hbgv(RT=float64, MT=C64);

#scope_file
_hbgv_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hbgv_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hbgv_pbstf :: (
    $$uplo: u8,
    $$n: int,
    $$kb: int,
    bb: *$BB/ArrayType2D,
    $$ldbb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && BB.data_type == MT;
} {
    #if MT == C32 {
        cpbstf(uplo, n, kb, bb, ldbb, info);
    } else {
        zpbstf(uplo, n, kb, bb, ldbb, info);
    }
}

_hbgv_hbgst :: (
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
    #if MT == C32 {
        chbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info);
    } else {
        zhbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info);
    }
}

_hbgv_hbtrd :: (
    $$vect: u8,
    $$uplo: u8,
    $$n: int,
    $$kd: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           Q.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        chbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
    } else {
        zhbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
    }
}

_hbgv_steqr :: (
    $$compz: u8,
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           D.data_type == RT &&
           E.data_type == RT &&
           Z.data_type == MT &&
           RWORK.data_type == RT;
} {
    #if MT == C32 {
        csteqr(compz, n, d, e, z, ldz, work, info);
    } else {
        zsteqr(compz, n, d, e, z, ldz, work, info);
    }
}

_hbgv_sterf :: (
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           E.data_type == RT;
} {
    #if RT == float32 {
        ssterf(n, d, e, info);
    } else {
        dsterf(n, d, e, info);
    }
}

_hbgv :: (
    $$jobz: u8,
    $$uplo: u8,
    $$n: int,
    $$ka: int,
    $$kb: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    bb: *$BB/ArrayType2D,
    $$ldbb: int,
    w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
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
           W.data_type == RT &&
           Z.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    zero_t :: #run _hbgv_complex_from_real(zero(RT), RT, MT);
    cone_t :: #run _hbgv_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHBGV" } else { "ZHBGV" };

    iinfo, inde, indwrk : int = ---;
    upper, wantz : bool = ---;
    vect : u8 = ---;

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
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
    } else if ldz < 1 || (wantz && ldz < n) {
        info.* = -12;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    _hbgv_pbstf(uplo, n, kb, bb, ldbb, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    inde = 1;
    indwrk = inde + n;
    _hbgv_hbgst(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, z, ldz, work, rwork, *iinfo, RT, MT);

    if wantz {
        vect = #char "U";
    } else {
        vect = #char "N";
    }

    e_view := _hbgv_subvector(rwork.*, inde, n);
    _hbgv_hbtrd(vect, uplo, n, ka, ab, ldab, w, *e_view, z, ldz, work, *iinfo, RT, MT);

    if !wantz {
        _hbgv_sterf(n, w, *e_view, info, RT);
    } else {
        rwork_tail2 := _hbgv_subvector(rwork.*, indwrk, 3*n - indwrk + 1);
        _hbgv_steqr(jobz, n, w, *e_view, z, ldz, *rwork_tail2, info, RT, MT);
    }
}
