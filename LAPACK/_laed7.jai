// ? LAPACK *laed7()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaed7/claed7/dlaed7/zlaed7
// *LAED7 computes the updated eigensystem of a diagonal matrix after
// modification by a rank-one symmetric matrix.

slaed7 :: #bake_arguments _laed7_real(RT=float32);
dlaed7 :: #bake_arguments _laed7_real(RT=float64);
claed7 :: #bake_arguments _laed7_complex(RT=float32, MT=C32);
zlaed7 :: #bake_arguments _laed7_complex(RT=float64, MT=C64);

#scope_file
_laed7_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_laed7_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _laed7_raw_slice(x.data + (offset - 1), count);
}

_laed7_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _laed7_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laed7_gemm_real :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int, alpha: $RT,
    a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int, beta: $RT,
    c: *$C/ArrayType2D, $$ldc: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           B.data_type == RT &&
           C.data_type == RT;
} {
    #if RT == float32 {
        sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_laed7_clacrm :: (
    $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int,
    c: *$C/ArrayType2D, $$ldc: int, rwork: *$RWORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == RT &&
           C.data_type == MT &&
           RWORK.data_type == RT;
} {
    #if MT == C32 {
        clacrm(m, n, a, lda, b, ldb, c, ldc, rwork);
    } else {
        zlacrm(m, n, a, lda, b, ldb, c, ldc, rwork);
    }
}

_laed7_lamrg :: ($$n1: int, $$n2: int, a: $A/ArrayType1D, $$dtrd1: int, $$dtrd2: int, index: *$INDEX/ArrayType1D, $RT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           INDEX.data_type == int;
} {
    #if RT == float32 {
        slamrg(n1, n2, a, dtrd1, dtrd2, index);
    } else {
        dlamrg(n1, n2, a, dtrd1, dtrd2, index);
    }
}

_laed7_laed9 :: (
    $$k: int, $$kstart: int, $$kstop: int, $$n: int,
    d: *$D/ArrayType1D, q: *$Q/ArrayType2D, $$ldq: int, $$rho: $RT,
    dlamda: *$DLAMDA/ArrayType1D, w: *$W/ArrayType1D, s: *$S/ArrayType1D, $$lds: int,
    info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           Q.data_type == RT &&
           DLAMDA.data_type == RT &&
           W.data_type == RT &&
           S.data_type == RT;
} {
    #if RT == float32 {
        slaed9(k, kstart, kstop, n, d, q, ldq, rho, dlamda, w, s, lds, info);
    } else {
        dlaed9(k, kstart, kstop, n, d, q, ldq, rho, dlamda, w, s, lds, info);
    }
}

_laed7_laeda :: (
    $$n: int, $$tlvls: int, $$curlvl: int, $$curpbm: int,
    prmptr: *$PRMPTR/ArrayType1D, perm: *$PERM/ArrayType1D, givptr: *$GIVPTR/ArrayType1D,
    givcol: *$GIVCOL/ArrayType2D, givnum: *$GIVNUM/ArrayType2D, qstore: *$QSTORE/ArrayType1D,
    qptr: *$QPTR/ArrayType1D, z: *$Z/ArrayType1D, ztemp: *$ZTEMP/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           PRMPTR.data_type == int &&
           PERM.data_type == int &&
           GIVPTR.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           QSTORE.data_type == RT &&
           QPTR.data_type == int &&
           Z.data_type == RT &&
           ZTEMP.data_type == RT;
} {
    #if RT == float32 {
        slaeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, qstore, qptr, z, ztemp, info);
    } else {
        dlaeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, qstore, qptr, z, ztemp, info);
    }
}

_laed7_real :: (
    $$icompq: int,
    $$n: int,
    $$qsiz: int,
    $$tlvls: int,
    $$curlvl: int,
    $$curpbm: int,
    d: *$D/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    indxq: *$INDXQ/ArrayType1D,
    $$rho: $RT,
    $$cutpnt: int,
    qstore: *$QSTORE/ArrayType1D,
    qptr: *$QPTR/ArrayType1D,
    prmptr: *$PRMPTR/ArrayType1D,
    perm: *$PERM/ArrayType1D,
    givptr: *$GIVPTR/ArrayType1D,
    givcol: *$GIVCOL/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           Q.data_type == RT &&
           INDXQ.data_type == int &&
           QSTORE.data_type == RT &&
           QPTR.data_type == int &&
           PRMPTR.data_type == int &&
           PERM.data_type == int &&
           GIVPTR.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { "SLAED7" } else { "DLAED7" };

    coltpy, curr, i, idlmda, indx, indxc, indxp, iq2, is, iw, iz, k, ldq2, n1, n2, ptr : int = ---;

    info.* = 0;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if icompq == 1 && qsiz < n {
        info.* = -3;
    } else if ldq < max(1, n) {
        info.* = -9;
    } else if min(1, n) > cutpnt || n < cutpnt {
        info.* = -12;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    if icompq == 1 {
        ldq2 = qsiz;
    } else {
        ldq2 = n;
    }

    iz = 1;
    idlmda = iz + n;
    iw = idlmda + n;
    iq2 = iw + n;
    is = iq2 + n*ldq2;

    indx = 1;
    indxc = indx + n;
    coltpy = indxc + n;
    indxp = coltpy + n;

    ptr = 1 + (1 << tlvls);
    for i:1..curlvl - 1 {
        ptr += 1 << (tlvls - i);
    }
    curr = ptr + curpbm;

    z_vec := _laed7_subvector(work.*, iz, n);
    ztemp_vec := _laed7_subvector(work.*, iz + n, n);
    _laed7_laeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, qstore, qptr, *z_vec, *ztemp_vec, info, RT);

    if curlvl == tlvls {
        fset(qptr, curr, 1);
        fset(prmptr, curr, 1);
        fset(givptr, curr, 1);
    }

    dlambda_vec := _laed7_subvector(work.*, idlmda, n);
    q2_raw := _laed7_raw_slice(work.data + (iq2 - 1), n*ldq2);
    q2_mat := array2d_view(q2_raw, ldq2, row_major=false);
    w_vec := _laed7_subvector(work.*, iw, n);
    perm_curr := _laed7_subvector(perm.*, fget(prmptr.*, curr), max(0, perm.count - fget(prmptr.*, curr) + 1));
    givcol_curr_raw := _laed7_raw_slice(givcol.data + (fget(givptr.*, curr) - 1)*2, max(0, givcol.count - (fget(givptr.*, curr) - 1)*2));
    givcol_curr := array2d_view(givcol_curr_raw, 2, row_major=false);
    givnum_curr_raw := _laed7_raw_slice(givnum.data + (fget(givptr.*, curr) - 1)*2, max(0, givnum.count - (fget(givptr.*, curr) - 1)*2));
    givnum_curr := array2d_view(givnum_curr_raw, 2, row_major=false);
    indxp_vec := _laed7_subvector(iwork.*, indxp, n);
    indx_vec := _laed7_subvector(iwork.*, indx, n);

    _laed8_real(icompq, &k, n, qsiz, d, q, ldq, indxq, &rho, cutpnt, *z_vec, *dlambda_vec, *q2_mat, ldq2, *w_vec,
                *perm_curr, givptr.data + curr, *givcol_curr, *givnum_curr, *indxp_vec, *indx_vec, info, RT);
    fset(prmptr, curr + 1, fget(prmptr.*, curr) + n);
    fset(givptr, curr + 1, fget(givptr.*, curr + 1) + fget(givptr.*, curr));

    if k != 0 {
        s_vec := _laed7_subvector(work.*, is, k*k);
        qstore_sub := _laed7_subvector(qstore.*, fget(qptr.*, curr), max(0, qstore.count - fget(qptr.*, curr) + 1));
        _laed7_laed9(k, 1, k, n, d, *array2d_view(_laed7_raw_slice(work.data + (is - 1), k*k), k, row_major=false), k, rho,
                     *dlambda_vec, *w_vec, *qstore_sub, k, info, RT);
        if info.* != 0 then return;

        if icompq == 1 {
            q2k_raw := _laed7_raw_slice(work.data + (iq2 - 1), k*ldq2);
            q2k := array2d_view(q2k_raw, ldq2, row_major=false);
            qstore_k_raw := _laed7_raw_slice(qstore.data + (fget(qptr.*, curr) - 1), k*k);
            qstore_k := array2d_view(qstore_k_raw, k, row_major=false);
            qk_raw := _laed7_matrix_storage(q.*, ldq, 1, 1, qsiz, k);
            qk := array2d_view(qk_raw, ldq, row_major=false);
            _laed7_gemm_real(#char "N", #char "N", qsiz, k, k, one_r, q2k.*, ldq2, qstore_k.*, k, zero_r, *qk, ldq, RT);
        }
        fset(qptr, curr + 1, fget(qptr.*, curr) + k*k);

        n1 = k;
        n2 = n - k;
        _laed7_lamrg(n1, n2, d, 1, -1, indxq, RT);
    } else {
        fset(qptr, curr + 1, fget(qptr.*, curr));
        for i:1..n {
            fset(indxq, i, i);
        }
    }
}

_laed7_complex :: (
    $$n: int,
    $$cutpnt: int,
    $$qsiz: int,
    $$tlvls: int,
    $$curlvl: int,
    $$curpbm: int,
    d: *$D/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    $$rho: $RT,
    indxq: *$INDXQ/ArrayType1D,
    qstore: *$QSTORE/ArrayType1D,
    qptr: *$QPTR/ArrayType1D,
    prmptr: *$PRMPTR/ArrayType1D,
    perm: *$PERM/ArrayType1D,
    givptr: *$GIVPTR/ArrayType1D,
    givcol: *$GIVCOL/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
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
           Q.data_type == MT &&
           INDXQ.data_type == int &&
           QSTORE.data_type == RT &&
           QPTR.data_type == int &&
           PRMPTR.data_type == int &&
           PERM.data_type == int &&
           GIVPTR.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    coltpy, curr, i, idlmda, indx, indxc, indxp, iq, iw, iz, k, n1, n2, ptr : int = ---;
    routine_name :: #run if MT == C32 { "CLAED7" } else { "ZLAED7" };

    info.* = 0;
    if n < 0 {
        info.* = -1;
    } else if min(1, n) > cutpnt || n < cutpnt {
        info.* = -2;
    } else if qsiz < n {
        info.* = -3;
    } else if ldq < max(1, n) {
        info.* = -9;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    iz = 1;
    idlmda = iz + n;
    iw = idlmda + n;
    iq = iw + n;

    indx = 1;
    indxc = indx + n;
    coltpy = indxc + n;
    indxp = coltpy + n;

    ptr = 1 + (1 << tlvls);
    for i:1..curlvl - 1 {
        ptr += 1 << (tlvls - i);
    }
    curr = ptr + curpbm;

    z_vec := _laed7_subvector(rwork.*, iz, n);
    ztemp_vec := _laed7_subvector(rwork.*, iz + n, n);
    _laed7_laeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, qstore, qptr, *z_vec, *ztemp_vec, info, RT);

    if curlvl == tlvls {
        fset(qptr, curr, 1);
        fset(prmptr, curr, 1);
        fset(givptr, curr, 1);
    }

    dlambda_vec := _laed7_subvector(rwork.*, idlmda, n);
    w_vec := _laed7_subvector(rwork.*, iw, n);
    work_mat_raw := _laed7_raw_slice(work.data, qsiz*n);
    work_mat := array2d_view(work_mat_raw, qsiz, row_major=false);
    indxp_vec := _laed7_subvector(iwork.*, indxp, n);
    indx_vec := _laed7_subvector(iwork.*, indx, n);
    perm_curr := _laed7_subvector(perm.*, fget(prmptr.*, curr), max(0, perm.count - fget(prmptr.*, curr) + 1));
    givcol_curr_raw := _laed7_raw_slice(givcol.data + (fget(givptr.*, curr) - 1)*2, max(0, givcol.count - (fget(givptr.*, curr) - 1)*2));
    givcol_curr := array2d_view(givcol_curr_raw, 2, row_major=false);
    givnum_curr_raw := _laed7_raw_slice(givnum.data + (fget(givptr.*, curr) - 1)*2, max(0, givnum.count - (fget(givptr.*, curr) - 1)*2));
    givnum_curr := array2d_view(givnum_curr_raw, 2, row_major=false);

    _laed8_complex(&k, n, qsiz, q, ldq, d, &rho, cutpnt, *z_vec, *dlambda_vec, *work_mat, qsiz, *w_vec,
                   *indxp_vec, *indx_vec, indxq, *perm_curr, givptr.data + curr, *givcol_curr, *givnum_curr, info, RT, MT);
    fset(prmptr, curr + 1, fget(prmptr.*, curr) + n);
    fset(givptr, curr + 1, fget(givptr.*, curr + 1) + fget(givptr.*, curr));

    if k != 0 {
        qstore_sub := _laed7_subvector(qstore.*, fget(qptr.*, curr), max(0, qstore.count - fget(qptr.*, curr) + 1));
        _laed7_laed9(k, 1, k, n, d, *array2d_view(_laed7_raw_slice(rwork.data + (iq - 1), k*k), k, row_major=false), k, rho,
                     *dlambda_vec, *w_vec, *qstore_sub, k, info, RT);
        _laed7_clacrm(qsiz, k, work, qsiz, *array2d_view(_laed7_raw_slice(qstore.data + (fget(qptr.*, curr) - 1), k*k), k, row_major=false), k,
                      q, ldq, *_laed7_subvector(rwork.*, iq, max(0, rwork.count - iq + 1)), RT, MT);
        fset(qptr, curr + 1, fget(qptr.*, curr) + k*k);
        if info.* != 0 then return;

        n1 = k;
        n2 = n - k;
        _laed7_lamrg(n1, n2, d, 1, -1, indxq, RT);
    } else {
        fset(qptr, curr + 1, fget(qptr.*, curr));
        for i:1..n {
            fset(indxq, i, i);
        }
    }
}
