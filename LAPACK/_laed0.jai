// ? LAPACK *laed0()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaed0/claed0/dlaed0/zlaed0
// *LAED0 computes all eigenvalues and corresponding eigenvectors of an
// unreduced symmetric tridiagonal matrix using the divide and conquer
// method.

slaed0 :: #bake_arguments _laed0_real(RT=float32);
dlaed0 :: #bake_arguments _laed0_real(RT=float64);
claed0 :: #bake_arguments _laed0_complex(RT=float32, MT=C32);
zlaed0 :: #bake_arguments _laed0_complex(RT=float64, MT=C64);

#scope_file
_laed0_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_laed0_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _laed0_raw_slice(x.data + (offset - 1), count);
}

_laed0_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _laed0_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laed0_ilaenv :: ($$ispec: int, $$name: string, $RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    return ilaenv(ispec, name, " ", 0, 0, 0, 0);
}

_laed0_real_abs_e :: inline (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    return abs(x);
}

_laed0_copy :: ($$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $T: Type) #modify {
    return (T == float32 || T == float64 || T == C32 || T == C64) &&
           X.data_type == T &&
           Y.data_type == T;
} {
    #if T == float32 {
        scopy(n, x, incx, y, incy);
    } else if T == float64 {
        dcopy(n, x, incx, y, incy);
    } else if T == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_laed0_gemm_real :: (
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

_laed0_lacpy_real :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           B.data_type == RT;
} {
    #if RT == float32 {
        slacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        dlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_laed0_steqr :: (
    $$compz: u8,
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           E.data_type == RT &&
           Z.data_type == RT &&
           WORK.data_type == RT;
} {
    #if RT == float32 {
        ssteqr(compz, n, d, e, z, ldz, work, info);
    } else {
        dsteqr(compz, n, d, e, z, ldz, work, info);
    }
}

_laed0_laed1 :: (
    $$n: int,
    d: *$D/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    indxq: *$INDXQ/ArrayType1D,
    $$rho: $RT,
    $$cutpnt: int,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           Q.data_type == RT &&
           INDXQ.data_type == int &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    #if RT == float32 {
        slaed1(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info);
    } else {
        dlaed1(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info);
    }
}

_laed0_clacrm :: (
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    rwork: *$RWORK/ArrayType1D,
    $RT: Type,
    $MT: Type
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

_laed0_real :: (
    $$icompq: int,
    $$qsiz: int,
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    qstore: *$QSTORE/ArrayType2D,
    $$ldqs: int,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           E.data_type == RT &&
           Q.data_type == RT &&
           QSTORE.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    routine_name :: #run if RT == float32 { "SLAED0" } else { "DLAED0" };

    curlvl, curprb, curr, i, igivcl, igivnm, igivpt, indxq, iperm, iprmpt, iq, iqptr, iwrem,
    j, k, lgn, matsiz, msd2, smlsiz, smm1, spm1, spm2, submat, subpbs, tlvls : int = ---;

    info.* = 0;
    if icompq < 0 || icompq > 2 {
        info.* = -1;
    } else if icompq == 1 && qsiz < max(0, n) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ldq < max(1, n) {
        info.* = -7;
    } else if ldqs < max(1, n) {
        info.* = -9;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    smlsiz = _laed0_ilaenv(9, routine_name, RT);

    fset(iwork, 1, n);
    subpbs = 1;
    tlvls = 0;
    while fget(iwork.*, subpbs) > smlsiz {
        for <j:1..subpbs {
            fset(iwork, 2*j, (fget(iwork.*, j) + 1)/2);
            fset(iwork, 2*j - 1, fget(iwork.*, j)/2);
        }
        tlvls += 1;
        subpbs *= 2;
    }

    for j:2..subpbs {
        fset(iwork, j, fget(iwork.*, j) + fget(iwork.*, j - 1));
    }

    spm1 = subpbs - 1;
    for i:1..spm1 {
        submat = fget(iwork.*, i) + 1;
        smm1 = submat - 1;
        fset(d, smm1, fget(d.*, smm1) - _laed0_real_abs_e(fget(e.*, smm1), RT));
        fset(d, submat, fget(d.*, submat) - _laed0_real_abs_e(fget(e.*, smm1), RT));
    }

    indxq = 4*n + 3;
    if icompq != 2 {
        lgn = 0;
        k = 1;
        while k < n {
            k *= 2;
            lgn += 1;
        }
        if k < n {
            lgn += 1;
        }

        iprmpt = indxq + n + 1;
        iperm = iprmpt + n*lgn;
        iqptr = iperm + n*lgn;
        igivpt = iqptr + n + 2;
        igivcl = igivpt + n*lgn;

        igivnm = 1;
        iq = igivnm + 2*n*lgn;
        iwrem = iq + n*n + 1;

        for i:0..subpbs {
            fset(iwork, iprmpt + i, 1);
            fset(iwork, igivpt + i, 1);
        }
        fset(iwork, iqptr, 1);
    } else {
        iprmpt = 0; iperm = 0; iqptr = 0; igivpt = 0; igivcl = 0;
        igivnm = 0; iq = 0; iwrem = 0;
    }

    curr = 0;
    for i:0..spm1 {
        if i == 0 {
            submat = 1;
            matsiz = fget(iwork.*, 1);
        } else {
            submat = fget(iwork.*, i) + 1;
            matsiz = fget(iwork.*, i + 1) - fget(iwork.*, i);
        }

        d_sub := _laed0_subvector(d.*, submat, matsiz);
        e_sub := _laed0_subvector(e.*, submat, matsiz);

        if icompq == 2 {
            q_sub_raw := _laed0_matrix_storage(q.*, ldq, submat, submat, matsiz, matsiz);
            q_sub := array2d_view(q_sub_raw, ldq, row_major=false);
            _laed0_steqr(#char "I", matsiz, *d_sub, *e_sub, *q_sub, ldq, work, info, RT);
            if info.* != 0 {
                info.* = submat*(n + 1) + submat + matsiz - 1;
                return;
            }
        } else {
            ll := iq - 1 + fget(iwork.*, iqptr + curr);
            qwork_raw := _laed0_raw_slice(work.data + (ll - 1), matsiz*matsiz);
            qwork := array2d_view(qwork_raw, matsiz, row_major=false);
            _laed0_steqr(#char "I", matsiz, *d_sub, *e_sub, *qwork, matsiz, work, info, RT);
            if info.* != 0 {
                info.* = submat*(n + 1) + submat + matsiz - 1;
                return;
            }

            if icompq == 1 {
                q_dense_raw := _laed0_matrix_storage(q.*, ldq, 1, submat, qsiz, matsiz);
                q_dense := array2d_view(q_dense_raw, ldq, row_major=false);
                qstore_sub_raw := _laed0_matrix_storage(qstore.*, ldqs, 1, submat, qsiz, matsiz);
                qstore_sub := array2d_view(qstore_sub_raw, ldqs, row_major=false);
                _laed0_gemm_real(#char "N", #char "N", qsiz, matsiz, matsiz, one_r, q_dense.*, ldq, qwork.*, matsiz, zero_r, *qstore_sub, ldqs, RT);
            }

            fset(iwork, iqptr + curr + 1, fget(iwork.*, iqptr + curr) + matsiz*matsiz);
            curr += 1;
        }

        k = 1;
        for j:submat..fget(iwork.*, i + 1) {
            fset(iwork, indxq + j, k);
            k += 1;
        }
    }

    curlvl = 1;
    while subpbs > 1 {
        spm2 = subpbs - 2;
        curprb = 0;
        step := 2;
        for i:0..spm2 {
            if mod(i, step) != 0 then continue;

            if i == 0 {
                submat = 1;
                matsiz = fget(iwork.*, 2);
                msd2 = fget(iwork.*, 1);
                curprb = 0;
            } else {
                submat = fget(iwork.*, i) + 1;
                matsiz = fget(iwork.*, i + 2) - fget(iwork.*, i);
                msd2 = matsiz/2;
                curprb += 1;
            }

            d_sub := _laed0_subvector(d.*, submat, matsiz);
            indxq_sub := _laed0_subvector(iwork.*, indxq + submat, matsiz);
            iwork_subpbs := _laed0_subvector(iwork.*, subpbs + 1, max(0, n - subpbs));

            if icompq == 2 {
                q_sub_raw := _laed0_matrix_storage(q.*, ldq, submat, submat, matsiz, matsiz);
                q_sub := array2d_view(q_sub_raw, ldq, row_major=false);
                rho := fget(e.*, submat + msd2 - 1);
                work_sub := _laed0_subvector(work.*, 1, max(0, work.count));
                _laed0_laed1(matsiz, *d_sub, *q_sub, ldq, *indxq_sub, rho, msd2, work, *iwork_subpbs, info, RT);
            } else {
                qstore_sub_raw := _laed0_matrix_storage(qstore.*, ldqs, 1, submat, qsiz, matsiz);
                qstore_sub := array2d_view(qstore_sub_raw, ldqs, row_major=false);
                qstore_vec := _laed0_raw_slice(qstore_sub.data, ifx qsiz <= 0 || matsiz <= 0 then 0 else (matsiz - 1)*ldqs + qsiz);
                iwork_iqptr := _laed0_subvector(iwork.*, iqptr, max(0, iwork.count - iqptr + 1));
                iwork_iprmpt := _laed0_subvector(iwork.*, iprmpt, max(0, iwork.count - iprmpt + 1));
                iwork_iperm := _laed0_subvector(iwork.*, iperm, max(0, iwork.count - iperm + 1));
                iwork_igivpt := _laed0_subvector(iwork.*, igivpt, max(0, iwork.count - igivpt + 1));
                givcol_raw := _laed0_raw_slice(iwork.data + (igivcl - 1), max(0, iwork.count - igivcl + 1));
                givcol_mat := array2d_view(givcol_raw, 2, row_major=false);
                givnum_raw := _laed0_raw_slice(work.data + (igivnm - 1), max(0, work.count - igivnm + 1));
                givnum_mat := array2d_view(givnum_raw, 2, row_major=false);
                work_iq := _laed0_subvector(work.*, iq, max(0, work.count - iq + 1));
                rho := fget(e.*, submat + msd2 - 1);
                _laed7_real(icompq, matsiz, qsiz, tlvls, curlvl, curprb, *d_sub, *qstore_sub, ldqs, *indxq_sub, rho, msd2,
                            *qstore_vec, *iwork_iqptr, *iwork_iprmpt, *iwork_iperm, *iwork_igivpt, *givcol_mat, *givnum_mat,
                            *work_iq, *iwork_subpbs, info, RT);
            }

            if info.* != 0 {
                info.* = submat*(n + 1) + submat + matsiz - 1;
                return;
            }
            fset(iwork, i/2 + 1, fget(iwork.*, i + 2));
        }
        subpbs /= 2;
        curlvl += 1;
    }

    if icompq == 1 {
        for i:1..n {
            j = fget(iwork.*, indxq + i);
            fset(work, i, fget(d.*, j));
            qstore_col := _laed0_matrix_column(qstore.*, ldqs, 1, j, qsiz);
            q_col := _laed0_matrix_column(q.*, ldq, 1, i, qsiz);
            _laed0_copy(qsiz, qstore_col, 1, *q_col, 1, RT);
        }
        _laed0_copy(n, _laed0_subvector(work.*, 1, n), 1, d, 1, RT);
    } else if icompq == 2 {
        for i:1..n {
            j = fget(iwork.*, indxq + i);
            fset(work, i, fget(d.*, j));
            row_copy_dst := _laed0_subvector(work.*, n*i + 1, n);
            _laed0_copy(n, _laed0_matrix_column(q.*, ldq, 1, j, n), 1, *row_copy_dst, 1, RT);
        }
        _laed0_copy(n, _laed0_subvector(work.*, 1, n), 1, d, 1, RT);
        work_mat_raw := _laed0_raw_slice(work.data + n, n*n);
        work_mat := array2d_view(work_mat_raw, n, row_major=false);
        _laed0_lacpy_real(#char "A", n, n, work_mat.*, n, q, ldq, RT);
    } else {
        for i:1..n {
            j = fget(iwork.*, indxq + i);
            fset(work, i, fget(d.*, j));
        }
        _laed0_copy(n, _laed0_subvector(work.*, 1, n), 1, d, 1, RT);
    }
}

_laed0_complex :: (
    $$qsiz: int,
    $$n: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    qstore: *$QSTORE/ArrayType2D,
    $$ldqs: int,
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
           Q.data_type == MT &&
           QSTORE.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CLAED0" } else { "ZLAED0" };

    curlvl, curprb, curr, i, igivcl, igivnm, igivpt, indxq, iperm, iprmpt, iq, iqptr, iwrem,
    j, k, lgn, ll, matsiz, msd2, smlsiz, smm1, spm1, spm2, submat, subpbs, tlvls : int = ---;
    info.* = 0;
    if qsiz < max(0, n) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if ldq < max(1, n) {
        info.* = -6;
    } else if ldqs < max(1, n) {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    smlsiz = _laed0_ilaenv(9, routine_name, RT);

    fset(iwork, 1, n);
    subpbs = 1;
    tlvls = 0;
    while fget(iwork.*, subpbs) > smlsiz {
        for <j:1..subpbs {
            fset(iwork, 2*j, (fget(iwork.*, j) + 1)/2);
            fset(iwork, 2*j - 1, fget(iwork.*, j)/2);
        }
        tlvls += 1;
        subpbs *= 2;
    }

    for j:2..subpbs {
        fset(iwork, j, fget(iwork.*, j) + fget(iwork.*, j - 1));
    }

    spm1 = subpbs - 1;
    for i:1..spm1 {
        submat = fget(iwork.*, i) + 1;
        smm1 = submat - 1;
        fset(d, smm1, fget(d.*, smm1) - abs(fget(e.*, smm1)));
        fset(d, submat, fget(d.*, submat) - abs(fget(e.*, smm1)));
    }

    indxq = 4*n + 3;
    lgn = 0;
    k = 1;
    while k < n {
        k *= 2;
        lgn += 1;
    }
    if k < n {
        lgn += 1;
    }
    iprmpt = indxq + n + 1;
    iperm = iprmpt + n*lgn;
    iqptr = iperm + n*lgn;
    igivpt = iqptr + n + 2;
    igivcl = igivpt + n*lgn;

    igivnm = 1;
    iq = igivnm + 2*n*lgn;
    iwrem = iq + n*n + 1;

    for i:0..subpbs {
        fset(iwork, iprmpt + i, 1);
        fset(iwork, igivpt + i, 1);
    }
    fset(iwork, iqptr, 1);

    curr = 0;
    for i:0..spm1 {
        if i == 0 {
            submat = 1;
            matsiz = fget(iwork.*, 1);
        } else {
            submat = fget(iwork.*, i) + 1;
            matsiz = fget(iwork.*, i + 1) - fget(iwork.*, i);
        }

        d_sub := _laed0_subvector(d.*, submat, matsiz);
        e_sub := _laed0_subvector(e.*, submat, matsiz);
        ll = iq - 1 + fget(iwork.*, iqptr + curr);
        eigvec_raw := _laed0_raw_slice(rwork.data + (ll - 1), matsiz*matsiz);
        eigvec := array2d_view(eigvec_raw, matsiz, row_major=false);
        _laed0_steqr(#char "I", matsiz, *d_sub, *e_sub, *eigvec, matsiz, rwork, info, RT);

        q_sub_raw := _laed0_matrix_storage(q.*, ldq, 1, submat, qsiz, matsiz);
        q_sub := array2d_view(q_sub_raw, ldq, row_major=false);
        qstore_sub_raw := _laed0_matrix_storage(qstore.*, ldqs, 1, submat, qsiz, matsiz);
        qstore_sub := array2d_view(qstore_sub_raw, ldqs, row_major=false);
        rwork_iwrem := _laed0_subvector(rwork.*, iwrem, max(0, rwork.count - iwrem + 1));
        _laed0_clacrm(qsiz, matsiz, q_sub.*, ldq, eigvec.*, matsiz, *qstore_sub, ldqs, rwork_iwrem, RT, MT);
        fset(iwork, iqptr + curr + 1, fget(iwork.*, iqptr + curr) + matsiz*matsiz);
        curr += 1;
        if info.* > 0 {
            info.* = submat*(n + 1) + submat + matsiz - 1;
            return;
        }

        k = 1;
        for j:submat..fget(iwork.*, i + 1) {
            fset(iwork, indxq + j, k);
            k += 1;
        }
    }

    curlvl = 1;
    while subpbs > 1 {
        spm2 = subpbs - 2;
        curprb = 0;
        step := 2;
        for i:0..spm2 {
            if mod(i, step) != 0 then continue;

            if i == 0 {
                submat = 1;
                matsiz = fget(iwork.*, 2);
                msd2 = fget(iwork.*, 1);
                curprb = 0;
            } else {
                submat = fget(iwork.*, i) + 1;
                matsiz = fget(iwork.*, i + 2) - fget(iwork.*, i);
                msd2 = matsiz/2;
                curprb += 1;
            }

            d_sub := _laed0_subvector(d.*, submat, matsiz);
            qstore_sub_raw := _laed0_matrix_storage(qstore.*, ldqs, 1, submat, qsiz, matsiz);
            qstore_sub := array2d_view(qstore_sub_raw, ldqs, row_major=false);
            indxq_sub := _laed0_subvector(iwork.*, indxq + submat, matsiz);
            qstore_real := _laed0_subvector(rwork.*, iq, max(0, rwork.count - iq + 1));
            qptr_sub := _laed0_subvector(iwork.*, iqptr, max(0, iwork.count - iqptr + 1));
            prmptr_sub := _laed0_subvector(iwork.*, iprmpt, max(0, iwork.count - iprmpt + 1));
            perm_sub := _laed0_subvector(iwork.*, iperm, max(0, iwork.count - iperm + 1));
            givptr_sub := _laed0_subvector(iwork.*, igivpt, max(0, iwork.count - igivpt + 1));
            givcol_raw := _laed0_raw_slice(iwork.data + (igivcl - 1), max(0, iwork.count - igivcl + 1));
            givcol_sub := array2d_view(givcol_raw, 2, row_major=false);
            givnum_raw := _laed0_raw_slice(rwork.data + (igivnm - 1), max(0, rwork.count - igivnm + 1));
            givnum_sub := array2d_view(givnum_raw, 2, row_major=false);
            q_work := _laed0_raw_slice(q.data + (submat - 1)*ldq, ifx qsiz <= 0 || matsiz <= 0 then 0 else (matsiz - 1)*ldq + qsiz);
            rwork_iwrem := _laed0_subvector(rwork.*, iwrem, max(0, rwork.count - iwrem + 1));
            rho := fget(e.*, submat + msd2 - 1);
            iwork_subpbs := _laed0_subvector(iwork.*, subpbs + 1, max(0, iwork.count - subpbs));

            _laed7_complex(matsiz, msd2, qsiz, tlvls, curlvl, curprb, *d_sub, *qstore_sub, ldqs, rho, *indxq_sub,
                           *qstore_real, *qptr_sub, *prmptr_sub, *perm_sub, *givptr_sub, *givcol_sub, *givnum_sub,
                           *q_work, *rwork_iwrem, *iwork_subpbs, info, RT, MT);
            if info.* > 0 {
                info.* = submat*(n + 1) + submat + matsiz - 1;
                return;
            }
            fset(iwork, i/2 + 1, fget(iwork.*, i + 2));
        }
        subpbs /= 2;
        curlvl += 1;
    }

    for i:1..n {
        j = fget(iwork.*, indxq + i);
        fset(rwork, i, fget(d.*, j));
        qstore_col := _laed0_matrix_column(qstore.*, ldqs, 1, j, qsiz);
        q_col := _laed0_matrix_column(q.*, ldq, 1, i, qsiz);
        _laed0_copy(qsiz, qstore_col, 1, *q_col, 1, MT);
    }
    _laed0_copy(n, _laed0_subvector(rwork.*, 1, n), 1, d, 1, RT);
}
