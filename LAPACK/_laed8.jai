// ? LAPACK *laed8()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaed8/claed8/dlaed8/zlaed8
// *LAED8 merges eigenvalues and deflates the secular equation.

slaed8 :: #bake_arguments _laed8_real(RT=float32);
dlaed8 :: #bake_arguments _laed8_real(RT=float64);
claed8 :: #bake_arguments _laed8_complex(RT=float32, MT=C32);
zlaed8 :: #bake_arguments _laed8_complex(RT=float64, MT=C64);

#scope_file
_laed8_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_laed8_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _laed8_raw_slice(x.data + (offset - 1), count);
}

_laed8_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _laed8_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laed8_copy :: ($$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $T: Type) #modify {
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

_laed8_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        slacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == float64 {
        dlacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_laed8_scal_real :: ($$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        sscal(n, alpha, x, incx);
    } else {
        dscal(n, alpha, x, incx);
    }
}

_laed8_iamax_1based :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx) + 1;
    } else {
        return idamax(n, x, incx) + 1;
    }
}

_laed8_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "E");
    } else {
        return dlamch(#char "E");
    }
}

_laed8_lapy2 :: (x: $RT, y: RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slapy2(x, y);
    } else {
        return dlapy2(x, y);
    }
}

_laed8_lamrg :: ($$n1: int, $$n2: int, a: $A/ArrayType1D, $$dtrd1: int, $$dtrd2: int, index: *$INDEX/ArrayType1D, $RT: Type) #modify {
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

_laed8_rot :: ($$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $$c: $RT, $$s: RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else if MT == float64 {
        drot(n, x, incx, y, incy, c, s);
    } else if MT == C32 {
        csrot(n, x, incx, y, incy, c, s);
    } else {
        zdrot(n, x, incx, y, incy, c, s);
    }
}

_laed8_real :: (
    $$icompq: int,
    k: *int,
    $$n: int,
    $$qsiz: int,
    d: *$D/ArrayType1D,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    indxq: *$INDXQ/ArrayType1D,
    rho: *RT,
    $$cutpnt: int,
    z: *$Z/ArrayType1D,
    dlambda: *$DLAMBDA/ArrayType1D,
    q2: *$Q2/ArrayType2D,
    $$ldq2: int,
    w: *$W/ArrayType1D,
    perm: *$PERM/ArrayType1D,
    givptr: *int,
    givcol: *$GIVCOL/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    indxp: *$INDXP/ArrayType1D,
    indx: *$INDX/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT &&
           Q.data_type == RT &&
           INDXQ.data_type == int &&
           Z.data_type == RT &&
           DLAMBDA.data_type == RT &&
           Q2.data_type == RT &&
           W.data_type == RT &&
           PERM.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT &&
           INDXP.data_type == int &&
           INDX.data_type == int;
} {
    mone_r  :: #run -one(RT);
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    two_r   :: #run ncast(RT, 2);
    eight_r :: #run ncast(RT, 8);
    routine_name :: #run if RT == float32 { "SLAED8" } else { "DLAED8" };

    c, eps, s, t, tau, tol : RT = ---;
    found_nondeflated : bool = false;
    i, imax, j, jlam, jmax, jp, k2, n1, n1p1, n2 : int = ---;

    info.* = 0;
    if icompq < 0 || icompq > 1 {
        info.* = -1;
    } else if n < 0 {
        info.* = -3;
    } else if icompq == 1 && qsiz < n {
        info.* = -4;
    } else if ldq < max(1, n) {
        info.* = -7;
    } else if cutpnt < min(1, n) || cutpnt > n {
        info.* = -10;
    } else if ldq2 < max(1, n) {
        info.* = -14;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    givptr.* = 0;
    if n == 0 then return;

    n1 = cutpnt;
    n2 = n - n1;
    n1p1 = n1 + 1;

    if rho.* < zero_r {
        z_tail := _laed8_subvector(z.*, n1p1, n2);
        _laed8_scal_real(n2, mone_r, *z_tail, 1, RT);
    }

    t = one_r / sqrt(two_r);
    for j:1..n {
        fset(indx, j, j);
    }
    _laed8_scal_real(n, t, z, 1, RT);
    rho.* = abs(two_r*rho.*);

    for i:cutpnt + 1..n {
        fset(indxq, i, fget(indxq.*, i) + cutpnt);
    }
    for i:1..n {
        fset(dlambda, i, fget(d.*, fget(indxq.*, i)));
        fset(w, i, fget(z.*, fget(indxq.*, i)));
    }
    _laed8_lamrg(n1, n2, dlambda, 1, 1, indx, RT);
    for i:1..n {
        fset(d, i, fget(dlambda.*, fget(indx.*, i)));
        fset(z, i, fget(w.*, fget(indx.*, i)));
    }

    imax = _laed8_iamax_1based(n, z.*, 1, RT);
    jmax = _laed8_iamax_1based(n, d.*, 1, RT);
    eps = _laed8_lamch_epsilon(RT);
    tol = eight_r*eps*abs(fget(d.*, jmax));

    if rho.* * abs(fget(z.*, imax)) <= tol {
        k.* = 0;
        if icompq == 0 {
            for j:1..n {
                fset(perm, j, fget(indxq.*, fget(indx.*, j)));
            }
        } else {
            for j:1..n {
                fset(perm, j, fget(indxq.*, fget(indx.*, j)));
                q_src := _laed8_matrix_column(q.*, ldq, 1, fget(perm.*, j), qsiz);
                q_dst := _laed8_matrix_column(q2.*, ldq2, 1, j, qsiz);
                _laed8_copy(qsiz, q_src, 1, *q_dst, 1, RT);
            }
            _laed8_lacpy(#char "A", qsiz, n, q2.*, ldq2, q, ldq, RT, RT);
        }
        return;
    }

    k.* = 0;
    k2 = n + 1;
    for j:1..n {
        if rho.* * abs(fget(z.*, j)) <= tol {
            k2 -= 1;
            fset(indxp, k2, j);
        } else {
            jlam = j;
            found_nondeflated = true;
            break;
        }
    }

    if found_nondeflated {
        for j:jlam + 1..n {
            if rho.* * abs(fget(z.*, j)) <= tol {
                k2 -= 1;
                fset(indxp, k2, j);
            } else {
                s = fget(z.*, jlam);
                c = fget(z.*, j);
                tau = _laed8_lapy2(c, s, RT);
                t = fget(d.*, j) - fget(d.*, jlam);
                c /= tau;
                s = -s / tau;
                if abs(t*c*s) <= tol {
                    fset(z, j, tau);
                    fset(z, jlam, zero_r);

                    givptr.* += 1;
                    fset(givcol, 1, givptr.*, fget(indxq.*, fget(indx.*, jlam)));
                    fset(givcol, 2, givptr.*, fget(indxq.*, fget(indx.*, j)));
                    fset(givnum, 1, givptr.*, c);
                    fset(givnum, 2, givptr.*, s);
                    if icompq == 1 {
                        q_jlam := _laed8_matrix_column(q.*, ldq, 1, fget(indxq.*, fget(indx.*, jlam)), qsiz);
                        q_j := _laed8_matrix_column(q.*, ldq, 1, fget(indxq.*, fget(indx.*, j)), qsiz);
                        _laed8_rot(qsiz, *q_jlam, 1, *q_j, 1, c, s, RT, RT);
                    }
                    t = fget(d.*, jlam)*c*c + fget(d.*, j)*s*s;
                    fset(d, j, fget(d.*, jlam)*s*s + fget(d.*, j)*c*c);
                    fset(d, jlam, t);
                    k2 -= 1;
                    i = 1;
                    while k2 + i <= n {
                        if fget(d.*, jlam) < fget(d.*, fget(indxp.*, k2 + i)) {
                            fset(indxp, k2 + i - 1, fget(indxp.*, k2 + i));
                            fset(indxp, k2 + i, jlam);
                            i += 1;
                        } else {
                            fset(indxp, k2 + i - 1, jlam);
                            break;
                        }
                    }
                    if k2 + i > n {
                        fset(indxp, k2 + i - 1, jlam);
                    }
                    jlam = j;
                } else {
                    k.* += 1;
                    fset(w, k.*, fget(z.*, jlam));
                    fset(dlambda, k.*, fget(d.*, jlam));
                    fset(indxp, k.*, jlam);
                    jlam = j;
                }
            }
        }
        k.* += 1;
        fset(w, k.*, fget(z.*, jlam));
        fset(dlambda, k.*, fget(d.*, jlam));
        fset(indxp, k.*, jlam);
    }

    for j:1..n {
        jp = fget(indxp.*, j);
        fset(dlambda, j, fget(d.*, jp));
        fset(perm, j, fget(indxq.*, fget(indx.*, jp)));
        if icompq == 1 {
            q_src := _laed8_matrix_column(q.*, ldq, 1, fget(perm.*, j), qsiz);
            q_dst := _laed8_matrix_column(q2.*, ldq2, 1, j, qsiz);
            _laed8_copy(qsiz, q_src, 1, *q_dst, 1, RT);
        }
    }

    if k.* < n {
        _laed8_copy(n - k.*, _laed8_subvector(dlambda.*, k.* + 1, n - k.*), 1, *_laed8_subvector(d.*, k.* + 1, n - k.*), 1, RT);
        if icompq == 1 {
            q2_tail_raw := _laed8_raw_slice(q2.data + k.*ldq2, ifx qsiz <= 0 || n - k.* <= 0 then 0 else (n - k.* - 1)*ldq2 + qsiz);
            q2_tail := array2d_view(q2_tail_raw, ldq2, row_major=false);
            q_tail_raw := _laed8_raw_slice(q.data + k.*ldq, ifx qsiz <= 0 || n - k.* <= 0 then 0 else (n - k.* - 1)*ldq + qsiz);
            q_tail := array2d_view(q_tail_raw, ldq, row_major=false);
            _laed8_lacpy(#char "A", qsiz, n - k.*, q2_tail.*, ldq2, *q_tail, ldq, RT, RT);
        }
    }
}

_laed8_complex :: (
    k: *int,
    $$n: int,
    $$qsiz: int,
    q: *$Q/ArrayType2D,
    $$ldq: int,
    d: *$D/ArrayType1D,
    rho: *RT,
    $$cutpnt: int,
    z: *$Z/ArrayType1D,
    dlambda: *$DLAMBDA/ArrayType1D,
    q2: *$Q2/ArrayType2D,
    $$ldq2: int,
    w: *$W/ArrayType1D,
    indxp: *$INDXP/ArrayType1D,
    indx: *$INDX/ArrayType1D,
    indxq: *$INDXQ/ArrayType1D,
    perm: *$PERM/ArrayType1D,
    givptr: *int,
    givcol: *$GIVCOL/ArrayType2D,
    givnum: *$GIVNUM/ArrayType2D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           Q.data_type == MT &&
           D.data_type == RT &&
           Z.data_type == RT &&
           DLAMBDA.data_type == RT &&
           Q2.data_type == MT &&
           W.data_type == RT &&
           INDXP.data_type == int &&
           INDX.data_type == int &&
           INDXQ.data_type == int &&
           PERM.data_type == int &&
           GIVCOL.data_type == int &&
           GIVNUM.data_type == RT;
} {
    mone_r  :: #run -one(RT);
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    two_r   :: #run ncast(RT, 2);
    eight_r :: #run ncast(RT, 8);
    routine_name :: #run if MT == C32 { "CLAED8" } else { "ZLAED8" };

    c, eps, s, t, tau, tol : RT = ---;
    found_nondeflated : bool = false;
    i, imax, j, jlam, jmax, jp, k2, n1, n1p1, n2 : int = ---;

    info.* = 0;
    if n < 0 {
        info.* = -2;
    } else if qsiz < n {
        info.* = -3;
    } else if ldq < max(1, n) {
        info.* = -5;
    } else if cutpnt < min(1, n) || cutpnt > n {
        info.* = -8;
    } else if ldq2 < max(1, n) {
        info.* = -12;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    givptr.* = 0;
    if n == 0 then return;

    n1 = cutpnt;
    n2 = n - n1;
    n1p1 = n1 + 1;

    if rho.* < zero_r {
        z_tail := _laed8_subvector(z.*, n1p1, n2);
        _laed8_scal_real(n2, mone_r, *z_tail, 1, RT);
    }

    t = one_r / sqrt(two_r);
    for j:1..n {
        fset(indx, j, j);
    }
    _laed8_scal_real(n, t, z, 1, RT);
    rho.* = abs(two_r*rho.*);

    for i:cutpnt + 1..n {
        fset(indxq, i, fget(indxq.*, i) + cutpnt);
    }
    for i:1..n {
        fset(dlambda, i, fget(d.*, fget(indxq.*, i)));
        fset(w, i, fget(z.*, fget(indxq.*, i)));
    }
    _laed8_lamrg(n1, n2, dlambda, 1, 1, indx, RT);
    for i:1..n {
        fset(d, i, fget(dlambda.*, fget(indx.*, i)));
        fset(z, i, fget(w.*, fget(indx.*, i)));
    }

    imax = _laed8_iamax_1based(n, z.*, 1, RT);
    jmax = _laed8_iamax_1based(n, d.*, 1, RT);
    eps = _laed8_lamch_epsilon(RT);
    tol = eight_r*eps*abs(fget(d.*, jmax));

    if rho.* * abs(fget(z.*, imax)) <= tol {
        k.* = 0;
        for j:1..n {
            fset(perm, j, fget(indxq.*, fget(indx.*, j)));
            q_src := _laed8_matrix_column(q.*, ldq, 1, fget(perm.*, j), qsiz);
            q_dst := _laed8_matrix_column(q2.*, ldq2, 1, j, qsiz);
            _laed8_copy(qsiz, q_src, 1, *q_dst, 1, MT);
        }
        _laed8_lacpy(#char "A", qsiz, n, q2.*, ldq2, q, ldq, RT, MT);
        return;
    }

    k.* = 0;
    k2 = n + 1;
    for j:1..n {
        if rho.* * abs(fget(z.*, j)) <= tol {
            k2 -= 1;
            fset(indxp, k2, j);
        } else {
            jlam = j;
            found_nondeflated = true;
            break;
        }
    }

    if found_nondeflated {
        for j:jlam + 1..n {
            if rho.* * abs(fget(z.*, j)) <= tol {
                k2 -= 1;
                fset(indxp, k2, j);
            } else {
                s = fget(z.*, jlam);
                c = fget(z.*, j);
                tau = _laed8_lapy2(c, s, RT);
                t = fget(d.*, j) - fget(d.*, jlam);
                c /= tau;
                s = -s / tau;
                if abs(t*c*s) <= tol {
                    fset(z, j, tau);
                    fset(z, jlam, zero_r);

                    givptr.* += 1;
                    fset(givcol, 1, givptr.*, fget(indxq.*, fget(indx.*, jlam)));
                    fset(givcol, 2, givptr.*, fget(indxq.*, fget(indx.*, j)));
                    fset(givnum, 1, givptr.*, c);
                    fset(givnum, 2, givptr.*, s);
                    q_jlam := _laed8_matrix_column(q.*, ldq, 1, fget(indxq.*, fget(indx.*, jlam)), qsiz);
                    q_j := _laed8_matrix_column(q.*, ldq, 1, fget(indxq.*, fget(indx.*, j)), qsiz);
                    _laed8_rot(qsiz, *q_jlam, 1, *q_j, 1, c, s, RT, MT);

                    t = fget(d.*, jlam)*c*c + fget(d.*, j)*s*s;
                    fset(d, j, fget(d.*, jlam)*s*s + fget(d.*, j)*c*c);
                    fset(d, jlam, t);
                    k2 -= 1;
                    i = 1;
                    while k2 + i <= n {
                        if fget(d.*, jlam) < fget(d.*, fget(indxp.*, k2 + i)) {
                            fset(indxp, k2 + i - 1, fget(indxp.*, k2 + i));
                            fset(indxp, k2 + i, jlam);
                            i += 1;
                        } else {
                            fset(indxp, k2 + i - 1, jlam);
                            break;
                        }
                    }
                    if k2 + i > n {
                        fset(indxp, k2 + i - 1, jlam);
                    }
                    jlam = j;
                } else {
                    k.* += 1;
                    fset(w, k.*, fget(z.*, jlam));
                    fset(dlambda, k.*, fget(d.*, jlam));
                    fset(indxp, k.*, jlam);
                    jlam = j;
                }
            }
        }
        k.* += 1;
        fset(w, k.*, fget(z.*, jlam));
        fset(dlambda, k.*, fget(d.*, jlam));
        fset(indxp, k.*, jlam);
    }

    for j:1..n {
        jp = fget(indxp.*, j);
        fset(dlambda, j, fget(d.*, jp));
        fset(perm, j, fget(indxq.*, fget(indx.*, jp)));
        q_src := _laed8_matrix_column(q.*, ldq, 1, fget(perm.*, j), qsiz);
        q_dst := _laed8_matrix_column(q2.*, ldq2, 1, j, qsiz);
        _laed8_copy(qsiz, q_src, 1, *q_dst, 1, MT);
    }

    if k.* < n {
        _laed8_copy(n - k.*, _laed8_subvector(dlambda.*, k.* + 1, n - k.*), 1, *_laed8_subvector(d.*, k.* + 1, n - k.*), 1, RT);
        q2_tail_raw := _laed8_raw_slice(q2.data + k.*ldq2, ifx qsiz <= 0 || n - k.* <= 0 then 0 else (n - k.* - 1)*ldq2 + qsiz);
        q2_tail := array2d_view(q2_tail_raw, ldq2, row_major=false);
        q_tail_raw := _laed8_raw_slice(q.data + k.*ldq, ifx qsiz <= 0 || n - k.* <= 0 then 0 else (n - k.* - 1)*ldq + qsiz);
        q_tail := array2d_view(q_tail_raw, ldq, row_major=false);
        _laed8_lacpy(#char "A", qsiz, n - k.*, q2_tail.*, ldq2, *q_tail, ldq, RT, MT);
    }
}
