// ? LAPACK *hsein()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chsein/zhsein
// *HSEIN uses inverse iteration to find specified right and/or left
// eigenvectors of a complex upper Hessenberg matrix.

chsein :: #bake_arguments _hsein(RT=float32, MT=C32);
zhsein :: #bake_arguments _hsein(RT=float64, MT=C64);

#scope_file
_hsein_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hsein_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hsein_raw_slice(x.data + offset - 1, count);
}

_hsein_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hsein_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hsein_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _hsein_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hsein_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hsein_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_hsein_lamch :: ($$cmach: u8, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(cmach); } else { return dlamch(cmach); }
}

_hsein_lanhs :: (
    $$norm: u8, $$n: int, h: $H/ArrayType2D, $$ldh: int, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT &&
           WORK.data_type == RT;
} {
    #if MT == C32 { return clanhs(norm, n, h, ldh, work); } else { return zlanhs(norm, n, h, ldh, work); }
}

_hsein_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sisnan(x); } else { return disnan(x); }
}

_hsein_laein :: (
    $$rightv: bool, $$noinit: bool, $$n: int, h: $H/ArrayType2D, $$ldh: int, w: $MT,
    v: *$V/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int, rwork: *$RWORK/ArrayType1D,
    eps3: $RT, smlnum: RT, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT &&
           V.data_type == MT &&
           B.data_type == MT &&
           RWORK.data_type == RT;
} {
    #if MT == C32 {
        claein(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info);
    } else {
        zlaein(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info);
    }
}

_hsein :: (
    $$side: u8, $$eigsrc: u8, $$initv: u8, select: $SELECT/ArrayType1D, $$n: int,
    h: $H/ArrayType2D, $$ldh: int, w: *$W/ArrayType1D, vl: *$VL/ArrayType2D, $$ldvl: int,
    vr: *$VR/ArrayType2D, $$ldvr: int, $$mm: int, m: *int, work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D, ifaill: *$IFAILL/ArrayType1D, ifailr: *$IFAILR/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           SELECT.data_type == bool &&
           H.data_type == MT &&
           W.data_type == MT &&
           VL.data_type == MT &&
           VR.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IFAILL.data_type == int &&
           IFAILR.data_type == int;
} {
    zero_t :: #run zero(MT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == C32 { "CHSEIN" } else { "ZHSEIN" };

    bothv, fromqr, leftv, noinit, rightv : bool = ---;
    i, iinfo, k, kl, kln, kr, ks, ldwork : int = ---;
    eps3, hnorm, smlnum, ulp, unfl : RT = ---;
    wk : MT = ---;

    bothv = lsame(side, #char "B");
    rightv = lsame(side, #char "R") || bothv;
    leftv = lsame(side, #char "L") || bothv;
    fromqr = lsame(eigsrc, #char "Q");
    noinit = lsame(initv, #char "N");

    m.* = 0;
    for k:1..n {
        if fget(select, k) {
            m.* += 1;
        }
    }

    info.* = 0;
    if !(rightv || leftv) {
        info.* = -1;
    } else if !(fromqr || lsame(eigsrc, #char "N")) {
        info.* = -2;
    } else if !(noinit || lsame(initv, #char "U")) {
        info.* = -3;
    } else if n < 0 {
        info.* = -5;
    } else if ldh < max(1, n) {
        info.* = -7;
    } else if ldvl < 1 || (leftv && ldvl < n) {
        info.* = -10;
    } else if ldvr < 1 || (rightv && ldvr < n) {
        info.* = -12;
    } else if mm < m.* {
        info.* = -13;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    unfl = _hsein_lamch(#char "S", RT);
    ulp = _hsein_lamch(#char "P", RT);
    smlnum = unfl*(ncast(RT, n) / ulp);

    ldwork = n;
    work_storage := _hsein_subvector(work.*, 1, n*n);
    work_mat := array2d_view(work_storage, ldwork, row_major=false);

    kl = 1;
    kln = 0;
    if fromqr {
        kr = 0;
    } else {
        kr = n;
    }
    ks = 1;

    for k:1..n {
        if !fget(select, k) continue;

        if fromqr {
            i = k;
            while i >= kl + 1 {
                if fget(h, i, i - 1) == zero_t {
                    break;
                }
                i -= 1;
            }
            kl = i;

            if k > kr {
                i = k;
                while i <= n - 1 {
                    if fget(h, i + 1, i) == zero_t {
                        break;
                    }
                    i += 1;
                }
                kr = i;
            }
        }

        if kl != kln {
            kln = kl;

            hsub_raw := _hsein_matrix_storage(h, ldh, kl, kl, kr - kl + 1, kr - kl + 1);
            hsub := array2d_view(hsub_raw, ldh, row_major=false);
            hnorm = _hsein_lanhs(#char "I", kr - kl + 1, hsub, ldh, rwork, RT, MT);
            if _hsein_isnan(hnorm, RT) {
                info.* = -6;
                return;
            } else if hnorm > zero_r {
                eps3 = hnorm*ulp;
            } else {
                eps3 = smlnum;
            }
        }

        wk = fget(w.*, k);
        while true {
            i = k - 1;
            while i >= kl {
                if fget(select, i) && _hsein_abs1(fget(w.*, i) - wk, RT, MT) < eps3 {
                    wk += _hsein_complex_from_real(eps3, RT, MT);
                    break;
                }
                i -= 1;
            }
            if i < kl then break;
        }
        fset(w, k, wk);

        if leftv {
            hleft_raw := _hsein_matrix_storage(h, ldh, kl, kl, n - kl + 1, n - kl + 1);
            hleft := array2d_view(hleft_raw, ldh, row_major=false);
            vl_col := _hsein_matrix_column(vl.*, ldvl, kl, ks, n - kl + 1);
            _hsein_laein(false, noinit, n - kl + 1, hleft, ldh, wk, *vl_col, *work_mat, ldwork, rwork, eps3, smlnum, *iinfo, RT, MT);
            if iinfo > 0 {
                info.* += 1;
                fset(ifaill, ks, k);
            } else {
                fset(ifaill, ks, 0);
            }
            for i:1..kl - 1 {
                fset(vl, i, ks, zero_t);
            }
        }

        if rightv {
            hright_raw := _hsein_matrix_storage(h, ldh, 1, 1, kr, kr);
            hright := array2d_view(hright_raw, ldh, row_major=false);
            vr_col := _hsein_matrix_column(vr.*, ldvr, 1, ks, kr);
            _hsein_laein(true, noinit, kr, hright, ldh, wk, *vr_col, *work_mat, ldwork, rwork, eps3, smlnum, *iinfo, RT, MT);
            if iinfo > 0 {
                info.* += 1;
                fset(ifailr, ks, k);
            } else {
                fset(ifailr, ks, 0);
            }
            for i:kr + 1..n {
                fset(vr, i, ks, zero_t);
            }
        }

        ks += 1;
    }
}
