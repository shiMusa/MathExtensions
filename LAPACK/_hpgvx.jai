// ? LAPACK *hpgvx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpgvx/zhpgvx
// *HPGVX computes selected eigenvalues and, optionally, eigenvectors of a
// complex generalized Hermitian-definite eigenproblem in packed storage.

chpgvx :: #bake_arguments _hpgvx(RT=float32, MT=C32);
zhpgvx :: #bake_arguments _hpgvx(RT=float64, MT=C64);

#scope_file
_hpgvx_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld;
    view.count = count;
    return view;
}

_hpgvx_pptrf :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT;
} {
    #if MT == C32 { cpptrf(uplo, n, ap, info); } else { zpptrf(uplo, n, ap, info); }
}

_hpgvx_hpevx :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, vl: $RT, vu: RT,
    $$il: int, $$iu: int, abstol: RT, m: *int, w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && W.data_type == RT && Z.data_type == MT && WORK.data_type == MT &&
           RWORK.data_type == RT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    #if MT == C32 { chpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info); }
    else { zhpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info); }
}

_hpgvx_tpmv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, ap: $AP/ArrayType1D, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { ctpmv(uplo, trans, diag, n, ap, x, incx); } else { ztpmv(uplo, trans, diag, n, ap, x, incx); }
}

_hpgvx_tpsv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, ap: $AP/ArrayType1D, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { ctpsv(uplo, trans, diag, n, ap, x, incx); } else { ztpsv(uplo, trans, diag, n, ap, x, incx); }
}

_hpgvx :: (
    $$itype: int, $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, bp: *$BP/ArrayType1D,
    vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT, m: *int, w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && BP.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    alleig, indeig, upper, valeig, wantz : bool = ---;
    j : int = ---;
    trans : u8 = ---;
    routine_name :: #run if MT == C32 { "CHPGVX" } else { "ZHPGVX" };

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");
    alleig = lsame(range, #char "A");
    valeig = lsame(range, #char "V");
    indeig = lsame(range, #char "I");

    info.* = 0;
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !(wantz || lsame(jobz, #char "N")) {
        info.* = -2;
    } else if !(alleig || valeig || indeig) {
        info.* = -3;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -4;
    } else if n < 0 {
        info.* = -5;
    } else {
        if valeig {
            if n > 0 && vu <= vl {
                info.* = -9;
            }
        } else if indeig {
            if il < 1 {
                info.* = -10;
            } else if iu < min(n, il) || iu > n {
                info.* = -11;
            }
        }
    }
    if info.* == 0 {
        if ldz < 1 || (wantz && ldz < n) {
            info.* = -16;
        }
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    _hpgvx_pptrf(uplo, n, bp, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    _hpgst(itype, uplo, n, ap, bp.*, info, RT, MT);
    _hpgvx_hpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info, RT, MT);

    if wantz {
        if info.* > 0 {
            m.* = info.* - 1;
        }

        if itype == 1 || itype == 2 {
            if upper { trans = #char "N"; } else { trans = #char "C"; }
            for j:1..m.* {
                zcol := _hpgvx_matrix_column(z.*, ldz, j, n);
                _hpgvx_tpsv(uplo, trans, #char "N", n, bp.*, *zcol, 1, MT);
            }
        } else if itype == 3 {
            if upper { trans = #char "C"; } else { trans = #char "N"; }
            for j:1..m.* {
                zcol := _hpgvx_matrix_column(z.*, ldz, j, n);
                _hpgvx_tpmv(uplo, trans, #char "N", n, bp.*, *zcol, 1, MT);
            }
        }
    }
}
