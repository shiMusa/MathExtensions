// ? LAPACK *hpgv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpgv/zhpgv
// *HPGV computes all eigenvalues and, optionally, eigenvectors of a
// complex generalized Hermitian-definite eigenproblem in packed storage.

chpgv :: #bake_arguments _hpgv(RT=float32, MT=C32);
zhpgv :: #bake_arguments _hpgv(RT=float64, MT=C64);

#scope_file
_hpgv_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld;
    view.count = count;
    return view;
}

_hpgv_pptrf :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT;
} {
    #if MT == C32 { cpptrf(uplo, n, ap, info); } else { zpptrf(uplo, n, ap, info); }
}

_hpgv_hpev :: (
    $$jobz: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { chpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info); }
    else { zhpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info); }
}

_hpgv_tpmv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, ap: $AP/ArrayType1D, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { ctpmv(uplo, trans, diag, n, ap, x, incx); } else { ztpmv(uplo, trans, diag, n, ap, x, incx); }
}

_hpgv_tpsv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, ap: $AP/ArrayType1D, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { ctpsv(uplo, trans, diag, n, ap, x, incx); } else { ztpsv(uplo, trans, diag, n, ap, x, incx); }
}

_hpgv :: (
    $$itype: int, $$jobz: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, bp: *$BP/ArrayType1D,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && BP.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    upper, wantz : bool = ---;
    j, neig : int = ---;
    trans : u8 = ---;
    routine_name :: #run if MT == C32 { "CHPGV " } else { "ZHPGV " };

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");

    info.* = 0;
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !(wantz || lsame(jobz, #char "N")) {
        info.* = -2;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if ldz < 1 || (wantz && ldz < n) {
        info.* = -9;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    _hpgv_pptrf(uplo, n, bp, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    _hpgst(itype, uplo, n, ap, bp.*, info, RT, MT);
    _hpgv_hpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info, RT, MT);

    if wantz {
        neig = n;
        if info.* > 0 { neig = info.* - 1; }

        if itype == 1 || itype == 2 {
            if upper { trans = #char "N"; } else { trans = #char "C"; }
            for j:1..neig {
                zcol := _hpgv_matrix_column(z.*, ldz, j, n);
                _hpgv_tpsv(uplo, trans, #char "N", n, bp.*, *zcol, 1, MT);
            }
        } else if itype == 3 {
            if upper { trans = #char "C"; } else { trans = #char "N"; }
            for j:1..neig {
                zcol := _hpgv_matrix_column(z.*, ldz, j, n);
                _hpgv_tpmv(uplo, trans, #char "N", n, bp.*, *zcol, 1, MT);
            }
        }
    }
}
