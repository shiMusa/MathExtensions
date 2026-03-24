// ? LAPACK *la_gbrcond_c()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cla_gbrcond_c/zla_gbrcond_c
// *LA_GBRCOND_C computes the infinity norm condition number of
// op(A) * inv(diag(C)) for general banded matrices.

cla_gbrcond_c :: #bake_arguments _la_gbrcond_c(RT=float32, MT=C32);
zla_gbrcond_c :: #bake_arguments _la_gbrcond_c(RT=float64, MT=C64);

#scope_file
_la_gbrcond_c_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_la_gbrcond_c_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_la_gbrcond_c_lacn2 :: (
    $$n: int, v: *$V/ArrayType1D, x: *$X/ArrayType1D, est: *$RT, kase: *int, isave: *[3]int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           V.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 {
        clacn2(n, v, x, est, kase, isave);
    } else {
        zlacn2(n, v, x, est, kase, isave);
    }
}

_la_gbrcond_c_gbtrs :: (
    $$trans: u8, $$n: int, $$kl: int, $$ku: int, afb: $AFB/ArrayType2D, $$ldafb: int,
    ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D, $$ldb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           AFB.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 {
        cgbtrs(trans, n, kl, ku, 1, afb, ldafb, ipiv, work, ldb, info);
    } else {
        zgbtrs(trans, n, kl, ku, 1, afb, ldafb, ipiv, work, ldb, info);
    }
}

_la_gbrcond_c :: (
    $$trans: u8, $$n: int, $$kl: int, $$ku: int, ab: $AB/ArrayType2D, $$ldab: int,
    afb: $AFB/ArrayType2D, $$ldafb: int, ipiv: $IPIV/ArrayType1D, c: $C/ArrayType1D,
    $$capply: bool, info: *int, work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AB.data_type == MT &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           C.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CLA_GBRCOND_C" } else { "ZLA_GBRCOND_C" };

    notrans : bool = ---;
    kase, i, j, kd, ke : int = ---;
    ainvnm, anorm, tmp : RT = ---;
    isave: [3] int;

    result : RT = 0;

    info.* = 0;
    notrans = lsame(trans, #char "N");
    if !notrans && !lsame(trans, #char "T") && !lsame(trans, #char "C") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if kl < 0 || kl > n - 1 {
        info.* = -3;
    } else if ku < 0 || ku > n - 1 {
        info.* = -4;
    } else if ldab < kl + ku + 1 {
        info.* = -6;
    } else if ldafb < 2*kl + ku + 1 {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return result;
    }

    anorm = 0;
    kd = ku + 1;
    ke = kl + 1;
    if notrans {
        for i:1..n {
            tmp = 0;
            if capply {
                for j:max(i - kl, 1)..min(i + ku, n) {
                    tmp += _la_gbrcond_c_abs1(fget(ab, kd + i - j, j), RT, MT) / fget(c, j);
                }
            } else {
                for j:max(i - kl, 1)..min(i + ku, n) {
                    tmp += _la_gbrcond_c_abs1(fget(ab, kd + i - j, j), RT, MT);
                }
            }
            fset(rwork, i, tmp);
            anorm = max(anorm, tmp);
        }
    } else {
        for i:1..n {
            tmp = 0;
            if capply {
                for j:max(i - kl, 1)..min(i + ku, n) {
                    tmp += _la_gbrcond_c_abs1(fget(ab, ke - i + j, i), RT, MT) / fget(c, j);
                }
            } else {
                for j:max(i - kl, 1)..min(i + ku, n) {
                    tmp += _la_gbrcond_c_abs1(fget(ab, ke - i + j, i), RT, MT);
                }
            }
            fset(rwork, i, tmp);
            anorm = max(anorm, tmp);
        }
    }

    if n == 0 {
        return 1;
    } else if anorm == 0 {
        return result;
    }

    ainvnm = 0;
    kase = 0;
    while true {
        work_v := _la_gbrcond_c_subvector(work.*, n + 1, n);
        work_x := _la_gbrcond_c_subvector(work.*, 1, n);
        _la_gbrcond_c_lacn2(n, *work_v, *work_x, *ainvnm, *kase, *isave, RT, MT);
        if kase == 0 then break;

        if kase == 2 {
            for i:1..n {
                fset(work, i, fget(work.*, i) * fget(rwork, i));
            }

            if notrans {
                _la_gbrcond_c_gbtrs(#char "N", n, kl, ku, afb, ldafb, ipiv, work, n, info, MT);
            } else {
                _la_gbrcond_c_gbtrs(#char "C", n, kl, ku, afb, ldafb, ipiv, work, n, info, MT);
            }

            if capply {
                for i:1..n {
                    fset(work, i, fget(work.*, i) * fget(c, i));
                }
            }
        } else {
            if capply {
                for i:1..n {
                    fset(work, i, fget(work.*, i) * fget(c, i));
                }
            }

            if notrans {
                _la_gbrcond_c_gbtrs(#char "C", n, kl, ku, afb, ldafb, ipiv, work, n, info, MT);
            } else {
                _la_gbrcond_c_gbtrs(#char "N", n, kl, ku, afb, ldafb, ipiv, work, n, info, MT);
            }

            for i:1..n {
                fset(work, i, fget(work.*, i) * fget(rwork, i));
            }
        }
    }

    if ainvnm != 0 {
        result = 1 / ainvnm;
    }
    return result;
}
