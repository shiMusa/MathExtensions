// ? LAPACK *la_porcond_x()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cla_porcond_x/zla_porcond_x
// *LA_PORCOND_X computes the infinity norm condition number of
// op(A) * diag(X) for Hermitian positive-definite matrices.

cla_porcond_x :: #bake_arguments _la_porcond_x(RT=float32, MT=C32);
zla_porcond_x :: #bake_arguments _la_porcond_x(RT=float64, MT=C64);

#scope_file
_la_porcond_x_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_la_porcond_x_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_la_porcond_x_lacn2 :: (
    $$n: int, v: *$V/ArrayType1D, x: *$X/ArrayType1D, est: *$RT, kase: *int, isave: *[3]int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           V.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { clacn2(n, v, x, est, kase, isave); } else { zlacn2(n, v, x, est, kase, isave); }
}

_la_porcond_x_potrs :: (
    $$uplo: u8, $$n: int, af: $AF/ArrayType2D, $$ldaf: int,
    work: *$WORK/ArrayType1D, $$ldb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AF.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 { cpotrs(uplo, n, 1, af, ldaf, work, ldb, info); } else { zpotrs(uplo, n, 1, af, ldaf, work, ldb, info); }
}

_la_porcond_x :: (
    $$uplo: u8, $$n: int, a: $A/ArrayType2D, $$lda: int,
    af: $AF/ArrayType2D, $$ldaf: int, x: $X/ArrayType1D,
    info: *int, work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AF.data_type == MT &&
           X.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 { "CLA_PORCOND_X" } else { "ZLA_PORCOND_X" };

    up, upper : bool = ---;
    kase, i, j : int = ---;
    ainvnm, anorm, tmp : RT = ---;
    isave: [3] int;

    result : RT = 0;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !upper && !lsame(uplo, #char "L") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if ldaf < max(1, n) {
        info.* = -6;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return result;
    }
    up = false;
    if lsame(uplo, #char "U") { up = true; }

    anorm = 0;
    if up {
        for i:1..n {
            tmp = 0;
            for j:1..i {
                tmp += _la_porcond_x_abs1(fget(a, j, i) * fget(x, j), RT, MT);
            }
            for j:i + 1..n {
                tmp += _la_porcond_x_abs1(fget(a, i, j) * fget(x, j), RT, MT);
            }
            fset(rwork, i, tmp);
            anorm = max(anorm, tmp);
        }
    } else {
        for i:1..n {
            tmp = 0;
            for j:1..i {
                tmp += _la_porcond_x_abs1(fget(a, i, j) * fget(x, j), RT, MT);
            }
            for j:i + 1..n {
                tmp += _la_porcond_x_abs1(fget(a, j, i) * fget(x, j), RT, MT);
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
        work_v := _la_porcond_x_subvector(work.*, n + 1, n);
        work_x := _la_porcond_x_subvector(work.*, 1, n);
        _la_porcond_x_lacn2(n, *work_v, *work_x, *ainvnm, *kase, *isave, RT, MT);
        if kase == 0 then break;

        if kase == 2 {
            for i:1..n {
                fset(work, i, fget(work.*, i) * fget(rwork, i));
            }

            if up {
                _la_porcond_x_potrs(#char "U", n, af, ldaf, work, n, info, MT);
            } else {
                _la_porcond_x_potrs(#char "L", n, af, ldaf, work, n, info, MT);
            }

            for i:1..n {
                fset(work, i, fget(work.*, i) / fget(x, i));
            }
        } else {
            for i:1..n {
                fset(work, i, fget(work.*, i) / fget(x, i));
            }

            if up {
                _la_porcond_x_potrs(#char "U", n, af, ldaf, work, n, info, MT);
            } else {
                _la_porcond_x_potrs(#char "L", n, af, ldaf, work, n, info, MT);
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
