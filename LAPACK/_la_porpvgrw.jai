// ? LAPACK *la_porpvgrw()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cla_porpvgrw/zla_porpvgrw
// *LA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U)
// for a Hermitian positive-definite matrix.

cla_porpvgrw :: #bake_arguments _la_porpvgrw(RT=float32, MT=C32);
zla_porpvgrw :: #bake_arguments _la_porpvgrw(RT=float64, MT=C64);

#scope_file
_la_porpvgrw_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_la_porpvgrw :: (
    $$uplo: u8, $$ncols: int, a: $A/ArrayType2D, $$lda: int,
    af: $AF/ArrayType2D, $$ldaf: int, work: *$WORK/ArrayType1D,
    $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AF.data_type == MT &&
           WORK.data_type == RT;
} {
    i, j : int = ---;
    amax, rpvgrw, umax : RT = ---;
    upper : bool = ---;

    upper = lsame(uplo, #char "U");
    rpvgrw = 1;
    for i:1..2*ncols {
        fset(work, i, 0);
    }

    if upper {
        for j:1..ncols {
            for i:1..j {
                fset(work, ncols + j, max(_la_porpvgrw_abs1(fget(a, i, j), RT, MT), fget(work.*, ncols + j)));
            }
        }
    } else {
        for j:1..ncols {
            for i:j..ncols {
                fset(work, ncols + j, max(_la_porpvgrw_abs1(fget(a, i, j), RT, MT), fget(work.*, ncols + j)));
            }
        }
    }

    if upper {
        for j:1..ncols {
            for i:1..j {
                fset(work, j, max(_la_porpvgrw_abs1(fget(af, i, j), RT, MT), fget(work.*, j)));
            }
        }
    } else {
        for j:1..ncols {
            for i:j..ncols {
                fset(work, j, max(_la_porpvgrw_abs1(fget(af, i, j), RT, MT), fget(work.*, j)));
            }
        }
    }

    for i:1..ncols {
        umax = fget(work.*, i);
        amax = fget(work.*, ncols + i);
        if umax != 0 {
            rpvgrw = min(amax / umax, rpvgrw);
        }
    }

    return rpvgrw;
}
