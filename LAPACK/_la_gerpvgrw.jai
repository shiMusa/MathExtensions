// ? LAPACK *la_gerpvgrw()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sla_gerpvgrw/cla_gerpvgrw/dla_gerpvgrw/zla_gerpvgrw
// *LA_GERPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U)
// for a general matrix.

sla_gerpvgrw :: #bake_arguments _la_gerpvgrw(RT=float32, MT=float32);
dla_gerpvgrw :: #bake_arguments _la_gerpvgrw(RT=float64, MT=float64);
cla_gerpvgrw :: #bake_arguments _la_gerpvgrw(RT=float32, MT=C32);
zla_gerpvgrw :: #bake_arguments _la_gerpvgrw(RT=float64, MT=C64);

#scope_file
_la_gerpvgrw_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_la_gerpvgrw :: (
    $$n: int, $$ncols: int, a: $A/ArrayType2D, $$lda: int,
    af: $AF/ArrayType2D, $$ldaf: int, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           AF.data_type == MT;
} {
    i, j : int = ---;
    amax, rpvgrw, umax : RT = ---;

    rpvgrw = 1;
    for j:1..ncols {
        amax = 0;
        umax = 0;
        for i:1..n {
            amax = max(_la_gerpvgrw_abs1(fget(a, i, j), RT, MT), amax);
        }
        for i:1..j {
            umax = max(_la_gerpvgrw_abs1(fget(af, i, j), RT, MT), umax);
        }
        if umax != 0 {
            rpvgrw = min(amax / umax, rpvgrw);
        }
    }
    return rpvgrw;
}
