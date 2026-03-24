// ? LAPACK *la_gbrpvgrw()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sla_gbrpvgrw/cla_gbrpvgrw/dla_gbrpvgrw/zla_gbrpvgrw
// *LA_GBRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U)
// for a general banded matrix.

sla_gbrpvgrw :: #bake_arguments _la_gbrpvgrw(RT=float32, MT=float32);
dla_gbrpvgrw :: #bake_arguments _la_gbrpvgrw(RT=float64, MT=float64);
cla_gbrpvgrw :: #bake_arguments _la_gbrpvgrw(RT=float32, MT=C32);
zla_gbrpvgrw :: #bake_arguments _la_gbrpvgrw(RT=float64, MT=C64);

#scope_file
_la_gbrpvgrw_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_la_gbrpvgrw :: (
    $$n: int, $$kl: int, $$ku: int, $$ncols: int, ab: $AB/ArrayType2D, $$ldab: int,
    afb: $AFB/ArrayType2D, $$ldafb: int, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT;
} {
    i, j, kd : int = ---;
    amax, rpvgrw, umax : RT = ---;

    rpvgrw = 1;
    kd = ku + 1;
    for j:1..ncols {
        amax = 0;
        umax = 0;
        for i:max(j - ku, 1)..min(j + kl, n) {
            amax = max(_la_gbrpvgrw_abs1(fget(ab, kd + i - j, j), RT, MT), amax);
        }
        for i:max(j - ku, 1)..j {
            umax = max(_la_gbrpvgrw_abs1(fget(afb, kd + i - j, j), RT, MT), umax);
        }
        if umax != 0 {
            rpvgrw = min(amax / umax, rpvgrw);
        }
    }
    return rpvgrw;
}
