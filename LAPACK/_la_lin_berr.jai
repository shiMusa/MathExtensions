// ? LAPACK *la_lin_berr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sla_lin_berr/dla_lin_berr/cla_lin_berr/zla_lin_berr
// *LA_LIN_BERR computes a component-wise relative backward error.

sla_lin_berr :: #bake_arguments _la_lin_berr(RT=float32, MT=float32);
dla_lin_berr :: #bake_arguments _la_lin_berr(RT=float64, MT=float64);
cla_lin_berr :: #bake_arguments _la_lin_berr(RT=float32, MT=C32);
zla_lin_berr :: #bake_arguments _la_lin_berr(RT=float64, MT=C64);

#scope_file
_la_lin_berr_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_la_lin_berr_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "S"); } else { return dlamch(#char "S"); }
}

_la_lin_berr :: (
    $$n: int, $$nz: int, $$nrhs: int,
    res: $RES/ArrayType1D, ayb: $AYB/ArrayType1D, berr: *$BERR/ArrayType1D,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           RES.data_type == MT &&
           AYB.data_type == RT &&
           BERR.data_type == RT;
} {
    i, j : int = ---;
    safe1, tmp : RT = ---;

    safe1 = _la_lin_berr_lamch_safe_minimum(RT);
    safe1 = ncast(RT, nz + 1) * safe1;

    for j:1..nrhs {
        fset(berr, j, 0);
        for i:1..n {
            idx := i + (j - 1)*n;
            if fget(ayb, idx) != 0 {
                tmp = (safe1 + _la_lin_berr_abs1(fget(res, idx), RT, MT)) / fget(ayb, idx);
                fset(berr, j, max(fget(berr.*, j), tmp));
            }
        }
    }
}
