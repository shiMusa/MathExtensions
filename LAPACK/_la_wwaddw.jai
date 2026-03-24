// ? LAPACK *la_wwaddw()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sla_wwaddw/dla_wwaddw/cla_wwaddw/zla_wwaddw
// *LA_WWADDW adds a vector into a doubled-single vector.

sla_wwaddw :: #bake_arguments _la_wwaddw(RT=float32, MT=float32);
dla_wwaddw :: #bake_arguments _la_wwaddw(RT=float64, MT=float64);
cla_wwaddw :: #bake_arguments _la_wwaddw(RT=float32, MT=C32);
zla_wwaddw :: #bake_arguments _la_wwaddw(RT=float64, MT=C64);

#scope_file
_la_wwaddw :: (
    $$n: int,
    x: *$X/ArrayType1D,
    y: *$Y/ArrayType1D,
    w: $W/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT &&
           W.data_type == MT;
} {
    s : MT = ---;
    i : int = ---;

    for i:1..n {
        s = fget(x.*, i) + fget(w, i);
        s = (s + s) - s;
        fset(y, i, ((fget(x.*, i) - s) + fget(w, i)) + fget(y.*, i));
        fset(x, i, s);
    }
}
