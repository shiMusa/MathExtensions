// ? LAPACK *larscl2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarscl2/dlarscl2/clarscl2/zlarscl2
// *LARSCL2 performs reciprocal diagonal scaling on a matrix.

slarscl2 :: #bake_arguments _larscl2(RT=float32, MT=float32);
dlarscl2 :: #bake_arguments _larscl2(RT=float64, MT=float64);
clarscl2 :: #bake_arguments _larscl2(RT=float32, MT=C32);
zlarscl2 :: #bake_arguments _larscl2(RT=float64, MT=C64);

_larscl2 :: (
    $$m: int, $$n: int, d: $D/ArrayType1D, x: *$X/ArrayType2D, $$ldx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           D.data_type == RT && X.data_type == MT;
} {
    i, j : int = ---;
    for j:1..n {
        for i:1..m {
            fset(x, i, j, fget(x.*, i, j) / fget(d, i));
        }
    }
}
