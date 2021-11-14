

// ? LAPACK *scal
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_gab2d569bbfe6356599c69c057ecc2b3f7.html#gab2d569bbfe6356599c69c057ecc2b3f7

_scal :: ($$n: int, ca: $T, cx: *$X/VectorType, $$incx: int) {
    if n <= 0 then return;
    if incx == 1 {
        for i: 0..n-1 {
            set(cx,i, ca*get(cx,i));
        }
    } else {
        ix := 0;
        if incx < 0 then ix = (-n+1)*incx;
        for i: 0..n-1 {
            set(cx,ix, ca*get(cx,ix));
            ix += incx;
        }
    }
}