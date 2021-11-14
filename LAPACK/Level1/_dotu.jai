

// ? LAPACK *dotc
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_gadd72f1b633553acc77250e32bc704a78.html#gadd72f1b633553acc77250e32bc704a78


_dotu :: ($$n: int, cx: $X/VectorType, $$incx: int, cy: $Y/VectorType, $$incy: int) -> X.T {
    #if CHECKS {
        #assert(X.T == Y.T);
    }

    res := zero(X.T);
    if n < 0 then return res;;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            res += get(cx,i)*get(cy,i);
        }
    } else {
        ix := 1;
        iy := 1;
        if incx < 0 then ix = (-n+1)*incx;
        if incy < 0 then iy = (-n+1)*incy;
        for i: 0..n-1 {
            res += get(cx,i) * get(cy,i);
            ix += incx;
            iy += incy;
        }
    }
    return res;
}