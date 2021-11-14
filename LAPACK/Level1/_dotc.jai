
// ? LAPACK *dotc
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga0f02c96fa8498b4aa6b144deee725c0d.html#ga0f02c96fa8498b4aa6b144deee725c0d


_dotc :: ($$n: int, cx: $X/VectorType, $$incx: int, cy: $Y/VectorType, $$incy: int) -> X.T {
    #if CHECKS {
        #assert(X.T == Y.T);
    }

    res := zero(X.T);
    if n < 0 then return res;;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            res += conjugate(get(cx,i))*get(cy,i);
        }
    } else {
        ix := 1;
        iy := 1;
        if incx < 0 then ix = (-n+1)*incx;
        if incy < 0 then iy = (-n+1)*incy;
        for i: 0..n-1 {
            res += conjugate(get(cx,i)) * get(cy,i);
            ix += incx;
            iy += incy;
        }
    }
    return res;
}