
// ? LAPACK *axpy()
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga9605cb98791e2038fd89aaef63a31be1.html#ga9605cb98791e2038fd89aaef63a31be1

_axpy :: ($$n: int, ca: $T, cx: $X/VectorType, $$incx: int, cy: *$Y/VectorType, $$incy: int) {
    #if CHECKS {
        #assert(T == X.T && T == Y.T);
    }

    if n <= 0 then return;
    if ca == zero(T) then return;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            set(cy, i, get(cy, i) + ca * get(cx,i));
        }
    } else {
        ix := 1;
        iy := 1;
        if incx < 0 then ix = (-n+1)*incx;
        if incy < 0 then iy = (-n+1)*incy;
        for i: 0..n-1 {
            set(cy,i, get(cy,i) + ca*get(cx,i));
            ix += incx;
            iy += incy;
        }
    }
}