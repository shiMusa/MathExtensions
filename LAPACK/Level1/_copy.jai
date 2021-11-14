
// ? LAPACK *copy
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga6113a670d3df40b1b081af52af8c29e1.html#ga6113a670d3df40b1b081af52af8c29e1

_copy :: (n: int, cx: $X/VectorType, incx: int, cy: *$Y/VectorType, incy: int) {
    #if CHECKS {
        #assert(X.T == Y.T);
    }

    if n <= 0 then return;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            set(cy,i,get(cx,i));
        }
    } else {
        ix := 1;
        iy := 1;
        if incx < 0 then ix = (-n+1)*incx + 1;
        if incy < 0 then iy = (-n+1)*incy + 1;
        for i: 0..n-1 {
            set(cy,i, get(cx,i));
            ix += incx;
            iy += incy;
        }
    }
}