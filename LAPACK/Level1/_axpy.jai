
// ? LAPACK *axpy()
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga9605cb98791e2038fd89aaef63a31be1.html#ga9605cb98791e2038fd89aaef63a31be1

caxpy :: _axpy;
zaxpy :: _axpy;
daxpy :: _axpy;
saxpy :: _axpy;
_axpy :: ($$n: int, ca: $T, cx: $X/ArrayType1D, $$incx: int, cy: *$Y/ArrayType1D, $$incy: int) {
    #assert(is_minor_number_type(Y.data_type, T));
    #assert(is_minor_number_type(Y.data_type, X.data_type));

    if n <= 0 then return;
    if ca == zero(T) then return;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            println("set(%,%) = %", <<cy, i, get(cy,i) + ca*get(cx,i));
            set(cy,i, get(cy,i) + ca * get(cx,i));
            println("<<cy = %", <<cy);
        }
    } else {
        ix := 0;
        iy := 0;
        if incx < 0 then ix = (-n+1)*incx;
        if incy < 0 then iy = (-n+1)*incy;
        for i: 0..n-1 {
            set(cy, iy, get(cy,iy) + ca * get(cx,ix));
            ix += incx;
            iy += incy;
        }
    }
}

#scope_module
test_lapack_axpy :: () {
    println_push("_axpy", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        x := array1d(int.[1,2,3,4,5,6]);
        y := array1d(int.[2,3,4,5,6,7]);
        println("y = %", y);
        _axpy(6, -2, x, 1, *y, 1);
        println("y = %", y);
        assert(y == array1d(0,-1,-2,-3,-4,-5));
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        x := array1d(int.[1,2,3,4,5,6]);
        y := array1d(int.[2,3,4,5,6,7]);
        println("y = %", y);
        _axpy(3, -2, x, 2, *y, 1);
        println("y = %", y);
        assert(y == array1d(0,-3,-6,5,6,7));
    }
}