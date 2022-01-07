
// ? LAPACK *copy
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga6113a670d3df40b1b081af52af8c29e1.html#ga6113a670d3df40b1b081af52af8c29e1
ccopy :: _copy;
zcopy :: _copy;
dcopy :: _copy;
scopy :: _copy;
_copy :: ($$n: int, cx: $X/ArrayType1D, $$incx: int, cy: *$Y/ArrayType1D, $$incy: int) {
    #run assert(X.data_type == Y.data_type);

    if n <= 0 then return;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            set(cy, i, get(cx, i));
        }
    } else {
        ix := 0;
        iy := 0;
        if incx < 0 then ix = (n-1)*(-incx);
        if incy < 0 then iy = (n-1)*(-incy);
        for i: 0..n-1 {
            set(cy, iy, get(cx, ix));
            ix += incx;
            iy += incy;
        }
    }
}


#scope_module
test_lapack_copy :: () {
    println_push("_copy", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        x := array1d(6,  1,2,3,4,5,6);
        y := array1d(6,  2,3,4,5,6,7);
        println("y = %", y);
        _copy(6, x, 1, *y, 1);
        println("y = %", y);
        assert(y == x);
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        x := array1d(6,  1,2,3,4,5,6);
        y := array1d(6,  2,3,4,5,6,7);
        println("y = %", y);
        _copy(3, x, 2, *y, 1);
        println("y = %", y);
        assert(y == array1d(6,  1,3,5,5,6,7));
    }
}