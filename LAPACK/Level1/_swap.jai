

// ? LAPACK *swap
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_gaa50b533bbd2eceac1f59dbc780d7a182.html#gaa50b533bbd2eceac1f59dbc780d7a182

sswap :: _swap;
dswap :: _swap;
cswap :: _swap;
zswap :: _swap;
_swap :: ($$n: int, cx: $X/Array1D, $$incx: int, cy: $Y/Array1D, $$incy: int) {
    if n <= 0 then return;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            tmp := get(cx,i);
            set(cx,i, get(cy,i));
            set(cy,i, tmp);
        }
    } else {
        ix := 0;
        iy := 0;
        if incx < 0 then ix = (-n+1)*incx;
        if incy < 0 then iy = (-n+1)*incy;
        for i: 0..n-1 {
            tmp := get(cx,ix);
            set(cx,ix, get(cy,iy));
            set(cy,iy, tmp);
            ix += incx;
            iy += incy;
        }
    }
}

#scope_module
test_lapack_swap :: () {
    println_push("_swap", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        x := int.[1,2,3,4,5,6];
        y := int.[2,3,4,5,6,7];
        println("y = %", y);
        _swap(6, x, 1, y, 1);
        println("y = %", y);
        assert(array_equals(y, int.[1,2,3,4,5,6]));
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        x := int.[1,2,3,4,5,6];
        y := int.[2,3,4,5,6,7];
        println("y = %", y);
        _swap(3, x, 2, y, 1);
        println("y = %", y);
        assert(array_equals(y, int.[1,3,5,5,6,7]));
    }
}