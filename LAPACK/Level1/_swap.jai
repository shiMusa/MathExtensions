

// ? LAPACK *swap
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_gaa50b533bbd2eceac1f59dbc780d7a182.html#gaa50b533bbd2eceac1f59dbc780d7a182

sswap :: _swap;
dswap :: _swap;
cswap :: _swap;
zswap :: _swap;
_swap :: ($$n: int, cx: *$X/ArrayType1D, $$incx: int, cy: *$Y/ArrayType1D, $$incy: int) {
    if n <= 0 then return;
    if incx == 1 && incy == 1 {
        for i: 1..n {
            tmp := fget(cx.*, i);
            fset(cx, i, fget(cy.*, i));
            fset(cy, i, tmp);
        }
    } else {
        ix := 1;
        iy := 1;
        if incx < 0 then ix = (-n+1)*incx + 1;
        if incy < 0 then iy = (-n+1)*incy + 1;
        for i: 1..n {
            tmp := fget(cx.*, ix);
            fset(cx, ix, fget(cy.*, iy));
            fset(cy, iy, tmp);
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
        x := array1d(6,  1,2,3,4,5,6);
        y := array1d(6,  2,3,4,5,6,7);
        println("y = %", y);
        _swap(6, *x, 1, *y, 1);
        println("y = %", y);
        assert(y == array1d(6,  1,2,3,4,5,6));
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        x := array1d(6,  1,2,3,4,5,6);
        y := array1d(6,  2,3,4,5,6,7);
        println("y = %", y);
        _swap(3, *x, 2, *y, 1);
        println("y = %", y);
        assert(y == array1d(6,  1,3,5,5,6,7));
    }
}
