
// ? LAPACK *dotc
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga0f02c96fa8498b4aa6b144deee725c0d.html#ga0f02c96fa8498b4aa6b144deee725c0d

cdotc :: _dotc;
zdotc :: _dotc;
ddot :: _dotc;
sdot :: _dotc;
_dotc :: ($$n: int, cx: $X/ArrayType1D, $$incx: int, cy: $Y/ArrayType1D, $$incy: int) -> UpCast(X.data_type,Y.data_type).T {
    T :: UpCast(X.data_type, Y.data_type).T;
    res := zero(T);
    if n < 0 then return res;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            res += conjugate(get(cx,i))*get(cy,i);
        }
    } else {
        ix := 0;
        iy := 0;
        if incx < 0 then ix = (-n+1)*incx;
        if incy < 0 then iy = (-n+1)*incy;
        for i: 0..n-1 {
            res += conjugate(get(cx,ix)) * get(cy,iy);
            ix += incx;
            iy += incy;
        }
    }
    return res;
}


#scope_module
test_lapack_dotc :: () {
    println_push("_dotc", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        x := array1d(1,2,3,4,5,6);
        y := array1d(2,3,4,5,6,7);
        println("y = %", y);
        res := _dotc(6, x, 1, y, 1);
        println("res = %", res);
        assert(res == 2+6+12+20+30+42);
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        x := array1d(1,2,3,4,5,6);
        y := array1d(2,3,4,5,6,7);
        println("y = %", y);
        res := _dotc(3, x, 2, y, 1);
        println("res = %", res);
        assert(res == 2+9+20);
    }

    {
        println_push("test 3", color = .FG_YELLOW);
        x := array1d(R, I, 2*R, 2*I);
        y := array1d(2*R, 2*R, 2*I, 2*I);
        res := _dotc(4, x, 1, y, 1);
        println("res = %", res);
        assert(res == 2*R-2*I+4*I+4*R);
    }
}