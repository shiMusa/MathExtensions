
// ? LAPACK *dotc
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga0f02c96fa8498b4aa6b144deee725c0d.html#ga0f02c96fa8498b4aa6b144deee725c0d

cdotu :: _dotu;
zdotu :: _dotu;
_dotu :: inline ($$n: int, cx: $X/ArrayType1D, $$incx: int, cy: $Y/ArrayType1D, $$incy: int) -> UpCast(X.data_type, Y.data_type).T {
    return _dotu_T(n,cx,incx,cy,incy,UpCast(X.data_type, Y.data_type).T);
}
_dotu_T :: ($$n: int, cx: $X/ArrayType1D, $$incx: int, cy: $Y/ArrayType1D, $$incy: int, $T: Type) -> T {
    res := zero(T);
    if n < 0 then return res;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            res += ncast(T,get(cx,i)) * ncast(T,get(cy,i));
        }
    } else {
        ix := 0;
        iy := 0;
        if incx < 0 then ix = (-n+1)*incx;
        if incy < 0 then iy = (-n+1)*incy;
        for i: 0..n-1 {
            res += ncast(T,get(cx,ix)) * ncast(T,get(cy,iy));
            ix += incx;
            iy += incy;
        }
    }
    return res;
}



// ? http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga32d6fccb43cb13feafc23825f2685ba0.html#ga32d6fccb43cb13feafc23825f2685ba0
dsdot :: inline ($$n: int, cx: $X/ArrayType1D, $$incx: int, cy: $Y/ArrayType1D, $$incy: int) -> float64 {
    #assert(X.data_type == float32 && Y.data_type == float32);
    return _dotu_T(n,cx,incx,cy,incy,float64);
}

sdsdot :: inline ($$n: int, cx: $X/ArrayType1D, $$incx: int, cy: $Y/ArrayType1D, $$incy: int) -> float32 {
    #assert(X.data_type == float32 && Y.data_type == float32);
    return cast(float32)_dotu_T(n,cx,incx,cy,incy,float64);
}


#scope_module
test_lapack_dotu :: () {
    println_push("_dotu", color = .FG_GREEN);

    {
        println_push("test 1", color = .FG_YELLOW);
        x := int.[1,2,3,4,5,6];
        y := int.[2,3,4,5,6,7];
        println("y = %", y);
        res := _dotu(6, x, 1, y, 1);
        println("res = %", res);
        assert(res == 2+6+12+20+30+42);
    }

    {
        println_push("test 2", color = .FG_YELLOW);
        x := int.[1,2,3,4,5,6];
        y := int.[2,3,4,5,6,7];
        println("y = %", y);
        res := _dotu(3, x, 2, y, 1);
        println("res = %", res);
        assert(res == 2+9+20);
    }

    {
        println_push("test 3", color = .FG_YELLOW);
        x := array(4, R, I, 2*R, 2*I);
        y := array(4, 2*R, 2*R, 2*I, 2*I);
        res := _dotu(4, x, 1, y, 1);
        println("res = %", res);
        assert(res == 2*R+2*I+4*I-4*R);
    }
}