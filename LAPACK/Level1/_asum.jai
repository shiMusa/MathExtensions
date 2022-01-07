
// ? LAPACK *asum
// ? http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga7372361a44de0649813750b3280c58cc.html#ga7372361a44de0649813750b3280c58cc

scasum :: sasum;
sasum :: #bake_arguments _asum(U=float32);
dzasum :: dasum;
dasum :: #bake_arguments _asum(U=float64);

_asum :: ($U: Type, $$n: int, cx: $A/ArrayType1D, $$incx: int) -> U {
    #run assert(U == get_number_inner_type(A.data_type));
    res := zero(U);
    // * tmp not necessary

    if n <= 0 || incx < 0 then return res;

    if incx == 1 {
        for i:0..n-1 {
            res += abs1(get(cx,i));
        }
    } else {
        ix := 0;
        for i:0..n-1 {
            res += abs1(get(cx,ix));
            ix += incx;
        }
    }

    return res;
}

#scope_module
test_lapack_asum :: () {
    println_push("_asum", color = .FG_GREEN);
    {
        println("test 1", color = .FG_YELLOW);
        x := array1d(4,  1,2,3,4);
        res := _asum(float64, 4, x, 1);
        assert(res == 10.0);
    }
    {
        println_push("test 2", color = .FG_YELLOW);
        x := array1d(3,  R,I,R+I);
        res := _asum(float64, 3, x, 1);
        println("res = %", res);
        assert( approx( res, 2.0+sqrt(2.0) ) );
    }
}