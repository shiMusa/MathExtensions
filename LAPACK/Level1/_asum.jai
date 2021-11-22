
// ? LAPACK *asum
// ? http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga7372361a44de0649813750b3280c58cc.html#ga7372361a44de0649813750b3280c58cc

scasum :: sasum;
sasum :: inline ($$n: int, x: $A/ArrayType1D, $$incx: int) -> float32 {
    return _asum(float32, n, x, incx);
}
dzasum :: dasum;
dasum :: inline ($$n: int, x: $A/ArrayType1D, $$incx: int) -> float64 {
    return _asum(float64, n, x, incx);
}
_asum :: ($U: Type, $$n: int, x: $A/ArrayType1D, $$incx: int) -> U {
    #assert(#run is_float_number(U) && is_intrinsic_number(U));
    sum := zero(U);
    // * tmp not necessary

    if n <= 0 then return sum; // * allow for negative stride as well - F. Fehse, 2021-11-16
    if incx == 1 {
        m := n % 6; // ! why blocks of 6 values?
        if m != 0 {
            for i: 0..m-1 {
                sum += abs(U, get(x,i));
            }
            if n < 6 then return sum;
        }
        ix := m+1;
        while ix < n {
            sum += abs(U, get(x,ix  )) + abs(U, get(x,ix+1)) + abs(U, get(x,ix+2)) 
                 + abs(U, get(x,ix+3)) + abs(U, get(x,ix+4)) + abs(U, get(x,ix+5));
            ix += 6;
        }
        return sum;
    } else {
        ix := 0;
        for i: 0..n-1 {
            sum += abs(U, get(x,ix));
            ix += incx;
        }
        return sum;
    }
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