

// ? LAPACK *scal
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_gab2d569bbfe6356599c69c057ecc2b3f7.html#gab2d569bbfe6356599c69c057ecc2b3f7

sscal  :: _scal;
dscal  :: _scal;
cscal  :: _scal;
csscal :: _scal;
zscal  :: _scal;
zdscal :: _scal;
_scal :: ($$n: int, ca: $U, cx: *$X/ArrayType1D, $$incx: int) #modify {
    if !is_intrinsic_number(U) {U = null; return;}
} {
    if n <= 0 then return;
    if incx == 1 {
        for i: 0..n-1 {
            set(cx,i, get(cx,i) * ca);
        }
    } else {
        ix := 0;
        if incx < 0 then ix = (-n+1)*incx;
        for i: 0..n-1 {
            set(cx,ix, get(cx,ix) * ca);
            ix += incx;
        }
    }
}

#scope_module
test_lapack_scal :: () {
    println_push("_scal", color = .FG_GREEN);

    {
        println("test 1", color = .FG_YELLOW);
        y := array1d(6,  1,2,3,4,5,6);
        _scal(6, 2, *y, 1);
        assert(y == array1d(6,  2,4,6,8,10,12));
    }

    {
        println("test 2", color = .FG_YELLOW);
        y := array1d(6,  1,2,3,4,5,6);
        _scal(3, -2, *y, 2);
        assert(y == array1d(6,  -2,2,-6,4,-10,6));
    }
}