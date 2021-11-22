
// ? LAPACK *srot
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_gaf660f9cda67820f314f36c9668272987.html#gaf660f9cda67820f314f36c9668272987

srot :: _rot;
drot :: _rot;
csrot :: _rot;
zdrot :: _rot;
_rot :: ($$n: int, cx: *$X/ArrayType1D, $$incx: int, cy: *$Y/ArrayType1D, $$incy: int, c:$U, s: U) #modify {
    if !is_intrinsic_float_number(U) {U = null; return;}
} {
    if n <= 0 then return;
    if incx == 1 && incy == 1 {
        for i: 0..n-1 {
            tmp :=    c * get(cx,i) + s * get(cy,i);
            set(cy,i, c * get(cy,i) - s * get(cx,i));
            set(cx,i, tmp);
        }
    } else {
        ix := 0;
        iy := 0;
        if incx < 0 then ix = (-n+1)*incx;
        if incy < 0 then iy = (-n+1)*incy;
        for i: 0..n-1 {
            tmp :=     c * get(cx,ix) + s * get(cy,iy);
            set(cy,iy, c * get(cy,iy) - s * get(cx,ix));
            set(cx,ix, tmp);
            ix += incx;
            iy += incy;
        }
    }
}

#scope_module
test_lapack_rot :: () {
    println_push("_rot", color = .FG_GREEN);

    println("TEST NOT YET IMPLEMENTED!", color = .FG_RED);
}