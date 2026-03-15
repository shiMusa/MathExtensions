
// ? LAPACK *srot
// ? http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_gaf660f9cda67820f314f36c9668272987.html#gaf660f9cda67820f314f36c9668272987

srot :: _rot;
drot :: _rot;
csrot :: _rot;
zdrot :: _rot;
_rot :: ($$n: int, cx: *$X/ArrayType1D, $$incx: int, cy: *$Y/ArrayType1D, $$incy: int, c:$U, s: U) #modify {
    return is_float_number(U);
} {
    if n <= 0 then return;
    if incx == 1 && incy == 1 {
        for i: 1..n {
            tmp :=    c * fget(cx.*, i) + s * fget(cy.*, i);
            fset(cy, i, c * fget(cy.*, i) - s * fget(cx.*, i));
            fset(cx, i, tmp);
        }
    } else {
        ix := 1;
        iy := 1;
        if incx < 0 then ix = (-n+1)*incx + 1;
        if incy < 0 then iy = (-n+1)*incy + 1;
        for i: 1..n {
            tmp :=     c * fget(cx.*, ix) + s * fget(cy.*, iy);
            fset(cy, iy, c * fget(cy.*, iy) - s * fget(cx.*, ix));
            fset(cx, ix, tmp);
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
