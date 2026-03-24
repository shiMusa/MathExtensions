// ? LAPACK *lacrt()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clacrt/zlacrt
// *LACRT performs a linear transformation of a pair of complex vectors.

clacrt :: #bake_arguments _lacrt(MT=C32);
zlacrt :: #bake_arguments _lacrt(MT=C64);

#scope_file
_lacrt :: (
    $$n: int,
    cx: *$CX/ArrayType1D,
    $$incx: int,
    cy: *$CY/ArrayType1D,
    $$incy: int,
    c: $MT,
    s: MT,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           CX.data_type == MT &&
           CY.data_type == MT;
} {
    ctemp : MT = ---;
    i, ix, iy : int = ---;

    if n <= 0 then return;
    if incx == 1 && incy == 1 {
        for i:1..n {
            ctemp = c*fget(cx.*, i) + s*fget(cy.*, i);
            fset(cy, i, c*fget(cy.*, i) - s*fget(cx.*, i));
            fset(cx, i, ctemp);
        }
        return;
    }

    ix = 1;
    iy = 1;
    if incx < 0 {
        ix = (-n + 1)*incx + 1;
    }
    if incy < 0 {
        iy = (-n + 1)*incy + 1;
    }

    for i:1..n {
        ctemp = c*fget(cx.*, ix) + s*fget(cy.*, iy);
        fset(cy, iy, c*fget(cy.*, iy) - s*fget(cx.*, ix));
        fset(cx, ix, ctemp);
        ix += incx;
        iy += incy;
    }
}
