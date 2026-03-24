// ? LAPACK *lacgv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clacgv/zlacgv
// *LACGV conjugates a complex vector.

clacgv :: #bake_arguments _lacgv(MT=C32);
zlacgv :: #bake_arguments _lacgv(MT=C64);

#scope_file
_lacgv :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    i, ioff : int = ---;

    if incx == 1 {
        for i:1..n {
            fset(x, i, conjugate(fget(x.*, i)));
        }
    } else {
        ioff = 1;
        if incx < 0 {
            ioff = 1 - (n - 1)*incx;
        }

        for i:1..n {
            fset(x, ioff, conjugate(fget(x.*, ioff)));
            ioff += incx;
        }
    }
}
