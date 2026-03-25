// ? LAPACK *lartv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slartv/dlartv/clartv/zlartv
// *LARTV applies a vector of plane rotations.

slartv :: #bake_arguments _lartv_real(RT=float32);
dlartv :: #bake_arguments _lartv_real(RT=float64);
clartv :: #bake_arguments _lartv_complex(RT=float32, MT=C32);
zlartv :: #bake_arguments _lartv_complex(RT=float64, MT=C64);

_lartv_real :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int,
    c: $C/ArrayType1D, s: $S/ArrayType1D, $$incc: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT && Y.data_type == RT && C.data_type == RT && S.data_type == RT;
} {
    ix, iy, ic : int = ---;
    xi, yi : RT = ---;
    ix = 1; iy = 1; ic = 1;
    for _:1..n {
        xi = fget(x.*, ix);
        yi = fget(y.*, iy);
        fset(x, ix, fget(c, ic)*xi + fget(s, ic)*yi);
        fset(y, iy, fget(c, ic)*yi - fget(s, ic)*xi);
        ix = ix + incx;
        iy = iy + incy;
        ic = ic + incc;
    }
}

_lartv_complex :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int,
    c: $C/ArrayType1D, s: $S/ArrayType1D, $$incc: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT && Y.data_type == MT && C.data_type == RT && S.data_type == MT;
} {
    ix, iy, ic : int = ---;
    xi, yi : MT = ---;
    ix = 1; iy = 1; ic = 1;
    for _:1..n {
        xi = fget(x.*, ix);
        yi = fget(y.*, iy);
        fset(x, ix, fget(c, ic)*xi + fget(s, ic)*yi);
        fset(y, iy, fget(c, ic)*yi - conjugate(fget(s, ic))*xi);
        ix = ix + incx;
        iy = iy + incy;
        ic = ic + incc;
    }
}
