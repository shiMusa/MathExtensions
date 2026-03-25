// ? LAPACK *lar2v()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slar2v/dlar2v/clar2v/zlar2v
// *LAR2V applies a vector of plane rotations from both sides to 2x2 blocks.

slar2v :: #bake_arguments _lar2v_real(RT=float32);
dlar2v :: #bake_arguments _lar2v_real(RT=float64);
clar2v :: #bake_arguments _lar2v_complex(RT=float32, MT=C32);
zlar2v :: #bake_arguments _lar2v_complex(RT=float64, MT=C64);

_lar2v_real :: (
    $$n: int, x: *$X/ArrayType1D, y: *$Y/ArrayType1D, z: *$Z/ArrayType1D,
    $$incx: int, c: $C/ArrayType1D, s: $S/ArrayType1D, $$incc: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT && Y.data_type == RT && Z.data_type == RT &&
           C.data_type == RT && S.data_type == RT;
} {
    ix, ic : int = ---;
    xi, yi, zi, ci, si, t1, t2, t3, t4, t5, t6 : RT = ---;

    ix = 1;
    ic = 1;
    for _:1..n {
        xi = fget(x.*, ix);
        yi = fget(y.*, ix);
        zi = fget(z.*, ix);
        ci = fget(c, ic);
        si = fget(s, ic);
        t1 = si*zi;
        t2 = ci*zi;
        t3 = t2 - si*xi;
        t4 = t2 + si*yi;
        t5 = ci*xi + t1;
        t6 = ci*yi - t1;
        fset(x, ix, ci*t5 + si*t4);
        fset(y, ix, ci*t6 - si*t3);
        fset(z, ix, ci*t4 - si*t5);
        ix = ix + incx;
        ic = ic + incc;
    }
}

_lar2v_complex :: (
    $$n: int, x: *$X/ArrayType1D, y: *$Y/ArrayType1D, z: *$Z/ArrayType1D,
    $$incx: int, c: $C/ArrayType1D, s: $S/ArrayType1D, $$incc: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT && Y.data_type == MT && Z.data_type == MT &&
           C.data_type == RT && S.data_type == MT;
} {
    ix, ic : int = ---;
    ci, sii, sir, t1i, t1r, t5, t6, xi, yi, zii, zir : RT = ---;
    zi, si, t2, t3, t4 : MT = ---;

    ix = 1;
    ic = 1;
    for _:1..n {
        xi = real(fget(x.*, ix));
        yi = real(fget(y.*, ix));
        zi = fget(z.*, ix);
        zir = real(zi);
        zii = imag(zi);
        ci = fget(c, ic);
        si = fget(s, ic);
        sir = real(si);
        sii = imag(si);
        t1r = sir*zir - sii*zii;
        t1i = sir*zii + sii*zir;
        t2 = ci*zi;
        t3 = t2 - conjugate(si)*xi;
        t4 = conjugate(t2) + si*yi;
        t5 = ci*xi + t1r;
        t6 = ci*yi - t1r;
        fset(x, ix, cpx(ci*t5 + (sir*real(t4) + sii*imag(t4)), zero(RT)));
        fset(y, ix, cpx(ci*t6 - (sir*real(t3) - sii*imag(t3)), zero(RT)));
        fset(z, ix, ci*t3 + conjugate(si)*cpx(t6, t1i));
        ix = ix + incx;
        ic = ic + incc;
    }
}
