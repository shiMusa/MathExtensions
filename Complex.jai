


I :: Complex(float64).{0.0, 1.0};
R :: Complex(float64).{1.0, 0.0};
CZERO :: Complex(float64).{0.0, 0.0};

C64_NAN :: Complex(float64).{Math.FLOAT64_NAN, Math.FLOAT64_NAN};




str :: (z: Complex($T)) -> string {
    if z.imag >= 0 {
        return sprint("% + %i", z.real, z.imag);
    } else {
        return sprint("% - %i", z.real, -z.imag);
    }
}




Complex :: struct (T: Type) {
    real: T;
    imag: T;
}

// Complex64 :: Complex(float64);
// Complex32 :: Complex(float32);

cpx :: inline ($$re: $T = 0.0, $$im: T = 0.0) -> Complex(T) {
    #if is_constant(re) && is_constant(im) {
        return Complex(T).{re, im};
    } else {
        r: Complex(T) = ---;
        r.real = re;
        r.imag = im;
        return r;
    }
}

operator + :: inline (a: Complex($T), b: Complex($U)) -> Complex(UpCast(T,U).T) {
    return cpx(a.real + b.real, a.imag + b.imag);
}

operator + :: inline (a: Complex($T), b: $U) -> Complex(UpCast(T,U).T) #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} #symmetric {
    return cpx(a.real + b, a.imag);
}

operator - :: inline (a: Complex($T), b: Complex($U)) -> Complex(UpCast(T,U).T) {
    return cpx(a.real - b.real, a.imag - b.imag);
}

operator - :: inline (a: Complex($T), b: $U) -> Complex(UpCast(T,U).T) #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} {
    return cpx(a.real - b, a.imag);
}

operator - :: inline (a: $U, b: Complex($T)) -> Complex(UpCast(T,U).T) #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} {
    return cpx(a - b.real, -b.imag);
}

operator - :: inline (a: Complex($T)) -> Complex(T) {
    return cpx(-a.real, -a.imag);
}

operator * :: inline (a: Complex($T), b: Complex($U)) -> Complex(UpCast(T,U).T) {
    return cpx(
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real
    );
}

operator / :: (a: Complex($T), b: Complex($U)) -> Complex(UpCast(T,U).T) {
    abs := b.real * b.real + b.imag * b.imag;
    return cpx(
        ( a.real * b.real + a.imag * b.imag )/abs,
        ( b.real * a.imag - a.real * b.imag )/abs
    );
}

operator * :: inline (a: Complex($T), v : $U) -> Complex(UpCast(T,U).T) #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} #symmetric {
    return cpx(a.real * v, a.imag * v);
}
// operator * :: inline (v: float64, a: Complex($T)) -> Complex($T) {
//     return cpx(a.real * v, a.imag * v);
// }

operator / :: inline (a: Complex($T), v : $U) -> Complex(UpCast(T,U).T) #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} {
    return cpx(a.real / v, a.imag / v);
}

operator / :: inline (v : $U, a: Complex($T)) -> Complex(UpCast(T,U).T) #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} {
    return cpx(v, 0.0)/a;
}

operator == :: (a: Complex($T), b: Complex($U)) -> bool {
    // TODO UpCast here?
    return a.real == b.real && a.imag == b.imag;
}

operator == :: (a: Complex($T), b: $U) -> bool #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} #symmetric {
    return a.real == b && a.imag == 0.0;
}






test_complex_generic :: () {
    println_push("Complex(T)");

    a := cpx(1.0, 0.0);
    b := cpx(0.0, 1.0);

    assert( cpx(0.0, 0.0) == cpx(-0.0, -0.0) );


    assert(b * b == -1.0 * R);
    assert(a * a == 1.0 * R);
    assert(a/b == -1.0 * I);
    assert(2.0*a+b == Complex(float64).{2.0, 1.0});
    assert(2.0*a-3*b == Complex(float64).{2.0, -3.0});
    assert(a != b);
    assert(abs_sq(a) == 1.0);
    assert(abs_sq(b) == 1.0);

    // print("%\n", abs(Complex(float64).{1.0, 1.0}));
    // print("%\n", Math.sqrt(2.0));


    {
        println_push("arg");
        println(
            "%, %, %, %, %, %, %, %",
            arg(Complex(float64).{ 1.0,  0.0}) / Math.PI * 180.0,
            arg(Complex(float64).{ 1.0,  1.0}) / Math.PI * 180.0,
            arg(Complex(float64).{ 0.0,  1.0}) / Math.PI * 180.0,
            arg(Complex(float64).{-1.0,  1.0}) / Math.PI * 180.0,
            arg(Complex(float64).{-1.0,  0.0}) / Math.PI * 180.0,
            arg(Complex(float64).{-1.0, -1.0}) / Math.PI * 180.0,
            arg(Complex(float64).{ 0.0, -1.0}) / Math.PI * 180.0,
            arg(Complex(float64).{ 1.0, -1.0}) / Math.PI * 180.0,
        );
    }

    println("%", str(2.0*R));

    println("% == %", str(3.0-4.0*I), str(Complex(float64).{3.0, -4.0}));
    assert(3.0 - 4.0 * I == Complex(float64).{3.0, -4.0});
}