


I      :: C64.{0.0, 1.0};
R      :: C64.{1.0, 0.0};
I64 :: I;
R64 :: R;
I32 :: C32.{0.0,1.0};
R32 :: C32.{1.0,0.0};
C_ZERO :: C64.{0.0, 0.0};
C_NAN  :: C64.{Math.FLOAT64_NAN, Math.FLOAT64_NAN};



C64 :: Complex(float64);
C32 :: Complex(float32);




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


c64 :: inline ($$re: $T = 0.0, $$im: $U = 0.0) -> C64 #modify {
    if !is_intrinsic_number(T) || !is_intrinsic_number(U) {T = null; U = null; return;}    
} {
    return cpx(ncast(float64, re), ncast(float64, im));
}

c32 :: inline ($$re: $T = 0.0, $$im: $U = 0.0) -> C32 #modify {
    if !is_intrinsic_number(T) || !is_intrinsic_number(U) {T = null; U = null; return;}    
} {
    return cpx(ncast(float32, re), ncast(float32, im));
}

cpx_r :: inline ($T: Type) -> Complex(T) {
    #run assert(is_intrinsic_number(T), "can only use intrinsic numbers to generate complex: Complex(T), T = % does not work.", T);
    return Complex(T).{#run one(T), #run zero(T)};
}

cpx_i :: inline ($T: Type) -> Complex(T) {
    #run assert(is_intrinsic_number(T), "can only use intrinsic numbers to generate complex: Complex(T), T = % does not work.", T);
    return Complex(T).{#run zero(T), #run one(T)};
}

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


cpx_inv :: inline (z: Complex($T)) -> Complex(T) {
    z_abs_sq := abs_sq(z);
    return cpx(z.real / z_abs_sq, -z.imag / z_abs_sq);
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

operator == :: inline (a: Complex($T), b: Complex($U)) -> bool {
    return a.real == b.real && a.imag == b.imag;
}

operator == :: inline (a: Complex($T), b: $U) -> bool #modify {
    if !is_intrinsic_number(U) { U = null; return; }
} #symmetric {
    return a.real == b && a.imag == 0.0;
}





#scope_module
test_complex :: () {
    println_push("Complex", color = .FG_WHITE);

    a := c64(1.0, 0.0);
    b := c64(0.0, 1.0);

    assert( c64(0.0, 0.0) == c64(-0.0, -0.0) );


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
        println_push("arg", color = .FG_GREEN);
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