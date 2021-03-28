


I :: Complex64.{0.0, 1.0};
R :: Complex64.{1.0, 0.0};
CZERO :: Complex64.{0.0, 0.0};

C64_NAN :: Complex64.{Math.FLOAT64_NAN, Math.FLOAT64_NAN};


Complex64 :: struct {
    real: float64;
    imag: float64;
}

c64 :: (re: float64 = 0.0, im: float64 = 0.0) -> Complex64 {
    r: Complex64 = ---;
    r.real = re;
    r.imag = im;
    return r;
}

operator + :: inline (a: Complex64, b: Complex64) -> Complex64 {
    return c64(a.real + b.real, a.imag + b.imag);
}

operator + :: inline (a: Complex64, b: float64) -> Complex64 #symmetric {
    return c64(a.real + b, a.imag);
}

operator - :: inline (a: Complex64, b: Complex64) -> Complex64 {
    return c64(a.real - b.real, a.imag - b.imag);
}

operator - :: inline (a: Complex64, b: float64) -> Complex64 {
    return c64(a.real - b, a.imag);
}

operator - :: inline (a: float64, b: Complex64) -> Complex64 {
    return c64(a - b.real, b.imag);
}

operator - :: inline (a: Complex64) -> Complex64 {
    return c64(-a.real, -a.imag);
}

operator * :: inline (a: Complex64, b: Complex64) -> Complex64 {
    return c64(
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real
    );
}

operator / :: (a: Complex64, b: Complex64) -> Complex64 {
    abs := b.real * b.real + b.imag * b.imag;
    return c64(
        ( a.real * b.real + a.imag * b.imag )/abs,
        ( b.real * a.imag - a.real * b.imag )/abs
    );
}

operator * :: inline (a: Complex64, v : float64) -> Complex64 #symmetric {
    return c64(a.real * v, a.imag * v);
}

operator / :: inline (a: Complex64, v : float64) -> Complex64 {
    return c64(a.real / v, a.imag / v);
}

operator / :: inline (v : float64, a: Complex64) -> Complex64 {
    return c64(v, 0.0)/a;
}

operator == :: (a: Complex64, b: Complex64) -> bool {
    return a.real == b.real && a.imag == b.imag;
}

operator != :: (a: Complex64, b: Complex64) -> bool {
    return a.real != b.real || a.imag != b.imag;
}



conjugate :: inline (a: Complex64) -> Complex64 {
    return c64(a.real, -a.imag);
}

abs_sq :: inline (a: Complex64) -> float64 {
    return a.real * a.real + a.imag * a.imag;
}

abs :: inline (a: Complex64) -> float64 {
    return Math.sqrt(a.real * a.real + a.imag * a.imag);
}

arg :: (a: Complex64) -> float64 {

    // prevent division by zero
    if a.real == 0 {
        return ifx a.imag >= 0 then 0.5 * Math.PI else 1.5 * Math.PI;
    }

    if a.real > 0 && a.imag >= 0 {
        return Math.atan(a.imag / a.real);
    } else if a.real < 0 {
        return Math.atan(a.imag / a.real) + Math.PI;
    } else {
        return Math.atan(a.imag / a.real) + 2.0 * Math.PI;
    }
}

phase :: inline (theta: float64) -> Complex64 {
    return c64(Math.cos(theta), Math.sin(theta));
}

polar :: inline (mag: float64, theta: float64) -> Complex64 {
    return c64(mag * Math.cos(theta), mag * Math.sin(theta));
}






test_complex :: () {
    print("\n\nComplex -------------------------------\n\n");

    a := c64(1.0, 0.0);
    b := c64(0.0, 1.0);


    assert(b * b == -1.0 * R);
    assert(a * a == 1.0 * R);
    assert(a/b == -1.0 * I);
    assert(2.0*a+b == Complex64.{2.0, 1.0});
    assert(2.0*a-3*b == Complex64.{2.0, -3.0});
    assert(a != b);
    assert(abs_sq(a) == 1.0);
    assert(abs_sq(b) == 1.0);

    print("%\n", abs(Complex64.{1.0, 1.0}));
    print("%\n", Math.sqrt(2.0));
    // assert(abs(Complex64.{1.0, 1.0}) - Math.sqrt(2.0));

    print(
        "%, %, %, %, %, %, %, %\n",
        arg(Complex64.{ 1.0,  0.0}) / Math.PI * 180.0,
        arg(Complex64.{ 1.0,  1.0}) / Math.PI * 180.0,
        arg(Complex64.{ 0.0,  1.0}) / Math.PI * 180.0,
        arg(Complex64.{-1.0,  1.0}) / Math.PI * 180.0,
        arg(Complex64.{-1.0,  0.0}) / Math.PI * 180.0,
        arg(Complex64.{-1.0, -1.0}) / Math.PI * 180.0,
        arg(Complex64.{ 0.0, -1.0}) / Math.PI * 180.0,
        arg(Complex64.{ 1.0, -1.0}) / Math.PI * 180.0,
    );


    print("%\n", 2.0*R);

}