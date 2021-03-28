


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.1, page 49
solve_quadratic :: (a: float64, b: float64, c: float64) -> (Complex64, Complex64) {
    sigma := max(Math.abs(a), Math.abs(b), Math.abs(c));
    if sigma == 0 {
        return (C64_NAN, C64_NAN);
    }

    sigma2 := cast(float64)next_power_of_two(cast(int)Math.ceil(sigma));

    a /= sigma;
    b /= sigma;
    c /= sigma;

    if Math.abs(a) == 0 {
        if Math.abs(b) == 0 then return (C64_NAN, C64_NAN);
        return (-c/b * R, C64_NAN);
    }

    b /= -2.0*a;
    c /= a;

    d := b*b - c;

    if d >= 0 {
        r1 := b + sign(b) * Math.sqrt(d);

        if Math.abs(r1) > 0 then return (r1 * R, c/r1 * R);
        return (r1 * R, 0.0 * R);
    }

    d = Math.sqrt(-d);

    return (c64(b, d), c64(b, -d));
}

// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.1, page 49
solve_quadratic_real :: (a: float64, b: float64, c: float64) -> (float64, float64) {
    sigma := max(Math.abs(a), Math.abs(b), Math.abs(c));
    if sigma == 0 {
        return (Math.FLOAT64_NAN, Math.FLOAT64_NAN);
    }

    sigma2 := cast(float64)next_power_of_two(cast(int)Math.ceil(sigma));

    a /= sigma;
    b /= sigma;
    c /= sigma;

    if Math.abs(a) == 0 {
        if Math.abs(b) == 0 then return (Math.FLOAT64_NAN, Math.FLOAT64_NAN);
        return (-c/b, Math.FLOAT64_NAN);
    }

    b /= -2.0*a;
    c /= a;

    d := b*b - c;

    if d >= 0 {
        r1 := b + sign(b) * Math.sqrt(d);

        if Math.abs(r1) > 0 then return (r1, c/r1);
        return (r1, 0.0);
    }

    // d = Math.sqrt(-d);
    
    // no real solutions
    return (Math.FLOAT64_NAN, Math.FLOAT64_NAN);
}



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.2 Horner's Rule, page 52
// a: list of prefactors from highest to lowest exponent, i.e.
// a[0] x^2 + a[1] x + a[2] x^0
polynom :: (x: float64, a: ..float64) -> float64 {
    p := a[0];
    for 1..a.count-1 {
        p = a[it] + x * p;
    }
    return p;
}


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.3 Synthetic Division, page 53
// calculates for polynom p = a[0] x^n + a[1] x^n-1 + ... + a[n]
// the division, such that p = (x-r) * q + d,
// where q is a polynom of degree n-1.
// returns (d, factors of q)
synthetic_division :: (r: float64, a: ..float64) -> (float64, []float64) {
    b : [..]float64;
    d : float64;
    array_resize(*b, a.count - 1, false);

    b[0] = a[0];
    for 1..a.count-2 {
        b[it] = a[it] + r * b[it - 1];
    }
    d = a[a.count-1] + r * b[a.count-2];

    return (d, b);
}





// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.4 Repeated Synthetic Division, page 54
// calculates synthetic division as before, but returns prefactors 
// for q and ∂p/∂x
// it is p(r) = b[n], ∂p(r)/∂x = c[n-1]
repeated_synthetic_division :: (r: float64, a: ..float64) -> ([]float64, []float64) {
    n := a.count-1;

    b : [..]float64;
    array_resize(*b, n, false);
    c : [..]float64;
    array_resize(*c, n-1, false);

    b[0] = a[0];
    c[0] = a[0];
    for 1..n-1 {
        b[it] = a[it] + r * b[it - 1];
        c[it] = b[it] + r * c[it - 1];
    }
    b[n] = a[n] + r * b[n-1];

    return (b, c);
}











test_polynomials :: () {
    print("\n\nPolynomials -------------------------------\n\n");
    r1 : float64;
    r2 : float64;
    r1, r2 = solve_quadratic_real(0.0, 0.0, 0.0); print("%, %\n", r1, r2); assert(Math.isnan(r1)); assert(Math.isnan(r2)); 
    r1, r2 = solve_quadratic_real(1.0, 0.0,-1.0); print("%, %\n", r1, r2); assert(r1 == -1); assert(r2 == 1); 
    r1, r2 = solve_quadratic_real(1.0, 0.0, 1.0); print("%, %\n", r1, r2); assert(Math.isnan(r1)); assert(Math.isnan(r2)); 
    r1, r2 = solve_quadratic_real(1.0, 0.0,-9.0); print("%, %\n", r1, r2); assert(r1 == -3); assert(r2 == 3); 
    r1, r2 = solve_quadratic_real(1.0, 6.0, 9.0); print("%, %\n", r1, r2); assert(r1 == -3); assert(r2 == -3); 
    r1, r2 = solve_quadratic_real(1.0, 2.0,-3.0); print("%, %\n", r1, r2); assert(r1 == -3); assert(r2 == 1); 

    print("\nDifficult quadratic problems (check manually)\n");
    r1, r2 = solve_quadratic_real(1.0e160, -2.0e160, 1.0e160); print("%, %\n", r1, r2); // result: 1, 1
    r1, r2 = solve_quadratic_real(1.0e-20, -1.0, 1.0); print("%, %\n", r1, r2); // result: 10^20, 1
    r1, r2 = solve_quadratic_real(1.0, -1.0, 1.0e-20); print("%, %\n", r1, r2); // result: 1, 10^-20
    r1, r2 = solve_quadratic_real(1.0e-10, -1.0, 1.0e-10); print("%, %\n", r1, r2); // result: 10^10, 1
    r1, r2 = solve_quadratic_real(1.0e-5, -1.0, 1.0e-5); print("%, %\n", r1, r2); // result: 2e5 - 0.5e-5, 0.5e-5 + 0.125e-15
    r1, r2 = solve_quadratic_real(1.0e-220, -1.0e100, 1.0e100); print("%, %\n", r1, r2); // result: 10e320, 1

    print("\nPolynom evaluation\n");

    // print("%\n", polynom(2.0, 3.0, -1.0, 2.0));
    assert(polynom(2.0, 2.0, -1.0, 3.0) == 9.0);

    d, b := synthetic_division(3.0, 2.0, -1.0, 3.0);
    print("%, %\n", d, b);
    assert(d == 18 && b[0] == 2 && b[1] == 5);
}