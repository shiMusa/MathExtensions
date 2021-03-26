


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

    print("Difficult problems\n");
    // print("%\n", 2.0e-64/1.0e-64);
    a : float64 = 2.0e-45;
    b : float64 = 1.0e-45;
    c : float64 = a/b;
    print("%\n", c);
    r1, r2 = solve_quadratic_real(1.0e160, -2.0e160, 1.0e160); print("%, %\n", r1, r2);
}