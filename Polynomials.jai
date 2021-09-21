


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.1, page 49
solve_quadratic :: (a: $T, b: T, c: T) -> ($U, U) #modify {
    if !is_number(T) {U = null; return;}
    if !is_complex_number(T) {
        U = get_complex_type_info_from_intrinsic_type_info(T);
    } else {
        U = T;
    }
} {
    X :: #run get_complex_inner_type(T);

    sigma := max(abs(a), abs(b), abs(c));
    if sigma == 0 {
        return (nan(U), nan(U));
    }

    sigma2 := next_power_of_two(cast(int)Math.ceil(sigma));

    a /= sigma;
    b /= sigma;
    c /= sigma;

    if abs(a) == 0 {
        if abs(b) == 0 then return (nan(U), nan(U));
        return (-c/b * R, nan(U));
    }

    b /= -2.0*a;
    c /= a;

    d := b*b - c;

    if d >= 0 {
        r1 := b + sign(b) * sqrt(d);

        if abs(r1) > 0 then return (r1 * R, c/r1 * R);
        return (r1 * R, 0.0 * R);
    }

    d = sqrt(-d);

    return (cpx(b, d), cpx(b, -d));
}

// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.1, page 49
solve_quadratic_real :: (a: $T, b: T, c: T) -> (T, T) #modify {
    if !is_number(T) {T = null; return;}
    if is_complex_number(T) {T = null; return;}   
} {
    sigma := max(Math.abs(a), Math.abs(b), Math.abs(c));
    if sigma == 0 {
        return (nan(T), nan(T));
    }

    sigma2 := next_power_of_two(cast(int)Math.ceil(sigma));

    a /= sigma;
    b /= sigma;
    c /= sigma;

    if Math.abs(a) == 0 {
        if Math.abs(b) == 0 then return (nan(T), nan(T));
        return (-c/b, nan(T));
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
    return (nan(T), nan(T));
}



// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.2 Horner's Rule, page 52
// a: list of prefactors from highest to lowest exponent, i.e.
// a[0] x^2 + a[1] x + a[2] x^0
polynom :: (x: $T, a: ..T) -> T #modify {
    if !is_number(T) {T = null; return;}    
} {
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
synthetic_division :: (r: $T, a: ..T) -> (T, []T) #modify {
    if !is_number(T) {T = null; return;}    
} {
    b : [..]T;
    d : T;
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
repeated_synthetic_division :: (r: $T, a: ..T) -> ([]T, []T) #modify {
    if !is_number(T) {T = null; return;}    
} {
    n := a.count-1;

    b : [..]T;
    array_resize(*b, n, false);
    c : [..]T;
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
    println_push("Polynomials", color = .FG_WHITE);

    r1 : float64;
    r2 : float64;
    r1, r2 = solve_quadratic_real(0.0, 0.0, 0.0); println("%, %", r1, r2); assert(Math.isnan(r1)); assert(Math.isnan(r2)); 
    r1, r2 = solve_quadratic_real(1.0, 0.0,-1.0); println("%, %", r1, r2); assert(r1 == 1.0);      assert(r2 == -1.0); 
    r1, r2 = solve_quadratic_real(1.0, 0.0, 1.0); println("%, %", r1, r2); assert(Math.isnan(r1)); assert(Math.isnan(r2)); 
    r1, r2 = solve_quadratic_real(1.0, 0.0,-9.0); println("%, %", r1, r2); assert(r1 == 3);        assert(r2 == -3); 
    r1, r2 = solve_quadratic_real(1.0, 6.0, 9.0); println("%, %", r1, r2); assert(r1 == -3);       assert(r2 == -3); 
    r1, r2 = solve_quadratic_real(1.0, 2.0,-3.0); println("%, %", r1, r2); assert(r1 == -3);       assert(r2 == 1); 

    {
        println_push("Difficult quadratic problems (check manually)", color = .FG_GREEN);
        r1, r2 = solve_quadratic_real(1.0e160, -2.0e160, 1.0e160);  println("%, %", r1, r2); // result: 1, 1
        r1, r2 = solve_quadratic_real(1.0e-20, -1.0, 1.0);          println("%, %", r1, r2); // result: 10^20, 1
        r1, r2 = solve_quadratic_real(1.0, -1.0, 1.0e-20);          println("%, %", r1, r2); // result: 1, 10^-20
        r1, r2 = solve_quadratic_real(1.0e-10, -1.0, 1.0e-10);      println("%, %", r1, r2); // result: 10^10, 1
        r1, r2 = solve_quadratic_real(1.0e-5, -1.0, 1.0e-5);        println("%, %", r1, r2); // result: 2e5 - 0.5e-5, 0.5e-5 + 0.125e-15
        r1, r2 = solve_quadratic_real(1.0e-220, -1.0e100, 1.0e100); println("%, %", r1, r2); // result: 10e320, 1
    }

    {
        println_push("Polynom evaluation", color = .FG_GREEN);

        // print("%\n", polynom(2.0, 3.0, -1.0, 2.0));
        assert(polynom(2.0, 2.0, -1.0, 3.0) == 9.0);

        d, b := synthetic_division(3.0, 2.0, -1.0, 3.0);
        println("%, %", d, b);
        assert(d == 18 && b[0] == 2 && b[1] == 5);
    }
}