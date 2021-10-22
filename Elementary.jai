






// ? a great source for quaternion functions is
// ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
// ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
// ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf

// ? There's also a Java implementation: https://ece.uwaterloo.ca/~dwharder/Java/








sign :: inline (x: $T) -> T #modify {
    if !is_number(T) {T = null; return;}    
} {
    #if #run is_quaternion_number(T) || is_complex_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        return x / abs(x);
    } else {
        return ifx x >= 0 then one(T) else -one(T);
    }
}


factorial :: (n: int, $T: Type) -> T {
    res := one(T);
    for 2..n res *= it;
    return res;
}



inv :: inline (n: $T) -> T #modify {
    if !is_number(T) {T = null; return;}
} {
    #if #run is_quaternion_number(T) {
        return qtn_inv(n);
    } else #if #run is_complex_number(T) {
        return #run cpx_inv(n);
    } else return 1/n;
}




// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 2.3.1.10 Binomial Coefficients, page 55
binomial :: inline (from: $T, choose: T) -> T {
    #if CHECKS {
        assert(choose <= from, "Out of range: from = %, choose = %", from, choose);
    }

    if choose == zero(T) || from == choose {
        return one(T);
    }

    return binomial(from - one(T), choose) + binomial(from - one(T), choose - one(T));
}



quat_r :: real;
real :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    U = get_number_inner_type(T);
} {
    #if #run is_complex_number(T) {
        return n.real;
    } else #if #run is_quaternion_number(T) {
        return n.r;    
    } else {
        return n;
    }
}

quat_i :: imag;
imag :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    U = get_number_inner_type(T);
} {
    #if #run is_complex_number(T) {
        return n.imag;
    } else #if #run is_quaternion_number(T) {
        return n.i;
    } else {
        return zero(T);
    }
}

quat_j :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    U = get_number_inner_type(T);
} {
    #if #run is_quaternion_number(T) {
        return n.j;
    } else {
        return zero(T);
    }
}

quat_k :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    U = get_number_inner_type(T);
} {
    #if #run is_quaternion_number(T) {
        return n.j;
    } else {
        return zero(T);
    }
}



inverse :: inline (n: $T) -> T #modify {
    if !is_number(T) {T = null; return;}
} {
    #if #run is_quaternion_number(n) {
        return qtn_inv(n);
    } else #if #run is_complex_number(n) {
        return #run cpx_inv(n);
    } else {
        return one(T)/n;
    }
}


conjugate :: inline (n: $T) -> T #modify {
    if !is_number(T) { T = null; return; }    
} {
    #if #run is_quaternion_number(T) {
        return qtn(n.r, -n.i, -n.j, -n.k);
    } else #if #run is_complex_number(T) {
        return cpx(n.real, -n.imag);
    } else {
        return n;
    }
}


// ? |z|²
abs_sq :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    U = get_number_inner_type(T);
} {
    #if #run is_quaternion_number(T) {
        return n.r*n.r + n.i*n.i + n.j*n.j + n.k*n.k;
    } else #if #run is_complex_number(T) {
        return n.real * n.real + n.imag * n.imag;
    } else {
        return n * n;
    }
}

abs :: inline (a: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    U = get_number_inner_type(T);
} {
    return Math.sqrt(abs_sq(a));
}


// ? given n = |n| exp(iθ), whis calculates θ.
// ? for quaternions: n = |n| exp(iθ eᵥ), where eᵥ = v/|v| with n = n.r + v
arg :: inline (a: $T) -> float64 #modify {
    if !is_number(T) {T = null; return;}    
} {
    #if #run is_quaternion_number(T) {
        // ? https://www.wikiwand.com/en/Quaternion
        return Math.acos(a.r/abs(a)); // ! branch cuts?
    } else #if #run is_complex_number(T) {
        // return carg(a);
        return Math.atan2(a.imag, a.real);
    } else 
        return 0.0;
}




// ? exp(z) where z ∈ ℂ,ℝ
exp :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://www.wikiwand.com/en/Quaternion
        v := qtn_vec_abs(n);
        return Math.exp(n.r) * (Math.cos(v) + qtn_vec(n)/v * Math.sin(v));
    } else #if #run is_complex_number(T) {
       return Math.exp(n.real) * phase(n.imag);
    } else {
        return Math.exp(n);
    }
}


// ? log(z) where z ∈ ℂ,ℝ
log :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://www.wikiwand.com/en/Quaternion
        // ?
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        if qtn_vec_is_zero(n) then return ncast(T, Math.log(abs(n.r)));
        
        q := abs(n);
        return ncast(T, Math.log(q) + sign(qtn_vec(n)) * Math.acos(n.r/q));

    } else #if #run is_complex_number(T) {
        return cpx(Math.log(abs(n)), Math.atan2(n.imag, n.real));
    } else {
        return cast(T)Math.log(cast(float64)n);
    }
}



// ? n^p, where n,p ∈ ℂ,ℝ
// ? n,p ∈ ℝ: build in pow(n,p)
// ? n ∈ ℂ, p ∈ ℝ: |n|^p exp(i θ p)
// ? n ∈ ℝ, p ∈ ℂ: pow(n,Re{p}) * exp(i ln(n) Im{p})
// ? n,p ∈ ℂ: |n|^p * exp(iθ p) referring back to above versions of exp and pow(ℝ,ℂ)
// ? 
// ? there might be better algorithms calculating this...
// ? but it's basically the same as in the Gnu Scientific Library (GSL)
// ? see complex/math.c:
// ?     gsl_complex gsl_complex_pow (gsl_complex a, gsl_complex b)
// ?
// ! the `type` thingy is kinda akward...
// TODO remove `type`
// ?
// ? Quaternion power: 
// ? "Taschenbuch der Mathematik", Bronstein et al, Verlag Harri Deutsch, 7. Auflage 2008
pow :: inline (n: $T, power: $U, $type: Type) -> type {
    #if #run is_quaternion_number(T) || is_quaternion_number(U) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        #if CHECKS { assert(n != zero(T), "cannot take the quaternion power of 0!"); }
        return exp(log(n)*power);
    } else #if #run is_complex_number(T) {
        #if #run is_complex_number(U) {
            return pow(abs(n), power, UpCast(T,U).T) * phase(arg(n) * power);
        } else 
            return Math.pow(abs(n), power) * phase(arg(n) * power);
    } else {
        #if #run is_complex_number(U) {
            return Math.pow(n, power.real) * phase(Math.log(n) * power.imag);
        } else
            return Math.pow(n, power);
    }
}


// ? in the hopes sqrt is faster than pow
// ? need to wait for compiler fix to get rid of `$type: Type`...
// TODO ^ I might be able to do that now!
sqrt :: inline(n: $T, $type: Type) -> type {
    #if #run is_quaternion_number(T) {
        return pow(n, 0.5, type);
    } else #if #run is_complex_number(T) {
        return Math.sqrt(abs(n)) * phase(arg(n) * 0.5);
    } else {
        #if #run is_complex_number(type) {
            return sqrt(cpx(n), type);
        } else
            return Math.sqrt(n);
    }
}




// ? exp(iθ) where θ ∈ ℂ,ℝ
// ? is also implemented for quaternions, where the real part of θ is the angle and the vector part is the direction
phase :: inline (theta: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    if is_non_real_number(T) {
        U = T;
    } else { // return complex number by default
        U = get_complex_type_info_from_intrinsic_type_info(T);
    }
} {
    #if #run is_quaternion_number(T) {
        return exp(theta.r * qtn_vec(theta));
    } else #if #run is_complex_number(T) {
        return Math.exp(-theta.imag) * phase(theta.real);
    } else {
        return cpx(Math.cos(cast(float64)theta), Math.sin(cast(float64)theta));
    }
}


// ? mag * exp(iθ) where m,θ ∈ ℂ,ℝ
// ! not sure if UpCast gives the right return value here...
polar :: inline (mag: $T, theta: $U) -> UpCast(T,U).T {
    return mag * phase(theta);
}




// ? https://www.wikiwand.com/en/Trigonometric_functions
// ? trigonometric functions
sin :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        if qtn_vec_is_zero(n) then return qtn(Math.sin(n.r));

        n_vec := qtn_vec(n);
        n_vec_sign := sign(n_vec);

        return ncast(T, -0.5 * n_vec_sign * (exp(n*n_vec_sign) - exp(-n*n_vec_sign)));

    } else #if #run is_complex_number(T) {
        return cpx( Math.sin(n.real)*rcosh(n.imag), Math.cos(n.real)*rsinh(n.imag));

    } else
        return Math.sin(n);
}

cos :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        if qtn_vec_is_zero(n) then return qtn(Math.cos(n.r));

        n_vec := qtn_vec(n);
        n_vec_sign := sign(n_vec);

        return ncast(T, 0.5 * (exp(n*n_vec_sign) + exp(-n*n_vec_sign)));

    } else #if #run is_complex_number(T) {
        return cpx( Math.cos(n.real)*rcosh(n.imag), -Math.sin(n.real)*rsinh(n.imag) );
    } else
        return Math.cos(n);
}

tan :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return sin(n)/cos(n);
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        r := cpx_r(X);
        tmp := exp(cpx(0.0, 2.0)*n);
        return -i * (tmp - r)/(tmp + r);
    } else 
        return Math.tan(n);
}

cot :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return cos(n)/sin(n);
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        r := cpx_r(X);
        tmp := exp(cpx(0.0, 2.0)*n);
        return i * (tmp + r)/(tmp - r);
    } else 
        return 1.0/Math.tan(n);
}

sec :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return inv(cos(n));
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        return sech(-i * n);
    } else {
        return 1.0/cos(n);
    }
}

csc :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return inv(sin(n));
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        return -i * csch(-i*n);
    } else {
        return 1.0/sin(n);        
    }
}


// ? https://www.wikiwand.com/en/Trigonometric_functions
// ? https://mathworld.wolfram.com/
// ? quaterions: see source at top of file
// ? inverse trigonometric functions
asin :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        nvs := sign(qtn_vec(n));
        return nvs * asinh(n * nvs);
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        r := cpx_r(X);
        return -i * log(i*n + sqrt(r - n*n, T));
    } else 
        return Math.asin(n);
}

acos :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return sign(qtn_vec(n)) * acosh(n);
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        r := cpx_r(X);
        return cpx(Math.PI * 0.5) + i * log(i*n + sqrt(r - n*n, T));
    } else 
        return Math.acos(n);
}

atan :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        nvs := sign(qtn_vec(n));
        return nvs * atanh(n * nvs);
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        r := cpx_r(X);
        return cpx(0.0, 0.5) * (log(r-i*n) - log(r + i*n));
    } else 
        return Math.atan(n);
}

acot :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://people.eng.unimelb.edu.au/shiltona/anion/
        return atan(inv(n)); 
        // return ncast(Math.PI*0.5, T) - atan(n); // ! not sure if this is correct!
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        return cpx(0.0, 0.5) * ( log((n - i)/n) - log((n + i)/n) );
    } else 
        return Math.PI*0.5 - Math.atan(n);
}

asec :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://people.eng.unimelb.edu.au/shiltona/anion/
        return acos(inv(n));
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        r := cpx_r(X);
        return cpx(Math.PI * 0.5) + i * log( sqrt(r - r/(n*n), T) + i/n );
    } else 
        return acos(1.0/n);
}

acsc :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://people.eng.unimelb.edu.au/shiltona/anion/
        return asin(inv(n));
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        i := cpx_i(X);
        r := cpx_r(X);
        return -i * log( sqrt(r - r/(n*n), T) + i/n );
    } else 
        return asin(1.0/n);
}



// https://www.wikiwand.com/en/Hyperbolic_functions
// hyperbolic functions
sinh :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return ncast(T, 0.5 * (exp(n) + exp(-n)));
    } else #if #run is_complex_number(T) {
        return cpx( rsinh(n.real)*Math.cos(n.imag), rcosh(n.real)*Math.sin(n.imag) );
    } else 
        return rsinh(n);
}

cosh :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return ncast(T, 0.5 * (exp(n) - exp(-n)));
    } else #if #run is_complex_number(T) {
        return cpx( rcosh(n.real)*Math.cos(n.imag), rsinh(n.real)*Math.sin(n.imag) );
    } else { 
        return rcosh(n);
    }
}

tanh :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return sinh(n)/cosh(n);
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        r := cpx_r(X);
        tmp := exp(2.0 * n);
        return (tmp - r)/(tmp + r);
    } else 
        return rtanh(n);
}

coth :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return cosh(n)/sinh(n);
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        r := cpx_r(X);
        tmp := exp(2.0 * n);
        return (tmp + r)/(tmp - r);
    } else 
        return rcoth(n);
}

sech :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return inv(cosh(n));
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        r := cpx_r(X);
        return 2.0 * exp(n) / (exp(2.0*n) + r);
    } else 
        return rsech(n);
}

csch :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return inv(sinh(n));
    } else #if #run is_complex_number(T) {
        X :: #run get_number_inner_type(T);
        r := cpx_r(X);
        return 2.0 * exp(n) / (exp(2.0*n) - r);
    } else 
        return rcsch(n);
}



// https://www.wikiwand.com/en/Hyperbolic_functions
// inverse hyperbolic functions
// from what I've seen, the same formulas work for quaternions, octonions as well
asinh :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return log(n + sqrt(n*n + one(T), T));
    } else 
        return rasinh(n);
}

acosh :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return log( n + sqrt(n*n - 1, T));
    } else 
        return racosh(n);
}

atanh :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        // ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
        // ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf
        // ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
        return ncast(T, 0.5 * ( log(one(T)+n) - log(one(T)-n) ));
    } else 
        return ratanh(n);
}

acoth :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://people.eng.unimelb.edu.au/shiltona/anion/
        return atanh(inv(n));
    } else #if #run is_complex_number(T) {
        return ncast(T, 0.5 * ( log(one(T)+one(T)/n) - log(one(T)-one(T)/n) ));
    } else 
        return racoth(n);
}

asech :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://people.eng.unimelb.edu.au/shiltona/anion/
        return acosh(inv(n));
    } else #if #run is_complex_number(T) || is_quaternion_number(T) {
        return log( sqrt(inv(n)-one(T), T) * sqrt(inv(n)+n, T) + inv(n) );
    } else 
        return rasech(n);
}

acsch :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://people.eng.unimelb.edu.au/shiltona/anion/
        return asinh(inv(n));
    } else #if #run is_complex_number(T) || is_quaternion_number(T) {
        return log( sqrt(one(T)+one(T)/(n*n), T) + inv(n) );
    } else 
        return racsch(n);
}








#scope_module



// ? https://www.wikiwand.com/en/Hyperbolic_functions
// ? hyperbolic functions
rsinh :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    return 0.5 * (Math.exp(x) - Math.exp(-x));
}

rcosh :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    return 0.5 * (Math.exp(x) + Math.exp(-x));
}

rtanh :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    tmp := Math.exp(2.0 * x);
    return (tmp - 1.0) / (tmp + 1.0);
}

rcoth :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    tmp := Math.exp(2.0 * x);
    return (tmp + 1.0) / (tmp - 1.0);
}

rsech :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    return 2.0 * Math.exp(x) / (Math.exp(2.0*x) + 1.0);
}

rcsch :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    return 2.0 * Math.exp(x) / (Math.exp(2.0*x) - 1.0);
}

// ? https://www.wikiwand.com/en/Hyperbolic_functions
// ? inverse hyperbolic functions
rasinh :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    return Math.log( x + Math.sqrt(x*x + 1.0) );
}

racosh :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    #if CHECKS {
        assert( x >= 1.0);
    }
    return Math.log( x + Math.sqrt(x*x - 1.0) );
}

ratanh :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    #if CHECKS {
        assert( Math.abs(x) < 1.0 );
    }
    return 0.5 * Math.log( (1.0+x)/(1.0-x) );
}

racoth :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    #if CHECKS {
        assert( Math.abs(x) > 1.0 );
    }
    return 0.5 * Math.log( (x+1.0)/(x-1.0) );
}

rasech :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    #if CHECKS {
        assert( x > 0.0 && x <= 1.0 );
    }
    return Math.log( (1.0 + Math.sqrt(1.0 - x*x))/x );
}

rasch :: inline (x: $T) -> T #modify {
    if !is_intrinsic_number(T) {T = null; return;}
} {
    #if CHECKS {
        assert( x != 0.0 );
    }
    return Math.log( 1.0/x + Math.sqrt(1.0/(x*x) + 1.0) );
}






test_elementary :: () {
    println_push("Elementary", color = .FG_WHITE);

    aa :: (a: $T, b: T) {
        assert(approx(a,b), "% != %", a, b);
    }

    {
        println_push("pow", color = .FG_GREEN);
        assert( pow(2.0, 2.0, float64) == 4);
        assert( approx(pow(I, 2.0, C64), -R) );

        println("%", str(pow(2.0, R+I, C64)));
        assert( approx(pow(2.0, R + I, C64), cpx(1.538477803, 1.277922553)) );

        assert( approx(
            pow(2.0 - 3.0*I, -0.5 + 2.5*I, C64),
            cpx(-5.219991933, -3.243591148)
        ));
    }

    {
        println_push("trigonometric functions real", color = .FG_GREEN);

        println("sin");
        aa( sin(2.0), 0.9092974268 );
        println("cos");
        aa( cos(2.0), -0.4161468365);
        println("tan");
        aa( tan(2.0), -2.185039863);
        println("cot");
        aa( cot(2.0), -0.4576575544);
        println("sec");
        aa( sec(2.0), -2.402997962);
        println("csc");
        aa( csc(2.0), 1.099750170);

        println("asin");
        println("acos");
        println("atan");
        aa( atan(2.0), 1.107148718);
        println("acot");
        aa( acot(2.0), 0.4636476090);
        println("asec");
        aa( asec(2.0), 1.047197551);
        println("acsc");
        aa( acsc(2.0), 0.5235987758);

        // TODO test hyperbolic functions as well...
    }

    {
        println_push("trigonometric functions complex", color = .FG_GREEN);

        println("sin");
        aa( sin(I)  ,  1.175201194 * I);
        println("cos");
        aa( cos(I)  , 1.543080635 * R);
        println("tan");
        aa( tan(I)  , 0.7615941560 * I);
        println("cot");
        aa( cot(I)  , -1.313035285 * I);
        println("sec");
        aa( sec(I)  , 0.6480542737 * R);
        println("csc");
        aa( csc(I)  , -0.8509181282 * I );

        println("asin");
        aa( asin(2.0*R), cpx(1.570796327, -1.316957897) );
        aa( asin(I)  , 0.8813735869 * I);
        println("acos");
        aa( acos(2.0*R), 1.316957897 * I);
        aa( acos(I)  , cpx(1.570796327, -0.8813735869));
        println("atan");
        aa( atan(2.0*I)  , cpx(1.570796327, 0.5493061443));
        println("acot");
        aa( acot(2.0*I)  , -0.5493061443 * I);
        println("asec");
        aa( asec(I)  , cpx(1.570796327, 0.8813735869));
        println("acsc");
        aa( acsc(I)  , -0.8813735869 * I );

        // TODO test hyperbolic functions as well...
    }

    {
        println_push("trigonometric functions quaternion", color = .FG_GREEN);

        q := q64(-3, 4, -1, 9);

        // ? the following results were calculated with Mathematica

        println("sin");
        aa( sin(q)  ,  q64(-1405.576937229209, -3984.23443158488, 996.0586078962197, -8964.527471065978));
        println("cos");
        aa( cos(q)  , q64(-9860.477195134747, 567.9388450199349, -141.98471125498367, 1277.8624012948533));
        println("tan");
        aa( tan(q)  , q64(0.0000000014082781322253716, 0.40406101586549537, -0.10101525396637381, 0.9091372856973644));
        println("cot");
        aa( cot(q)  , q64(0.0000000014082781458556475, -0.40406101977627334, 0.10101525494406831, -0.909137294496615));
        println("sec");
        aa( sec(q)  , q64(-0.00009939530553404795, -0.000005724921209012474, 0.000001431230302253118, -0.000012881072720278064));
        println("csc");
        aa( csc(q)  , q64(-0.000014168457338822158, 0.0000401617685069994, -0.000010040442126749847, 0.00009036397914074864));

        println("asin");
        aa( asin(q)  , -q64(-0.2929592549548952, 1.2249117226255266, -0.30622793065638154, 2.7560513759074343));
        println("acos");
        aa( acos(q)  , q64(1.8637555817497917, -1.2249117226255266, 0.30622793065638154, -2.7560513759074343));
        println("atan");
        aa( atan(q)  , -q64(-1.542524608322488, 0.037460519419390825, -0.009365129854847703, 0.08428616869362934));
        println("acot");
        aa( acot(q)  , -q64(-0.02827171847240869, -0.037460519419390825, 0.009365129854847703, -0.08428616869362934));
        println("asec");
        aa( asec(q)  , q64(1.5987180145786841, 0.03734456556349926, -0.009336141390874811, 0.08402527251787331));
        println("acsc");
        aa( acsc(q)  , -q64(-0.027921687783787664, -0.03734456556349926, 0.009336141390874811, -0.08402527251787331));

        // TODO test hyperbolic functions as well...
    }
}