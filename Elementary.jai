






// ? a great source for quaternion functions is
// ? "New Aspects on Elementary Functions in the Context of Quaternionic Analysis"
// ? in CUBO A Mathematical Journal, Vol.14, No 01, (93–110). March 2012
// ? https://pdfs.semanticscholar.org/8ef1/761d8402f56c8a6258257e6b166d4c38a314.pdf









// maybe return int?
sign :: inline (x: $T) -> T #modify {
    if !is_number(T) {T = null; return;}    
} {
    #if #run is_quaternion_number(T) || is_complex_number(T) {
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
        return cpx_inv(n);
    } else return 1/n;
}




// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// 2.3.1.10 Binomial Coefficients, page 55
binomial :: inline (from: $T, choose: T) -> T {
    #if CHECKS {
        assert(choose <= from, "Out of range: from = %, choose = %", from, choose);
    }

    // assert(from >= choose); // TODO disable for release build?!
    if choose == zero(T) || from == choose {
        return one(T);
    }

    return binomial(from - one(T), choose) + binomial(from - one(T), choose - one(T));
}



quat_r :: real;
real :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else if is_quaternion_number(T) {
        U = get_quaternion_inner_type(T);    
    } else {
        U = T;
    }
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
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else if is_quaternion_number(T) {
        U = get_quaternion_inner_type(T);    
    } else {
        U = T;
    }
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
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else if is_quaternion_number(T) {
        U = get_quaternion_inner_type(T);    
    } else {
        U = T;
    }
} {
    #if #run is_quaternion_number(T) {
        return n.j;
    } else {
        return zero(T);
    }
}

quat_k :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else if is_quaternion_number(T) {
        U = get_quaternion_inner_type(T);    
    } else {
        U = T;
    }
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
        return cpx_inv(n);
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


// |z|²
abs_sq :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else if is_quaternion_number(T) {
        U = get_quaternion_inner_type(T);    
    } else {
        U = T;
    }
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
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else if is_quaternion_number(T) {
        U = get_quaternion_inner_type(T);    
    } else {
        U = T;
    }
} {
    return Math.sqrt(abs_sq(a));
}


// given n = |n| exp(iθ), whis calculates θ.
// for quaternions: n = |n| exp(iθ eᵥ), where eᵥ = v/|v| with n = n.r + v
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




// exp(z) where z ∈ ℂ,ℝ
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


// log(z) where z ∈ ℂ,ℝ
log :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        // ? https://www.wikiwand.com/en/Quaternion
        q := abs(n);
        return Math.log(q) + sign(n) * Math.acos(n.r/q);
    } else #if #run is_complex_number(T) {
        return cpx(Math.log(abs(n)), Math.atan2(n.imag, n.real));
    } else {
        return cast(T)Math.log(cast(float64)n);
    }
}



// z^p, where z,p ∈ ℂ,ℝ
// z,p ∈ ℝ: build in pow(z,p)
// z ∈ ℂ, p ∈ ℝ: |z|^p exp(i θ p)
// z ∈ ℝ, p ∈ ℂ: pow(z,Re{p}) * exp(i ln(z) Im{p})
// z,p ∈ ℂ: |z|^p * exp(iθ p) referring back to above versions of exp and pow(ℝ,ℂ)
//
// there might be better algorithms calculating this...
// but it's basically the same as in the Gnu Scientific Library (GSL)
// see complex/math.c:
//     gsl_complex gsl_complex_pow (gsl_complex a, gsl_complex b)
//
// ! the `type` thingy is kinda akward...
pow :: inline (z: $T, power: $U, $type: Type) -> type {
    #if #run is_quaternion_number(T) || is_quaternion_number(U) {
        return exp(log(z)*power); // yeah, this works... but there might be a better way to do it!
    } else #if #run is_complex_number(T) {
        #if #run is_complex_number(U) {
            return pow(abs(z), power, UpCast(T,U).T) * phase(arg(z) * power);
        } else 
            return Math.pow(abs(z), power) * phase(arg(z) * power);
    } else {
        #if #run is_complex_number(U) {
            return Math.pow(z, power.real) * phase(Math.log(z) * power.imag);
        } else
            return Math.pow(z, power);
    }
}


// in the hopes sqrt is faster than pow
// need to wait for compiler fix to get rid of `$type: Type`...
// ^ I might be able to do that now!
sqrt :: inline(z: $T, $type: Type) -> type {
    #if #run is_quaternion_number(T) {
        return pow(z, 0.5, type);
    } else #if #run is_complex_number(T) {
        return Math.sqrt(abs(z)) * phase(arg(z) * 0.5);
    } else {
        #if #run is_complex_number(type) {
            return sqrt(cpx(z), type);
        } else
            return Math.sqrt(z);
    }
}




// exp(iθ) where θ ∈ ℂ,ℝ
// is also implemented for quaternions, where the real part of θ is the angle and the vector part is the direction
phase :: inline (theta: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    if is_quaternion_number(T) || is_complex_number(T) {
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


// mag * exp(iθ) where m,θ ∈ ℂ,ℝ
// ! not sure if UpCast gives the right return value here...
polar :: inline (mag: $T, theta: $U) -> UpCast(T,U).T {
    return mag * phase(theta);
}




// https://www.wikiwand.com/en/Trigonometric_functions
// trigonometric functions
sin :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return qtn_i(T) * sinh(-qtn_i(T) * n);
    } else #if #run is_complex_number(T) {
        return cpx( Math.sin(n.real)*rcosh(n.imag), Math.cos(n.real)*rsinh(n.imag));
    } else
        return Math.sin(n);
}

cos :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return cosh(-qtn_i(T) * n);
    } else #if #run is_complex_number(T) {
        return cpx( Math.cos(n.real)*rcosh(n.imag), -Math.sin(n.real)*rsinh(n.imag) );
    } else
        return Math.cos(n);
}

tan :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return qtn_i(T) * tanh(-qtn_i(T) * n);
    } else #if #run is_complex_number(T) {
        tmp := exp(cpx(0.0, 2.0)*n);
        return -I * (tmp - R)/(tmp + R);
    } else 
        return Math.tan(n);
}

cot :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return -qtn_i(T) * coth(-qtn_i(T) * n);
    } else #if #run is_complex_number(T) {
        tmp := exp(cpx(0.0, 2.0)*n);
        return I * (tmp + R)/(tmp - R);
    } else 
        return 1.0/Math.tan(n);
}

sec :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return sech(-qtn_i(T) * n);
    } else #if #run is_complex_number(T) {
        return sech(-I * n);
    } else {
        return 1.0/cos(n);
    }
}

csc :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return -qtn_i(T) * csch(-qtn_i(T)*n);
    } else #if #run is_complex_number(T) {
        return -I * csch(-I*n);
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
        return sign(qtn_vec(n)) * asinh(n * sgn(qtn_vec(n)));
    } else #if #run is_complex_number(T) {
        return -I * log(I*n + sqrt(R - n*n, T)); // TODO make I and R generic over type
    } else 
        return Math.asin(n);
}

acos :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return sgn(qtn_vec(n)) * acosh(n);
    } else #if #run is_complex_number(T) {
        return cpx(Math.PI * 0.5) + I * log(I*n + sqrt(R - n*n, T)); // TODO make I and R generic over type
    } else 
        return Math.acos(n);
}

atan :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return sign(qtn_vec(n)) * atanh(n * sign(qtn_vec(n)));
    } else #if #run is_complex_number(T) {
        return cpx(0.0, 0.5) * (log(R-I*n) - log(R + I*n)); // TODO make I and R generic over type
    } else 
        return Math.atan(n);
}

acot :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return ncast(Math.PI*0.5, T) - atan(n); // ! not sure if this is correct!
    } else #if #run is_complex_number(T) {
        return cpx(0.0, 0.5) * ( log((n - I)/n) - log((n + I)/n) ); // TODO make I and R generic over type
    } else 
        return Math.PI*0.5 - Math.atan(n);
}

asec :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return acos(inv(n)); // ! not sure if this is correct!
    } else #if #run is_complex_number(T) {
        return cpx(Math.PI * 0.5) + I * log( sqrt(R - R/(n*n), T) + I/n ); // TODO make I and R generic over type
    } else 
        return acos(1.0/n);
}

acsc :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return asin(inv(n)); // ! not sure if this is correct!
    } else #if #run is_complex_number(T) {
        return -I * log( sqrt(R - R/(n*n), T) + I/n ); // TODO make I and R generic over type
    } else 
        return asin(1.0/n);
}



// https://www.wikiwand.com/en/Hyperbolic_functions
// hyperbolic functions
sinh :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return 0.5 * (exp(n) + exp(-n));
    } else #if #run is_complex_number(T) {
        return cpx( rsinh(n.real)*Math.cos(n.imag), rcosh(n.real)*Math.sin(n.imag) );
    } else 
        return rsinh(n);
}

cosh :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return 0.5 * (exp(n) - exp(-n));
    } #if #run is_complex_number(T) {
        return cpx( rcosh(n.real)*Math.cos(n.imag), rsinh(n.real)*Math.sin(n.imag) );
    } else 
        return rcosh(n);
}

tanh :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return sinh(n)/cosh(n);
    } else #if #run is_complex_number(T) {
        tmp := exp(2.0 * n);
        return (tmp - R)/(tmp + R);
    } else 
        return rtanh(n);
}

coth :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return cosh(n)/sinh(n);
    } else #if #run is_complex_number(T) {
        tmp := exp(2.0 * n);
        return (tmp + R)/(tmp - R);
    } else 
        return rcoth(n);
}

sech :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return qtn_inv(cosh(n));
    } else #if #run is_complex_number(T) {
        return 2.0 * exp(n) / (exp(2.0*n) + R);
    } else 
        return rsech(n);
}

csch :: inline (n: $T) -> T {
    #if #run is_quaternion_number(T) {
        return qtn_inv(sinh(n));
    } else #if #run is_complex_number(T) {
        return 2.0 * exp(n) / (exp(2.0*n) - R);
    } else 
        return rcsch(n);
}



// https://www.wikiwand.com/en/Hyperbolic_functions
// inverse hyperbolic functions
// from what I've seen, the same formulas work for quaternions, octonions as well
asinh :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        return log(n + sqrt(one(T) + n*n, T));
    } else 
        return rasinh(n);
}

acosh :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        return log( n + sqrt(n - one(T), T) * sqrt(n + one(T), T) );
    } else 
        return racosh(n);
}

atanh :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        return ncast(0.5 * ( log(one(T)+n) - log(one(T)-n) ), T);
    } else 
        return ratanh(n);
}

acoth :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        return ncast(0.5 * ( log(one(T)+one(T)/n) - log(one(T)-one(T)/n) ), T);
    } else 
        return racoth(n);
}

asech :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        return log( sqrt(inv(n)-one(T), T) * sqrt(inv(n)+n, T) + inv(n) );
    } else 
        return rasech(n);
}

acsch :: inline (n: $T) -> T {
    #if #run is_complex_number(T) || is_quaternion_number(T) {
        return log( sqrt(one(T)+one(T)/(n*n), T) + inv(n) );
    } else 
        return racsch(n);
}








#scope_module



// https://www.wikiwand.com/en/Hyperbolic_functions
// hyperbolic functions
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

// https://www.wikiwand.com/en/Hyperbolic_functions
// inverse hyperbolic functions
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
        aa( sec(q)  , q64(-0.00009939530553404795, -0.000005724921209012474, 1.431230302253118, -0.000012881072720278064));
        println("csc");
        aa( csc(q)  , q64(-0.000014168457338822158, 0.0000401617685069994, -0.000010040442126749847, 0.00009036397914074864));

        println("asin");
        println("acos");
        println("atan");
        println("acot");
        println("asec");
        println("acsc");

        // TODO test hyperbolic functions as well...
    }
}