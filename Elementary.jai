// maybe return int?
sign :: inline (x: $T) -> T {
    #if #run is_complex_number(T) {
        if x.real == 0.0 then return R * sign(x.imag);
        return R * sign(x.real);
    } else {
        if x >= 0 then return one(T);
        return -one(T);
    }
}


factorial :: (n: int, $T: Type) -> T {
    res := one(T);
    for 2..n res *= it;
    return res;
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




real :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else {
        U = T;
    }
} {
    #if #run is_complex_number(T) {
        return n.real;
    } else {
        return n;
    }
}

imag :: inline (n: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else {
        U = T;
    }
} {
    #if #run is_complex_number(T) {
        return n.imag;
    } else {
        return n;
    }
}



conjugate :: inline (n: $T) -> T #modify {
    if !is_number(T) { T = null; return; }    
} {
    #if #run is_complex_number(T) {
        cpx(n.real, -n.imag);
    } else {
        return n;
    }
}


// |z|²
abs_sq :: inline (n: $T) -> $U #modify {
    if !is_number(T) {T = null; return;}
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else {
        U = T;
    } 
} {
    #if #run is_complex_number(T) {
        return n.real * n.real + n.imag * n.imag;
    } else {
        return n * n;
    }
}

abs :: inline (a: $T) -> $U #modify {
    if !is_number(T) {T = null; return;}
    if is_complex_number(T) {
        U = get_complex_inner_type(T);
    } else {
        U = T;
    }
} {
    return Math.sqrt(real(abs_sq(a)));
}


// given z = |z| exp(iθ), whis calculates θ.
arg :: inline (a: $T) -> float64 {
    #if #run is_complex_number(T) {
        // return carg(a);
        return Math.atan2(a.imag, a.real);
    } else 
        return 0.0;
}




// exp(z) where z ∈ ℂ,ℝ
exp :: inline (z: $T) -> T {
    #if !#run is_complex_number(T) {
        return Math.exp(z);
    } else {
       return Math.exp(z.real) * phase(z.imag);
    }
}


// log(z) where z ∈ ℂ,ℝ
log :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx(Math.log(abs(z)), Math.atan2(z.imag, z.real));
    } else {
        return cast(T)Math.log(cast(float64)z);
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
pow :: inline (z: $T, power: $U, $type: Type) -> type {
    #if #run is_complex_number(T) {
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
sqrt :: inline(z: $T, $type: Type) -> type {
    #if !#run is_complex_number(T) {
        #if #run !is_complex_number(type) {
            return Math.sqrt(z);
        } else
            return sqrt(cpx(z));
    } else {
        return Math.sqrt(abs(z)) * phase(arg(z) * 0.5);
    }
}




// exp(iθ) where θ ∈ ℂ,ℝ
phase :: inline (theta: $T) -> $U #modify {
    if !is_number(T) {U = null; return;}
    if is_complex_number(T) {
        U = T;
    } else {
        U = get_complex_type_info_from_intrinsic_type_info(T);
    }
} {
    #if #run !is_complex_number(T) {
        // return cphase(cast(float64) theta);
        return cpx(Math.cos(cast(float64)theta), Math.sin(cast(float64)theta));
    } else {
        return Math.exp(-theta.imag) * phase(theta.real);
    }
}


// mag * exp(iθ) where m,θ ∈ ℂ,ℝ
polar :: inline (mag: $T, theta: $U) -> UpCast(T,U).T {
    return mag * phase(theta);
}




// https://www.wikiwand.com/en/Trigonometric_functions
// trigonometric functions
sin :: inline (x: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx( Math.sin(x.real)*rcosh(x.imag), Math.cos(x.real)*rsinh(x.imag));
    } else
        return Math.sin(x);
}

cos :: inline (x: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx( Math.cos(x.real)*rcosh(x.imag), -Math.sin(x.real)*rsinh(x.imag) );
    } else
        return Math.cos(x);
}

tan :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        tmp := exp(cpx(0.0, 2.0)*z);
        return -I * (tmp - R)/(tmp + R);
    } else 
        return Math.tan(z);
}

cot :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        tmp := exp(cpx(0.0, 2.0)*z);
        return I * (tmp + R)/(tmp - R);
    } else 
        return 1.0/Math.tan(z);
}

sec :: inline (z: $T) -> T {
    return 1.0/cos(z);
}

csc :: inline (z: $T) -> T {
    return 1.0/sin(z);
}


// https://www.wikiwand.com/en/Trigonometric_functions
// https://mathworld.wolfram.com/
// inverse trigonometric functions
asin :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return -I * log(I*z + sqrt(R - z*z, T));
    } else 
        return Math.asin(z);
}

acos :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx(Math.PI * 0.5) + I * log(I*z + sqrt(R - z*z, T));
    } else 
        return Math.acos(z);
}

atan :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx(0.0, 0.5) * (log(R-I*z) - log(R + I*z));
    } else 
        return Math.atan(z);
}

acot :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx(0.0, 0.5) * ( log((z - I)/z) - log((z + I)/z) );
    } else 
        return Math.PI*0.5 - Math.atan(z);
}

asec :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx(Math.PI * 0.5) + I * log( sqrt(R - R/(z*z), T) + I/z );
    } else 
        return acos(1.0/z);
}

acsc :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return -I * log( sqrt(R - R/(z*z), T) + I/z );
    } else 
        return asin(1.0/z);
}



// https://www.wikiwand.com/en/Hyperbolic_functions
// hyperbolic functions
sinh :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx( rsinh(z.real)*Math.cos(z.imag), rcosh(z.real)*Math.sin(z.imag) );
    } else 
        return rsinh(z);
}

cosh :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx( rcosh(z.real)*Math.cos(z.imag), rsinh(z.real)*Math.sin(z.imag) );
    } else 
        return rcosh(z);
}

tanh :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        tmp := exp(2.0 * z);
        return (tmp - R)/(tmp + R);
    } else 
        return rtanh(z);
}

coth :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        tmp := exp(2.0 * z);
        return (tmp + R)/(tmp - R);
    } else 
        return rcoth(z);
}

sech :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return 2.0 * exp(z) / (exp(2.0*z) + R);
    } else 
        return rsech(z);
}
csch :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return 2.0 * exp(z) / (exp(2.0*z) - R);
    } else 
        return rcsch(z);
}



// https://www.wikiwand.com/en/Hyperbolic_functions
// inverse hyperbolic functions
asinh :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return log(z + sqrt(R + z*z, T));
    } else 
        return rasinh(z);
}

acosh :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return log( z + sqrt(z - R, T) * sqrt(z + R, T) );
    } else 
        return racosh(z);
}

atanh :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx(0.5) * ( log(R+z) - log(R-z) );
    } else 
        return ratanh(z);
}

acoth :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return cpx(0.5) * ( log(R+R/z) - log(R-R/z) );
    } else 
        return racoth(z);
}

asech :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return log( sqrt(R/z-R, T) * sqrt(R/z+z, T) + r/z );
    } else 
        return rasech(z);
}

acsch :: inline (z: $T) -> T {
    #if #run is_complex_number(T) {
        return log( sqrt(R+R/(z*z), T) + r/z );
    } else 
        return racsch(z);
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
    println_push("Elementary");

    aa :: (a: $T, b: T) {
        assert(approx(a,b), "% != %", a, b);
    }

    {
        println_push("pow");
        assert( pow(2.0, 2.0, float64) == 4);
        assert( approx(pow(I, 2.0, Complex(float64)), -R) );

        println("%", str(pow(2.0, R+I, Complex(float64))));
        assert( approx(pow(2.0, R + I, Complex(float64)), cpx(1.538477803, 1.277922553)) );

        assert( approx(
            pow(2.0 - 3.0*I, -0.5 + 2.5*I, Complex(float64)),
            cpx(-5.219991933, -3.243591148)
        ));
    }

    {
        println_push("trigonometric functions");

        println("sin");
        aa( sin(2.0), 0.9092974268 );
        aa( sin(I)  ,  1.175201194 * I);
        println("cos");
        aa( cos(2.0), -0.4161468365);
        aa( cos(I)  , 1.543080635 * R);
        println("tan");
        aa( tan(2.0), -2.185039863);
        aa( tan(I)  , 0.7615941560 * I);
        println("cot");
        aa( cot(2.0), -0.4576575544);
        aa( cot(I)  , -1.313035285 * I);
        println("sec");
        aa( sec(2.0), -2.402997962);
        aa( sec(I)  , 0.6480542737 * R);
        println("csc");
        aa( csc(2.0), 1.099750170);
        aa( csc(I)  , -0.8509181282 * I );

        println("asin");
        aa( asin(2.0*R), cpx(1.570796327, -1.316957897) );
        aa( asin(I)  , 0.8813735869 * I);
        println("acos");
        aa( acos(2.0*R), 1.316957897 * I);
        aa( acos(I)  , cpx(1.570796327, -0.8813735869));
        println("atan");
        aa( atan(2.0), 1.107148718);
        aa( atan(2.0*I)  , cpx(1.570796327, 0.5493061443));
        println("acot");
        aa( acot(2.0), 0.4636476090);
        aa( acot(2.0*I)  , -0.5493061443 * I);
        println("asec");
        aa( asec(2.0), 1.047197551);
        aa( asec(I)  , cpx(1.570796327, 0.8813735869));
        println("acsc");
        aa( acsc(2.0), 0.5235987758);
        aa( acsc(I)  , -0.8813735869 * I );

        // TODO test hyperbolic functions as well...
    }
}