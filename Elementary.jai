sign :: inline (x: $T) -> T {
    return 2 * ncast(T, x >= 1) - one(T);
}


factorial :: (n: $T) -> T {
    res := 1;
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




conjugate :: inline (a: $T) -> T {
    #if T == Complex64 {
        return c64(a.real, -a.imag);
    } else
        return a;
}


abs_sq :: inline (a: $T) -> float64 {
    #if T == Complex64 {
        return a.real * a.real + a.imag * a.imag;
    } else 
        return cast(float64)(a * a);
}

abs :: inline (a: $T) -> float64 {
    return Math.sqrt(abs_sq(a));
}


// given z = |z| exp(iθ), whis calculates θ.
arg :: inline (a: $T) -> float64 {
    #if T == Complex64 {
        // return carg(a);
        return Math.atan2(a.imag, a.real);
    } else 
        return 0.0;
}




// exp(z) where z ∈ ℂ,ℝ
exp :: inline (z: $T) -> T {
    #if T != Complex64 {
        return Math.exp(z);
    } else {
       return Math.exp(z.real) * phase(z.imag);
    }
}


// log(z) where z ∈ ℂ,ℝ
log :: inline (z: $T) -> T {
    #if T == Complex64 {
        return c64(Math.log(abs(z)), Math.atan2(z.imag, z.real));
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
    #if T == Complex64 {
        #if U == Complex64 {
            return pow(abs(z), power, Complex64) * phase(arg(z) * power);
        } else 
            return Math.pow(abs(z), power) * phase(arg(z) * power);
    } else {
        #if U == Complex64 {
            return Math.pow(z, power.real) * phase(Math.log(z) * power.imag);
        } else
            return Math.pow(z, power);
    }
}


// in the hopes sqrt is faster than pow
sqrt :: inline(z: $T, $type: Type) -> type {
    #if T != Complex64 {
        #if type != Complex64 {
            return Math.sqrt(z);
        } else
            return c64(Math.sqrt(z));
    } else {
        return Math.sqrt(abs(z)) * phase(arg(z) * 0.5);
    }
}




// exp(iθ) where θ ∈ ℂ,ℝ
phase :: inline (theta: $T) -> Complex64 {
    #if T != Complex64 {
        // return cphase(cast(float64) theta);
        return c64(Math.cos(cast(float64)theta), Math.sin(cast(float64)theta));
    } else {
        return Math.exp(-theta.imag) * phase(theta.real);
    }
}


// mag * exp(iθ) where m,θ ∈ ℂ,ℝ
polar :: inline (mag: $T, theta: $U) -> Complex64 {
    return mag * phase(theta);
}




// https://www.wikiwand.com/en/Trigonometric_functions
// trigonometric functions
sin :: inline (x: $T) -> T {
    #if T == Complex64 {
        return c64( Math.sin(x.real)*rcosh(x.imag), Math.cos(x.real)*rsinh(x.imag));
    } else
        return Math.sin(x);
}

cos :: inline (x: $T) -> T {
    #if T == Complex64 {
        return c64( Math.cos(x.real)*rcosh(x.imag), -Math.sin(x.real)*rsinh(x.imag) );
    } else
        return Math.cos(x);
}

tan :: inline (z: $T) -> T {
    #if T == Complex64 {
        tmp := exp(c64(0.0, 2.0)*z);
        return -I * (tmp - R)/(tmp + R);
    } else 
        return Math.tan(z);
}

cot :: inline (z: $T) -> T {
    #if T == Complex64 {
        tmp := exp(c64(0.0, 2.0)*z);
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
    #if T == Complex64 {
        return -I * log(I*z + sqrt(R - z*z));
    } else 
        return Math.asin(z);
}

acos :: inline (z: $T) -> T {
    #if T == Complex64 {
        return c64(Math.PI * 0.5) + I * log(I*z + sqrt(R - z*z));
    } else 
        return Math.acos(z);
}

atan :: inline (z: $T) -> T {
    #if T == Complex64 {
        return c64(0.0, 0.5) * (log(R-I*z) - log(R + I*z));
    } else 
        return Math.atan(z);
}

acot :: inline (z: $T) -> T {
    #if T == Complex64 {
        return c64(0.0, 0.5) * ( log((z - I)/z) - log((z + I)/z) );
    } else 
        return Math.PI*0.5 - Math.atan(z);
}

asec :: inline (z: $T) -> T {
    #if T == Complex64 {
        return c64(Math.PI * 0.5) + I * log( sqrt(R - R/(z*z)) + I/z );
    } else 
        return acos(1.0/z);
}

acsc :: inline (z: $T) -> T {
    #if T == Complex64 {
        return -I * log( sqrt(R - R/(z*z)) + I/z );
    } else 
        return asin(1.0/z);
}



// https://www.wikiwand.com/en/Hyperbolic_functions
// hyperbolic functions
sinh :: inline (z: $T) -> T {
    #if T == Complex64 {
        return c64( rsinh(z.real)*Math.cos(z.imag), rcosh(z.real)*Math.sin(z.imag) );
    } else 
        return rsinh(z);
}

cosh :: inline (z: $T) -> T {
    #if T == Complex64 {
        return c64( rcosh(z.real)*Math.cos(z.imag), rsinh(z.real)*Math.sin(z.imag) );
    } else 
        return rcosh(z);
}

tanh :: inline (z: $T) -> T {
    #if T == Complex64 {
        tmp := exp(2.0 * z);
        return (tmp - R)/(tmp + R);
    } else 
        return rtanh(z);
}

coth :: inline (z: $T) -> T {
    #if T == Complex64 {
        tmp := exp(2.0 * z);
        return (tmp + R)/(tmp - R);
    } else 
        return rcoth(z);
}

sech :: inline (z: $T) -> T {
    #if T == Complex64 {
        return 2.0 * exp(z) / (exp(2.0*z) + R);
    } else 
        return rsech(z);
}
csch :: inline (z: $T) -> T {
    #if T == Complex64 {
        return 2.0 * exp(z) / (exp(2.0*z) - R);
    } else 
        return rcsch(z);
}



// https://www.wikiwand.com/en/Hyperbolic_functions
// inverse hyperbolic functions









#scope_module



// https://www.wikiwand.com/en/Hyperbolic_functions
// hyperbolic functions
rsinh :: inline (x: float64) -> float64 {
    return 0.5 * (Math.exp(x) - Math.exp(-x));
}

rcosh :: inline (x: float64) -> float64 {
    return 0.5 * (Math.exp(x) + Math.exp(-x));
}

rtanh :: inline (x: float64) -> float64 {
    tmp := Math.exp(2.0 * x);
    return (tmp - 1.0) / (tmp + 1.0);
}

rcoth :: inline (x: float64) -> float64 {
    tmp := Math.exp(2.0 * x);
    return (tmp + 1.0) / (tmp - 1.0);
}

rsech :: inline (x: float64) -> float64 {
    return 2.0 * Math.exp(x) / (Math.exp(2.0*x) + 1.0);
}

rcsch :: inline (x: float64) -> float64 {
    return 2.0 * Math.exp(x) / (Math.exp(2.0*x) - 1.0);
}

// https://www.wikiwand.com/en/Hyperbolic_functions
// inverse hyperbolic functions
rasinh :: inline (x: float64) -> float64 {
    return Math.log( x + Math.sqrt(x*x + 1.0) );
}

racosh :: inline (x: float64) -> float64 {
    #if CHECKS {
        assert( x >= 1.0);
    }
    return Math.log( x + Math.sqrt(x*x - 1.0) );
}

ratanh :: inline (x: float64) -> float64 {
    #if CHECKS {
        assert( Math.abs(x) < 1.0 );
    }
    return 0.5 * Math.log( (1.0+x)/(1.0-x) );
}

racoth :: inline (x: float64) -> float64 {
    #if CHECKS {
        assert( Math.abs(x) > 1.0 );
    }
    return 0.5 * Math.log( (x+1.0)/(x-1.0) );
}

rasech :: inline (x: float64) -> float64 {
    #if CHECKS {
        assert( x > 0.0 && x <= 1.0 );
    }
    return Math.log( (1.0 + Math.sqrt(1.0 - x*x))/x );
}

rasch :: inline (x: float64) -> float64 {
    #if CHECKS {
        assert( x != 0.0 );
    }
    return Math.log( 1.0/x + Math.sqrt(1.0/(x*x) + 1.0) );
}






test_elementary :: () {
    print("\n\nElementary -------------------------------\n\n");

    {
        assert( pow(2.0, 2.0, float64) == 4);
        assert( approx(pow(I, 2.0, Complex64), -R) );

        print("%\n", pow(2.0, R+I, Complex64));
        assert( approx(pow(2.0, R + I, Complex64), c64(1.538477803, 1.277922553)) );

        assert( approx(
            pow(2.0 - 3.0*I, -0.5 + 2.5*I, Complex64),
            c64(-5.219991933, -3.243591148)
        ));
    }
}