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
pow :: inline (z: $T, power: $U, $type: Type) -> type {
    #if T == Complex64 {
        #if U == Complex64 {
            return pow(abs(z), power) * phase(arg(z) * power);
        } else 
            return Math.pow(abs(z), power) * phase(arg(z) * power);
    } else {
        #if U == Complex64 {
            return Math.pow(z, power.real) * phase(Math.log(z) * power.imag);
        } else
            return Math.pow(z, power);
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