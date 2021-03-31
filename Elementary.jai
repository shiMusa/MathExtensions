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



arg :: inline (a: $T) -> float64 {
    #if T == Complex64 {
        // return carg(a);
        return Math.atan2(a.imag, a.real);
    } else 
        return 0.0;
}





exp :: inline (e: $T) -> T {
    #if T != Complex64 {
        return Math.exp(e);
    } else {
        return Math.exp(e.real) * phase(e.imag);
    }
}


log :: inline (z: $T) -> T {
    #if T == Complex64 {
        return c64(Math.log(abs(z)), Math.atan2(z.imag, z.real));
    } else {
        return cast(T)Math.log(cast(float64)z);
    }
}





phase :: inline (theta: $T) -> Complex64 {
    #if T != Complex64 {
        // return cphase(cast(float64) theta);
        return c64(Math.cos(cast(float64)theta), Math.sin(cast(float64)theta));
    } else {
        return Math.exp(-theta.imag) * phase(theta.real);
    }
}


polar :: inline (mag: $T, theta: $U) -> Complex64 {
    return mag * phase(theta);
}