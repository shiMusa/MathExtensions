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






// conjugate :: inline (a: Complex64) -> Complex64 {
//     return c64(a.real, -a.imag);
// }

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


// abs_sq :: inline (a: Complex64) -> float64 {
//     return a.real * a.real + a.imag * a.imag;
// }

// abs :: inline (a: Complex64) -> float64 {
//     return Math.sqrt(a.real * a.real + a.imag * a.imag);
// }

arg :: inline (a: $T) -> float64 {
    #if T == Complex64 {
        return carg(a);
    } else 
        return 0.0;
}

carg :: (a: Complex64) -> float64 {

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




exp :: inline (e: $T) -> T {
    #if T != Complex64 {
        return Math.exp(e);
    } else {
        return Math.exp(e.real) * cphase(e.imag);
    }
}

phase :: inline (theta: $T) -> Complex64 {
    #if T != Complex64 {
        return cphase(cast(float64) theta);
    } else {
        return Math.exp(-theta.imag) * cphase(theta.real);
    }
}

cphase :: inline (theta: float64) -> Complex64 {
    return c64(Math.cos(theta), Math.sin(theta));
}


polar :: inline (mag: $T, theta: $U) -> Complex64 {
    return mag * phase(theta);
}

// polar :: inline (mag: float64, theta: float64) -> Complex64 {
//     return c64(mag * Math.cos(theta), mag * Math.sin(theta));
// }