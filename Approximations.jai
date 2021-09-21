











/*

These approximations use Newtons method which requires an initial guess.

The initial guess for sqrt can be given by simply dividing the exponent by 2,
or in other words, bit shift the exponent by 1.
Not sure how to do that atm in Jai...


*/









// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.1.4, page 50
approx_inverse :: (val: $T, guess: T, rtol: T, $show_steps: bool = false) -> T #modify {
    if T.type != .FLOAT {T = null; return;}    
} {
    x : T = guess + rtol + rtol;
    xnew : T = guess;
    i := 0;

    #if show_steps {
        println_push("approx_inverse (rtol)");
    }

    while Math.abs(xnew - x) > rtol {
        x = xnew;
        xnew = x * (2.0 - val * x);
        #if show_steps {
            i += 1;
            println("%: x = %, rel err = %", i, x, Math.abs(xnew - x));
        }
    }

    return x;
}

// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.1.4, page 50
approx_inverse :: (val: $T, guess: T, iterations: int, $show_steps: bool = false) -> T #modify {
    if T.type != .FLOAT {T = null; return;}    
} {
    x : T = guess;

    #if show_steps {
        println_push("approx_inverse (iterations)");
    }

    for 1..iterations {
        x = x * (2.0 - val * x);
        #if show_steps {
            println("%: x = %", it, x);
        }
    }

    return x;
}




// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.1.5, page 51
// TODO: bit shift exponent (divid by 2) for initial guess
approx_sqrt :: (val: $T, guess: T, iterations: int, $show_steps: bool = false) -> T #modify {
    if T.type != .FLOAT {T = null; return;}    
} {
    x : T = guess;

    #if show_steps {
        println_push("approx_sqrt");
    }

    for 1..iterations {
        x = 0.5 * (x + val/x);
        #if show_steps {
            println("%: x = %", it, x);
        }
    }

    return x;
}


// ugh, messy stuff

// approx_sqrt :: inline (val: float64, iterations: int, $show_steps: bool = false) -> float64 {

//     // test : s64 = 0b0_11111111111_1000000000000000000000000000000000000000000000000000;
//     test : s64 = 0b0_01111111111_0000000000000000000000000000000000000000000000000000;
//     print("test %\n", <<cast(*float64)(*test));

//     i : s64 = <<cast(*s64)(*val);
//     // mask out the exponent only
//     mask : s64 = 0b0_11111111111_0000000000000000000000000000000000000000000000000000;
//     print("%\n", <<cast(*float64)(*i));
//     i = mask & i;
//     print("%\n", <<cast(*float64)(*i));
//     i = i >> 1;
//     print("%\n", <<cast(*float64)(*i));
//     // i = mask & i;
//     guess : float64 = <<cast(*float64)(*i);

//     print("initial %, guess %\n", val, guess);

//     return approx_sqrt(val, guess, iterations, show_steps);
// }


// Scientific Computing, Vol I: Linear and nonlinear equations
// Texts in computational science and engineering 18
// Springer
// Algorithm 2.3.1.6, page 51
approx_inv_sqrt :: (val: $T, guess: T, iterations: int, $show_steps: bool = false) -> T #modify {
    if T.type != .FLOAT {T = null; return;}    
} {
    x : T = guess;

    #if show_steps {
        println_push("approx_inv_sqrt");
    }

    for 1..iterations {
        x = 0.5 * x * (3.0 - val * x * x);
        #if show_steps {
            println("%: x = %", it, x);
        }
    }

    return x;
}





test_approximations :: () {
    println_push("Approximations", color = C_FG_WHITE);

    println("% = %", 1.0/42.0, approx_inverse(42.0, 0.01, 10));
    println("% = %", 1.0/42.0, approx_inverse(42.0, 0.01, 1.0e-12, true));
    println("err = %", Math.abs(approx_inverse(42.0, 0.01, 1.0e-12) - 1.0/42.0));


    println("% = %", Math.sqrt(42.0), approx_sqrt(42.0, 6.0, 3, true));
    // println("% = %\n", Math.sqrt(42.0), approx_sqrt(42, 3, true));
    println("% = %", 1.0/Math.sqrt(42.0), approx_inv_sqrt(42.0, 0.1, 3, true));
}