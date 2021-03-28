

// Copied over from Hash_Table.jai
next_power_of_two :: inline (x : int) -> int {
    assert(x != 0);
    p := 1;
    while x > p p += p;
    return p;
}

sign :: inline (x: $T) -> T {
    return 2 * cast(T)(x >= 1) - 1;
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
    // assert(from >= choose); // TODO disable for release build?!
    if choose == 0 || from == choose {
        return 1;
    }

    return binomial(from-1, choose) + binomial(from-1, choose-1);
}



test_utils :: () {
    print("\n\nUtils -------------------------------\n\n");
    assert(sign(2.0) == 1.0);
    assert(sign(-2.0) == -1.0);
    assert(sign(-2) == -1);

    assert(binomial(42, 42) == 1);
    assert(binomial(42, 0) == 1);
    assert(binomial(5, 3) == 10);
    assert(binomial(11, 6) == 462);
    // print("%\n", binomial(67, 30));
}