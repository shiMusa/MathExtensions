






array_view_range :: inline (arr: []$T, from: s64, to: s64) -> []T {
    #if CHECKS {
        assert(to < arr.count, "Index out of bounds: to = %, arr.count = %", to, arr.count);
        assert(from < to, "Negative range: from = %, to = %", from, to);
    }
    res := arr;
    res.data += from;
    res.count = (to - from);
    return res;
}

array_view_count :: inline (arr: []$T, from: s64, count: s64) -> []T {
    #if CHECKS {
        assert(from + count <= arr.count, "Index out of bounds: from+count = %, arr.count = %", from+count, arr.count);
        assert(count > 0, "Negative/zero range: count = %", count);
    }

    res := arr;
    res.data += from;
    res.count = count;
    return res;
}





f64 :: inline (a: $T) -> float64 {
    return cast(float64) a;
}






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
    #if CHECKS {
        assert(choose <= from, "Out of range: from = %, choose = %", from, choose);
    }

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

    {
        a := float64.[0,1,2,3,4,5,6];
        v := array_view_range(a, 2, 4);
        print("%\n", a);
        print("%\n", v);
    }
}