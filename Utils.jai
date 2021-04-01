


#scope_module




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





f64 :: inline (num: $T) -> float64 {
    return ncast(float64, num);
}




Super :: struct (T1: Type, T2: Type) {
    T :: #run super_type(T1, T2);
}


super_type :: inline ($a: Type, $b: Type) -> Type {
    #if a == Complex64 || b == Complex64 {
        return Complex64;
    } else if a == float64 || b == float64 {
        return float64;
    } else if a == float32 || b == float32 {
        return float32;
    } else if a == s64 || b == s64 {
        return s64;
    } else
        return s32;
    
    // if a == s32 || b == s32 {
    //     return s32;
    // } else if a == s16 || b == s16 {
    //     return s16;
    // } else if a == s8 || 
}


ncast :: inline ($T: Type, num: $N) -> T {
    #if T == Complex64 {
        #if N == Complex64 {
            return num;
        } else {
            return c64(cast(float64)num);
            // res : Complex64;
            // res.real = cast(float64) num;
            // res.imag = 0.0;
            // return res;
        }
    } else {
        #if N == Complex64 {
            return cast(T) num.real;
        } else
            return cast(T) num;
    }
}

zero :: inline ($T: Type) -> T {
    return ncast(T, 0);
}

one :: inline ($T: Type) -> T {
    return ncast(T, 1);
}



DELTA :: 1.0e-6;
approx :: (a: $T, b: T) -> bool {
    // print("approx(%, %): abs(1 - b/a) = %\n", a, b, abs(one(T) - b/a));
    return abs(one(T) - b/a) < DELTA;
}




// Copied over from Hash_Table.jai
next_power_of_two :: inline (x : int) -> int {
    assert(x != 0);
    p := 1;
    while x > p p += p;
    return p;
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

    {
        assert(Complex64 == super_type(Complex64, u8));
        assert(float64 == super_type(float32, float64));
    }
}