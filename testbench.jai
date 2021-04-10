#import "Basic";

CHECKS :: true;

is_a_scalar :: inline (T: Type) -> bool {
    return T == Complex64
        || T == float64
        || T == float32
        || T == float
        || T == int
        || T == s64
        || T == s32
        || T == s16
        || T == u64
        || T == u32
        || T == u16
        || T == u8
        ;
}

Complex64 :: struct {
    real: float64;
    imag: float64;
}

c64 :: (re: float64 = 0.0, im: float64 = 0.0) -> Complex64 {
    r: Complex64 = ---;
    r.real = re;
    r.imag = im;
    return r;
}


ncast :: inline ($T: Type, num: $N) -> T {
    #if CHECKS {
        assert(is_a_scalar(T) && is_a_scalar(N), "cannot cast between non-scalars! T = %, N = %\n", T, N);
    }

    print("ncast T = %, num = %, N = %\n", T, num, N);

    #if T == N {
        return num;
    } else #if T == Complex64 {
        return c64(cast(float64)num);
    } else #if N != Complex64 {
        // print("num = %, N = %\n", num, N);
        return cast,no_check(T) num; // TODO other types?
    } else {
        assert(N != Complex64, "cannot cast complex to real number!");
        return zero(T);
    }
}

zero :: inline ($T: Type) -> T {
    return ncast(T, 0);
}

one :: inline ($T: Type) -> T {
    return ncast(T, 1);
}

main :: () {
    print("%\n", ncast(Complex64, 1));
    print("%\n", ncast(float64, 1));
    print("%\n", ncast(s64, 1));
    print("%\n", ncast(u8, 1));
    print("%\n", ncast(u8, c64(1.0)));
}