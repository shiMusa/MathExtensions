
// ? LAPACK *abs
// ? http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga53914429b84cb315da483e71e27ed5c1.html#ga53914429b84cb315da483e71e27ed5c1

dcabs1 :: inline (z: $T) -> float64 #modify {
    if !is_number(T) {T = null; return;}
} {
    return inline abs1(float64, z);
}
scabs1 :: inline (z: $T) -> float32 #modify {
    if !is_number(T) { T = null; return; }
} {
    return inline abs1(float32, z);
}

_abs1 :: abs1;