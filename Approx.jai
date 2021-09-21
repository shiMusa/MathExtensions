DELTA :: 1.0e-6;

// approx :: inline (a: $A/VectorType, b: $B/VectorType) -> bool {
//     return norm_2(a - b) < DELTA;
// }

// approx :: inline (a: $A/MatrixType, b: $B/MatrixType) -> bool {
//     return norm_frobenius(a - b) < DELTA;
// }

approx :: inline (a: $T, b: $U) -> bool #modify {
    if !is_number(T) || !is_number(U) {
        U = null; // reject non-scalars
        return;
    }
} {
    return abs(a - b) < DELTA;
}