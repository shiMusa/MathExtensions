DELTA : float64 : 1.0e-6;

// ? L2 norm
approx :: inline (a: $A/VectorType, b: $B/VectorType) -> bool {
    #if CHECKS {
        assert(dim(a) == dim(b));
    }

    sum : float64 = 0.0;
    for a {
        sum += abs(it - get(b, it_index));
    }
    return sqrt(sum, float64) < DELTA;
}


// ? Frobenius norm
approx :: inline (a: $A/MatrixType, b: $B/MatrixType) -> bool {
    #if CHECKS {
        assert(rows(a) == rows(b) && cols(a) == cols(b));
    }
    
    sum : float64 = 0.0;
    for a {
        sum += abs_sq(it - get(b, it_row, it_column));
    }
    return sqrt(sum, float64) < DELTA;
}

approx :: inline (a: $T, b: $U) -> bool #modify {
    return is_number(cast(*Type_Info)T) && is_number(cast(*Type_Info)U);
} {
    return abs(a - b) < DELTA;
}