// ? LAPACK *lag2z()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clag2z
// *LAG2Z converts a complex single precision matrix to a complex
// double precision matrix.

clag2z :: (
    $$m: int,
    $$n: int,
    sa: $SA/ArrayType2D,
    $$ldsa: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    info: *int
) #modify {
    return SA.data_type == C32 &&
           A.data_type == C64;
} {
    info.* = 0;
    for j:1..n {
        for i:1..m {
            v := fget(sa, i, j);
            fset(a, i, j, cpx(ncast(float64, real(v)), ncast(float64, imag(v))));
        }
    }
}
